#
# Copyright 2015 by Justin MacCallum, Alberto Perez, Ken Dill
# All rights reserved
#


class SlaveReplicaExchangeRunner(object):
    '''
    This class coordinates running replica exchange on the slaves.

    '''

    def __init__(self, step, max_steps):
        self._step = step
        self._max_steps = max_steps

    @classmethod
    def from_master(cls, master):
        '''
        Initialize a new slave from a master.

        :param master: a :class:`meld.remd.master_runner.
                                 MasterReplicaExchangeRunner`
                       to serve as a template
        :return: a :class:`SlaveReplicaExchangeRunner`

        '''
        new_slave = cls(master.step, master.max_steps)
        return new_slave

    @property
    def step(self):
        return self._step

    @property
    def max_steps(self):
        return self._max_steps

    def run(self, communicator, system_runner):
        '''
        Continue running slave jobs until done.

        :param communicator: a communicator object for talking to the master
        :param system_runner: a system_runner object for actually running the
                              simulations

        '''
        my_alpha = None

        # we always minimize when we first start, either on the first
        # stage or the first stage after a restart
        minimize = True
        while self._step <= self._max_steps:
            # update simulation conditions
            new_alpha = communicator.receive_alpha_from_master()
            state = communicator.receive_state_from_master()

            my_alpha = new_alpha
            system_runner.prepare_for_timestep(my_alpha, self._step)
            state.alpha = my_alpha

            # do one round of simulation
            if minimize:
                state = system_runner.minimize_then_run(state)
                minimize = False    # we don't need to minimize again
            else:
                state = system_runner.run(state)

            # compute energies
            states = communicator.exchange_states_for_energy_calc(state)

            energies = []
            for state in states:
                energy = system_runner.get_energy(state)
                energies.append(energy)
            communicator.send_energies_to_master(energies)

            self._step += 1
