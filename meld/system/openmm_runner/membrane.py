#!/usr/bin/env python
# encoding: utf-8

from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import *
from sys import stdout
from simtk.openmm.app import forcefield as ff
from xml.etree import ElementTree as ET

def _keep_force(root, force_name):
    forces = root.find('Forces')
    force = forces.findall('./Force[@type="{}"]'.format(force_name))
    print force_name
    print forces.findall('./Force')
    if force:
        assert len(force) == 1
        force = force[0]
        return force
    else:
        return None

def _extract_force(root, force_name):
    forces = root.find('Forces')
    force = forces.findall('./Force[@type="{}"]'.format(force_name))
    if force:
        assert len(force) == 1
        force = force[0]
        forces.remove(force)
        return force
    else:
        return None


def _get_membrane_force(gb_charges, gb_radii, gb_scale, solventDielectric=78.5, soluteDielectric=1, thickness=3., SA=None, reference_force=None):
    custom = mm.CustomGBForce()
    custom.setNonbondedMethod(reference_force.getNonbondedMethod())
    custom.setCutoffDistance(reference_force.getCutoffDistance())
    custom.addPerParticleParameter("q");
    custom.addPerParticleParameter("radius");
    custom.addPerParticleParameter("scale");
    custom.addGlobalParameter("thickness", thickness);
    custom.addGlobalParameter("solventDielectric", solventDielectric);
    custom.addGlobalParameter("soluteDielectric", soluteDielectric);
    custom.addComputedValue("Imol", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                             "U=r+sr2;"
                             "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                             "L=max(or1, D);"
                             "D=abs(r-sr2);"
                             "sr2 = scale2*or2;"
                             "or1 = radius1-0.009; or2 = radius2-0.009", mm.CustomGBForce.ParticlePairNoExclusions)
    custom.addComputedValue("Imem", "(1/radius+2*log(2)/thickness)/(1+exp(7.2*(abs(z)+radius-0.5*thickness)))",  mm.CustomGBForce.SingleParticle)
    custom.addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                             "psi=max(Imol,Imem)*or; or=radius-0.009",  mm.CustomGBForce.SingleParticle)
    custom.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B", mm.CustomGBForce.SingleParticle)
    custom.addEnergyTerm("-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                          "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", mm.CustomGBForce.ParticlePairNoExclusions)
    
    for c, r, s in zip(gb_charges, gb_radii, gb_scale):
        custom.addParticle([c, r, s])

    return custom


def add_membrane(system, sigma_min=0.151):
    system_xml = mm.XmlSerializer.serializeSystem(system)
    fo = open('nomembrane.dat','w')
    fo.write(system_xml)
    fo.close()

    # grab and remove the non-bonded force
    root = ET.fromstring(system_xml)
    nb_string = _keep_force(root, 'NonbondedForce')
    gb_string = _extract_force(root, 'GBSAOBCForce')

    nb_force = mm.XmlSerializer.deserialize(ET.tostring(nb_string))

    n_particles = nb_force.getNumParticles()

    # extract the non-bonded parameters
    new_system = mm.XmlSerializer.deserializeSystem(ET.tostring(root))

    # extract the gb_parameters
    if gb_string is not None:
        gb_force = mm.XmlSerializer.deserialize(ET.tostring(gb_string))
        gb_params = [gb_force.getParticleParameters(i) for i in range(n_particles)]
        gb_charge = [p[0] for p in gb_params]
        gb_radius = [p[1] for p in gb_params]
        gb_scale = [p[2] for p in gb_params]
        membrane_force = _get_membrane_force(gb_charge, gb_radius, gb_scale,SA='ACE',reference_force=gb_force)
        new_system.addForce(membrane_force)
    new_xml = mm.XmlSerializer.serializeSystem(new_system)
    fo = open('membrane.dat','w')
    fo.write(new_xml)
    fo.close()
    return new_system
