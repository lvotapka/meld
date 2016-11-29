#!/usr/bin/env python
# encoding: utf-8

from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import *
from sys import stdout
from simtk.openmm.app import forcefield as ff
from xml.etree import ElementTree as ET
from simtk.openmm.openmm import Continuous1DFunction
from simtk.openmm.openmm import Continuous2DFunction
import numpy as np

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

def _get_membrane_force_eastman(gb_charges, gb_ors, go_srs, gb_alphas, gb_betas, gb_gammas, gb_radindices, solventDielectric=78.5, soluteDielectric=1, thickness=5., offset=0.0195141, SA=None, reference_force=None):
    custom = mm.CustomGBForce()
    custom.setNonbondedMethod(reference_force.getNonbondedMethod())
    custom.setCutoffDistance(reference_force.getCutoffDistance())
    custom.addPerParticleParameter("q");
    custom.addPerParticleParameter("or");
    custom.addPerParticleParameter("sr");
    custom.addPerParticleParameter("alpha");
    custom.addPerParticleParameter("beta");
    custom.addPerParticleParameter("gamma");
    custom.addPerParticleParameter("radindex");
    custom.addGlobalParameter("memb_thickness", thickness);
    custom.addGlobalParameter("solventDielectric", solventDielectric);
    custom.addGlobalParameter("soluteDielectric", soluteDielectric);
    custom.addGlobalParameter("offset", offset);
    custom.addComputedValue("Imol", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                             "U=r+sr2;"
                             "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                             "L=max(or1, D);"
                             "D=abs(r-sr2);", mm.CustomGBForce.ParticlePairNoExclusions)
                             #"sr2 = scale2*or2;"
                             #"or1 = radius1-0.009; or2 = radius2-0.009", mm.CustomGBForce.ParticlePairNoExclusions)
    custom.addComputedValue("", mm.CustomGBForce.SingleParticle)
    custom.addComputedValue("Imem", "(1/radius+2*log(2)/thickness)/(1+exp(7.2*(abs(z)+radius-0.5*thickness)))",  mm.CustomGBForce.SingleParticle)
    custom.addComputedValue("B", "1/(1/or-tanh(alpha*psi-beta*psi^2+gamma*psi^3)/radius);"
                             "psi=max(Imol,Imem)*or; or=radius-0.009",  mm.CustomGBForce.SingleParticle)
    custom.addComputedValue("ratio", "0.580127*soluteDielectric/solventDielectric",mm.CustomGBForce.SingleParticle)
    custom.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B", mm.CustomGBForce.SingleParticle)
    custom.addEnergyTerm("-138.935456*(1/soluteDielectric-1/solventDielectric)*(1/(1+ratio))*q1*q2/f;"
                          "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", mm.CustomGBForce.ParticlePairNoExclusions)
    
    for c, ors, srs, a, b, g, ri in zip(gb_charges, gb_ors, go_srs, gb_alphas, gb_betas, gb_gammas, gb_radindices):
        custom.addParticle([c, ors, srs, a, b, g, ri])

    return custom

def make_membrane_2d_func(thickness, xmin, xmax, xres, ymin, ymax, yres): # TODO: remove
    x_counter = 0
    y_counter = 0
    values = []
    half = thickness / 2.0
    for z in np.arange(xmin, xmax, xres):
        y_counter = 0
        for r_old in np.arange(ymin, ymax, yres):
            if abs(z) < half - r_old: # then the atom is fully within the membrane
                d1 = abs(z + half)
                d2 = abs(z - half)
                r_new = 4.0/(1.0/d1 + 1.0/d2)
      
            elif (abs(z) >= half - r_old) and (abs(z) < half + r_old): # the atom is partially embedded in the memb
                d2 = abs(z) + half
                z0 = half - abs(z)
                r_new = 4.0/((2.0*r_old - z0)/r_old**2 + 1.0/d2) #4.0/((r_old - z0)/r_old**2 + 1.0/d2 + 1.0/r_old)
      
            elif abs(z) >= half + r_old: # the atom is fully outside the membrane
                d1 = abs(z) - half
                d2 = abs(z) + half
                r_new = 1.0 / (1.0/r_old - 0.25*(1.0/d1 - 1.0/d2))
                
            else:
                raise Exception, "Should never get here! z = %f" % z
            
            values.append(r_new)
            y_counter += 1
        x_counter += 1
    
    return x_counter, y_counter, values

def _get_membrane_force(gb_charges, gb_radii, gb_scale, solventDielectric=78.5, membraneDielectric=10.0, thickness=2., SA=None, reference_force=None):
    custom = mm.CustomGBForce()
    custom.setNonbondedMethod(reference_force.getNonbondedMethod())
    custom.setCutoffDistance(reference_force.getCutoffDistance())
    print "DEBUG thickness:", thickness
    custom.addPerParticleParameter("q");
    custom.addPerParticleParameter("or");
    custom.addPerParticleParameter("sr");
    custom.addGlobalParameter("thickness", thickness);
    custom.addGlobalParameter("solventDielectric", solventDielectric);
    custom.addGlobalParameter("membraneDielectric", membraneDielectric);
    custom.addComputedValue("Imol", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                             "U=r+sr2;"
                             "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                             "L=max(or1, D);"
                             "D=abs(r-sr2);", mm.CustomGBForce.ParticlePairNoExclusions)
    #custom.addComputedValue("Imem", "(1/radius+2*log(2)/thickness)/(1+exp(7.2*(abs(z)+radius-0.5*thickness)));"
    #                         "radius=or+0.009",  mm.CustomGBForce.SingleParticle)
    custom.addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                             "psi=Imol*or; radius=or+0.009",  mm.CustomGBForce.SingleParticle)
    #custom.addComputedValue("B_old", "(1-step(B-half))*B + step(B-half)*half;"
    #                         "half=thickness/2.0", mm.CustomGBForce.SingleParticle)
    custom.addComputedValue("Bnew", "condition1+condition2+condition3+condition4+condition5+condition6;"
                             "condition1=(1-step(r_edge+half))/abs((1.0/B-0.25*(1.0/sqrt((abs(z)-half)^2+0.0001)-1.0/(sqrt((abs(z)+half)^2+0.0001))))+0.0001);"
                            "condition2=step(r_edge+half)*(1-step(r_edge-half))*(1-step(l_edge+half))*4.0/sqrt(((2.0*B-half-z)/(B^2)+1.0/sqrt((half-z)^2+0.0001))^2+0.0001);"
                            "condition3=step(r_edge+half)*(1-step(r_edge-half))*step(l_edge+half)*(1-step(l_edge-half))*2.0*(half-z)*(z+half)/half;"
                            "condition4=(1-step(l_edge+half))*step(r_edge-half)*B^2/sqrt((B - 0.5*half)^2+0.0001);"
                            "condition5=step(r_edge-half)*step(l_edge+half)*(1-step(l_edge-half))*4.0/sqrt(((2.0*B - half + z)/B^2 + 1.0/sqrt((z+half)^2+0.0001))^2+0.0001);"
                            "condition6=step(l_edge-half)/sqrt((1.0/B - 0.25*(1.0/sqrt((z-half)^2+0.0001) - 1.0/sqrt((z+half)^2+0.0001)))^2+0.0001);"
                             "half=thickness/2.0;l_edge=z-B;r_edge=z+B; diel_ratio=membraneDielectric/solventDielectric", mm.CustomGBForce.SingleParticle)
    #custom.addTabulatedFunction("get_new_Born_rad", Continuous2DFunction(xsize, ysize, memb_values, xmin, xmax, ymin, ymax))
    #custom.addComputedValue("Bnew", "B+(thickness-B)*exp(-z^2/(2*half^2));"
    #                        "half=thickness/2.0", mm.CustomGBForce.SingleParticle)
    #custom.addComputedValue("Bnew", "get_new_Born_rad(z,B)", mm.CustomGBForce.SingleParticle)
    custom.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/Bnew)^6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/Bnew;"
                          "radius=or+0.009", mm.CustomGBForce.SingleParticle)
    custom.addEnergyTerm("-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                          "f=sqrt(r^2+Bnew1*Bnew2*exp(-r^2/(4*Bnew1*Bnew2)))", mm.CustomGBForce.ParticlePairNoExclusions)
    
    for c, r, s in zip(gb_charges, gb_radii, gb_scale):
        custom.addParticle([c, r, s])

    return custom

def _get_membrane_force_old(gb_charges, gb_radii, gb_scale, solventDielectric=78.5, soluteDielectric=1, thickness=5., SA=None, reference_force=None):
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


def add_membrane(system, sigma_min=0.151, implicit_solvent="obc"):
    system_xml = mm.XmlSerializer.serializeSystem(system)
    fo = open('nomembrane.dat','w')
    fo.write(system_xml)
    fo.close()

    # grab and remove the non-bonded force
    print "grab and remove the non-bonded force"
    root = ET.fromstring(system_xml)
    nb_string = _keep_force(root, 'NonbondedForce')
    #gb_string = _extract_force(root, 'GBSAOBCForce')
    gb_string = _extract_force(root, 'CustomGBForce')
    
    assert implicit_solvent == 'obc', "Only gbNeck implemented currently. %s has not yet been included as a MELD membrane model." % implicit_solvent
    
    print "nb_string:", nb_string
    print "gb_string:", gb_string # this is None

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
