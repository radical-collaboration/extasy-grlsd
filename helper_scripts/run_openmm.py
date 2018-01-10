import os
import mdtraj
import mdtraj.reporters

from datetime import datetime
import argparse

from simtk.openmm.app import * 
from simtk.openmm import * 
from simtk.unit import *
from sys import stdout


parser = argparse.ArgumentParser()
parser.add_argument('--gro',dest='grofile_name',required=True,type=str)
parser.add_argument('--out',dest='output_grofile_name',required=True,type=str)
args = parser.parse_args()

file=mdtraj.load(grofile_name)

print("num of structures:",len(pdb))
for i in range(len(pdb)):
	topology = pdb.topology.to_openmm()
	#implicit forcefield
	forcefield = ForceField('amber99sbildn.xml', 'amber99_obc.xml')
	temp=300
	long_step=True
	#long_step=False
	if long_step:
	  system = forcefield.createSystem(topology,nonbondedMethod=CutoffNonPeriodic,constraints=AllBonds, hydrogenMass=4*amu)
	  dt=0.005*picoseconds
	  integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, dt)
	else:
	  system = forcefield.createSystem(topology, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds)
	  dt=0.002*picoseconds
	  integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, dt)
	simulation = Simulation(topology, system, integrator)
	simulation.context.setPositions(pdb[i].xyz[0])
	simulation.context.setVelocitiesToTemperature(temp*kelvin)
	#simulation.minimizeEnergy() 
	#simulation.reporters.append(PDBReporter('output.pdb', 1000)) 
	simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
	potentialEnergy=True, temperature=True)) 
	steps=10000
	start=datetime.now()
	simulation.step(steps)
	end = datetime.now()
	elapsed = end -start
	time=elapsed.seconds + elapsed.microseconds*1e-6
	print('Integrated %d steps in %g seconds' % (steps, time))
	print('%g ns/day' % (dt*steps*86400/time).value_in_unit(nanoseconds))
	state = simulation.context.getState(getPositions=True, getVelocities=True,getEnergy=True)
	pbv = state.getPeriodicBoxVectors(asNumpy=True)
	vel = state.getVelocities(asNumpy=True)
	pos = state.getPositions(asNumpy=True)
	print(state.getPotentialEnergy(), state.getKineticEnergy())
	PDBFile.writeFile(simulation.topology, pos, open('Ala12_output_concat.pdb', 'a'))


file=mdtraj.load('Ala12_output_concat.pdb')
file.save_gro(output_grofile_name)



