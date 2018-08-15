
from __future__ import absolute_import, print_function

import os
import argparse
from sys import stdout, exit
import socket
import numpy as np
import mdtraj as md

import time, random

import simtk.unit as u
from simtk.openmm import Platform, XmlSerializer
from simtk.openmm.app import PDBFile, Simulation, DCDReporter, StateDataReporter

from datetime import datetime
import argparse
import glob


def get_xml(xml_file):
    # TODO file access control
    attempt = 0
    retries = 500
    if not xml_file.endswith('.xml'):
        raise IOError("{} must end in '.xml' for reading as XML file".format(xml_file))
    while True:
        try:
            with open(xml_file) as f:
                xml = f.read()
                cereal = XmlSerializer.deserialize(xml)
            return xml, cereal
        except ValueError as e:
            if attempt < retries:
                attempt += 1
                time.sleep(5*random.random())
            else:
                raise e


def get_platform(platform_name):
    if platform_name == 'fastest':
        platform = None
    else:
        # TODO file access control
        attempt = 0
        retries = 500
        while True:
            try:
                platform = Platform.getPlatformByName(platform_name)
                return platform
            except IndexError as e:
                if attempt < retries:
                    attempt += 1
                    time.sleep(5*random.random())
                else:
                    raise e


def get_pdbfile(topology_pdb):
    # TODO file access control
    attempt = 0
    retries = 500
    if not topology_pdb.endswith('.pdb'):
        raise IOError("{} must end in '.pdb' for reading as PDB file".format(topology_pdb))
    while True:
        try:
            pdb = PDBFile(topology_pdb)
            return pdb
        except IndexError as e:
            if attempt < retries:
                attempt += 1
                time.sleep(5*random.random())
            else:
                raise e


def read_input(platform, pdbfile, system, integrator):
    return_order = ['get_platform', 'get_pdbfile',
                    'get_system', 'get_integrator']
    funcs = {'get_platform':   [get_platform, platform],
             'get_pdbfile':    [get_pdbfile, pdbfile],
             'get_system':     [get_xml, system],
             'get_integrator': [get_xml, integrator]}
    kfuncs = list(funcs.keys())
    random.shuffle(kfuncs)
    returns = dict()
    while kfuncs:
        op_name = kfuncs.pop(0)
        func, arg = funcs[op_name]
        returns.update({op_name: func(arg)})
    return [returns[nxt] for nxt in return_order]



platformNames = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print(platformNames)


parser = argparse.ArgumentParser()
parser.add_argument('--save_traj',dest='save_traj',required=True,type=str)
parser.add_argument('--md_steps',dest='md_steps',required=True,type=int)
parser.add_argument('--trajstride',dest='trajstride', required=False, type=int, default=10)
parser.add_argument('--idxstart',dest='idxstart',required=True,type=int)
parser.add_argument('--idxend',dest='idxend',required=True,type=int)
parser.add_argument('--iter',dest='iter',required=True,type=int)
parser.add_argument('--path',dest='path',required=True,type=str)
args = parser.parse_args()
#grofile_name='start2.gro

#pdb=mdtraj.load(grofile_name)
#save_traj=False
#save_traj='True'
print("num of structures:",str(args.idxend-args.idxstart))
print("found num:", str(len(glob.glob(args.path+'/iter'+str(args.iter)+'_input*.pdb')))) 
for i in range(args.idxstart,args.idxend):
  print("pdb "+str(i)+" from "+ str(args.idxend))
  a_topology_pdb = args.path+'/iter'+str(args.iter)+'_input'+str(i)+'.pdb' 
  a_platform = 'fastest'
  properties = None
  a_system_xml = 'system-5.xml'
  a_integrator_xml = 'integrator-5.xml'
  platform, pdb, (system_xml, system), (integrator_xml, integrator) \
   = read_input(a_platform, a_topology_pdb, a_system_xml, a_integrator_xml)
  simulation = Simulation( pdb.topology,  system,integrator, platform, properties)
  if i == args.idxstart:
    print('# platform used:', simulation.context.getPlatform().getName())
    print('# platforms available')
    for no_platform in range(Platform.getNumPlatforms()):
        # noinspection PyCallByClass,PyTypeChecker
        print('(%d) %s' % (no_platform, Platform.getPlatform(no_platform).getName()))
    print(os.environ)
    print(Platform.getPluginLoadFailures())
    print(Platform.getDefaultPluginsDirectory())
  simulation.context.setPositions(pdb.positions)
  pbv = system.getDefaultPeriodicBoxVectors()
  simulation.context.setPeriodicBoxVectors(*pbv)
  # set velocities to temperature in integrator
  temperature = integrator.getTemperature()
  dt = integrator.getStepSize()
  simulation.context.setVelocitiesToTemperature(temperature)
  simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
  potentialEnergy=True, temperature=True)) 
  if args.save_traj=='True':
    simulation.reporters.append(DCDReporter(args.path+'/iter'+str(args.iter)+'_traj'+str(i)+'.dcd', args.trajstride)) 
  steps=args.md_steps #1000=2sec each, 10000=20sec
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
  PDBFile.writeFile(simulation.topology, pos, open(args.path+'/iter'+str(args.iter)+'_out'+str(i)+'.pdb', 'a'))
  del simulation, integrator, system
