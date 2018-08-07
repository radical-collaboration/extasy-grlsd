__author__ = 'vivek'

'''
Purpose :   This file is used to split the large input coordinate file into
            smaller coordinate files that can are executed on by each compute unit.
            It uses the gro.py file to get the specifics of the system in question.

Arguments : num_tasks = number of compute units
            grofile_name = name of the coordinate file

'''


import os
import sys
import shutil
import argparse
import mdtraj

if __name__ == '__main__':
    #num_tasks = int(sys.argv[1])
    #grofile_name = str(sys.argv[2])
   
    curdir = os.path.dirname(os.path.abspath(__file__))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', dest='path', required=True, type=str)
    parser.add_argument('--gro',dest='grofile_name',required=True,type=str)
    parser.add_argument('--clone',dest='clone',required=True,type=int)
    args = parser.parse_args()
    file=mdtraj.load(args.path+'/'+args.grofile_name)
    file.save_pdb(args.path + '/tmp_in.pdb')

    print 'Cloning replicas'
    print 'cp '+args.path+'/tmp_in.pdb '+args.path+'/iter0_inputxxx.gro' 
    for i in range(args.clone):
        os.system('cp '+args.path+'/tmp_in.pdb '+args.path+'/iter0_input'+str(i)+'.pdb')

    print 'Finished cloning'
