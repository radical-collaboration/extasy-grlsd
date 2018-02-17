__author__ = 'vivek'

'''
Purpose :   This file is used to split the large input coordinate file into
            smaller coordinate files that can are executed on by each compute unit.
            It uses the gro.py file to get the specifics of the system in question.

Arguments : num_tasks = number of compute units
            grofile_name = name of the coordinate file

'''


import gro
import os
import sys
import shutil
import argparse

if __name__ == '__main__':
    #num_tasks = int(sys.argv[1])
    #grofile_name = str(sys.argv[2])
   
    curdir = os.path.dirname(os.path.abspath(__file__))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', dest='num_tasks', required=True, type=int)
    parser.add_argument('-gro',dest='grofile_name',required=True,type=str)
    parser.add_argument('--clone',dest='clone',required=False,type=int)
    args = parser.parse_args()
    
    if args.clone is not None:
      print 'Cloning replicas'
      for i in range(args.clone):
        os.system('cat '+args.grofile_name+' >> input-cloned_replicas.gro')
      print 'Prepare grofiles for splitting'
      grofile_obj = gro.GroFile(os.path.dirname(os.path.abspath(__file__)) + '/input-cloned_replicas.gro')
    else:
      print 'Prepare grofiles for splitting'
      grofile_obj = gro.GroFile(os.path.dirname(os.path.abspath(__file__)) + '/out.gro')

    if grofile_obj.nruns<args.num_tasks:
        print "WARNING: number of runs should be greater or equal to the number of tasks."
        #sys.exit(1)

    nruns_per_task = [grofile_obj.nruns/args.num_tasks for _ in xrange(args.num_tasks)]
    nextraruns=grofile_obj.nruns%args.num_tasks

    for idx in xrange(nextraruns):
        nruns_per_task[idx] += 1

    if os.path.isdir('%s/temp'%curdir) is True:
        shutil.rmtree('%s/temp' % curdir)
    os.mkdir('%s/temp'%curdir)

    with open(grofile_obj.filename, 'r') as grofile:
        for idx in xrange(args.num_tasks):
            start_grofile_name = curdir + '/temp/start%s.gro'%idx
            with open(start_grofile_name, 'w') as start_grofile:
                nlines_per_task = nruns_per_task[idx]*grofile_obj.nlines_per_run
                for jdx in xrange(nlines_per_task):
                    line=grofile.readline()
                    if line:
                        line = line.replace("\n", "")
                        print >> start_grofile, line
                    else:
                        break

    print 'Finished splitting'
