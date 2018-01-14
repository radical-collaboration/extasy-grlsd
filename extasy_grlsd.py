#!/usr/bin/env python

__author__ = 'Vivek <vivek.balasubramanian@rutgers.edu>'
__copyright__ = 'Copyright 2017, http://radical.rutgers.edu'
__license__ = 'MIT'
__use_case_name__ = 'Gromacs + LSDMap simulation-analysis using EnTK 0.6'


from radical.entk import Pipeline, Stage, Task, AppManager, ResourceManager
import argparse
import os
import glob
import sys
import imp
import json
import traceback

os.environ['RADICAL_PILOT_DBURL'] = 'mongodb://eh22:a3Qv*zs0@ds141209.mlab.com:41209/clementigroup'
#os.environ['RADICAL_PILOT_DBURL'] = 'mongodb://entk:entk@ds033196.mlab.com:33196/extasy_grlsd'
os.environ['RADICAL_ENTK_VERBOSE'] = 'INFO'

def create_workflow(Kconfig):

    # User settings
    ENSEMBLE_SIZE = int(Kconfig.num_CUs)          # Number of ensemble members
    TOTAL_ITERS = int(Kconfig.num_iterations)   # Number of iterations to run current trial

    wf = Pipeline()

    # ------------------------------------------------------------------------------------------------------------------
    '''
    pre_proc_stage :
            
        Purpose : Transfers files, Split the input file into smaller files to be used by each of the
            gromacs instances in the first iteration.

        Arguments :     
            inputfile = file to be split
            numCUs    = number of simulation instances/ number of smaller files
    '''
    pre_proc_stage = Stage()
    pre_proc_task = Task()
    pre_proc_task.pre_exec = ['module load bwpy',
                              'export iter=-1']
    pre_proc_task.executable = ['python']
    pre_proc_task.arguments = [ 'spliter.py',
                                Kconfig.num_CUs,
                                os.path.basename(Kconfig.md_input_file)
                            ]
    pre_proc_task.copy_input_data = ['$SHARED/%s' % os.path.basename(Kconfig.md_input_file),
                                     '$SHARED/spliter.py',
                                     '$SHARED/gro.py'
                                     ]
    pre_proc_task_ref = '$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, pre_proc_stage.uid, pre_proc_task.uid)

    pre_proc_stage.add_tasks(pre_proc_task)
    wf.add_stages(pre_proc_stage)
    # ------------------------------------------------------------------------------------------------------------------

    cur_iter = 0
    while(cur_iter < TOTAL_ITERS):

        # --------------------------------------------------------------------------------------------------------------
        # sim_stage:
        #     Purpose:  In iter=1, use the input files from pre_loop, else use the outputs of the analysis stage in the
        #               previous iteration. Run gromacs on each of the smaller files. Parameter files and executables
        #                are input from pre_loop. There are 'numCUs' number of instances of gromacs per iteration.
        #     Arguments :
        #           grompp    = gromacs parameters filename
        #           topol     = topology filename

        sim_stage = Stage()
        sim_task_ref = list()
        for sim_num in range(ENSEMBLE_SIZE):

            sim_task = Task()
            sim_task.pre_exec = [   'module load bwpy', 
                                     'export PYTHONPATH="/u/sciteam/hruska/local/lib/python2.7/site-packages:/u/sciteam/hruska/local:/u/sciteam/hruska/local/lib/python:$PYTHONPATH"', 
                                     'export PATH=/u/sciteam/hruska/local/bin:$PATH',
                                    'export iter=%s' % cur_iter]
            sim_task.executable = ['/sw/bw/bwpy/0.3.0/python-single/usr/bin/python']
            sim_task.cores = 16
            sim_task.arguments = ['run_openmm.py',
                                  '--gro', 'start.gro',
                                  '--out', 'out.gro', '>', 'md.log']
            sim_task.link_input_data = ['$SHARED/run_openmm.py > run_openmm.py']

            if Kconfig.ndx_file is not None:
                sim_task.link_input_data.append('$SHARED/{0}'.format(os.path.basename(Kconfig.ndx_file)))

            if (cur_iter == 0):
                sim_task.link_input_data.append('%s/temp/start%s.gro > start.gro' % (pre_proc_task_ref, sim_num))
            else:
                sim_task.link_input_data.append('%s/temp/start%s.gro > start.gro' % (post_ana_task_ref, sim_num))

            sim_task_ref.append('$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, sim_stage.uid, sim_task.uid))
            sim_stage.add_tasks(sim_task)

        wf.add_stages(sim_stage)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # pre_ana_task:
        #     Purpose:   The output of each gromacs instance in the simulation stage is a small coordinate file.
        #                 Concatenate such files from each of the gromacs instances to form a larger file.
        #     Arguments:
        #             numCUs = number of simulation instances / number of small files to be concatenated

        pre_ana_stage = Stage()
        pre_ana_task = Task()
        pre_ana_task.pre_exec = [    'module load bwpy', 
                                     'export PYTHONPATH="/u/sciteam/hruska/local/lib/python2.7/site-packages:/u/sciteam/hruska/local:/u/sciteam/hruska/local/lib/python:$PYTHONPATH"', 
                                     'export PATH=/u/sciteam/hruska/local/bin:$PATH',
                                    'export iter=%s' % cur_iter]
        pre_ana_task.executable = ['/sw/bw/bwpy/0.3.0/python-single/usr/bin/python']
        pre_ana_task.arguments = ['pre_analyze_openmm.py']

        pre_ana_task.link_input_data = ['$SHARED/pre_analyze_openmm.py > pre_analyze_openmm.py']
        
        for sim_num in range(ENSEMBLE_SIZE):
            pre_ana_task.link_input_data += ['%s/out.gro > out%s.gro' % (sim_task_ref[sim_num], sim_num)]

        pre_ana_task.copy_output_data = ['tmpha.gro > $SHARED/iter_%s/tmpha.gro' % cur_iter,
                                         'tmp.gro > $SHARED/iter_%s/tmp.gro' % cur_iter]
                                         #'tmp.gro > resource://iter_%s/tmp.gro' % cur_iter

        pre_ana_stage.add_tasks(pre_ana_task)
        wf.add_stages(pre_ana_stage)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # lsdmap:
        #     Purpose: Perform LSDMap on the large coordinate file to generate weights and eigen values.
        #     Arguments:
        #             config = name of the config file to be used during LSDMap

        ana_stage = Stage()
        ana_task = Task()
        ana_task.pre_exec = [   'module load bwpy',
                                'module load platform-mpi',
                                'export PYTHONPATH=/u/sciteam/balasubr/.local/lib/python2.7/site-packages:$PYTHONPATH',
                                'export PATH=/u/sciteam/balasubr/.local/bin:$PATH',
                                'source /u/sciteam/balasubr/ve-extasy/bin/activate',
                                'export iter=%s' % cur_iter
                                ]
        ana_task.executable = ['lsdmap'] #/u/sciteam/hruska/local/bin/lsdmap
        ana_task.arguments = ['-f', os.path.basename(Kconfig.lsdm_config_file),
                              '-c', 'tmpha.gro',
                              '-n', 'out.nn',
                              '-w', 'weight.w'
                              ]

        ana_task.cores = 1
        ana_task.link_input_data = ['$SHARED/{0} > {0}'.format(os.path.basename(Kconfig.lsdm_config_file)),
                                    '$SHARED/iter_%s/tmpha.gro > tmpha.gro' % cur_iter]
        ana_task.copy_output_data = ['tmpha.ev > $SHARED/iter_%s/tmpha.ev' % cur_iter,
                                     'tmpha.eg > $SHARED/iter_%s/tmpha.eg' % cur_iter,
                                     'lsdmap.log > output/iter%s/lsdmap.log'%cur_iter]
        if cur_iter > 0:
          ana_task.link_input_data += ['%s/weight.w > weight.w' % post_ana_task_ref]
          ana_task.copy_output_data += ['weight.w > $SHARED/iter_%s/weight.w' % cur_iter]

        if(cur_iter % Kconfig.nsave == 0):
            ana_task.download_output_data = ['lsdmap.log > output/iter%s/lsdmap.log'%cur_iter]


        ana_task_ref = '$Pipeline_%s_Stage_%s_Task_%s'%(wf.uid, ana_stage.uid, ana_task.uid)

        
        ana_stage.add_tasks(ana_task)
        wf.add_stages(ana_stage)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # post_lsdmap:
        #     Purpose:   Use the weights, eigen values generated in lsdmap along with other parameter files from pre_loop
        #                 to generate the new coordinate file to be used by the simulation_step in the next iteration.
        #     Arguments:
        #             num_runs              = number of configurations to be generated in the new coordinate file
        #             out                   = output filename
        #             cycle                 = iteration number
        #             max_dead_neighbors    = max dead neighbors to be considered
        #             max_alive_neighbors   = max alive neighbors to be considered
        #             numCUs                = number of simulation instances/ number of smaller files

        post_ana_stage = Stage()
        post_ana_task = Task()
        post_ana_task._name      = 'post_ana_task'
        post_ana_task.pre_exec = [   'module load bwpy', 
                                     'export PYTHONPATH="/u/sciteam/hruska/local/lib/python2.7/site-packages:/u/sciteam/hruska/local:/u/sciteam/hruska/local/lib/python:$PYTHONPATH"', 
                                     'export PATH=/u/sciteam/hruska/local/bin:$PATH',
                                    'export iter=%s' % cur_iter
                                ]
        post_ana_task.executable = ['/sw/bw/bwpy/0.3.0/python-single/usr/bin/python']
        post_ana_task.arguments = [ 'post_analyze.py',                                   
                                    Kconfig.num_runs,
                                    'tmpha.ev',
                                    'ncopies.nc',
                                    'tmp.gro',
                                    'out.nn',
                                    'weight.w',
                                    'out.gro',
                                    Kconfig.max_alive_neighbors,
                                    Kconfig.max_dead_neighbors,
                                    'input.gro',
                                    cur_iter,
                                    Kconfig.num_CUs]

        post_ana_task.link_input_data = ['$SHARED/post_analyze.py > post_analyze.py',
                                         '$SHARED/selection.py > selection.py',
                                         '$SHARED/selection-cluster.py > selection-cluster.py',
                                         '$SHARED/reweighting.py > reweighting.py',
                                         '$SHARED/spliter.py > spliter.py',
                                         '$SHARED/gro.py > gro.py',
                                         '$SHARED/iter_%s/tmp.gro > tmp.gro' % cur_iter,
                                         '$SHARED/iter_%s/tmpha.ev > tmpha.ev' % cur_iter,
                                         '$SHARED/iter_%s/tmpha.eg > tmpha.eg' % cur_iter,
                                         '$SHARED/iter_%s/out.nn > out.nn' % cur_iter,
                                         '$SHARED/input.gro > input.gro']

        #if cur_iter > 0:
        post_ana_task.link_input_data += ['%s/weight.w > weight_new.w' % ana_task_ref]

        if(cur_iter % Kconfig.nsave == 0):
            post_ana_task.download_output_data = ['out.gro > output/iter%s/out.gro' % cur_iter,
                                             'weight.w > output/iter%s/weight.w' % cur_iter,
                                             '$SHARED/iter_%s/tmp.gro > output/iter%s/tmp.gro' % (cur_iter,cur_iter) 
                                             ]

        post_ana_task.copy_output_data = ['out.nn > $SHARED/iter_%s/out.nn' % cur_iter,
                                     'ncopies.nc > $SHARED/iter_%s/out.nn' % cur_iter,
                                     'plot-scatter-cluster-10d.png > $SHARED/iter_%s/plot-scatter-cluster-10d.png' % cur_iter,
                                     'ncopies.nc > $SHARED/iter_%s/ncopies.nc' % cur_iter]
                                    # 'plot-scatter-cluster-10d.png > resource://iter_%s/plot-scatter-cluster-10d.png' % cur_iter,
        post_ana_task_ref = '$Pipeline_%s_Stage_%s_Task_%s'%(wf.uid, post_ana_stage.uid, post_ana_task.uid)

        post_ana_stage.add_tasks(post_ana_task)
        wf.add_stages(post_ana_stage)
        # --------------------------------------------------------------------------------------------------------------

        cur_iter += 1

    return wf


# ------------------------------------------------------------------------------
#
if __name__ == '__main__':

    try:

        parser = argparse.ArgumentParser()
        parser.add_argument('--RPconfig', help='link to Radical Pilot related configurations file')
        parser.add_argument('--Kconfig', help='link to Kernel configurations file')

        args = parser.parse_args()

        if args.RPconfig is None:
            parser.error('Please enter a RP configuration file')
            sys.exit(1)

        if args.Kconfig is None:
            parser.error('Please enter a Kernel configuration file')
            sys.exit(0)

        RPconfig = imp.load_source('RPconfig', args.RPconfig)
        Kconfig = imp.load_source('Kconfig', args.Kconfig)

        wf = create_workflow(Kconfig)

        # Create a dictionary describe four mandatory keys:
        # resource, walltime, cores and project
        res_dict = {

            'resource': RPconfig.REMOTE_HOST,
            'walltime': RPconfig.WALLTIME,
            'cores': RPconfig.PILOTSIZE,
            'project': RPconfig.ALLOCATION,
            #'queue': RPconfig.QUEUE,
            'access_schema': 'gsissh'
        }

        # Create Resource Manager object with the above resource description
        rman = ResourceManager(res_dict)
        # Data common to multiple tasks -- transferred only once to common staging area
        rman.shared_data = [Kconfig.md_input_file,
                            Kconfig.lsdm_config_file,
                            Kconfig.top_file,
                            Kconfig.mdp_file,
                            '%s/spliter.py' % Kconfig.helper_scripts,
                            '%s/gro.py' % Kconfig.helper_scripts,
                            '%s/run.py' % Kconfig.helper_scripts,
                            '%s/run_openmm.py' % Kconfig.helper_scripts,
                            '%s/pre_analyze.py' % Kconfig.helper_scripts,
                            '%s/pre_analyze_openmm.py' % Kconfig.helper_scripts,
                            '%s/post_analyze.py' % Kconfig.helper_scripts,
                            '%s/selection.py' % Kconfig.helper_scripts,
                            '%s/selection-cluster.py' % Kconfig.helper_scripts,
                            '%s/reweighting.py' % Kconfig.helper_scripts
                            ]

        if Kconfig.ndx_file is not None:
            rman.shared_data.append(Kconfig.ndx_file)

        # Create Application Manager
        appman = AppManager()
        # appman = AppManager(port=) # if using docker, specify port here.

        # Assign resource manager to the Application Manager
        appman.resource_manager = rman

        # Assign the workflow as a set of Pipelines to the Application Manager
        appman.assign_workflow(wf)

        # Run the Application Manager
        appman.run()

    except Exception as ex:

        print 'Error: {0}'.format(str(ex))
        print traceback.format_exc()
