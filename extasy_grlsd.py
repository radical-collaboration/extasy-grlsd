#!/usr/bin/env python

__author__ = 'Vivek <vivek.balasubramanian@rutgers.edu> Eugen <eh22@rice.edu>'
__copyright__ = 'Copyright 2018, http://radical.rutgers.edu,  http://clementiresearch.rice.edu'
__license__ = 'MIT'
__use_case_name__ = 'Gromacs + LSDMap simulation-analysis using EnTK'


from radical.entk import Pipeline, Stage, Task, AppManager, ResourceManager
import argparse
import os
import glob
import sys
import imp
import json
import traceback


def create_workflow(Kconfig):


    wf = Pipeline()

    # ------------------------------------------------------------------------------------------------------------------
    cur_iter = int(Kconfig.start_iter)#0
    #assumed of iteration non zero that files are in combined_path
    combined_path=str(Kconfig.remote_output_directory)  #'/u/sciteam/hruska/scratch/extasy-grlsd'
    if cur_iter==0:
    	restart_iter=0
    else:
    	restart_iter=cur_iter


    if cur_iter==0:
      pre_proc_stage = Stage()
      pre_proc_task = Task()
      pre_proc_task.pre_exec = ['module load bwpy', 'export tasks=pre_proc',  
                                'export iter=-1','export OMP_NUM_THREADS=1']
      pre_proc_task.executable = ['python']
      pre_proc_task.arguments = [ 'spliter.py','-n',
                                  Kconfig.num_parallel_MD_sim,'-gro',
                                  'input.gro','--clone',str(Kconfig.num_replicas)
                              ]
      pre_proc_task.copy_input_data = ['$SHARED/%s > %s/iter_%s/input.gro' % (os.path.basename(Kconfig.md_input_file),combined_path,cur_iter),
                                       '$SHARED/%s > input.gro' % os.path.basename(Kconfig.md_input_file),
                                       '$SHARED/spliter.py > spliter.py',
                                       '$SHARED/gro.py > gro.py']

                                       
      pre_proc_task_ref = '$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, pre_proc_stage.uid, pre_proc_task.uid)
      pre_proc_stage.add_tasks(pre_proc_task)
      wf.add_stages(pre_proc_stage)
      # ------------------------------------------------------------------------------------------------------------------
    else:
      pre_proc_stage = Stage()
      pre_proc_task = Task()
      pre_proc_task.pre_exec = ['module load bwpy',
                                'export tasks=pre_proc','export iter=-1', 'export OMP_NUM_THREADS=1' ]
      pre_proc_task.executable = ['python']
      pre_proc_task.arguments = [ 'spliter.py','-n',
                                  Kconfig.num_parallel_MD_sim,
                                  '-gro','input.gro'
                              ]
      pre_proc_task.copy_input_data = ['%s/iter_%s/out.gro > input.gro'  % (combined_path,cur_iter-1),
                                       '$SHARED/spliter.py > spliter.py',
                                       '$SHARED/gro.py > gro.py']
      pre_proc_task_ref = '$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, pre_proc_stage.uid, pre_proc_task.uid)
      pre_proc_stage.add_tasks(pre_proc_task)
      wf.add_stages(pre_proc_stage)
    
    while(cur_iter <  int(Kconfig.num_iterations)):

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
        for sim_num in range(min(int(Kconfig.num_parallel_MD_sim),int(Kconfig.num_replicas))):

            sim_task = Task()
            if Kconfig.use_gpus=='False':
              sim_task.executable = ['/sw/bw/bwpy/0.3.0/python-single/usr/bin/python']
              sim_task.pre_exec = [   'module load bwpy',
                                     'export PYTHONPATH="/u/sciteam/hruska/local/lib/python2.7/site-packages:/u/sciteam/hruska/local:/u/sciteam/hruska/local/lib/python:$PYTHONPATH"',
                                     'export PATH=/u/sciteam/hruska/local/bin:$PATH',
                                    'export iter=%s' % cur_iter]  
              sim_task.cores = int(Kconfig.num_CUs_per_MD_replica) #on bluewaters tasks on one node are executed concurently
            else:
              sim_task.executable = ['python']
              sim_task.pre_exec = [  'module swap PrgEnv-cray PrgEnv-gnu','module add bwpy','module add bwpy-mpi', 'module add fftw', 'module add cray-netcdf', 'module add cudatoolkit/7.5.18-1.0502.10743.2.1', 'module add cmake', 'module unload darshan, xalt','export CRAYPE_LINK_TYPE=dynamic', 'export CRAY_ADD_RPATH=yes', 'export FC=ftn', 'source /projects/sciteam/bamm/hruska/vpy2/bin/activate',
                                     'export tasks=md','export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1' ]
              sim_task.gpu_reqs = { 'processes': 1,
                                    'process_type': None,
                                    'threads_per_process': 1,
                                    'thread_type': None
                                }
              sim_task.cpu_reqs = { 'processes': 0, 
                                    'process_type': None, 
                                    'threads_per_process': 0, 
                                    'thread_type': None
                                  }
            sim_task.arguments = ['run_openmm.py',
                                  '--gro', 'start.gro',
                                  '--out', 'out.gro', '--md_steps',str(Kconfig.md_steps), '--save_traj', 'False','>', 'md.log']
            sim_task.link_input_data = ['$SHARED/%s > run_openmm.py' % (os.path.basename(Kconfig.md_run_file))]

            #if Kconfig.ndx_file is not None:
            #    sim_task.link_input_data.append('$SHARED/{0}'.format(os.path.basename(Kconfig.ndx_file)))
            if restart_iter==cur_iter:
            	sim_task.link_input_data.append('%s/temp/start%s.gro > start.gro' % (pre_proc_task_ref, sim_num))
            else:
                sim_task.link_input_data.append('%s/temp/start%s.gro > start.gro' % (post_ana_task_ref, sim_num))
            

            sim_task_ref.append('$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, sim_stage.uid, sim_task.uid))
            sim_stage.add_tasks(sim_task)

        wf.add_stages(sim_stage)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # pre_ana_task:
        #     Purpose:   The output of each gromacs instance in the simulaxftion stage is a small coordinate file.
        #                 Concatenate such files from each of the gromacs instances to form a larger file.
        #     Arguments:
        #             numCUs = number of simulation instances / number of small files to be concatenated

        pre_ana_stage = Stage()
        pre_ana_task = Task()
        pre_ana_task.pre_exec = [ 'module swap PrgEnv-cray PrgEnv-gnu','module add bwpy','module add bwpy-mpi', 'module add fftw', 'module add cray-netcdf', 'module add cudatoolkit/7.5.18-1.0502.10743.2.1', 'module add cmake', 'module unload darshan, xalt','export CRAYPE_LINK_TYPE=dynamic', 'export CRAY_ADD_RPATH=yes', 'export FC=ftn', 'source /projects/sciteam/bamm/hruska/vpy2/bin/activate', 'export tasks=pre_ana',
                                    'export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1' ]
        pre_ana_task.executable = ['python']
        pre_ana_task.arguments = ['pre_analyze_openmm.py']

        pre_ana_task.link_input_data = ['$SHARED/pre_analyze_openmm.py > pre_analyze_openmm.py']
        
        for sim_num in range(min(int(Kconfig.num_parallel_MD_sim),int(Kconfig.num_replicas))):
            pre_ana_task.link_input_data += ['%s/out.gro > out%s.gro' % (sim_task_ref[sim_num], sim_num)]

        pre_ana_task.copy_output_data = ['tmpha.gro > %s/iter_%s/tmpha.gro' % (combined_path,cur_iter),
                                         'tmp.gro > %s/iter_%s/tmp.gro' % (combined_path,cur_iter)]
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
        ana_task.pre_exec = [   
'module swap PrgEnv-cray PrgEnv-gnu','module add bwpy','module add bwpy-mpi', 'module add fftw', 'module add cray-netcdf', 'module add cudatoolkit/7.5.18-1.0502.10743.2.1', 'module add cmake', 'module unload darshan, xalt','export CRAYPE_LINK_TYPE=dynamic', 'export CRAY_ADD_RPATH=yes', 'export FC=ftn', 'source /projects/sciteam/bamm/hruska/vpy2/bin/activate',
 'export tasks=lsdmap',
'export iter=%s' % cur_iter,
'export OMP_NUM_THREADS=1' ]
        ana_task.executable = ['lsdmap'] #/u/sciteam/hruska/local/bin/lsdmap
        ana_task.arguments = ['-f', os.path.basename(Kconfig.lsdm_config_file),
                              '-c', 'tmpha.gro',
                              '-n', 'out.nn',
                              '-w', 'weight.w'
                              ]

        ana_task.cores = 1
        ana_task.link_input_data = ['$SHARED/{0} > {0}'.format(os.path.basename(Kconfig.lsdm_config_file)),
                                    '%s/iter_%s/tmpha.gro > tmpha.gro' % (combined_path,cur_iter)]
        ana_task.copy_output_data = ['tmpha.ev > $SHARED/iter_%s/tmpha.ev' % cur_iter,
                                     'tmpha.eg > $SHARED/iter_%s/tmpha.eg' % cur_iter,
                                     'lsdmap.log > output/iter_%s/lsdmap.log'%cur_iter,
                                     'tmpha.ev > %s/iter_%s/tmpha.ev' % (combined_path,cur_iter),
                                     'tmpha.eps > %s/iter_%s/tmpha.eps' % (combined_path,cur_iter),
                                     'tmpha.eg > %s/iter_%s/tmpha.eg' % (combined_path,cur_iter),
                                     'out.nn > %s/iter_%s/out.nn' % (combined_path,cur_iter),
                                     'lsdmap.log > %s/iter_%s/lsdmap.log' % (combined_path,cur_iter)
                                     ]
        if cur_iter > 0:
            ana_task.link_input_data += ['%s/iter_%s/weight_out.w > weight.w' % (combined_path,cur_iter-1)]

        if(cur_iter % Kconfig.nsave == 0):
            ana_task.download_output_data = ['lsdmap.log > output/iter_%s/lsdmap.log' % cur_iter]

        ana_task_ref = '$Pipeline_%s_Stage_%s_Task_%s'%(wf.uid, ana_stage.uid, ana_task.uid)

        
        ana_stage.add_tasks(ana_task)
        wf.add_stages(ana_stage)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # post_lsdmap:
        #     Purpose:   Use the weights, eigen values generated in lsdmap along with other parameter files from pre_loop
        #                 to generate the new coordinate file to be used by the simulation_step in the next iteration.
        #     Arguments:
        #             num_replicas              = number of configurations to be generated in the new coordinate file
        #             out                   = output filename
        #             cycle                 = iteration number
        #             max_dead_neighbors    = max dead neighbors to be considered
        #             max_alive_neighbors   = max alive neighbors to be considered
        #             numCUs                = number of simulation instances/ number of smaller files

        post_ana_stage = Stage()
        post_ana_task = Task()
        post_ana_task._name      = 'post_ana_task'
        if Kconfig.restarts == 'clustering':
          post_ana_task.pre_exec = [ 'module swap PrgEnv-cray PrgEnv-gnu','module add bwpy/0.3.0','module add bwpy-mpi', 'module add fftw', 'module add cray-netcdf', 'module add cudatoolkit/7.5.18-1.0502.10743.2.1', 'module add cmake', 'module unload darshan, xalt','export CRAYPE_LINK_TYPE=dynamic', 'export CRAY_ADD_RPATH=yes', 'export FC=ftn', 'source /projects/sciteam/bamm/hruska/vpy2/bin/activate', 
'export tasks=post_ana',
                                    'export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1'   ]
          post_ana_task.executable = ['python']
          post_ana_task.arguments = [ 'post_analyze.py',                                   
                                    Kconfig.num_replicas,
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
                                    Kconfig.num_parallel_MD_sim,
                                    'weight_out.w',
                                    'tmpha.eg']

          post_ana_task.link_input_data = ['$SHARED/post_analyze.py > post_analyze.py',
                                         '$SHARED/selection.py > selection.py',
                                         '$SHARED/selection-cluster.py > selection-cluster.py',
                                         '$SHARED/reweighting.py > reweighting.py',
                                         '$SHARED/spliter.py > spliter.py',
                                         '$SHARED/gro.py > gro.py',
                                         '%s/iter_%s/weight_out.w > weight.w' % (combined_path,cur_iter-1),
                                         '%s/iter_%s/tmp.gro > tmp.gro' % (combined_path,cur_iter),
                                         '%s/iter_%s/tmpha.ev > tmpha.ev' % (combined_path,cur_iter),
                                         '%s/iter_%s/tmpha.eg > tmpha.eg' % (combined_path,cur_iter),
                                         '%s/iter_%s/out.nn > out.nn' % (combined_path,cur_iter)]


          if(cur_iter % Kconfig.nsave == 0):
               post_ana_task.download_output_data = ['out.gro > output/iter_%s/out.gro' % cur_iter,
                                             'weight_out.w > output/iter_%s/weight_out.w' % cur_iter,
                                             'plot-scatter-cluster-10d.png > output/iter_%s/plot-scatter-cluster-10d.png' % (cur_iter),
                                             'ncopies.nc > output/iter_%s/ncopies.nc' % (cur_iter),
                                             '%s/iter_%s/tmp.gro > output/iter_%s/tmp.gro' % (combined_path,cur_iter,cur_iter) 
                                             ]

          post_ana_task.copy_output_data = ['ncopies.nc > %s/iter_%s/ncopies.nc' % (combined_path,cur_iter),
                                     'weight_out.w > %s/iter_%s/weight_out.w' % (combined_path,cur_iter),
                                     'out.gro > %s/iter_%s/out.gro' % (combined_path,cur_iter),
                                     'plot-scatter-cluster-10d.png > %s/iter_%s/plot-scatter-cluster-10d.png' % (combined_path,cur_iter),
                                     'plot-scatter-cluster-10d-counts.png > %s/iter_%s/plot-scatter-cluster-10d-counts.png' % (combined_path,cur_iter),
                                     'plot-scatter-cluster-10d-ncopiess.png > %s/iter_%s/plot-scatter-cluster-10d-ncopiess.png' % (combined_path,cur_iter)]

        post_ana_task_ref = '$Pipeline_%s_Stage_%s_Task_%s'%(wf.uid, post_ana_stage.uid, post_ana_task.uid)

        post_ana_stage.add_tasks(post_ana_task)
        wf.add_stages(post_ana_stage)
        # --------------------------------------------------------------------------------------------------------------

        cur_iter += 1
        Kconfig.start_iter=str(cur_iter)

    return wf


# ------------------------------------------------------------------------------
#
if __name__ == '__main__':

    try:

        parser = argparse.ArgumentParser()
        parser.add_argument('--Kconfig', help='link to Kernel configurations file')
        #parser.add_argument('--port', dest="port", help='port for RabbitMQ server', default=5672, type=int)
        args = parser.parse_args()
        
        if args.Kconfig is None:
            parser.error('Please enter a Kernel configuration file')
            sys.exit(0)

        Kconfig = imp.load_source('Kconfig', args.Kconfig)

        wf = create_workflow(Kconfig)

        # Create a dictionary describe four mandatory keys:
        # resource, walltime, cores and project
        if Kconfig.use_gpus=='False':
          res_dict = {
            'resource': Kconfig.REMOTE_HOST,
            'walltime': Kconfig.WALLTIME,
            'cores': Kconfig.PILOTSIZE,
            'project': Kconfig.ALLOCATION,
            'queue': Kconfig.QUEUE,
            'access_schema': 'gsissh'
          }
        elif Kconfig.use_gpus=='True':
          print "using gpus"
          res_dict = {
            'resource': Kconfig.REMOTE_HOST,
            'walltime': Kconfig.WALLTIME,
            'cores': Kconfig.PILOTSIZE,
            'cpus': Kconfig.PILOTSIZE,
            'cpu_processes': Kconfig.num_CUs_per_MD_replica,#PILOTSIZE,
            'gpus': Kconfig.PILOTSIZE/16,
            'project': Kconfig.ALLOCATION,
            'queue': Kconfig.QUEUE,
            'access_schema': 'gsissh'
          }	  
        else:
          print("use_gpus not recognized")
          
        print res_dict
        # Create Resource Manager object with the above resource description
        rman = ResourceManager(res_dict)
        # Data common to multiple tasks -- transferred only once to common staging area
        rman.shared_data = [Kconfig.md_input_file,
                            Kconfig.md_run_file,
                            Kconfig.lsdm_config_file,
                            '%s/spliter.py' % Kconfig.helper_scripts,
                            '%s/gro.py' % Kconfig.helper_scripts,
                            #'%s/run.py' % Kconfig.helper_scripts,
                            #'%s/run_openmm.py' % Kconfig.helper_scripts,
                            #'%s/pre_analyze.py' % Kconfig.helper_scripts,
                            '%s/pre_analyze_openmm.py' % Kconfig.helper_scripts,
                            '%s/post_analyze.py' % Kconfig.helper_scripts,
                            #'%s/selection.py' % Kconfig.helper_scripts,
                            '%s/selection-cluster.py' % Kconfig.helper_scripts,
                            '%s/reweighting.py' % Kconfig.helper_scripts
                            ]

        #if Kconfig.ndx_file is not None:
        #    rman.shared_data.append(Kconfig.ndx_file)

        # Create Application Manager, only one extasy script on one rabbit-mq server now
        appman = AppManager(hostname='two.radical-project.org', port=33134)#port=args.port)
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
