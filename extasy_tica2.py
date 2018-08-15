#!/usr/bin/env python

__author__ = 'Vivek <vivek.balasubramanian@rutgers.edu> Eugen <eh22@rice.edu>'
__copyright__ = 'Copyright 2018, http://radical.rutgers.edu,  http://clementiresearch.rice.edu'
__license__ = 'MIT'
__use_case_name__ = 'Gromacs + LSDMap simulation-analysis using EnTK'


from radical.entk import Pipeline, Stage, Task, AppManager
import argparse
import os
import glob
import sys
import imp
import json
import traceback
import time
import socket
print(socket.gethostname())

def create_workflow(Kconfig,args):


    wf = Pipeline()

    # ------------------------------------------------------------------------------------------------------------------
    cur_iter = int(Kconfig.start_iter)#0
    #assumed of iteration non zero that files are in combined_path
    if str(socket.gethostname())=='giotto.rice.edu':
      combined_path=str(Kconfig.remote_output_directory)+'-giotto' 
    else:
      combined_path=str(Kconfig.remote_output_directory)  #'/u/sciteam/hruska/scratch/extasy-tica'
    num_parallel=int(Kconfig.NODESIZE)
    num_replicas=int(Kconfig.num_replicas)
    #if cur_iter==0:
    #	restart_iter=0
    #else:
    #	restart_iter=cur_iter


    if cur_iter==0:
      pre_proc_stage = Stage()
      pre_proc_task = Task()
      pre_proc_task.pre_exec = ['export tasks=pre_proc_task','export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1']
      pre_proc_task.executable = ['mv']
      pre_proc_task.arguments = [ combined_path, combined_path + time.strftime("%Y-%m-%d-%H-%M") ]
      pre_proc_task_ref = '$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, pre_proc_stage.uid, pre_proc_task.uid)
      pre_proc_stage.add_tasks(pre_proc_task)
      wf.add_stages(pre_proc_stage)
      pre_proc_stage2 = Stage()
      pre_proc_task2 = Task()
      pre_proc_task2.pre_exec = ['export tasks=pre_proc_task','export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1']
      pre_proc_task2.executable = ['ls']
      pre_proc_task2.arguments = ['-l']
      pre_proc_task2.copy_input_data = ['$SHARED/%s > %s/%s' % (args.Kconfig,combined_path, args.Kconfig),
                                     '$SHARED/run-tica-msm.py > %s/run-tica-msm.py' % combined_path,
                                     '$SHARED/%s > %s/%s' % (Kconfig.md_run_file,combined_path,Kconfig.md_run_file)
                                       ]


      pre_proc_task_ref2 = '$Pipeline_%s_Stage_%s_Task_%s' % (wf.uid, pre_proc_stage2.uid, pre_proc_task2.uid)
      pre_proc_stage2.add_tasks(pre_proc_task2)
      wf.add_stages(pre_proc_stage2)
 
      # ------------------------------------------------------------------------------------------------------------------
    
    while(cur_iter <  int(Kconfig.num_iterations)):

        # --------------------------------------------------------------------------------------------------------------
        # sim_stage:
        #     Purpose:  In iter=1, use the input files from pre_loop, else use the outputs of the analysis stage in the
        #               previous iteration. Run gromacs on each of the smaller files. Parameter files and executables
        #                are input from pre_loop. There arei 'numCUs' number of instances of gromacs per iteration.
        #     Arguments :
        #           grompp    = gromacs parameters filename
        #           topol     = topology filename

        sim_stage = Stage()
        sim_task_ref = list()
        def_rep_per_thread=int(num_replicas/num_parallel)+1
        num_allocated_rep=0
        num_used_threads=0
        while(num_allocated_rep<num_replicas):
          if (num_used_threads==num_parallel):
             print("ALLERT tried use more gpus than allocated")
          if ((num_replicas-num_allocated_rep)>def_rep_per_thread):
             use_replicas=def_rep_per_thread
          else:
             use_replicas=(num_replicas-num_allocated_rep)
          sim_task = Task()
          sim_task.executable = ['python']
          
	  pre_exec_arr = [  'module unload PrgEnv-cray', 'module load PrgEnv-gnu','module unload bwpy','module load bwpy','module add bwpy-mpi', 'module add fftw', 'module add cray-netcdf', 'module add cudatoolkit/7.5.18-1.0502.10743.2.1', 'module add cmake', 'module unload darshan xalt','export CRAYPE_LINK_TYPE=dynamic', 'export CRAY_ADD_RPATH=yes', 'export FC=ftn', 'source /projects/sciteam/bamm/hruska/vpy2/bin/activate',
                                     'export tasks=md','export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1' ]
          #if cur_iter==0 and num_allocated_rep==0:
          #  pre_exec_arr = pre_exec_arr + [ 'mv %s']
          sim_task.pre_exec = pre_exec_arr
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
                                  '--trajstride', str(Kconfig.trajstride), '--idxstart',str(num_allocated_rep), '--idxend',str((num_allocated_rep+use_replicas)),
                                  '--path',combined_path,'--iter',str(cur_iter),
                                  '--md_steps',str(Kconfig.md_steps), '--save_traj', 'True','>', 'md.log']
          if Kconfig.md_use_xml=='yes':
            link_arr=['$SHARED/%s > run_openmm.py' % (os.path.basename(Kconfig.md_run_file)),
                      '$SHARED/system-5.xml > system-5.xml',
                      '$SHARED/integrator-5.xml > integrator-5.xml']            
          else:
            link_arr=['$SHARED/%s > run_openmm.py' % (os.path.basename(Kconfig.msd_run_file))]
          copy_arr=[]
          if cur_iter==0:
            for idx in range(num_allocated_rep, num_allocated_rep+use_replicas):
              copy_arr=copy_arr+['$SHARED/%s > %s/iter0_input%s.pdb' % (Kconfig.md_input_file, combined_path, idx)]           
            #if num_allocated_rep==0:
            #  copy_arr=copy_arr + ['$SHARED/%s > %s/%s' % (args.Kconfig,combined_path, args.Kconfig),
            #                         '$SHARED/run-tica-msm.py > %s/run-tica-msm.py' % combined_path,
            #                         '$SHARED/%s > %s/%s' % (Kconfig.md_run_file,combined_path,Kconfig.md_run_file)
            #                           ]

    
          #if cur_iter==0 and num_allocated_rep==0:
          #   copy_arr = copy_arr +['$SHARED/%s > %s/%s' % (args.Kconfig, combined_path, args.Kconfig)]
          sim_task.link_input_data = link_arr #+ copy_arr
          sim_task.copy_input_data = copy_arr
          if str(Kconfig.strategy)=='extend':
            copy_out=[]
            for idx in range(num_allocated_rep, num_allocated_rep+use_replicas):
              #copy_arr=copy_arr+['$SHARED/%s > iter0_input%s.pdb' % (Kconfig.md_input_file, idx)]
              copy_out=copy_out+['%s/iter%s_out%s.pdb > %s/iter%s_input%s.pdb' % (combined_path, cur_iter, idx, combined_path, (cur_iter+1), idx)]           
            sim_task.copy_output_data = copy_out  
            #if Kconfig.ndx_file is not None:
            #    sim_task.link_input_data.append('$SHARED/{0}'.format(os.path.basename(Kconfig.ndx_file)))
            
          num_allocated_rep=num_allocated_rep+use_replicas
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
        if str(Kconfig.strategy)!='extend':
          ana_stage = Stage()
          ana_task = Task()
          ana_task.pre_exec = [ 'module unload PrgEnv-cray','module load PrgEnv-gnu','module unload bwpy','module load bwpy/0.3.0','module add bwpy-mpi', 'module add fftw', 'module add cray-netcdf', 'module add cudatoolkit/7.5.18-1.0502.10743.2.1', 'module add cmake', 'module unload darshan xalt','export CRAYPE_LINK_TYPE=dynamic', 'export CRAY_ADD_RPATH=yes', 'export FC=ftn', 'source /projects/sciteam/bamm/hruska/vpy2/bin/activate', 'export tasks=tica_msm_ana',
                                   'export PYEMMA_NJOBS=1', 'export iter=%s' % cur_iter, 'export OMP_NUM_THREADS=1' ]
          ana_task.executable = ['python']
          ana_task.arguments = ['run-tica-msm.py', '--path',combined_path,'--n_select', str(num_replicas),'--cur_iter',str(cur_iter), '--Kconfig', str(args.Kconfig), '>', 'analyse.log']

          ana_task.cpu_reqs = { 'processes': 1,
                                    'process_type': None,
                                    'threads_per_process': 1,
                                    'thread_type': None
                                  }

          ana_task.link_input_data = ['$SHARED/run-tica-msm.py > run-tica-msm.py', '$SHARED/%s > %s'%(args.Kconfig,args.Kconfig)]
          
          #for sim_num in range(min(int(Kconfig.num_parallel_MD_sim),int(Kconfig.num_replicas))):
          ana_task.copy_output_data = ['analyse.log > %s/iter%s_analyse.log' % (combined_path, cur_iter)]

          #ana_task.copy_output_data = ['tmpha.gro > %s/iter_%s/tmpha.gro' % (combined_path,cur_iter),
           #                              'tmp.gro > %s/iter_%s/tmp.gro' % (combined_path,cur_iter)]
                                         #'tmp.gro > resource://iter_%s/tmp.gro' % cur_iter

          ana_task_ref = '$Pipeline_%s_Stage_%s_Task_%s'%(wf.uid, ana_stage.uid, ana_task.uid)
          ana_stage.add_tasks(ana_task)
          wf.add_stages(ana_stage)
        # --------------------------------------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------------------------------------
        # lsdmap:
        #     Purpose: Perform LSDMap on the large coordinate file to generate weights and eigen values.
        #     Arguments:
        #             config = name of the config file to be used during LSDMap
          
          #if(cur_iter % Kconfig.nsave == 0):
          #     post_ana_task.download_output_data = ['out.gro > output/iter_%s/out.gro' % cur_iter,
          #                                   'weight_out.w > output/iter_%s/weight_out.w' % cur_iter,
          #                                   'plot-scatter-cluster-10d.png > output/iter_%s/plot-scatter-cluster-10d.png' % (cur_iter),
          #                                   'ncopies.nc > output/iter_%s/ncopies.nc' % (cur_iter),
          #                                   '%s/iter_%s/tmp.gro > output/iter_%s/tmp.gro' % (combined_path,cur_iter,cur_iter) 
          #                                   ]

          #post_ana_task.copy_output_data = ['ncopies.nc > %s/iter_%s/ncopies.nc' % (combined_path,cur_iter),
          #                           'weight_out.w > %s/iter_%s/weight_out.w' % (combined_path,cur_iter),
          #                           'out.gro > %s/iter_%s/out.gro' % (combined_path,cur_iter),
          #                           'plot-scatter-cluster-10d.png > %s/iter_%s/plot-scatter-cluster-10d.png' % (combined_path,cur_iter),
          #                           'plot-scatter-cluster-10d-counts.png > %s/iter_%s/plot-scatter-cluster-10d-counts.png' % (combined_path,cur_iter),
          #                           'plot-scatter-cluster-10d-ncopiess.png > %s/iter_%s/plot-scatter-cluster-10d-ncopiess.png' % (combined_path,cur_iter)]

          #post_ana_task_ref = '$Pipeline_%s_Stage_%s_Task_%s'%(wf.uid, post_ana_stage.uid, post_ana_task.uid)

          #post_ana_stage.add_tasks(post_ana_task)
          #wf.add_stages(post_ana_stage)
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
        combined_path=str(Kconfig.remote_output_directory)
        wf = create_workflow(Kconfig, args)

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
            #'cores': Kconfig.PILOTSIZE,
            'cpus': Kconfig.NODESIZE*Kconfig.CPUs_per_NODE,
            #'cpu_processes': Kconfig.num_CUs_per_MD_replica,#PILOTSIZE,
            'gpus': Kconfig.NODESIZE,
            'project': Kconfig.ALLOCATION,
            'queue': Kconfig.QUEUE,
            'access_schema': 'gsissh'
          }	  
        else:
          print("use_gpus not recognized")
          
        print res_dict
        # Create Resource Manager object with the above resource description
        #rman = ResourceManager(res_dict)
        # Data common to multiple tasks -- transferred only once to common staging area
        shared_data_all = [args.Kconfig
                           ]
        if Kconfig.md_use_xml=='yes':
          shared_data_all=shared_data_all+['%s/system-5.xml' % Kconfig.md_dir,
                                           '%s/integrator-5.xml' % Kconfig.md_dir,
                                           Kconfig.md_run_dir+Kconfig.md_run_file,
                                           Kconfig.md_dir+Kconfig.md_input_file,
                                          '%s/run-tica-msm.py' % Kconfig.helper_scripts]
        else:
          shared_data_all=shared_data_all+[Kconfig.md_dir+Kconfig.md_input_file,
                                           Kconfig.md_dir+Kconfig.md_run_file,
                                          '%s/run-tica-msm.py' % Kconfig.helper_scripts]
        print "shared_data_all", shared_data_all 
       #if Kconfig.ndx_file is not None:
        #    rman.shared_data.append(Kconfig.ndx_file)

        # Create Application Manager, only one extasy script on one rabbit-mq server now
        appman = AppManager(hostname='two.radical-project.org', port=33134)#port=args.port)
        # appman = AppManager(port=) # if using docker, specify port here.
        appman.resource_desc = res_dict
        appman.shared_data = shared_data_all

        # Assign resource manager to the Application Manager
        #appman.resource_manager = rman

        # Assign the workflow as a set of Pipelines to the Application Manager
        appman.workflow = set([wf])

        # Run the Application Manager
        appman.run()

    except Exception as ex:

        print 'Error: {0}'.format(str(ex))
        print traceback.format_exc()
