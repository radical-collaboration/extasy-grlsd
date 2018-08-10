import sys, os
print(os.path.dirname(sys.executable))
import os
import argparse
import numpy as np
import pyemma
import matplotlib
matplotlib.use('Agg')
import glob
import mdtraj as md
import imp
from matplotlib.pyplot import *
from pyemma import plots
matplotlib.rcParams.update({'font.size': 14})
print("pyemma version",pyemma.__version__)

class Runticamsm(object):
    """
    A class used to perform the TICA and MSM
    """

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="run TICA, MSM")

        # required options
        parser.add_argument("--path",
           type=str,
           dest="path",
           required=True,
           help="Full path of all files input, trajs, out")

        parser.add_argument("--cur_iter",
           type=int,
           dest="cur_iter",
           required=True)

        parser.add_argument("--n_select",
           type=int,
           dest="n_select",
           required=True)
        
	parser.add_argument("--Kconfig",
           type=str,
           dest="Kconfig",
           required=True)

        return parser

    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()
        
	
	#parser = argparse.ArgumentParser()
        #parser.add_argument('--Kconfig', help='link to Kernel configurations file')
        #parser.add_argument('--port', dest="port", help='port for RabbitMQ server', default=5672, type=int)
        #args = parser.parse_args()

        Kconfig = imp.load_source('Kconfig', args.Kconfig)


        pdb_file=glob.glob(args.path+'/iter*_input*.pdb')[0]
        #pdb_file=glob.glob('iter*_input*.pdb')[0]
        #traj_files=glob.glob(args.path+'/iter*_traj*.dcd')
        p_cont=True
        p_iter=0
        traj_files=[]
        iter_arr=[]
        while(p_cont):
           traj_files_tmp=glob.glob(args.path+'/iter'+str(p_iter)+'_traj*.dcd')
           traj_files.sort()
           if len(traj_files_tmp)==0:
             p_cont=False
           else:
             print("iter", str(p_iter), " # files", str(len(traj_files_tmp))) 
             traj_files=traj_files+traj_files_tmp
             iter_arr=iter_arr+[p_iter]*len(traj_files_tmp)
             p_iter=p_iter+1

        p_iter_max=p_iter-1
        iter_arr=np.array(iter_arr)
        #traj_files=glob.glob('iter*_traj*.dcd')
        traj_files.sort()
        topfile = md.load(pdb_file)
        featurizer = pyemma.coordinates.featurizer(topfile)
        featurizer.add_residue_mindist(residue_pairs='all', scheme='closest-heavy')
        featurizer.add_backbone_torsions(cossin=True)
        featurizer.dimension()

        inp = pyemma.coordinates.source(traj_files, featurizer)
        #inp.get_output()

        tica_lag=Kconfig.tica_lag#1
        tica_dim=10
        tica_stride=1

        tica_obj = pyemma.coordinates.tica(inp, lag=tica_lag, dim=tica_dim, kinetic_map=True, stride=tica_stride)

        y = tica_obj.get_output()
        #y[0].shape

        msm_states=Kconfig.msm_states
        msm_stride=1
        msm_lag=Kconfig.msm_lag#1
        cl = pyemma.coordinates.cluster_kmeans(data=y, k=msm_states, max_iter=50, stride=msm_stride)
        
        m = pyemma.msm.estimate_markov_model(cl.dtrajs, msm_lag)



        print("n atoms",topfile.n_atoms)
        print("n frames total",inp.n_frames_total())
        print("n trajs",inp.number_of_trajectories())
        print(" traj lengths", inp.trajectory_lengths())
        print(" input dimension",inp.dimension())

        print("TICA eigenvalues", tica_obj.eigenvalues)
        #print(tica_obj.eigenvectors)

        print("MSM eigenvalues",m.eigenvalues(10))
        #print(m.eigenvectors_left(10))
        #print(m.eigenvectors_right(10))
        print("MSM P connected",m.P)  #only connected

        print("MSM clustercenters",cl.clustercenters)
        
        print("TICA timescales",tica_obj.timescales)
        print("MSM timescales", m.timescales(10))
        print("MSM stat", m.stationary_distribution)
        print("MSM active set", m.active_set)
        print('fraction of states used = ', m.active_state_fraction)
        print('fraction of counts used = ', m.active_count_fraction)


        c = m.count_matrix_full
        s =  np.sum(c, axis=1)
        print("count matrix sums",s)

	if Kconfig.strategy=='cmicro':
          if 0 not in s:
              q = 1.0 / s

        elif Kconfig.strategy=='cmacro':
           q = s

	else:
 	  print('ERROR strategy not recognized')

        n_states=c.shape[0]

        dtrajs = [ t for t in cl.dtrajs ]

        print("msm dtrajs", dtrajs)
        #get frame_list 


        frame_state_list = {n: [] for n in range(n_states)}
        for nn, dt in enumerate(dtrajs):
            for mm, state in enumerate(dt):
                    frame_state_list[state].append((nn,mm))

        for k in range(n_states):
         if len(frame_state_list[k]) == 0:
            print('removing state '+str(k)+' no frames')
            q[k] = 0.0

                    # and normalize the remaining one
        q /= np.sum(q)


        n_pick=int(args.n_select)#100

        state_picks = np.random.choice(np.arange(len(q)), size=n_pick, p=q)


        picks = [
            frame_state_list[state][np.random.randint(0,
            len(frame_state_list[state]))]
            for state in state_picks
            ]




        traj_select = [traj_files[pick[0]] for pick in picks]
        frame_select = [pick[1]*tica_stride*msm_stride for pick in picks]
        print('traj_select',traj_select)
        print('frame_select',traj_select)

        text_file = open(args.path + "/traj_select.txt", "w")
        for idx in range(n_pick):
          text_file.write(traj_select[idx]+' to iter '+str(args.cur_iter)+' idx '+str(idx)+' \n')

        text_file.close()


        # write new input files from frames


        for idx in range(n_pick):
          tmp =md.load(args.path+'/iter0_input0.pdb')
          files = md.load(traj_select[idx], top=args.path+'/iter0_input0.pdb')
          tmp.xyz[0,:,:]=files.xyz[frame_select[idx],:,:]
          tmp.save_pdb(args.path+'/iter'+str(args.cur_iter+1)+'_input'+str(idx)+'.pdb')

        #rg rmsd

        original_file = md.load(args.path+'/iter0_input0.pdb')
        out_files=glob.glob(args.path+'/iter*_out*.pdb')
        out_files.sort()

        q_thres=0.5
        rg_arr=[]
        rmsd_arr=[]
        q_arr=[]
        for file in out_files:
          file2 = md.load(file)
          rg_arr.append(md.compute_rg(file2)[0])
          rmsd_arr.append(md.rmsd(file2,original_file)[0])
          dist_arr=md.compute_contacts(file2)[0][0]
          q_arr.append(dist_arr[dist_arr<q_thres].shape[0])


        rg_arr=np.array(rg_arr)
        print("rg values", rg_arr.min(), rg_arr.max(), rg_arr)
        rmsd_arr=np.array(rmsd_arr)
        print("rmsd values", rmsd_arr.min(), rmsd_arr.max(), rmsd_arr)

        q_arr=np.array(q_arr)
        print("Q values", q_arr.min(), q_arr.max(), q_arr)


        ########################################

        tica0=np.array([])
        tica1=np.array([])
        for i in range(len(y)):
          tica0=np.append(tica0,y[i][:,0])
          tica1=np.append(tica1,y[i][:,1])


        clf()
        fig=figure()
        ax = fig.add_subplot(111)
        ax.scatter(np.arange(tica_obj.timescales.shape[0]),tica_obj.timescales)
        ax.set_ylabel('TICA Timescales (steps)')
        ax.set_xlabel('# TICA eigenvector')
        ax.set_yscale('log')
        savefig(args.path+'/plot_iter9_tica_timescales.png', bbox_inches='tight', dpi=200)

        cumvar = np.cumsum(tica_obj.timescales)
        cumvar /= cumvar[-1]
        clf()
        plot(cumvar, linewidth=2)
        for thres in [0.5,0.8,0.95]:
          threshold_index=np.argwhere(cumvar > thres)[0][0]
          print thres, threshold_index
          vlines(threshold_index, 0.0, 1.0, linewidth=2)
          hlines(thres, 0, cumvar.shape[0], linewidth=2)

        xlabel('Eigenvalue Number', fontsize = 16)
        ylabel('cumulative kinetic content', fontsize = 16)
        savefig(args.path+'/plot_iter9_tica_cumulative_kinetic_content.png', bbox_inches='tight', dpi=200)

        msm_timescales=m.timescales(100)
        clf()
        fig=figure()
        ax = fig.add_subplot(111)
        ax.scatter(np.arange(msm_timescales.shape[0]),msm_timescales*tica_stride)
        ax.set_ylabel('MSM Timescales (steps)')
        ax.set_xlabel('# MSM eigenvector')
        ax.set_yscale('log')
        savefig(args.path+'/plot_iter9_msm_timescales.png', bbox_inches='tight', dpi=200)

        cumvar = np.cumsum(m.timescales(100))
        cumvar /= cumvar[-1]
        clf()
        plot(cumvar, linewidth=2)
        for thres in [0.5,0.8,0.95]:
          threshold_index=np.argwhere(cumvar > thres)[0][0]
          print thres, threshold_index
          vlines(threshold_index, 0.0, 1.0, linewidth=2)
          hlines(thres, 0, cumvar.shape[0], linewidth=2)

        xlabel('Eigenvalue Number', fontsize = 16)
        ylabel('cumulative kinetic content', fontsize = 16)
        savefig(args.path+'/plot_iter9_msm_cumulative_kinetic_content.png', bbox_inches='tight', dpi=200)





        clf()
        xlabel("TICA ev0")
        ylabel("TICA ev1")
        cp = scatter(tica0, tica1, s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='MSM states')
        savefig(args.path+'/plot_iter9_tica_evs.png', bbox_inches='tight', dpi=200)


        clf()
        xlabel("TICA ev0")
        ylabel("TICA ev1")
        fig, ax = plots.plot_free_energy(tica0, tica1,cmap='Spectral')
        savefig(args.path+'/plot_iter9_tica_evs2.png', bbox_inches='tight', dpi=200)

        clf()
        xlabel("TICA ev0")
        ylabel("TICA ev1")
        fig, ax = plots.plot_free_energy(tica0, tica1,cmap='Spectral')
        cp = scatter(cl.clustercenters[:,0], cl.clustercenters[:,1], s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='MSM state centers')
        legend()
        savefig(args.path+'/plot_iter9_tica_evs3_centers.png', bbox_inches='tight', dpi=200)

        #plot msm ev
        clf()
        xlabel("MSM ev1")
        ylabel("MSM ev2")
        cp = scatter(m.eigenvectors_right(10)[:,1], m.eigenvectors_right(10)[:,2], s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='MSM states')
        savefig(args.path+'/plot_iter9_msm_evs.png', bbox_inches='tight', dpi=200)

        #plot msm ev
        clf()
        xlabel("MSM ev1")
        ylabel("MSM ev2")
        fig, ax = plots.plot_free_energy(m.eigenvectors_right(10)[:,1], m.eigenvectors_right(10)[:,2], cmap='Spectral', weights=m.stationary_distribution)
        savefig(args.path+'/plot_iter9_msm_evs2.png', bbox_inches='tight', dpi=200)

        clf()
        xlabel("RMSD")
        ylabel("Rg")
        cp = scatter(rmsd_arr, rg_arr, s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='MSM states')
        savefig(args.path+'/plot_iter9_rgrmsd.png', bbox_inches='tight', dpi=200)

        #plot msm ev
        clf()
        xlabel("RMSD")
        ylabel("Rg")
        fig, ax = plots.plot_free_energy(rmsd_arr, rg_arr, cmap='Spectral')
        savefig(args.path+'/plot_iter9_rgrmsd2.png', bbox_inches='tight', dpi=200)

        clf()
        xlabel("Q")
        ylabel("Rg")
        cp = scatter(q_arr, rg_arr, s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='MSM states')
        savefig(args.path+'/plot_iter9_qrg.png', bbox_inches='tight', dpi=200)

        #Q 1d free energy
        clf()
        z, x = np.histogram(q_arr, bins=10)
        F = -np.log(z)
        F=F-F.min()
        plot(x[1:], F)
        xlabel('Q', fontsize = 15)
        ylabel('Free Energy [kT]', fontsize =15)
        savefig(args.path+'/plot_iter9_free_energy_q.png', bbox_inches='tight', dpi=200)

        #MSM 1d free energy
        clf()
        n_step=int(m.P.shape[0]/10)
        bins=np.sort(m.eigenvectors_right(10)[:,1])[::n_step]
        bins=np.append(bins,np.sort(m.eigenvectors_right(10)[:,1])[-1])
        z, x = np.histogram(m.eigenvectors_right(10)[:,1], weights=m.stationary_distribution, density=True, bins=bins)
        F = -np.log(z)
        F=F-F.min()
        plot(x[1:], F)
        xlabel('MSM ev1', fontsize = 15)
        ylabel('Free Energy [kT]', fontsize =15)
        savefig(args.path+'/plot_iter9_msm_free_energy.png', bbox_inches='tight', dpi=200)


        #which tica frames seleted


        tica0_sel=np.array([])
        tica1_sel=np.array([])
        for i in range(n_pick):
          tica0_sel=np.append(tica0_sel,y[picks[i][0]][frame_select[i],0])
          tica1_sel=np.append(tica1_sel,y[picks[i][0]][frame_select[i],1])

        clf()
        xlabel("TICA ev0")
        ylabel("TICA ev1")
        cp = scatter(tica0, tica1, s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='all frames')
        cp = scatter(tica0_sel, tica1_sel, s=10, c='red', marker='o', linewidth=0.,cmap='jet', label='selected')
        legend()
        savefig(args.path+'/plot_iter9_tica_evs4_selected.png', bbox_inches='tight', dpi=200)





        #m.ck_test

        ck=m.cktest(2)

        clf()
        pyemma.plots.plot_cktest(ck, diag=True, figsize=(7,7), layout=(2,2), padding_top=0.1, y01=False, padding_between=0.3, dt=0.1, units='ns')
        savefig(args.path+'/plot_iter9_msm_cktest.png')

        its = pyemma.msm.its(dtrajs, nits=10)
        clf()
        pyemma.plots.plot_implied_timescales(its, ylog=False, units='steps', linewidth=2)
        #xlim(0, 40); ylim(0, 120);
        savefig(args.path+'/plot_iter9_msm_its.png')



        #which msm states selected
        #warning m only connected, c full -selected
        #m.active_set
        #state_picks
        #msm_states
        p_picks_active=[]
        for i in state_picks:
          if i in m.active_set:
            p_picks_active.append(np.argwhere(i==m.active_set)[0][0])

        p_picks_active=np.unique(np.array(p_picks_active)).astype(int)
          



        clf()
        xlabel("MSM ev1")
        ylabel("MSM ev2")
        cp = scatter(m.eigenvectors_right(10)[:,1], m.eigenvectors_right(10)[:,2], s=10, c='blue', marker='o', linewidth=0.,cmap='jet', label='MSM states')
        cp = scatter(m.eigenvectors_right(10)[p_picks_active,1], m.eigenvectors_right(10)[p_picks_active,2], s=10, c='red', marker='o', linewidth=0.,cmap='jet', label='selected')
        legend(loc='center left', bbox_to_anchor=(1, 0.5))
        savefig(args.path+'/plot_iter9_msm_evs_4_select.png', bbox_inches='tight', dpi=200)




        p_states=np.array([])
        p_unique=[]
        for p_iter in range(p_iter_max,-1,-1):
            p_arr=np.argwhere(iter_arr==p_iter)
            for i in p_arr:
              #print i[0]
              p_states=np.append(p_states,dtrajs[i[0]])
            p_states=np.unique(p_states).astype(int)
            p_unique.append(p_states.shape[0])

        p_unique=np.array(p_unique)

        clf()
        fig=figure()
        ax = fig.add_subplot(111)
        ax.scatter(np.arange(p_unique.shape[0]),p_unique)
        ax.set_ylabel('# of current msm states explored')
        ax.set_xlabel('iteration')
        #ax.set_yscale('log')
        savefig(args.path+'/plot_iter9_strategy.png', bbox_inches='tight', dpi=200)


        clf()
        xlabel("TICA ev0")
        ylabel("TICA ev1")

        for p_iter in range(p_iter_max,-1,-1):
            p_arr=np.argwhere(iter_arr==p_iter)
            tica0=np.array([])
            tica1=np.array([])
            for i in p_arr:
              #print i[0]
              tica0=np.append(tica0,y[i[0]][:,0])
              tica1=np.append(tica1,y[i[0]][:,1])
            cp = scatter(tica0, tica1, s=10, marker='o', linewidth=0.,cmap='jet', label='iter '+str(p_iter))

        legend(loc='center left', bbox_to_anchor=(1, 0.5))
        savefig(args.path+'/plot_iter9_tica_evs5_iters.png', bbox_inches='tight', dpi=200)


        clf()
        xlabel("MSM ev1")
        ylabel("MSM ev2")
        for p_iter in range(p_iter_max,-1,-1):
            p_arr=np.argwhere(iter_arr==p_iter)
            p_states=np.array([])
            for i in p_arr:
              #print i[0]
              p_states=np.append(p_states,dtrajs[i[0]])
            p_states=np.unique(p_states).astype(int)
            p_states_active=[]
            for i in p_states:
              if i in m.active_set:
                p_states_active.append(np.argwhere(i==m.active_set)[0][0])
            p_states_active=np.unique(np.array(p_states_active)).astype(int)
            cp = scatter(m.eigenvectors_right(10)[p_states_active,1], m.eigenvectors_right(10)[p_states_active,2], s=10,  marker='o', linewidth=0., cmap='spectral', label='iter '+str(p_iter))

        legend(loc='center left', bbox_to_anchor=(1, 0.5))
        savefig(args.path+'/plot_iter9_msm_evs_3_iter.png', bbox_inches='tight', dpi=200)



if __name__ == '__main__':
    Runticamsm().run()
