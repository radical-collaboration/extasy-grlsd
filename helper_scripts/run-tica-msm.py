import sys, os
print(os.path.dirname(sys.executable))
import os
import argparse
import numpy as np
import pyemma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import mdtraj as md


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

        return parser

    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()
        
        pdb_file=glob.glob(args.path+'/iter*_input*.pdb')[0]
        #pdb_file=glob.glob('iter*_input*.pdb')[0]
        traj_files=glob.glob(args.path+'/iter*_traj*.dcd')
        #traj_files=glob.glob('iter*_traj*.dcd')
        traj_files.sort()
        topfile = md.load(pdb_file)
        featurizer = pyemma.coordinates.featurizer(topfile)
        featurizer.add_residue_mindist(residue_pairs='all', scheme='closest-heavy')
        featurizer.add_backbone_torsions(cossin=True)
        featurizer.dimension()

        inp = pyemma.coordinates.source(traj_files, featurizer)
        #inp.get_output()

        tica_lag=1
        tica_dim=10
        tica_stride=1

        tica_obj = pyemma.coordinates.tica(inp, lag=tica_lag, dim=tica_dim, kinetic_map=True, stride=tica_stride)

        y = tica_obj.get_output()
        #y[0].shape

        msm_states=100
        msm_stride=1
        msm_lag=1
        cl = pyemma.coordinates.cluster_kmeans(data=y, k=msm_states, max_iter=50, stride=msm_stride)
        
        m = pyemma.msm.estimate_markov_model(cl.dtrajs, msm_lag)

        #plot msm ev
        plt.clf()
        plt.xlabel("MSM ev1")
        plt.ylabel("MSM ev2")
        cp = plt.scatter(m.eigenvectors_right(10)[:,1], m.eigenvectors_right(10)[:,2], s=10, c='blue', marker='o', linewidth=0.,cmap='jet')
        plt.savefig(args.path+'/plot_iter0_msm_evs.png', bbox_inches='tight', dpi=200)



        topfile.n_atoms
        inp.n_frames_total()
        inp.number_of_trajectories()
        inp.trajectory_lengths()
        inp.dimension()

        tica_obj.eigenvalues
        tica_obj.eigenvectors

        m.eigenvalues(10)
        m.eigenvectors_left(10)
        m.eigenvectors_right(10)
        m.P  #only connected

        cl.clustercenters

        c = m.count_matrix_full
        s =  np.sum(c, axis=1)
        print(s)
        if 0 not in s:
            q = 1.0 / s

        n_states=c.shape[0]

        dtrajs = [ t for t in cl.dtrajs ]


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



if __name__ == '__main__':
    Runticamsm().run()
