import sys, os
print(os.path.dirname(sys.executable))
import os
import argparse
import numpy as np
import pyemma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class SelectionStep(object):
    """
    SelectionStep()

    A class used to perform the selection step of Extended DM-d-MD 
    """

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="select n new configurations uniformily along the first and second DCs..")
        parser.add_argument('npoints', metavar='n', type=int, help='number of configurations to be selected')

        # required options
        parser.add_argument("-eg",
           type=str,
           dest="egfile",
           required=True,
           help="File containing the eigenvalues: eg")

        parser.add_argument("-ev",
           type=str,
           dest="evfile",
           required=False,
           help="File containing the eigenvectors (DC's) of each configuration (input): ev")

        parser.add_argument("-o",
           type=str,
           dest="ncfile",
           required=False,
           help='File containing a single column with the number of copies to be made for each configuration (output, opt.): nc')

        return parser

   
    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()

        #evfile = reader.open(args.evfile)
        #evs = evfile.readlines()
        
        evs = np.loadtxt(args.evfile)#"tmpha.ev")
        scale_time=False
        scale_eg=True
        if scale_time or scale_eg:
          egs = np.loadtxt(args.egfile)#"tmpha.eg")
        num_points = evs.shape[0]
        num_evs = evs.shape[1]
        #rescale
        for i in range(1,min(5,num_evs)):
          if scale_time:
            evs[:,i]=np.sqrt(-1./np.log(np.abs(egs[i])))*((evs[:,i]-evs[:,i].min())/(evs[:,i].max()-evs[:,i].min())-0.5)
          elif scale_eg:
            evs[:,i]=np.abs(egs[i])*((evs[:,i]-evs[:,i].min())/(evs[:,i].max()-evs[:,i].min())-0.5)
          else:
            evs[:,i]=2.*((evs[:,i]-evs[:,i].min())/(evs[:,i].max()-evs[:,i].min())-0.5)
            #print(i, evs[:,i].max()-evs[:,i].min())

        cluster_obj=pyemma.coordinates.cluster_kmeans(data = evs[:,1:], max_iter =20)# changing max_iter fails

        assigns=cluster_obj.get_output()[0][:,0]
        assigns.shape

        counts=np.zeros(assigns.max()+1)
        for i in range(assigns.max()+1):
          counts[i]=len(np.where(assigns==i)[0])

        assign_counts=counts[assigns]

        plt.clf()
        plt.xlabel("DC1")
        plt.ylabel("DC2")
        cp = plt.scatter(evs[:,1], evs[:,2], s=10, c=assigns, marker='o', linewidth=0.,cmap='jet')
        plt.colorbar(label="cluster id")
        plt.savefig('plot-scatter-cluster-10d.png', bbox_inches='tight', dpi=200)

        plt.clf()
        plt.xlabel("DC1")
        plt.ylabel("DC2")
        cp = plt.scatter(evs[:,1], evs[:,2], s=10, c=assign_counts,marker='o', linewidth=0.,cmap='jet')
        plt.colorbar(label="counts")
        plt.savefig('plot-scatter-cluster-10d-counts.png', bbox_inches='tight', dpi=200)

        n_states=args.npoints
        inv_values = 1.0 / counts
        p = inv_values / np.sum(inv_values)
        cluster_picks=np.random.choice(np.arange(assigns.max()+1), p = p, size=n_states)

        ncopiess = np.zeros(num_points, dtype='int')
        structures_pick=np.zeros(n_states, dtype='int')
        for i, cluster_pick in enumerate(cluster_picks):
          #print i, cluster_pick
          structures_pick[i]=np.random.choice(np.where(assigns==cluster_pick)[0])
          ncopiess[int(structures_pick[i])] +=1


        if args.ncfile is None:
          args.ncfile = os.path.splitext(args.evfile)[0] + '.nc'

        np.savetxt(args.ncfile, ncopiess, fmt='%i') 

        print "number of copies stored in %s" %args.ncfile

        plt.clf()
        plt.xlabel("DC1")
        plt.ylabel("DC2")
        cp = plt.scatter(evs[structures_pick,1], evs[structures_pick,2], s=10, marker='o', linewidth=0.,cmap='jet')
        plt.savefig('plot-scatter-cluster-10d-ncopiess.png', bbox_inches='tight', dpi=200)




if __name__ == '__main__':
    SelectionStep().run()
