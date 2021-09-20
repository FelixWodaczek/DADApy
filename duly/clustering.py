import multiprocessing
import time
import numpy as np
import scipy as sp

from duly.cython_ import cython_clustering as cf
from duly.density_estimation import DensityEstimation

cores = multiprocessing.cpu_count()


class Clustering(DensityEstimation):
    """This class contains various density-based clustering algorithms.

    Inherits from the DensityEstimation class.

    Attributes:

        N_clusters (int) : number of clusters found
        cluster_assignment (list(int)): cluster assignment. A list of length N containg the cluster assignment of each point
        as an integer from 0 to N_clusters-1.
        cluster_centers (int): Indices of the centroids of each cluster (density peak)
        cluster_indices (list(list(int))): a list of lists. Each sublist contains the indices belonging to the
        corresponding cluster.
        log_den_bord (array(float)): an array of dimensions N_clusters x N_clusters containg the estimated log density
        of the the saddle point between each couple of peaks.
        log_den_bord_err (array(float)): an array of dimensions N_clusters x N_clusters containg the estimated error
        on the log density of the saddle point between each couple of peaks.

    """

    def __init__(
        self, coordinates=None, distances=None, maxk=None, verbose=False, njobs=cores
    ):
        super().__init__(
            coordinates=coordinates,
            distances=distances,
            maxk=maxk,
            verbose=verbose,
            njobs=njobs,
        )

        self.cluster_indices = None
        self.N_clusters = None
        self.cluster_assignment = None
        self.cluster_centers = None
        self.log_den_bord_err = None
        self.log_den_bord = None

        self.delta = None  # Minimum distance from an element with higher density
        self.ref = None  # Index of the nearest element with higher density

    def compute_clustering(self, Z=1.65, halo=False):
        assert self.log_den is not None, "Compute density before clustering"
        if self.verb:
            print("Clustering started")

        # Make all values of log_den positives (this is important to help convergence)
        log_den_min = np.min(self.log_den)

        log_den_c = self.log_den
        log_den_c = log_den_c - log_den_min + 1

        # Putative modes of the PDF as preliminary clusters

        N = self.distances.shape[0]
        g = log_den_c - self.log_den_err
        # centers are point of max density  (max(g) ) within their optimal neighbolog_denod (defined by kstar)
        seci = time.time()

        out = cf._compute_clustering(
            Z,
            halo,
            self.kstar,
            self.dist_indices.astype(int),
            self.maxk,
            self.verb,
            self.log_den_err,
            log_den_min,
            log_den_c,
            g,
            N,
        )

        secf = time.time()

        self.cluster_indices = out[0]
        self.N_clusters = out[1]
        self.cluster_assignment = out[2]
        self.cluster_centers = out[3]
        out_bord = out[4]
        log_den_min = out[5]
        self.log_den_bord_err = out[6]

        self.log_den_bord = out_bord + log_den_min - 1 - np.log(N)

        if self.verb:
            print("Clustering finished, {} clusters found".format(self.N_clusters))
            print("total time is, {}".format(secf - seci))

    def compute_DecGraph(self):
        assert self.log_den is not None, "Compute density before"
        assert self.X is not None
        self.delta = np.zeros(self.N)
        self.ref = np.zeros(self.N, dtype="int")
        tt = np.arange(self.N)
        imax = []
        ncalls = 0
        for i in range(self.N):
            ll = tt[((self.log_den > self.log_den[i]) & (tt != i))]
            if ll.shape[0] > 0:
                a1, a2, a3 = np.intersect1d(
                    self.dist_indices[i, :], ll, return_indices=True, assume_unique=True
                )
                if a1.shape[0] > 0:
                    aa = np.min(a2)
                    self.delta[i] = self.distances[i, aa]
                    self.ref[i] = self.dist_indices[i, aa]
                else:
                    ncalls = ncalls + 1
                    dd = self.X[((self.log_den > self.log_den[i]) & (tt != i))]
                    ds = np.transpose(
                        sp.spatial.distance.cdist([np.transpose(self.X[i, :])], dd)
                    )
                    j = np.argmin(ds)
                    self.ref[i] = ll[j]
                    self.delta[i] = ds[j]
            else:
                self.delta[i] = -100.0
                imax.append(i)
        self.delta[imax] = 1.05 * np.max(self.delta)
        print("Number of points for which self.delta needed call to cdist=", ncalls)

    def compute_cluster_DP(self, dens_cut=0.0, delta_cut=0.0, halo=False):
        assert self.delta is not None
        ordered = np.argsort(-self.log_den)
        self.cluster_assignment = np.zeros(self.N, dtype="int")
        tt = np.arange(self.N)
        center_label = np.zeros(self.N, dtype="int")
        ncluster = -1
        for i in range(self.N):
            j = ordered[i]
            if (self.log_den[j] > dens_cut) & (self.delta[j] > delta_cut):
                ncluster = ncluster + 1
                self.cluster_assignment[j] = ncluster
                center_label[j] = ncluster
            else:
                self.cluster_assignment[j] = self.cluster_assignment[self.ref[j]]
                center_label[j] = -1
        self.centers = tt[(center_label != -1)]
        if halo:
            bord = np.zeros(self.N, dtype="int")
            halo = np.copy(self.cluster_assignment)

            for i in range(self.N):
                for j in self.dist_indices[i, :][(self.distances[i, :] <= self.dc[i])]:
                    if self.cluster_assignment[i] != self.cluster_assignment[j]:
                        bord[i] = 1
            halo_cutoff = np.zeros(ncluster + 1)
            halo_cutoff[:] = np.min(self.log_den) - 1
            for j in range(ncluster + 1):
                td = self.log_den[((bord == 1) & (self.cluster_assignment == j))]
                if td.size != 0:
                    halo_cutoff[j] = np.max(td)
            halo[tt[(self.log_den < halo_cutoff[self.cluster_assignment])]] = -1
            self.cluster_assignment = halo

    def compute_clustering_pure_python(self, Z=1.65, halo=False):
        assert self.log_den is not None, "Compute density before clustering"
        if self.verb:
            print("Clustering started")

        # Make all values of log_den positives (this is important to help convergence)
        log_den_min = np.min(self.log_den)
        log_den_c = self.log_den
        log_den_c = log_den_c - log_den_min + 1

        # Putative modes of the PDF as preliminary clusters
        sec = time.time()
        N = self.distances.shape[0]
        g = log_den_c - self.log_den_err
        centers = []
        for i in range(N):
            t = 0
            for j in range(1, self.kstar[i] + 1):
                if g[i] < g[self.dist_indices[i, j]]:
                    t = 1
                    break
            if t == 0:
                centers.append(i)
        for i in centers:
            l, m = np.where(self.dist_indices == i)
            for j in range(l.shape[0]):
                if (g[l[j]] > g[i]) & (m[j] <= self.kstar[l[j]]):
                    centers.remove(i)
                    break
        cluster_init = []
        for j in range(N):
            cluster_init.append(-1)
            Nclus = len(centers)
        if self.verb:
            print("Number of clusters before multimodality test=", Nclus)

        for i in centers:
            cluster_init[i] = centers.index(i)
        sortg = np.argsort(-g)  # Get the rank of the elements in the g vector
        # sorted in descendent order.
        # Perform preliminar assignation to clusters
        for j in range(N):
            ele = sortg[j]
            nn = 0
            while cluster_init[ele] == -1:
                nn = nn + 1
                cluster_init[ele] = cluster_init[self.dist_indices[ele, nn]]
        clstruct = []  # useful list of points in the clusters
        for i in range(Nclus):
            x1 = []
            for j in range(N):
                if cluster_init[j] == i:
                    x1.append(j)
            clstruct.append(x1)
        sec2 = time.time()
        if self.verb:
            print(
                "{0:0.2f} seconds clustering before multimodality test".format(
                    sec2 - sec
                )
            )
        log_den_bord = np.zeros((Nclus, Nclus), dtype=float)
        log_den_bord_err = np.zeros((Nclus, Nclus), dtype=float)
        Point_bord = np.zeros((Nclus, Nclus), dtype=int)
        # Find border points between putative clusters
        sec = time.time()
        for i in range(Nclus):
            for j in range(Nclus):
                Point_bord[i][j] = -1
        for c in range(Nclus):
            for p1 in clstruct[c]:
                for k in range(1, self.kstar[p1] + 1):
                    p2 = self.dist_indices[p1, k]
                    pp = -1
                    if cluster_init[p2] != c:
                        pp = p2
                        cp = cluster_init[pp]
                        break
                if pp != -1:
                    for k in range(1, self.maxk):
                        po = self.dist_indices[pp, k]
                        if po == p1:
                            break
                        if cluster_init[po] == c:
                            pp = -1
                            break
                if pp != -1:
                    if g[p1] > log_den_bord[c][cp]:
                        log_den_bord[c][cp] = g[p1]
                        log_den_bord[cp][c] = g[p1]
                        Point_bord[cp][c] = p1
                        Point_bord[c][cp] = p1
                    # if (g[pp]>log_den_bord[c][cp]):
                    #     log_den_bord[c][cp]=g[pp]
                    #     log_den_bord[cp][c]=g[pp]
                    #     Point_bord[cp][c]=pp
                    #     Point_bord[c][cp]=pp
        for i in range(Nclus - 1):
            for j in range(i + 1, Nclus):
                if Point_bord[i][j] != -1:
                    log_den_bord[i][j] = log_den_c[Point_bord[i][j]]
                    log_den_bord[j][i] = log_den_c[Point_bord[j][i]]
                    log_den_bord_err[i][j] = self.log_den_err[Point_bord[i][j]]
                    log_den_bord_err[j][i] = self.log_den_err[Point_bord[j][i]]
        for i in range(Nclus):
            log_den_bord[i][i] = -1.0
            log_den_bord_err[i][i] = 0.0
        sec2 = time.time()
        if self.verb:
            print("{0:0.2f} seconds identifying the borders".format(sec2 - sec))
        check = 1
        clsurv = []
        sec = time.time()
        for i in range(Nclus):
            clsurv.append(1)
        while check == 1:
            pos = []
            ipos = []
            jpos = []
            check = 0
            for i in range(Nclus - 1):
                for j in range(i + 1, Nclus):
                    a1 = log_den_c[centers[i]] - log_den_bord[i][j]
                    a2 = log_den_c[centers[j]] - log_den_bord[i][j]
                    e1 = Z * (self.log_den_err[centers[i]] + log_den_bord_err[i][j])
                    e2 = Z * (self.log_den_err[centers[j]] + log_den_bord_err[i][j])
                    if a1 < e1 or a2 < e2:
                        check = 1
                        pos.append(log_den_bord[i][j])
                        ipos.append(i)
                        jpos.append(j)
            if check == 1:
                barriers = pos.index(max(pos))
                imod = ipos[barriers]
                jmod = jpos[barriers]
                if log_den_c[centers[imod]] < log_den_c[centers[jmod]]:
                    tmp = jmod
                    jmod = imod
                    imod = tmp
                clsurv[jmod] = 0
                log_den_bord[imod][jmod] = -1.0
                log_den_bord[jmod][imod] = -1.0
                log_den_bord_err[imod][jmod] = 0.0
                log_den_bord_err[jmod][imod] = 0.0
                clstruct[imod].extend(clstruct[jmod])
                clstruct[jmod] = []
                for i in range(Nclus):
                    if i != imod and i != jmod:
                        if log_den_bord[imod][i] < log_den_bord[jmod][i]:
                            log_den_bord[imod][i] = log_den_bord[jmod][i]
                            log_den_bord[i][imod] = log_den_bord[imod][i]
                            log_den_bord_err[imod][i] = log_den_bord_err[jmod][i]
                            log_den_bord_err[i][imod] = log_den_bord_err[imod][i]
                        log_den_bord[jmod][i] = -1
                        log_den_bord[i][jmod] = log_den_bord[jmod][i]
                        log_den_bord_err[jmod][i] = 0
                        log_den_bord_err[i][jmod] = log_den_bord_err[jmod][i]
        sec2 = time.time()
        if self.verb:
            print("{0:0.2f} seconds with multimodality test".format(sec2 - sec))
        N_clusters = 0
        cluster_indices = []
        cluster_centers = []
        nnum = []
        for j in range(Nclus):
            nnum.append(-1)
            if clsurv[j] == 1:
                nnum[j] = N_clusters
                N_clusters = N_clusters + 1
                cluster_indices.append(clstruct[j])
                cluster_centers.append(centers[j])
        log_den_bord_m = np.zeros((N_clusters, N_clusters), dtype=float)
        log_den_bord_err = np.zeros((N_clusters, N_clusters), dtype=float)
        for j in range(Nclus):
            if clsurv[j] == 1:
                jj = nnum[j]
                for k in range(Nclus):
                    if clsurv[k] == 1:
                        kk = nnum[k]
                        log_den_bord_m[jj][kk] = log_den_bord[j][k]
                        log_den_bord_err[jj][kk] = log_den_bord_err[j][k]
        Last_cls = np.empty(N, dtype=int)
        for j in range(N_clusters):
            for k in cluster_indices[j]:
                Last_cls[k] = j
        Last_cls_halo = np.copy(Last_cls)
        nh = 0
        for j in range(N_clusters):
            log_den_halo = max(log_den_bord_m[j])
            for k in cluster_indices[j]:
                if log_den_c[k] < log_den_halo:
                    nh = nh + 1
                    Last_cls_halo[k] = -1
        if halo:
            cluster_assignment = Last_cls_halo
        else:
            cluster_assignment = Last_cls
        out_bord = np.copy(log_den_bord_m)
        self.cluster_indices = cluster_indices
        self.N_clusters = N_clusters
        self.cluster_assignment = cluster_assignment
        self.cluster_centers = cluster_centers
        self.log_den_bord = (
            out_bord + log_den_min - 1 - np.log(N)
        )  # remove wrong normalisation introduced earlier
        self.log_den_bord_err = log_den_bord_err
        if self.verb:
            print("Clustering finished, {} clusters found".format(self.N_clusters))
