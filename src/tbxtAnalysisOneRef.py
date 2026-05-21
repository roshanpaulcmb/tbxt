#### IMPORTS ####
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, pca
import mdtraj as mdt
from scipy.ndimage import gaussian_filter
import os
import time
import umap
import hdbscan
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA

#### FUNCTIONS ####

def setup(args):
    data = { }
    data["pdb"] = args.pdb
    data["dcd"] = args.dcd
    data["outputName"] = args.outputName
    return data

def makeU(data):
    print("\nMaking universe..")
    data["u"] = mda.Universe(data["pdb"], data["dcd"])
    data["protein"] = data["u"].select_atoms("protein")
    
    # manual PCA requires mdt.load
    # data["t"] = mdt.load(files["dcd"], top = files["pdb"])
    # data["t"].superpose(data["t"], 0) # align with 0th frame
    return data

def calcRmsd(data, verbose = True):
    print("\nCalculating RMSD..")
    rFirst = rms.RMSD(data["protein"],
                 data["protein"],
                 select = "name CA",
                 ref_frame = 0).run(verbose = verbose)
    data["rFirstM"] = rFirst.results.rmsd.T
    
    rLast = rms.RMSD(data["protein"],
                 data["protein"],
                 select = "name CA",
                 ref_frame = -1).run(verbose = verbose)
    data["rLastM"] = rLast.results.rmsd.T
    
    # Roshan was checking largest RMSD jump
    # maxDelta = 0
    # jumpI = 0
    # for i in range(1, len(data["rLastM"][2])):
    #     print(f'{data["rLastM"][2][i]} - {data["rLastM"][2][i - 1]}')
    #     delta = abs( data["rLastM"][2][i] - data["rLastM"][2][i - 1] )
    #     if delta > maxDelta:
    #         maxDelta = delta
    #         jumpI = i
    # print(i)

    # Visualization
    fig, ax = plt.subplots(figsize=(8,6))

    inferno_cmap = colormaps["inferno"]

    ax.plot(data["rFirstM"][1], data["rFirstM"][2], color=inferno_cmap(0.3), label="First Frame Reference")
    ax.plot(data["rLastM"][1],  data["rLastM"][2],  color=inferno_cmap(0.7), label="Last Frame Reference")

    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD ($\AA$)")
    ax.set_title("RMSD Over Trajectory")
    ax.legend()
    plt.savefig(data["outputName"] + "Rmsd.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    return data

def calcRmsf(data, ref_frame = 0, select = "name CA", verbose = True):
    print("\nCalculating RMSF..")
    ## Calculate the RMSF for the selection
    data["select"] = data["protein"].select_atoms("name CA")
    
    r = rms.RMSF(data["select"]).run(verbose = verbose)
    
    data["rmsfSelect"] = r.results.rmsf
    np.save(f'rmsf{select}.npy', data["rmsfSelect"]) #alpha carbon only rmsf

    plt.plot(data["select"].resids, data["rmsfSelect"])
    plt.xlabel('Residue number')
    plt.ylabel('RMSF ($\AA$)')
    plt.savefig(data["outputName"] + "Rmsf.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    r = rms.RMSF(data["protein"]).run(verbose = verbose)
    data["rmsfProtein"] = r.results.rmsf
    data["protein"].atoms.tempfactors = data["rmsfProtein"]
    np.save('rmsfProtein.npy', data["rmsfProtein"])
    data["protein"].atoms.write('rmsfProtein.pdb') # pdb with tempfactors set for every atom
    
    return data

# BACKUP PCA function if mdtraj fails
def runPca(data, verbose=True):
    print("\nRunning PCA..")

    # Select CA atoms
    data["ca"] = data["u"].select_atoms("name CA")
    n_frames = len(data["u"].trajectory)
    n_atoms = len(data["ca"])

    # Allocate flattened coordinate array: [frames, 3*N_CA]
    coords = np.zeros((n_frames, n_atoms * 3), dtype=np.float32)

    # Fill coordinate matrix
    for i, ts in enumerate(data["u"].trajectory):
        coords[i, :] = data["ca"].positions.reshape(-1)

    # Center the coordinates
    coords_centered = coords - coords.mean(axis=0)

    # scikit-learn PCA
    pca = PCA()

    # Store results
    data["pc"] = pca.fit_transform(coords_centered) # principal components per frame
    data["pcaModel"] = pca # eigenvectors, eigenvalues available via pca.components_ and pca.explained_variance_

    if verbose:
        print(" Finished PCA:")
        print(" Frames:", data["pc"].shape[0])
        print(" Components:", data["pc"].shape[1])

    return data

def applyPCA(data1, data2, verbose=True):
    print("\nApplying PCA from data1 to data2..")
    
    # Select CA atoms from data2
    data2["ca"] = data2["u"].select_atoms("name CA")
    n_frames = len(data2["u"].trajectory)
    n_atoms = len(data2["ca"])

    # Allocate flattened coordinate array: [frames, 3*N_CA]
    coords = np.zeros((n_frames, n_atoms * 3), dtype=np.float32)

    # Fill coordinate matrix
    for i, ts in enumerate(data2["u"].trajectory):
        coords[i, :] = data2["ca"].positions.reshape(-1)

    # Use data1's pre-trained PCA model to transform data2
    # .transform() automatically centers using data1's mean_
    data2["pc"] = data1["pcaModel"].transform(coords)
    data2["pcaModel"] = data1["pcaModel"]  # Reference same model

    if verbose:
        print(" Finished applying PCA:")
        print(" Frames:", data2["pc"].shape[0])
        print(" Components:", data2["pc"].shape[1])

    return data2
    

def calcSfe(data, bins = 200, sigma = 2, temperature = 310):
    print("\nCalculating surface free energy..")
    kB = 1.380649E-23 # J/K
    T = 273.15 + temperature # K
    
    states, xedges, yedges = np.histogram2d(data["pc"][:,0],
                                       data["pc"][:,1],
                                       bins=bins)
    prob = states / np.sum(states) # probability of a state (bin)
    F = -kB*T*np.log(prob, where=(prob>0)) # absolute free energy
    deltaF = F - np.nanmin(F) # relative free energy, sets min to 0

    data["deltaFsmooth"] = gaussian_filter(deltaF, sigma = sigma)
    
    # Visualization
    plt.contourf(data["deltaFsmooth"].T, cmap='inferno')
    plt.xlabel(f"PC1 ({100*data['pcaModel'].explained_variance_ratio_[0]:.2f}% Explained Variance)")
    plt.ylabel(f"PC2 ({100*data['pcaModel'].explained_variance_ratio_[1]:.2f}% Explained Variance)")
    plt.colorbar(label="ΔF (kcal/mol)")
    plt.title("Surface Free Energy Plot")
    plt.savefig(data["outputName"] + "SfePCA.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    return data

def runUmap(data, n_neighbors=70, min_dist=0.8, n_components=3):
    print("\nRunning UMAP..")
    start = time.perf_counter()
    
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric='euclidean',
        random_state=42
    )
    
    data["umap"] = reducer.fit_transform(data["pc"][:, :10]) # using first 10 pcs
    
    end = time.perf_counter()
    print(f"Finished finding contacts in {end - start} seconds")
    
    # Visualization
    plt.figure(figsize=(8,6))
    plt.scatter(data["umap"][:,0], data["umap"][:,1], 
                c = np.arange(len(data["umap"])), 
                cmap="inferno", s=3)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.colorbar(label="Frame")
    plt.title("UMAP of PCA Projection")
    plt.savefig(data["outputName"] + f"umapN{n_neighbors}D{min_dist}.png", 
                dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    return data

def runUmapSfe(data, bins = 200, temperature = 310, sigma = 2):
    print("\nRunning UMAP SFE..")
    kB = 1.380649E-23 # J/K
    T = 273.15 + temperature # K
    
    states, xedges, yedges = np.histogram2d(data["umap"][:,0],
                                       data["umap"][:,1],
                                       bins=bins)
    prob = states / np.sum(states) # probability of a state (bin)
    F = -kB*T*np.log(prob, where=(prob>0)) # absolute free energy
    deltaF = F - np.nanmin(F) # relative free energy, sets min to 0

    data["deltaFsmooth"] = gaussian_filter(deltaF, sigma = sigma)
    
    # Visualization
    
    plt.contourf(data["deltaFsmooth"].T, cmap='inferno')
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.colorbar(label="ΔF (kcal/mol)")
    plt.title("Surface Free Energy UMAP of PCA Projection")
    plt.savefig(data["outputName"] + "SfeUMAP.png", 
                dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    return data

def clusterUmap(data, min_cluster_size=500):
    print("\nClustering UMAP with no noise..")

    # HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=1,
        cluster_selection_epsilon=0.01
    )
    labels = clusterer.fit_predict(data["umap"])

    # Assign noise points (-1) to nearest cluster
    noise_idx = np.where(labels == -1)[0]
    if len(noise_idx) > 0:
        unique_labels = np.unique(labels[labels != -1])
        centroids = np.array([data["umap"][labels == lbl].mean(axis=0) for lbl in unique_labels])

        for i in noise_idx:
            nearest = np.argmin(cdist([data["umap"][i]], centroids))
            labels[i] = unique_labels[nearest]

    data["cluster"] = labels

    # Visualization
    plt.figure(figsize=(8,6))
    unique_clusters = np.unique(labels)
    cmap = plt.get_cmap("tab10")

    for idx, cluster_id in enumerate(unique_clusters):
        cluster_points = data["umap"][labels == cluster_id]
        plt.scatter(
            cluster_points[:,0], 
            cluster_points[:,1], 
            color=cmap(idx % 10), 
            s=4, 
            label=f"Cluster {cluster_id}"
        )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Clustered UMAP")
    plt.legend(markerscale=3, fontsize=10)
    plt.savefig(data["outputName"] + "umapClustered.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()

    return data

### RUN ####

def runAll(args1, args2=None):
    """
    If args2 is None, run single analysis on args1.
    If args2 is provided, compute PCA on data1 (from args1) and apply to data2 (from args2).
    """
    data1 = setup(args1)
    data1 = makeU(data1)
    data1 = calcRmsd(data1)
    data1 = calcRmsf(data1)
    data1 = runPca(data1)
    data1 = calcSfe(data1)
    data1 = runUmap(data1)
    data1 = clusterUmap(data1)
    data1 = runUmapSfe(data1)

    # If a second dataset is provided, apply data1's PCA to it
    if args2 is not None:
        data2 = setup(args2)
        data2 = makeU(data2)
        data2 = applyPCA(data1, data2)  # Apply data1's PCA model to data2
        data2 = calcSfe(data2)  # Use data1's PCA model (via data2["pcaModel"])
        data2 = runUmap(data2)
        data2 = clusterUmap(data2)
        data2 = runUmapSfe(data2)

    return None

if __name__ == "__main__":
    #### ARGPARSE ####
    parser = argparse.ArgumentParser(description = "TBXT Analysis with optional reference PCA")
    parser.add_argument("--pdb", required = True, help = "PDB file (reference for PCA)")
    parser.add_argument("--dcd", required = True, help = "DCD file (reference for PCA)")
    parser.add_argument("--pdb2", default=None, help = "Second PDB file (optional, to apply reference PCA)")
    parser.add_argument("--dcd2", default=None, help = "Second DCD file (optional, to apply reference PCA)")
    parser.add_argument("--outputName", required = True,
                        type = str, 
                        help="Base string name for files. Do not append .png .pdb or .dcd")
    parser.add_argument("--outputName2", default=None,
                        type = str,
                        help="Base string name for second dataset files (optional)")
    
    args = parser.parse_args()
    
    # Change working directory if you need to
    os.chdir("/volumes/rpaul1tb/pitt/ssdTbxt2035/unstitchedReplicas")
    
    # If both pdb2 and dcd2 are provided, set up args2
    if args.pdb2 and args.dcd2:
        args2 = argparse.Namespace(
            pdb=args.pdb2,
            dcd=args.dcd2,
            outputName=args.outputName2 if args.outputName2 else args.outputName + "_data2"
        )
        runAll(args, args2)
    else:
        runAll(args)


# Manual PCA
# def runPca(data, verbose = True):
#     print("\nRunning PCA..")
#     atoms = data["t"].topology.select("name CA")
#     coords = data["t"].xyz[:, atoms, :].reshape(len(data["t"]), -1) 
#     # flatten to 2d matrix for PCA
#     # flattened rows = frames or samples: [x1 y1 z1 x2 y2 z2 ... xn yn zn]
#     # columns = features
    
#     nFeatures = coords.shape[1]
#     coordsCentered = coords - coords.mean(axis=0) # centered matrix for PCA
#     # feature x sample @ sample x feature / nFeatures
#     cov = (coordsCentered.T @ coordsCentered) / nFeatures
    
#     eigvals, eigvecs = np.linalg.eigh(cov)
#     idx = np.argsort(eigvals)[::-1]
#     eigvals = eigvals[idx]
#     eigvecs = eigvecs[:, idx]
    
#     data["pc"] = coordsCentered@eigvecs # [frames x features] @ [features @ features]
#     return data