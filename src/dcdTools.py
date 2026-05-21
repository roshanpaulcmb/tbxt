import argparse

import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDFile
import MDAnalysis.transformations as transformations
from MDAnalysis.analysis import align

import mdtraj as mdt
import os

#### ARGPARSE ####

parser = argparse.ArgumentParser(description = "Tools to edit DCD")
parser.add_argument("--pdb", required = True, help = "PDB file")
parser.add_argument("--dcd", required = True, help = "DCD file")
parser.add_argument("--outputName", required = True,
                    type = str, default = "system",
                    help="Base string name for ALL files. Do not append .pdb or .dcd")
parser.add_argument("--dcdStitch", help = "Addable DCD file for stitching")

args = parser.parse_args()

#### FUNCTIONS ####

def setup(pdb = args.pdb, dcd = args.dcd):
    print("\nSetting up..")
    data = { }
    data["pdb"] = pdb
    data["dcd"] = dcd
    
    data["u"] = mda.Universe(data["pdb"], data["dcd"])
    return data

# https://userguide.mdanalysis.org/stable/reading_and_writing.html
def stitch(data, dcdStitch = args.dcdStitch, 
           outputName = args.outputName + "Stitched",
           save = True):
    print("\nStitching .dcd files")
    data["u"] = mda.Universe(data["pdb"], data["dcd"], dcdStitch)
    # Note for self. I merge the trajectories, then add align to ref = 0
    # Do I need to align the second dcd to the last frame of the first dcd?
    
    data["protein"] = data["u"].select_atoms("protein")

    # Write the stitched trajectory to a new DCD file
    if save:
        data["protein"].write(f'{outputName}.pdb')
        data["protein"].write(f'{outputName}.dcd', frames = "all")

    return data

def alignU(data, select = "backbone", verbose = True):
    print("\nAligning mdanalysis universe..")    
    
    aligner = align.AlignTraj(data["u"], 
                              data["u"],
                              select = select,
                              filename = args.outputName + 'Aligned.dcd').run(verbose = verbose)
    
    # Write over original dcd, keep aligned.dcd
    data["dcd"] = args.outputName + 'Aligned.dcd'
    data["u"] = mda.Universe(data["pdb"], data["dcd"])
    
    return data

# https://userguide.mdanalysis.org/stable/reading_and_writing.htmlv
def stripWater(data, 
               outputName = args.outputName + "Stripped",
               verbose = True,
               save = True):
    print("\nStripping water from .pdb and .dcd files..")
    data["protein"] = data["u"].select_atoms("protein")
    
    if save:
        data["protein"].write(f'{outputName}.pdb')
        data["protein"].write(f'{outputName}.dcd', frames = "all")
    
    return data

def downsample(data,
               start = 0, end = -1, step = 10,
               outputName = args.outputName + "Downsampled",
               save = True):
    print("\nWriting downsized .dcd file..")
    data["protein"] = data["u"].select_atoms("protein")

    # Write the downsampled trajectory to a new DCD file
    if save:
        data["protein"].write(f'{outputName}.dcd', 
                              frames = data["u"].trajectory[start:end:step])
        data["dcd"] = f'{outputName}.dcd'
        data["u"] = mda.Universe(data["pdb"], data["dcd"])
    else:
        data["dcd"] = data["u"].trajectory[start:end:step]
        data["u"] = mda.Universe(data["pdb"], data["dcd"])
    return data

def runAll():
    # use whatever functions you need to here
    data = setup()
    if args.dcdStitch is not None:
        data = stitch(data)
    data = downsample(data)
    data = alignU(data)
    data = stripWater(data)
    
    
    return None


if __name__ == "__main__":
    runAll()
