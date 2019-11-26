# MAINMASTseg
MAIHNMASTseg is a segmentation program for EM maps with symmetry.

Copyright (C) 2019 Genki Terashi, Daisuke Kihara, and Purdue University.

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing)

Contact: Daisuke Kihara (dkihara@purdue.edu)

Cite : Terashi, Genki, and Daisuke Kihara. "MAINMASTseg: Automated map segmentation method for cryo-EM density maps with symmetry" in submission. 


## Pre-required software

UCSF Chimera {EM map preparation & visualization} : https://www.cgl.ucsf.edu/chimera/  
Pymol{for visualiztion} : https://pymol.org/2/  
OpenMP library

## Installation for MacOS

1. Install Xcode  
2. Install Homebrew  

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

3. Install OpenMP  

brew install libomp

4. Complile code  

git clone https://github.com/kiharalab/MAINMASTseg.git<br>
cd MAINMASTseg<br>
rm MAINMASTseg *.o #remove compiled files <br>
make -f MakefileMacOS #For mac<br>

## Tutorial
[Example1](./Example1) contains all input files and result files.  
Please check http://kiharalab.org/mainmast_seg/Tutorials.html

### Input Data
#### Density Map (mrc format)
In [Example1](./Example1), there is a emd-0093.mrc

EMD-0093 has C4 symmetry. The size is (220*1.34,220*1.34,220*1.34). The center is (147.4,147.4,147.4).  
In Chimera, type:   
1. open 6gyn.pdb as #0
2. Type the following command in Chimera command line

sym #0 group c4 center 147.4,147.4,147.4
3. Save all pdbs as symmetry.pdb 4. Extract the rotation matrix from the saved pdb file as:

grep BIOMT[1-3] symmetry.pdb > MTX.txt 

### Rotation Matrix file (PDB format)

### Segmentation

### Visualization




