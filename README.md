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

