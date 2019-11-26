#HOW TO GENERATE ROTATION MATRIX
#generate rotation matrix by "sym" command in UCSF Chimera GUI command:
1. open 6gyn.pdb
2. Command (the center coordinat is 220*1.34*0.50=147.4)
sym #1 group c4 center 147.4,147.4,147.4
3. Save all pdbs as symmetry.pdb
4. Extract the rotation matrix from the saved pdb file.
grep BIOMT[1-3] symmetry.pdb > MTX.txt


#or use phenix.map_symmetry:
phenix.map_symmetry MAP_m4A.mrc symmetry=C4
../conv_ncs.pl symmetry_from_map.ncs_spec > MTX_fromPhenix.txt



#Example the density threshold=0.7
../MainmastSeg -i MAP_m4A.mrc -Y MTX.txt -c 8 -t 0.7 -M > test.cif 

or 

../MainmastSeg -i MAP_m4A.mrc -Y MTX_phenix.txt -c 8 -t 0.7 -M > test.cif

#visuaize MST
../bondtreeCIF.pl test.cif > a.txt
pymol -u a.txt

