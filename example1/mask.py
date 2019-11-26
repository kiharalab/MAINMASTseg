
from chimera import runCommand as rc

rc("open ./emd-0093.mrc")
rc("volume #0 level 0.7 step 1")
rc("open ./MAP.mrc")
rc("volume #1 step 1")
rc("sop hideDust #0 size 100 metric volume")
rc("mask #1 #0 extend 3 fullmap true")
rc("volume #2 save MAP_m4A.mrc")
rc("stop now")


