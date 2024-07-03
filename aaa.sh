#!/bin/bash

python ./source/mainmfcc.py  mfccModule.inp  ./test/1zdd4.pdb > 1.log
rm -rf EE-GMFCC_calculFile
python ./source/mainmfcc.py  mfccModule.inp  ./test/G.pdb    > 2.log
rm -rf EE-GMFCC_calculFile
python ./source/mainmfcc.py  mfccModule.inp  ./test/1le1.pdb > 3.log
rm -rf EE-GMFCC_calculFile
python ./source/mainmfcc.py  mfccModule.inp  ./test/2xl1.pdb > 4.log
rm -rf EE-GMFCC_calculFile
python ./source/mainmfcc.py  mfccModule.inp  ./test/1wn8.pdb > 5.log
