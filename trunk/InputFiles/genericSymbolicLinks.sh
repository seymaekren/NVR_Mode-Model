#!/bin/sh
#genericSymbolicLinks.sh: uses command line arguments
ln -s answerkey.m.$1 answerkey.m
ln -s myinput.m.$1 myinput.m
ln -s NOES.txt.$1 NOES.txt
ln -s order.m.$1 order.m
ln -s TOCSY.m.$1 TOCSY.m
ln -s MySHIFTS.7.model6.$1  SHIFTS/MySHIFTS.7.model6
ln -s MySHIFTX.7.model6.$1  SHIFTX/MySHIFTX.7.model6
ln -s MySHIFTX2.7.model6.$1 SHIFTX2/MySHIFTX2.7.model6
ln -s C-H_vectors.m.$1 C-H_vectors.m
ln -s C-H_medium1.m.$1 C-H_medium1.m
ln -s N-H_vectors.m.$1 N-H_vectors.m
ln -s N-H_medium1.m.$1 N-H_medium1.m
