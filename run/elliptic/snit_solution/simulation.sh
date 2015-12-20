#!/bin/bash

Kn=$(python -c "
import numpy, sys
numpy.savetxt(sys.stdout, numpy.logspace(-4,-1,19)[2:-2], fmt='%.6e')
")

problem=$(basename $(pwd))
problem="../_$problem"

rm -rf $problem
mkdir -p $problem
for kn in $Kn; do
    case=$(printf "%.4f" $kn)
    mkdir $problem/$case
    cp -r * $problem/$case/
    (
        echo "Simulate for Kn=$kn"
        cd $problem/$case
        ./Allrun asym $kn
    )
done
