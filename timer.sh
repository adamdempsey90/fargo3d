#!/bin/bash

PROCS=(1 2 3 4 5 6 7 8 9 10)
for i in ${PROCS[*]} 
do 
    mpirun -np $i fargo3d_748x1382 -o "nx=128,ny=128,ninterm=1,dt=1.0,ntot=0" in/accretion.par &> f.out
    echo "np="$i":    "`grep time f.out | awk {'print $5'}`
done

