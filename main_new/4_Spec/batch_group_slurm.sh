#!/bin/bash
NCLUSTER=1
for i in `seq 1 1 $NCLUSTER`;
do
{
	sbatch --job-name=Spec_$i 444.sh $i
}
done
