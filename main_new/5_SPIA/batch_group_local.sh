#!/bin/bash
NCLUSTER=1
for i in `seq 1 1 $NCLUSTER`; 
do
{
        matlab -nodesktop -nosplash -nodisplay -r "NCL=$i;SPIA;quit" > log_$i 2>&1 &
}
done
wait
