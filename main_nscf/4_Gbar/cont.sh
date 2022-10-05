#!/bin/bash
while :
do
	EXIST=`ps -ef | grep chenhr | grep MATLAB | grep -v "grep" | wc -l`
	if [ $EXIST -ne 0 ]; then
		sleep 600
	else
                nohup matlab < Gbar_gamma_sub.m > out.file2 2>&1 & 
		break
	fi
done

wait
