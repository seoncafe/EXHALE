#!/bin/bash
# Run LaRT spherical-illumination Ly-alpha RT for hd209.
exec < /dev/null 2>&1
trap "" HUP

EXEC=/nfs/lart4/kiseon/LaRT/combine/LaRT_v2.00/LaRT_calcP.x

HOSTS=lart4,lart3,lart2
host_file=/tmp/host_file_$RANDOM
for host in $(echo $HOSTS | tr "," "\n"); do
   [[ $host = "mocafe" ]] && num=88 || num=72
   echo $host:$num >> $host_file
done

echo "Running $EXEC on $HOSTS"
echo "   machinefile $host_file:"; cat $host_file
mpirun -machinefile $host_file $EXEC hd209_lya.in
