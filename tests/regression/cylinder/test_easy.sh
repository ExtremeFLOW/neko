#!/bin/sh
export DYLD_FALLBACK_LIBRARY_PATH=$JSON_INSTALL/lib
echo "Easy test for small cylinder (restart and deformations)"
mpirun -np $1 ../../../src/neko cylinder_part1.case > log1
mpirun -np $1 ../../../src/neko cylinder_part2.case > log2
echo "Initial run"
cat log1
echo "Restarted run"
cat log2
ref1=ref1_${2}.log
ref2=ref2_${2}.log
#Check that first residual is the same, allwoed to differ on last digit...
awk '/Iteration/{getline;print(substr($2, 1, length($2)-6))}' log1 | head -1 > l1
awk '/Iteration/{getline;print(substr($2, 1, length($2)-6))}' ${ref1} | head -1 > r1
diff l1 r1 > res
#Check that we do same number of residuals
awk '/Iteration/{getline;print($1)}' log1 > l1
awk '/Iteration/{getline;print($1)}' ${ref1} > r1
diff l1 r1 >> res
#Check that residual is same after restart
awk '/Pressure/ {i+=1; if(i== 8){getline;getline;print($1,$2)}}' log1 > l1
awk '/Pressure/ {i+=1; if(i== 2){getline;getline;print($1,$2)}}' log2 > l2
diff l1 l2 >> res

if  [ -s res ]
then
   echo 'Differences compared to reference logs'
   cat res
   echo "Did not pass sanity check" >> /dev/stderr
   exit 1
else
   echo "Passed sanity check"
fi

