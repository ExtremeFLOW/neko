#!/bin/sh
export DYLD_FALLBACK_LIBRARY_PATH=$JSON_INSTALL/lib
echo "Test that small cylinder (restart and deformations)"
mpirun -np $1 ../../neko cylinder_part1.case > log1
mpirun -np $1 ../../neko cylinder_part2.case > log2
echo "Initial run"
cat log1
echo "Restarted run"
cat log2
ref1=ref1_${2}.log
ref2=ref2_${2}.log
awk 'c&&c--{if(c !=17 && NF > 2){print}};/100.00%/{c=19}' log1 > l1
awk 'c&&c--{if(c !=17 && NF > 2){print}};/100.00%/{c=19}' log2 > l2
awk 'c&&c--{if(c !=17 && NF > 2){print}};/100.00%/{c=19}' ${ref1} > r1
awk 'c&&c--{if(c !=17 && NF > 2){print}};/100.00%/{c=19}' ${ref2} > r2
diff l1 r1 > res
diff l2 r1 >> res
diff l1 r2 >> res
diff l2 r2 >> res
if  [ -s res ]
then
   echo 'Differences compared to reference logs'
   cat res
   echo "Did not pass sanity check" >> /dev/stderr
   exit 1
else
   echo "Passed sanity check"
fi

