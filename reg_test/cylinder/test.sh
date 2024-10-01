#!/bin/sh
export DYLD_FALLBACK_LIBRARY_PATH=$JSON_INSTALL/lib
mpirun -np $1 neko cylinder_part1.case > log1
mpirun -np $1 neko cylinder_part2.case > log2
awk 'c&&c--{if(c !=17){print}};/100.00%/{c=19}' log1 > l1
awk 'c&&c--{if(c !=17){print}};/100.00%/{c=19}' log2 > l2
awk 'c&&c--{if(c !=17){print}};/100.00%/{c=19}' ref2.log > r2
awk 'c&&c--{if(c !=17){print}};/100.00%/{c=19}' ref1.log > r1
echo 'Differences compared Martins laptop'
diff l1 r1
diff l2 r1
diff l1 r2
diff l2 r2
