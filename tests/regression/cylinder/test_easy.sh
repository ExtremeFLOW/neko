#!/bin/sh
export DYLD_FALLBACK_LIBRARY_PATH=$JSON_INSTALL/lib
echo "Easy test for small cylinder (restart and deformations)"
mpirun -np $1 ../../../src/neko cylinder_part1.case >log1
mpirun -np $1 ../../../src/neko cylinder_part2.case >log2
echo "Initial run"
cat log1
echo "Restarted run"
cat log2
ref1=ref1_${2}.log
ref2=ref2_${2}.log

# Check that first residual is the same, allowed to differ on last digit...
grep -E '\d+\s+\| Pressure' log1 | head -1 >l1
grep -E '\d+\s+\| Pressure' ${ref1} | head -1 >r1
diff l1 r1 >res
# Check that we do same number of residuals
grep -E '\d+\s+\| Pressure' log1 | wc -l >l1
grep -E '\d+\s+\| Pressure' ${ref1} | wc -l >r1
diff l1 r1 >>res
# Check that residual is same after restart
grep -E ' 7\s+\| Pressure' log1 | sed 's/^[^|]*|/|/' >l1
grep -E ' 2\s+\| Pressure' log2 | sed 's/^[^|]*|/|/' >l2
diff l1 l2 >res

if [ -s res ]; then
    echo 'Differences compared to reference logs'
    cat res
    echo "Did not pass sanity check" >>/dev/stderr
    exit 1
else
    echo "Passed sanity check"
fi
