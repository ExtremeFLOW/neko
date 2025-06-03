#!/bin/sh
export DYLD_FALLBACK_LIBRARY_PATH=$JSON_INSTALL/lib
echo "Test that small cylinder (restart and deformations)"
mpirun -np $1 ../../../src/neko cylinder_part1.case >log1
mpirun -np $1 ../../../src/neko cylinder_part2.case >log2
echo "Initial run"
cat log1
echo "Restarted run"
cat log2
ref1=ref1_${2}.log
ref2=ref2_${2}.log

# Check all residuals of last time step
# We look for the 11'th pressure iterations and 3 lines following it.
# Here we trim the first column, which is the time step number, since that is
# not relevant for the comparison.
grep -A3 '11 | Pressure' log1 | sed 's/^[^|]*|/|/' >l1
grep -A3 ' 6 | Pressure' log2 | sed 's/^[^|]*|/|/' >l2
grep -A3 '11 | Pressure' $ref1 | sed 's/^[^|]*|/|/' >r1
grep -A3 ' 6 | Pressure' $ref2 | sed 's/^[^|]*|/|/' >r2

# Compare that they are identical to what was done before
diff l1 r1 >res
diff l2 r1 >>res
diff l1 r2 >>res
diff l2 r2 >>res
if [ -s res ]; then
    echo 'Differences compared to reference logs'
    cat res
    echo "Did not pass sanity check" >>/dev/stderr
    exit 1
else
    echo 'Differences compared to reference logs:'
    cat res
    echo 'None.'
    echo "Passed sanity check"
fi
