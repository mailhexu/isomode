#!/bin/bash
cd /Users/hexu/projects/isomode/tests/YAlO3/phbst
# OpenMp Environment
export OMP_NUM_THREADS=1
mpirun -n 1 anaddb < /Users/hexu/projects/isomode/tests/YAlO3/phbst/run.files > /Users/hexu/projects/isomode/tests/YAlO3/phbst/run.log 2> /Users/hexu/projects/isomode/tests/YAlO3/phbst/run.err
