#!/bin/bash
cd /home/hexu/projects/DionJ/phonon
# OpenMp Environment
export OMP_NUM_THREADS=1
# Commands before execution
export PATH=/home/hexu/.local/abinit/8.4.1/bin/:$PATH

mpirun -n 1 anaddb < /home/hexu/projects/DionJ/phonon/run.files > /home/hexu/projects/DionJ/phonon/run.log 2> /home/hexu/projects/DionJ/phonon/run.err
