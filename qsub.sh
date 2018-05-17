#!/bin/bash

# -- SGE options (whose lines must begin with #$)

#$ -S /bin/bash		# Our jobscript is written for the bash shell
#$ -V			# Inherit environment settings (e.g., from loaded modulefiles)
#$ -cwd			# Run the job in the current directory

#$ -t 1-5		# Tell SGE that this is an array job, 
			# with "tasks" numbered from 1 to 5
# -- the commands to be executed (programs to be run) on a compute node:

./run_MVR2.sh $MATLAB_HOME "$PWD" $SGE_TASK_ID
