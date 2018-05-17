# Job array
Application support for user wanting to speed up video processing.

Code modifications to make the code suitable for use in a job array on the CSF.
In short, the code needs to be one function which is then called in a loop (the job array).

# Missing files
- The video files (input).
- Polyparticle video analysis code.

# Instructions
- log in to the CSF
- load the git module file `module load tools/gcc/git/2.8.2`
- clone this repo into your scratch space
	- `cd scratch`
	- `git clone https://github.com/UoMResearchIT/jobarray_markjohnston.git`
- move the toolbox files into the same directory as the cloned repo
- compile the MATLAB code
	- load the MATLAB module file `module load apps/binapps/matlab/R2017a`
	- `mcc -m MVR2.m -R -singleCompThread -a Polyparticle`
- write jobscript for the job array, including the following command
	- `./run_MVR2.sh $MATLAB_HOME "${PWD}" $SGE_TASK_ID`
