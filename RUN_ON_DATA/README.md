



This is a fully self contained analysis template for PC study using RECOSIM trees

all tests and jobs are pathed relative to repo, but data lists are hardcoded to directory on cluster

job_test.sh and local_test.sh are are for debugging
these test cases are linked to simple file lists, and produce output in the pc_Execution directory

job_launch.sh launches the full recosim list from minbias 2018, the data list is hardcoded in the shell script
this produces output also in pc_execution directory

The work flow should be as follows--
	For each iteration of the analysis, rather than make balooning modifications to histset/hists/enum etc, pull and work on a fresh version of the repo
	for useful or repeated calculations that should be adopted into further templates, add these tools to PCTools.h (on a fresh repo) and they will be added into the master branch for later analyses


For different datatier/tree types I will make different templates and replace recosim.C/h

to build type
make -f recosim.mk

to setup root locally 
. setenv.sh


to launch a job test or full job set
sbatch job_test.sh
sbatch job_launch.sh

to run jobs locally (interactive without submission)
./job_test.sh
./local_test.sh


REMEMBER TO CHANGE EMAIL NOTIFICATION SETTINGS IN JOB LAUNCH SCRIPTS
