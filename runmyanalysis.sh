#!/bin/sh
# runmyanalysis.sh

date

#CODEDIR=/home/gwwilson/pc
#CODEDIR=/home/j342a201/DPG/graham_repo/pc
CODEDIR=$PWD

INPUT=${1-batch_test}
MACRO=${2-runmacro.py}
NTHREADS=${3-1}

echo "runmyanalysis.sh Input Args"
echo ${INPUT}
echo ${MACRO}
echo ${NTHREADS}



INPUTLIST=${CODEDIR}/Lists/${INPUT}.list

#Execute this from execution directory, so that we can have several output files in parallel

#EXEDIR=/home/gwwilson/pc_ExecutionDirectory/${INPUT}
EXEDIR=$PWD/pc_ExecutionDirectory/${INPUT}

rm -r ${EXEDIR}
# Check if it exists. If not make it.
mkdir ${EXEDIR}

cd ${EXEDIR}
pwd

#input is: num files, numthreads,  yourdata.list
# Note any change to the thread count needs to also be in the job description file ..
python2 ${CODEDIR}/${MACRO} 0 ${NTHREADS} ${INPUTLIST} ${CODEDIR}

#Need to find some way of having several of these in parallel ...
#cp Outfile.root ${CODEDIR}/PC_${INPUT}.root

date

exit
