#!/bin/bash
# document this script to stdout (assumes redirection from caller)
cat $0

# receive my worker number
export WRK_NB=$1

# create worker-private subdirectory in $SNIC_TMP
export WRK_DIR=$SNIC_TMP/WRK_${WRK_NB}
mkdir $WRK_DIR

# create a variable to address the "job directory"
export JOB_DIR=$SLURM_SUBMIT_DIR/job_${WRK_NB}

# now copy the input data and program from there

cd $JOB_DIR

cp -p input.mat job.m $WRK_DIR

# change to the execution directory

cd $WRK_DIR

# run the program

# Start program
module load matlab/2020b

matlab -singleCompThread -nodesktop -nodisplay -nosplash -r "job; quit" > outfile.txt

# rescue the results back to job directory
# including any text files or images ending with png
cp -p results.mat *.txt *.png ${JOB_DIR} 

# clean up the local disk and remove the worker-private directory

cd $SNIC_TMP

rm -rf WRK_${WRK_NB}
