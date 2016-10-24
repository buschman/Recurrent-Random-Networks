#!/usr/bin/env bash

#name job, output to slurm file, use partition all, run for 300 minutes and use 4GB of ram
#SBATCH -J 'SimWM_RecRnd_MultiTrials'
#SBATCH -o SimWM_RecRnd_MultiTrials_%A_%a.out
#SBATCH -p all
#SBATCH -t 600
#SBATCH --mem-per-cpu=4096  
#SBATCH --array=0-287

# run the latest matlab
module load matlab/R2015b

cd '/jukebox/buschman/Projects/Models of Working Memory/Recurrent-Random Networks'

#run a paranoid version of matlab that crashes gracefully and records why just incase.
# in addition, xvfb-run also creates a virtual desktop so that you can plot easily
xvfb-run -d matlab -nosplash -nodisplay -nodesktop -r "try; spockSimWM_RecRnd_MultiTrials_Chunk($SLURM_ARRAY_TASK_ID, 20, 1); catch me; fprintf('%s / %s\n',me.identifier,me.message); for i = 1:length(me.stack), me.stack(i), end; end; exit"