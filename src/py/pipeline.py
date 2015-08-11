#!/usr/bin/python3

# slurm options
#SBATCH --ntasks=1
#SBATCH --job-name="pipeline"
#SBATCH --mail-type=ALL
#SBATCH --output=ruffus/pipeline.%j.out

# imports
import os
import re
import datetime
from ruffus import *
import ruffus.cmdline as cmdline
from subprocess import Popen, PIPE, communicate

# command-line options
parser = cmdline.get_argparse(description = 'Run LMD analysis pipeline.')
parser.add_argument('--email', '-e',
                        help ='Logon email address for JGI',
                        type = str,
                        dest = 'jgiLogon')
parser.add_argument('--password', '-p',
                        help ='JGI password',
                        type = str,
                        dest = 'jgiPassword')
options = parser.parse_args()
jgiLogon = options.jgiLogon
jgiPassword = options.jgiPassword

# parse SLURM job-id
if 'SLURM_JOBID' in os.environ:
    slurm_jobid = os.environ.get('SLURM_JOBID')
else:
    slurm_jobid = 'local-' + str(datetime.date.today()) 

# Custom job submission step. Dirty hack. Need to exit 0 at end of each script
# and use set -e for safety
def submit_job(jobScript, ntasks, cpus_per_task, job_name):

    '''
    Submit the job using salloc hack. When complete return job id and write output to file.
    '''
    # call salloc with Popen
    proc = Popen(['salloc', '--ntasks=' + ntasks, '--cpus-per-task=' + cpus_per_task,
    '--job-name=' + job_name, jobScript], stdout = PIPE, stderr = PIPE)
    # get stdout and stderr    
    out, err = proc.communicate()
    # parse stderr (salloc output) for job id
    jobRegex = re.compile(b'\d+')
    jobIdBytes = jobRegex.search(err).group(0)
    jobId = jobIdBytes.decode("utf-8")
    # write stdout to log file    
    logFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.log'
    f = open(logFile, 'wb')
    f.write(out)
    f.close()
    
    return(jobId)

# folder for ruffus files
@mkdir("ruffus")

#---------------------------------------------------------------
# download rice genomes
#
@originate(['ruffus/os.genome'], jgiLogon, jgiPassword)

def download_osativa_genome(outputFiles, jgiLogon, jgiPassword):
    # make sure logon and password were set
    try:
        jgiLogon
    except NameError:
        raise Exception("--email is required")
    try:
        jgiPassword
    except NameError:
        raise Exception("--password is required")
    jobScript = 'src/sh/downloadGenomes.sh'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'osGen'
    logFile = 'ruffus/' + job_name + '.' + slurm_jobid + '.out'
    f = open(logFile, 'w')
    call(['salloc', '--ntasks=' + ntasks, '--cpus-per-task=' + cpus_per_task,
    '--job-name=' + job_name, jobScript, "-e", jgiLogon, "-p", jgiPassword,
    "-l", logFile], stdout = f, stderr = STDOUT)
    f.close()

#---------------------------------------------------------------
# download tomato genome
#
@originate(['ruffus/sl.genome'], jgiLogon, jgiPassword)

#---------------------------------------------------------------
# download arabidopsis genome
#
@originate(['ruffus/at.genome'], jgiLogon, jgiPassword)



#---------------------------------------------------------------
# subsample some libraries for testing
#
@originate(['ruffus/osativa.reads'])

def subsample_os_reads(outputFiles):
    jobScript = 'src/sh/subsample.sh'
    ntasks = '3'
    cpus_per_task = '1'
    job_name = 'subsamp'
    logFile = '/tmp/' + job_name + '.ruffus.out'
    submit_job(jobScript, ntasks, cpus_per_task, job_name, logFile)

#---------------------------------------------------------------
# create STAR index
#
@transform(download_osativa_genome, suffix('.genomes'), ".index")

def create_star_index(input_files, output_files):
    jobScript = 'src/sh/starGenomeGenerate.sh'
    ntasks = '1'
    cpus_per_task = '3'
    job_name = 'stargg'
    logFile = '/tmp/' + job_name + '.ruffus.out'
    submit_job(jobScript, ntasks, cpus_per_task, job_name, logFile)

#---------------------------------------------------------------
# run some mapping jobs
#
@merge([subsample_os_reads, create_star_index], ".bamfiles")

def map_os_reads(input_files, output_files):
    jobScript = 'src/sh/starQuickmap.sh'
    ntasks = '1'
    cpus_per_task = '6'
    job_name = 'star'
    logFile = '/tmp/' + job_name + '.ruffus.out'
    submit_job(jobScript, ntasks, cpus_per_task, job_name, logFile)

# options for visualising
pipeline_printout()
pipeline_printout_graph("ruffus/flowchart." + slurm_jobid + ".pdf", "pdf")

# run the pipeline (disabled for now)
# cmdline.run(options, multithread = 8)
