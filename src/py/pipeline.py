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
from subprocess import Popen, PIPE

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

# time function
def print_now():
    now = datetime.datetime.now()
    return(now.strftime("%Y-%m-%d %H:%M"))

# Custom job submission step. Dirty hack. Need to exit 0 at end of each script
# and use set -e for safety
def submit_job(jobScript, ntasks, cpus_per_task, job_name):
    '''
    Submit the job using salloc hack. When complete return job id and write output to file.
    '''
    # call salloc as subprocess
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
    # mail output
    mail = Popen(['mail', '-s', "[Tom@SLURM] Pipeline step " + jobId + " finished",
                  'tom'], stdin = PIPE)
    mail.communicate(out)
    # check subprocess exit code
    assert proc.returncode == 0, 'Job failed with non-zero exit code'  
    return(jobId)

# Variant of submit_job for JGI jobs where email/password is required
def submit_download_job(jobScript, job_name, jgiLogon, jgiPassword):
    # make sure logon and password were set
    assert jgiLogon, "--email is required"
    assert jgiPassword, "--password is required"

    ntasks = '1'

    # call download script    
    proc = Popen(['salloc', '--ntasks=' + ntasks,'--job-name=' + job_name,
                  jobScript, "-e", jgiLogon, "-p", jgiPassword],
                  stdout = PIPE, stderr = PIPE)
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
    # mail output
    mail = Popen(['mail', '-s', "[Tom@SLURM] Pipeline step " + jobId + " finished",
                  'tom'], stdin = PIPE)
    mail.communicate(out)
    # check completion    
    assert proc.returncode == 0, "Job " + job_name + " failed with non-zero exit code"
    return(jobId)

# touch function for updating ruffus flag files
def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd, **kwargs)

#---------------------------------------------------------------
# download rice genomes
#
@originate(['ruffus/os.genome'], jgiLogon, jgiPassword)

def download_os_genome(outputFiles, jgiLogon, jgiPassword):
    jobScript = 'src/sh/downloadGenomes.sh'
    job_name = 'osGen'
    jobId = submit_download_job(jobScript, job_name, jgiLogon, jgiPassword)

    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# download tomato genome
#
@originate(['ruffus/sl.genome'], jgiLogon, jgiPassword)

def download_sl_genome(outputFiles, jgiLogon, jgiPassword):
    jobScript = 'src/sh/downloadSlGenome.sh'
    job_name = 'slGen'
    jobId = submit_download_job(jobScript, job_name, jgiLogon, jgiPassword)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# download arabidopsis genome
#
@originate(['ruffus/at.genome'], jgiLogon, jgiPassword)

def download_at_genome(outputFiles, jgiLogon, jgiPassword):
    jobScript = 'src/sh/downloadAtGenome.sh'
    job_name = 'atGen'
    jobId = submit_download_job(jobScript, job_name, jgiLogon, jgiPassword)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)
    
#---------------------------------------------------------------
# download tomato reads
#
@originate(['ruffus/sl.reads'])

def download_sl_reads(outputFiles):
    jobScript = 'src/sh/downloadSlReads.sh'
    ntasks = '4'
    cpus_per_task = '1'
    job_name = 'slReads'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# download arabidopsis reads
#
@originate(['ruffus/at.reads'])

def download_at_reads(outputFiles):
    jobScript = 'src/sh/downloadAtReads.sh'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'atReads'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)
    

    
# options for visualising
pipeline_printout()
pipeline_printout_graph("ruffus/flowchart." + slurm_jobid + ".pdf", "pdf")

# run the pipeline (disabled for now)
cmdline.run(options, multithread = 8)
