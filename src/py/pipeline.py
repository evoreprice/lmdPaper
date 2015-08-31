#!/usr/bin/python3

# slurm options
#SBATCH --ntasks=1
#SBATCH --job-name="pipeline"
#SBATCH --mail-type=ALL
#SBATCH --output=ruffus/pipeline.%j.log

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
    # write stderr & stdout to log file    
    logFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.log'
    f = open(logFile, 'wb')
    f.write(b"Standard output:\n")
    f.write(out)
    f.write(b"\nStandard error:\n")
    f.write(err)
    f.close()
    # mail output
    if proc.returncode != 0:
        subject = "[Tom@SLURM] Pipeline step " + jobId + " FAILED"
    else:
        subject = "[Tom@SLURM] Pipeline step " + jobId + " finished"
    mail = Popen(['mail', '-s', subject, 'tom'], stdin = PIPE)
    mail.communicate(out + err)
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
    # write stderr & stdout to log file    
    logFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.log'
    f = open(logFile, 'wb')
    f.write(b"Standard output:\n")
    f.write(out)
    f.write(b"\nStandard error:\n")
    f.write(err)
    f.close()
    # mail output
    if proc.returncode != 0:
        subject = "[Tom@SLURM] Pipeline step " + jobId + " FAILED"
    else:
        subject = "[Tom@SLURM] Pipeline step " + jobId + " finished"
    mail = Popen(['mail', '-s', subject, 'tom'], stdin = PIPE)
    mail.communicate(out + err)
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
# SETUP TASKS
#---------------------------------------------------------------

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
    ntasks = '2'
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

#---------------------------------------------------------------
# download TFDB
#
@originate(['ruffus/tfdb.data'])

def download_tfdb(outputFiles):
    jobScript = 'src/sh/downloadTfdb.sh'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'tfdb'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# download MADS peptides
#
@originate(['ruffus/madsPeptides.data'])

def get_mads_peptides(outputFiles):
    jobScript = 'src/R/getMADSsequences.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'madsPept'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# define rice reads
#
@originate(['ruffus/os.reads'])

def define_os_reads(outputFiles):
    pathToReads = 'data/reads/os'
    assert os.path.isdir(pathToReads), "Error: reads folder " + pathToReads + " missing"
    readFiles = os.listdir(pathToReads)
    print("[", print_now(), ": Using Oryza sativa reads in folder " + pathToReads + " ]")
    for fileName in readFiles:
        qualName = pathToReads + '/' + fileName
        assert os.path.isfile(qualName), "Error: read file " + qualName + " missing"
        print(qualName)
    touch(outputFiles)

#---------------------------------------------------------------
# define zhang in situ data
#
@originate(['ruffus/os.insitu'])
def define_review_inSituDB(outputFiles):
    pathToFile = 'data/zhangGenes.tab'
    assert os.path.isfile(pathToFile), "Error: in situ data file " + qualName + " missing"
    touch(outputFiles)

#---------------------------------------------------------------
# ANALYSIS TASKS
#---------------------------------------------------------------

#---------------------------------------------------------------
# trim rice reads
#
@transform(define_os_reads, suffix('.reads'), '.trimmedReads')

def trim_os_reads(inputFiles, outputFiles):
    jobScript = 'src/sh/cutadapt.sh'
    ntasks = '7'
    cpus_per_task = '1'
    job_name = 'cutadapt'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)
    
#---------------------------------------------------------------
# generate STAR index for OS
#
@transform(download_os_genome, suffix('.genome'), '.index')

def generate_os_index(inputFiles, outputFiles):
    jobScript = 'src/sh/starGenomeGenerate.sh'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'stargg'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# map rice reads
#
@merge([trim_os_reads, generate_os_index], output = 'ruffus/os.bamfiles')
           
def map_os_reads(inputFiles, outputFiles):
    jobScript = 'src/sh/starMappingTwoStep.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'star2s'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)
 
#---------------------------------------------------------------
# map tomato reads
#
@merge([download_sl_reads, download_sl_genome], output = "ruffus/sl.bamfiles")

def map_sl_reads(inputFiles, outputFiles):
    jobScript = 'src/sh/pipelineSl.sh'
    ntasks = '1'
    cpus_per_task = '4'
    job_name = 'slMap'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# map arabidopsis reads
#
@merge([download_at_reads, download_at_genome], output = "ruffus/at.bamfiles")

def map_at_reads(inputFiles, outputFiles):
    jobScript = 'src/sh/pipelineAt.sh'
    ntasks = '1'
    cpus_per_task = '4'
    job_name = 'atMap'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    

#---------------------------------------------------------------
# run DESeq2
#
@transform(map_os_reads, suffix(".bamfiles"), ".deseq2")

def run_deseq2_os(inputFiles, outputFiles):
    jobScript = 'src/R/runDeseq.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'DESeq'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    
    
#---------------------------------------------------------------
# calculate TPM
#
@merge([run_deseq2_os, map_os_reads, download_os_genome], 'ruffus/os.tpm')

def calculate_tpm(inputFiles, outputFiles):
    jobScript = 'src/R/calculateTpm.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'calcTpm'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    

#---------------------------------------------------------------
# shuffle GTF
#
@merge([download_os_genome, calculate_tpm], 'ruffus/os.shuffle')

def shuffle_gtf(inputFiles, outputFiles):
    jobScript = 'src/sh/shuffle.sh'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'shuffle'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    

#---------------------------------------------------------------
# count reads in shuffled gtf
#
@merge([map_os_reads, shuffle_gtf], 'ruffus/os.shufCounts')

def shuffled_counts(inputFiles, outputFiles):
    jobScript = 'src/sh/htseqShuffle.sh'
    ntasks = '7'
    cpus_per_task = '1'
    job_name = 'shfcount'        
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    

#---------------------------------------------------------------
# calculate intergenic tpms
#
@merge([shuffle_gtf, shuffled_counts, map_os_reads, run_deseq2_os, calculate_tpm], 'ruffus/os.intTpm')
def intergenic_tpm(inputFiles, outputFiles):
    jobScript = 'src/R/shuffledTpm.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'shufTpm'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    

#---------------------------------------------------------------
# calculate list of expressed genes
#
@merge([shuffle_gtf, intergenic_tpm, calculate_tpm], 'ruffus/os.expgen')
def detect_expressed_genes(inputFiles, outputFiles):
    jobScript = 'src/R/calculateCutoffs.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'expGen'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    
    
#---------------------------------------------------------------
# fuzzy c-means clustering
#
@merge([detect_expressed_genes, run_deseq2_os, download_os_genome], 'ruffus/os.mfuzz')
def clustering(inputFiles, outputFiles):
    jobScript = 'src/R/mfuzz.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'mfuzz'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)    

#---------------------------------------------------------------
# hypergeometric tests for enrichment of TFs
#
@merge([clustering, download_tfdb], 'ruffus/os.hypergeom')
def hypergeom(inputFiles, outputFiles):
    jobScript = 'src/R/tfHypergeom.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'tfHyper'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# count tRNA/rRNA in libraries
#
@merge([download_os_genome, map_os_reads], 'ruffus/os.rnaStats')
def rnaStats(inputFiles, outputFiles):
    jobScript = 'src/sh/countRTrna.sh'
    ntasks = '2'
    cpus_per_task = '1'
    job_name = 'rtCount'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# parse library stats
#
@merge([map_os_reads, rnaStats, run_deseq2_os, detect_expressed_genes], 'ruffus/os.libStats')
def libStats(inputFiles, outputFiles):
    jobScript = 'src/R/parseQuantStats.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'libStats'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)
    
#---------------------------------------------------------------
# Geneset enrichment analysis
#
@merge([download_tfdb, run_deseq2_os, detect_expressed_genes], 'ruffus/os.gsea')
def gsea(inputFiles, outputFiles):
    jobScript = 'src/R/gsea.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'gsea'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# Compare published TF data
#
@merge([detect_expressed_genes, define_review_inSituDB], 'ruffus/os.compare')
def compare_inSitus(inputFiles, outputFiles):
    jobScript = 'src/R/compareVsInSitu.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'compare'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    touch(outputFiles)

#---------------------------------------------------------------
# FIGURES AND TABLES
#---------------------------------------------------------------

#---------------------------------------------------------------
# library stats table
#
@merge(libStats, "ruffus/table.st_libStats")
def st_libStats(inputFiles, outputFiles):
    touch(outputFiles)

#---------------------------------------------------------------
# Hypergeometric tests table
#
@merge(hypergeom, "ruffus/table.t_hypergeom")
def t_hypergeom(inputFiles, outputFiles):
    touch(outputFiles)

#---------------------------------------------------------------
# PCA plot
#
@merge([run_deseq2_os, detect_expressed_genes], "ruffus/figure.sf_pca")
def sf_pca(inputFiles, outputFiles):
    touch(outputFiles)

#---------------------------------------------------------------
# Clustering figures
#
@merge(clustering, "ruffus/figure.f_mfuzzClusters")
def f_mfuzzClusters(inputFiles, outputFiles):
    touch(outputFiles)
@merge(clustering, "ruffus/figure.sf_mfuzzCentroids")
def sf_mfuzzCentroids(inputFiles, outputFiles):
    touch(outputFiles)
@merge(clustering, "ruffus/figure.sf_mfuzzPca")
def sf_mfuzzPca(inputFiles, outputFiles):
    touch(outputFiles)

#---------------------------------------------------------------
# GSEA plot
#
@merge(gsea, "ruffus/figure.f_gsea")
def f_gsea(inputFiles, outputFiles):
    touch(outputFiles)

#---------------------------------------------------------------
# in situ comparison plot
#
@merge(compare_inSitus, "ruffus/table.st_reviewInSitu")
def st_reviewInSitu(inputFiles, outputFiles):
    touch(outputFiles)





# options for visualising
pipeline_printout()
pipeline_printout_graph("ruffus/flowchart." + slurm_jobid + ".pdf", "pdf")

# run the pipeline (disabled for now)
cmdline.run(options, multithread = 8)
