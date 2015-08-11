#!/usr/bin/python3

import os
import subprocess
import re


def submit_job():
    '''
    Submit the job using salloc hack. Return job id and write output to file.
    '''
    # call salloc with Popen
    proc = subprocess.Popen(['salloc', './bof.sh'], stdout = PIPE, stderr = PIPE)
    # get stdout and stderr    
    out, err = proc.communicate()
    # parse stderr (salloc output) for job id
    jobRegex = re.compile(b'\d+')
    jobIdBytes = jobRegex.search(err).group(0)
    jobId = jobIdBytes.decode("utf-8")
    # write stdout to log file    
    logFile = "log." + jobId + ".out"
    f = open(logFile, 'wb')
    f.write(out)
    f.close()
    
    return(jobId)


submit_job()