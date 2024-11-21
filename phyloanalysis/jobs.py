#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 15:26:53 2018

@author: nw14
"""

import time
import os
import sys
import subprocess
import argparse
from math import floor
from math import ceil
#def get_lsf():
    
THIS_SCRIPT=os.path.realpath(__file__)
MAX_RETRIES=3

LSF_RUNLIMIT='TERM_RUNLIMIT'
LSF_MEMLIMIT='MEMLIMIT'
LSF_ABNORMAL_EXIT='Exited with exit code'
LSF_SUCCESS='Successfully completed.'
LSF_NO_LOG_FILE="NO_LSF_LOG_FILE"
LSF_UNKNOWN="UNKNOWN"

def callShell(cmd,b_fail_on_error=True):
    print(cmd,flush=True)
    cmd=[ '/bin/bash', '-c', 'set -eo pipefail; ' + cmd ]
    retcode=subprocess.call(cmd)
    return retcode;

        
def submit_as_job_array(jobs,job_prefix,job_dir,queue="normal",mem_gb=1,b_wait=True,max_chunks=200,ncores=1,retry_number=0,noreason_retry_number=0):
    N=len(jobs)
    
    runAsLsfJobArray="python {}".format(THIS_SCRIPT);  #/nfs/team78pc/nw14/mpn_code/runAsLsfJobArray.pl"
    mb=ceil(1000*mem_gb)
    mem_flag="-R\"select[mem>{mb}] rusage[mem={mb}] span[hosts=1]\" -M{mb} -n{ncores}".format(mb=mb,ncores=ncores);
    mem_flag2="-R\"select[mem>{mb}] rusage[mem={mb}]\" -M{mb}".format(mb=100);
    tt=time.time()
    #numeric_suffix=int(round((tt-floor(tt)) * 1e6))
    numeric_suffix=int(round(tt* 1e6))
    jobid="{}{}".format(job_prefix,numeric_suffix)
    lsf_log="{}/{}%I".format(job_dir,job_prefix)
    os.makedirs(job_dir,exist_ok=True)
    job_file="{}/{}.cmd".format(job_dir,job_prefix)
    done_file="{}.DONE".format(job_file)
    if os.path.exists(done_file):
        print("Removing ",done_file,file=sys.stderr,flush=True)
        os.remove(done_file)
    jf=open(job_file,mode='w')
    i=1
    done_files=[]
    for job in jobs:
        this_done_file="{}.{}.DONE".format(job_file,i)
        job_cmd="{} ; touch {}".format(job,this_done_file)
        print(job_cmd,file=jf)
        print(job_cmd,file=sys.stderr,flush=True)
        i=i+1
        done_files.append(this_done_file)
    jf.close()
    chunk=""
    if len(jobs)>max_chunks:
        chunk="%{}".format(max_chunks)
        
    cmd1="bsub -J\"{}[1-{}]{}\" -o {}.log -q {} {} {} {}".format(jobid,N,chunk,lsf_log,queue,mem_flag,runAsLsfJobArray,job_file)
    
    print(cmd1,flush=True)
    v1 = callShell(cmd1)
    failed_jobs=[]
    if b_wait:
        cmd2="bsub -w \"ended({})\" -q {} {} -o {}.done.log touch {}".format(jobid,queue,mem_flag2,job_file,done_file)
        print(cmd2,flush=True)
        v2 = callShell(cmd2)
        while not os.path.isfile(done_file):
            time.sleep(1)
        ##Give 1 minute extra for all log files to be written too...
        print("Jobs finished checking for success",flush=True)
        time.sleep(60)
        i=0
        memf=1
        newqueue=queue
        orig_lsf_log=[]
        for dfile in done_files :
            if not os.path.isfile(dfile):
                logfile="{}.log".format(lsf_log.replace("%I",str(i+1)))
                status=get_cause(logfile)
            #if status[0] != LSF_SUCCESS or (not os.path.isfile(dfile)):
                print("Could not find {}".format(dfile),flush=True)
                if status[0]==LSF_MEMLIMIT:
                    memf=2
                if status[0]==LSF_RUNLIMIT:
                    newqueue=next_queue(queue)
                print("CODE:{}:DETAILS:{}".format(*status),flush=True)
                failed_jobs.append(jobs[i])
                orig_lsf_log.append(logfile)
            i=i+1
        if len(failed_jobs)==0:
            return 0
        if len(failed_jobs)==len(jobs) and memf==1 and newqueue==queue:
            raise ValueError("All LSF jobs faileds")
        else:
            ##Only fail if remaining jobs are failing for non-obvious reasons and above MAX_RETRIES
            if memf==1 and newqueue==queue:
                noreason_retry_number=noreason_retry_number+1
            if noreason_retry_number < MAX_RETRIES:
                retry_number=retry_number+1
                print("Rerunning failed jobs...",file=sys.stderr,flush=True)
                ##Append to original LSF log
                for i in range(len(failed_jobs)):
                    with open(orig_lsf_log[i], "a") as orig_log:
                        orig_log.write("Retrying : see new log file {}\n".format("{}/{}v{}{}.log".format(job_dir,job_prefix,retry_number,i+1)))
                return submit_as_job_array(failed_jobs,"{}v{}".format(job_prefix,retry_number),
                                           job_dir,queue=newqueue,mem_gb=memf*mem_gb,
                                           b_wait=b_wait,max_chunks=max_chunks,
                                           retry_number=retry_number,
                                           noreason_retry_number=noreason_retry_number,ncores=ncores)
            else:
                raise ValueError("Too many retries")
            

def next_queue(queue):
    nextqueue={"normal":"long","long":"basement","basement":"basement"}
    if queue not in nextqueue:
        raise ValueError("Unsupported queue {} requested ".format(queue))
    else:
        return nextqueue[queue]

def run_job_array_element(job_file):
    job_index=int(os.environ['LSB_JOBINDEX'])-1
    print("JOB_INDEX=",job_index,flush=True);
    print("JOB_FILE=",job_file,flush=True)
    cmds=[entry for entry in open(job_file)]
    #print(cmds)
    cmd=cmds[job_index]
    print("CMD:",cmd,flush=True)
    status=callShell(cmd);
    if status != 0:
        raise ValueError("Failed executing jobindex={} : {}".format(job_index+1,cmd))
    
 
def get_cause(lsf_log):
    keywords={LSF_RUNLIMIT,LSF_MEMLIMIT,LSF_ABNORMAL_EXIT,LSF_SUCCESS}
    if not os.path.isfile(lsf_log):
        return (LSF_NO_LOG_FILE,"UNKNOWN")
    for x in open(lsf_log):
        for keyword in keywords:
            if keyword in x:
                return (keyword,x)
    return (LSF_UNKNOWN,"UNKNOWN")

        
if __name__ == "__main__": 
    if 'LSB_JOBINDEX' in os.environ:
        run_job_array_element(sys.argv[1])
    else:
         ##submit_as_job_array(jobs,job_prefix,job_dir,queue="normal",mem_gb=1,b_wait=True,max_chunks=200,retry_number=0,noreason_retry_number=0)
        parser = argparse.ArgumentParser(description="Kicks of LSF array job")
        parser.add_argument("--jobfile",help="File containing commands",required=True)
        parser.add_argument("--jobname",type=str,required=True)
        parser.add_argument("--jobdir",type=str,required=True)
        parser.add_argument("--queue",type=str,default="normal")
        parser.add_argument("--mem_gb",type=int,default=1)
        parser.add_argument("--ncores",type=int,default=1)
        parser.add_argument("--max_chunks",type=int,default=200)
        args=parser.parse_args()
        ##print("running test..")
        jobs=[line.rstrip() for line in open(args.jobfile)]
        submit_as_job_array(jobs,args.jobname,args.jobdir,queue=args.queue,mem_gb=args.mem_gb,max_chunks=args.max_chunks,ncores=args.ncores)

    
