#!/usr/bin/env python

__version__ = "htsqc v0.1.0"


__description__ = """
{version}
=================

htsqc pipeline verifying the content of HTS samples

Fabian Kilpert - July 22, 2015
email: kilpert@ie-freiburg.mpg.de

This software is distributed WITHOUT ANY WARRANTY!
--------------------------------------------------
    """.format(version=__version__)

import argparse
from collections import OrderedDict
import datetime
import gzip
import os
import os.path
from Queue import Queue
import random
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import textwrap
from threading import Thread
import time


#### PATHS ####################################################################
temp_dir = tempfile.gettempdir()
script_path = os.path.dirname(os.path.realpath(__file__)) + "/"
python_path=sys.executable


#### Development settings ######################################################

if socket.gethostname() == "pc305":

    ## Reading start options (if available)
    start_options_file = os.path.join(script_path, "htsqc", "start_options.txt")
    if os.path.isfile(start_options_file):
        start_options = []
        with open(start_options_file, "r") as f:
            for line in f.readlines():
                line = line.strip()
                if not line.startswith("#"):
                    if line:
                        start_options.extend(line.split())
        if start_options:
            print "Using start options:", " ".join(start_options)
            sys.argv = [sys.argv[0]] + start_options

    if "--featureCounts" in sys.argv:
        sys.argv2 = sys.argv
        sys.argv2[sys.argv.index( '--featureCounts')+1] = "'" + sys.argv[sys.argv2.index( '--featureCounts')+1] + "'"   # just to output all trim galore option as one string
        print " ".join(sys.argv2)   # output command line
    else:
        print " ".join(sys.argv)   # output command line



## Path defaults to all needed scripts and programs + versions ################
fastqc_path = "/package/FastQC-0.11.3/"; fastqc_ver = "FastQC-0.11.3"
# trim_galore_path = "/package/trim_galore_v0.4.0/"; trim_galore_ver = "TrimGalore-v0.4.0"
# cutadapt_activate = "source /package/cutadapt-1.8.1/bin/activate &&"; cutadapt_ver = "Cutadapt-1.8.1"
# rseqc_path = "/package/RSeQC-2.6.1/bin/"; rseqc_ver = "RSeQC-2.6.1"
# rseqc_activate = "source /package/RSeQC-2.6.1/bin/activate &&"
# bowtie2_path = "/package/bowtie2-2.2.3/"; bowtie2_ver = "Bowtie2-2.2.3"
# bowtie2_export = "export PATH={}:$PATH &&".format(bowtie2_path)
picardtools_path = "/package/picard-tools-1.121/"; picardtools_ver = "Picard-tools-1.1.21"
# tophat2_path = "/package/tophat-2.0.13.Linux_x86_64/"; tophat2_ver = "TopHat-2.0.13"
# feature_counts_path = "/package/subread-1.4.6-p2/bin/"; feature_counts_ver = "featureCounts (subread-1.4.6-p2)"
# htseq_count_path = "/package/HTSeq-0.6.1/bin/"; htseq_count_ver = "HTSeq-0.6.1"
R_path = "/package/R-3.2.0/bin/"; R_ver = "R-3.2.0"
samtools_path = "/package/samtools-1.2/"; samtools_ver = "Samtools-1.2"
samtools_export = "export PATH={}:$PATH &&".format(samtools_path)
# ucsctools_dir_path = "/package/UCSCtools/"
#hisat_path = "/package/hisat-0.1.5-beta/bin/hisat"
hisat_path = "/package/hisat-0.1.6-beta/bin/"; hisat_ver = "HISAT-0.1.6-beta"
R_libraries_export = "export R_LIBS_USER=/data/manke/repository/scripts/rna-seq-qc/R/x86_64-unknown-linux-gnu-library/3.2 &&"
# deseq2_ver = "DESeq2-1.8.1"
pandoc_path = "/usr/bin/"; pandoc_ver = "pandoc 1.12.3.1"

valid_genomes = ["mm10", "hg38", "dm6", "Zv9"]

## Different configurations for other physical machines
if socket.gethostname() == "pc305":
    fastqc_path = "/home/kilpert/Software/bin/"; fastqc_ver = "FastQC"
    # trim_galore_path = "/home/kilpert/Software/trim_galore/trim_galore_v0.4.0/"; trim_galore_ver = "TrimGalore-v0.4.0"
    # cutadapt_activate = ""; cutadapt_ver = "Cutadapt"
    # rseqc_path = "/home/kilpert/Software/RSeQC/RSeQC-2.6.1/scripts/"; rseqc_ver = "RSeQC-2.6.1"
    # rseqc_activate = "source activate pc305_RSeQC_2.6.1 &&"
    # bowtie2_path = ""; bowtie2_ver = "Bowtie2"
    # bowtie2_export = ""
    picardtools_path = "/home/kilpert/Software/picard-tools/picard-tools-1.115/"; picardtools_ver = "Picard-tools-1.115"
    # tophat2_path = ""; tophat2_ver = "TopHat-2"
    # feature_counts_path = ""; feature_counts_ver = "featureCounts"
    # htseq_count_path = ""; htseq_count_ver = "HTSeq"
    R_path = "/usr/bin/"; R_ver = "R-3.2.0"
    samtools_path = "/home/kilpert/Software/samtools/samtools-1.2/"; samtools_ver = "Samtools"
    samtools_export = ""
    # ucsctools_dir_path = ""
    hisat_path = "/home/kilpert/Software/hisat/hisat-0.1.5-beta/"; hisat_ver = "HISAT"
    R_libraries_export = "export R_LIBS_USER=/data/manke/repository/scripts/rna-seq-qc/R/x86_64-pc-linux-gnu-library/3.2 &&"
    # deseq2_ver = "DESeq2-1.8.1"
    pandoc_path = "/usr/bin/"; pandoc_ver = "pandoc 1.12.4.2"


#### DEFAULT VARIABLES #################################################################################################
is_error = False
default_threads = 3           # Threads per process
default_parallel = 1          # Parallel files
samtools_mem = 1
samtools_threads = 1

if socket.gethostname().startswith("deep"):
    temp_dir = "/data/extended"
    default_threads = 4
    default_parallel = 6
    samtools_mem = 2
    samtools_threads = 2

elif socket.gethostname().startswith("solserv"):
    default_threads = 8
    default_parallel = 1
    samtools_mem = 2
    samtools_threads = 2



#### COMMAND LINE PARSING ##############################################################################################

def parse_args():
    """
    Parse arguments from the command line.
    """
    parser = argparse.ArgumentParser(
    prog='htsqc.py',
    formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(__description__))

    ### Optional arguments
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument("-i", "--input-dir", dest="main_indir", required=True, help="Input dir containing (FASTQ)")
    parser.add_argument("-o", "--output-dir", dest="main_outdir", required=True, help="Output directory")
    parser.add_argument("--overwrite", dest="overwrite", action = "store_true", default=False, help="Overwrite results in existing folders!")
    parser.add_argument("-p", "--parallel", dest="parallel", metavar="INT", help="Number of files in parallel processing (default: {})".format(default_parallel), type=int, default=default_parallel)
    parser.add_argument("-t", "--threads", dest="threads", metavar="INT", help="Maximum number of threads for a single process (default: {})".format(default_threads), type=int, default=default_threads)
    parser.add_argument("--seed", dest="seed", metavar="INT", help="Random sampling seed", type=int, default=None)
    parser.add_argument("--fastq-downsample", dest="fastq_downsample", metavar="INT", help="Subsample first n fastq sequences (for testing only!)", type=int, default=1000000)
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose output")
    parser.add_argument("--non-random", dest="non_random", action="store_true", default=False, help="Deactivate random downsampling. Downsampling from head instead.")
    parser.add_argument("-g", "--genomes", dest="genomes", metavar="STR", help="Reference genomes used for mapping", type=str, default="all")

    ### Positional arguments (required!)
    args = parser.parse_args()

    ## Tools path for report
    args.script_path = script_path

    ### Add useful paths
    args.cwd = os.getcwd()
    args.now = datetime.datetime.now()

    args.indir = args.main_indir
    args.outdir = args.main_outdir

    ### Sanity checking
    ## Input dir
    ##print "indir:", args.indir
    args.indir = os.path.expanduser(args.indir)
    args.indir = os.path.realpath(args.indir)
    if (not os.path.isdir(args.indir)):
        print "Error: The specified input directory does not exist or is not a directory:\n", args.indir
        exit(1)

    ## Output dir
    args.outdir = os.path.join(args.cwd, os.path.expanduser(args.outdir))

    args.error = 0

    if args.genomes == "all":
        print "Mapping to all available reference genomes!"
        args.genomes = valid_genomes
    else:
        print "Mapping to user specified reference genomes only!"
        args.genomes = re.split(';|,', re.sub('"', '', args.genomes) )

        for g in args.genomes:
            if g not in valid_genomes:
                print "Valid mapping references are:", ",".join(valid_genomes)
                print "Error! User specified genome is not available:", g
                exit(1)

    return args


#### GENERIC FUNCTIONS #################################################################################################

class Qjob():
    def __init__(self, cmds=None, cwd=None, logfile=None, backcopy=True, shell=False, keep_temp=False):
        self.cmds = cmds
        self.cwd = cwd
        self.logfile = logfile
        self.backcopy = backcopy
        self.shell = shell
        self.keep_temp = keep_temp
    def run(self):
        pass


def parse_config(file_path):
    """
    Parse a configuration text file, e.g. *.cfg
    """
    options = {}
    for line in file(file_path):
        line = line.rstrip("\r\n")
        if not line.startswith("#"):
            if "=" in line:
                key, value = line.split('=', 1)
                value = value.strip()
                value = value.strip('"')
                value = value.strip("'")
                options[key] = value
    return options


def get_genomes(all):
    """
    Get reference data paths from config file
    """
    available_genomes = {}
    for genome in all:
        ref_cfg_file_path = os.path.join(script_path, "htsqc/{}.cfg".format(genome))
        if not os.path.isfile(ref_cfg_file_path):
            print "Error! Configuration file NOT found for {}: {}".format(genome, ref_cfg_file_path)
            exit(1)
        configs = parse_config(ref_cfg_file_path)
        available_genomes[genome] = configs
    return available_genomes


def check_for_paired_infiles(args, indir, ext, verbose=False):
    """
    Read input files from given dir that have the user specified extension.
    """
    infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(ext)])

    paired_infiles = OrderedDict()
    for infile in infiles:
        fname = os.path.basename(infile).replace(ext, "")
        m = re.match("^(\w+?)(_R1|_R2)*(_\d*)*$", fname)
        if m:
            bname = m.group(1)
            if bname not in paired_infiles:
                paired_infiles[bname] = [infile]
            else:
                paired_infiles[bname].append(infile)
        else:
            print "No match!"
    infiles = paired_infiles.values()

    for files in infiles:
        if len(files) > 2:
            print "Error! More than 2 files are sharing the same base name:"
            print " ".join(files)
            exit(1)

    ## Paired
    if len(infiles[0]) == 2:
        args.paired = True
        print "Paired end FASTQ files ({} pairs)".format(len(infiles))
        if verbose:
            for pair in infiles:
                print "  {}  {}".format(os.path.basename(pair[0]), os.path.basename(pair[1]))
    elif len(infiles[0]) == 1:
    ## Single
        args.paired = False
        infiles = map(lambda x: x[0], infiles)
        print "Single end FASTQ files ({})".format(len(infiles))
        if verbose:
            for infile in infiles:
                print "  {}".format(os.path.basename(infile))

    return infiles


def run_subprocess(cmd, cwd, td, shell=False, logfile=None, backcopy=True, verbose=False, keep_temp=False):
    """
    Run the subprocess.
    """
    try:
        if verbose:
            print "Temp dir:", td
            print cmd

        if shell == True:
            return subprocess.check_call("cd {} && ".format(td)+cmd, shell=True, executable='/bin/bash')

        else:
            if type(cmd) is str:
                cmd = cmd.split()

            # print "CMD:", cmd
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=td)
            try:
                if verbose:
                    print "PID:", p.pid
                (stdoutdata, stderrdata) = p.communicate()
                if logfile:
                    with open(logfile, "a+") as log:
                        if stderrdata:
                            log.write("{}\n".format(stderrdata))
                        if stdoutdata:
                            log.write("{}\n".format(stdoutdata))
                if stderrdata:
                    print stderrdata
                if stdoutdata:
                    print stdoutdata
                if p.returncode != 0:
                    print "Error! Command failed with return code:", p.returncode
                    print " ".join(cmd) + "\n"
            finally:
                # print "poll:", p.poll()
                if p.poll() == None:
                    print "Error! Trying to terminate PID {}".format(p.pid)
                    p.terminate()
                    time.sleep(20)
                    if p.poll() == None:
                        print "Error! Trying to kill PID {}!".format(p.pid)
                        p.kill()
                    exit(1)
            return p.wait()
    finally:
        if backcopy:
            for f in os.listdir(td):
                if verbose:
                    print "Backcopy:", os.path.join(td, f)
                if os.path.isfile(os.path.join(cwd,f)) or os.path.isdir(os.path.join(cwd,f)):     # delete existing files with the same name in cwd
                    os.remove(os.path.join(cwd,f))
                #shutil.rmtree(os.path.join(cwd,f), ignore_errors=True)
                shutil.move(os.path.join(td,f), cwd)


def queue_worker(q, verbose=False, rest=0.2):
    """
    Worker executing jobs (command lines) sent to the queue.
    """
    global is_error

    while True:
        job = q.get()
        if job.logfile:
            logfile = job.logfile
        else:
            logfile = None
        td = tempfile.mkdtemp(prefix="htsqc.", dir=temp_dir)
        try:
            for cmd in job.cmds:
                if job.logfile:
                    with open(job.logfile, "a+") as log:
                        log.write("{}\n\n".format(cmd))
                return_code = run_subprocess(cmd, job.cwd, td, shell=job.shell, logfile=logfile, backcopy=job.backcopy, verbose=verbose, keep_temp=job.keep_temp)
                if return_code:
                    is_error = True
                time.sleep(rest)
        finally:
            if os.path.isdir(td):
                if not job.keep_temp:
                    shutil.rmtree(td)
        q.task_done()


#### TOOLS #############################################################################################################

#### FASTQ DOWNSAMPLING ################################################################################################

def run_fastq_downsampling(args, q, indir, analysis_name="FASTQ_downsampling"):
    """
    Reduce the number of sequences in FASTQ file to n first sequences
    or random downsampling (with seed).
    """
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name.replace(" ", "_"))
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    print "Outdir:", outdir
    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])

        logfile = os.path.join(cwd, "LOG")
        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        if args.non_random:
            print "Downsampling from head of the files."
            for infile in infiles:
                jobs = ["zcat {} | head -n{} | gzip > {}".format(infile, 4*int(args.fastq_downsample), os.path.join(cwd, os.path.basename(infile)) ),]
                ## Note that there is a not misleading "gzip: stdout: Broken pipe" message. The output is fine though!!!
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True) )
        else:
            print "Random downsampling"
            if args.seed:
                for infile in infiles:
                    print "Using seed:", args.seed
                    jobs = ["{} {}htsqc/downsample_fastq.py -v -s {} -n {} {} {}".format(python_path, script_path, args.seed, args.fastq_downsample, infile, os.path.join(cwd, os.path.basename(infile)) ),]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True) )
            else:
                infiles = check_for_paired_infiles(args, indir, ".fastq.gz")
                if args.paired:
                    for pair in infiles:
                        bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))
                        jobs = ["bash {}htsqc/downsample_reservoir/downsample_reservoir_pe.sh {} {} {} {} {}".format(script_path,
                                    args.fastq_downsample,
                                    pair[0],
                                    pair[1],
                                    os.path.join(cwd,bname+"_R1.fastq.gz"),
                                    os.path.join(cwd,bname+"_R2.fastq.gz")),]
                        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True))
                else:
                    for infile in infiles:
                        bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                        jobs = ["bash {}htsqc/downsample_reservoir/downsample_reservoir_se.sh {} {} {}".format(script_path,
                                    args.fastq_downsample,
                                    infile,
                                    os.path.join(cwd,bname+".fastq.gz")),]
                        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True))
        q.join()
        if is_error:
            exit(is_error)
        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### FASTQC ############################################################################################################

def run_fastqc(args, q, indir, analysis_name="FastQC"):
    """
    Run FastQC on FASTQ files.
    """
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])

        logfile = os.path.join(cwd, "LOG")
        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        jobs = ["{}fastqc --version".format(fastqc_path)]
        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True))
        q.join()

        # FastQC
        for infile in infiles:
            jobs = ["{}fastqc --extract -o {} {}".format(fastqc_path, cwd, infile)]
            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True))
        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)



def run_hisat(args, q, indir, available_genomes):
    """
    Run HISAT mapping.
    """
    analysis_name = "HISAT"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        ## set library-type
        if args.paired:
            library_type = "--rna-strandness RF"
        else:
            library_type = ""

        hisat_opts = ""

        for genome in args.genomes:
            print "Start mapping to {}...".format(genome)

            if not os.path.isdir( os.path.join(cwd, genome) ):
                os.mkdir( os.path.join(cwd, genome) )

            ## PE
            if args.paired:
                for pair in infiles:
                    bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))

                    if not os.path.isdir( os.path.join(cwd, genome, bname) ):
                        os.mkdir( os.path.join(cwd, genome, bname) )

                    cmdl = "{}hisat {} -p {} -x {} {} -1 {} -2 {} --novel-splicesite-outfile {} --un-conc-gz {} --al-conc-gz {} --met-file {} 2> {} | {}samtools view -Sb - | {}samtools sort -@ {} -m {}G - {}"\
                                .format(hisat_path, hisat_opts, args.threads, available_genomes[genome]['hisat_index'], library_type, pair[0], pair[1],
                                        os.path.join(cwd, genome, bname+"/"+"splice_sites.txt"),
                                        os.path.join(cwd, genome, bname+"/"+"un-conc.fastq.gz"),        # --un-conc
                                        os.path.join(cwd, genome, bname+"/"+"al-conc.fastq.gz"),        # --al-conc
                                        os.path.join(cwd, genome, bname+"/"+"metrics.txt"),
                                        os.path.join(cwd, genome, bname+"/"+"align_summary.txt"),
                                        samtools_path,
                                        samtools_path, samtools_threads, samtools_mem,
                                        os.path.join(cwd, genome, bname+"/"+"accepted_hits"),
                                        )

                    jobs = [cmdl,
                            ## add command line to bam header
                            "cat <({}samtools view -H {}) <(echo '@PG\tCL:{}') > {}"\
                                .format(samtools_path, os.path.join(cwd, genome, bname+"/"+"accepted_hits.bam"), cmdl, os.path.join(cwd, genome, bname+"/"+"header.sam")),
                            "{}samtools reheader {} {} | {}samtools view -bS - > {}"\
                                .format(samtools_path, os.path.join(cwd, genome, bname+"/"+"header.sam"), os.path.join(cwd, genome, bname+"/"+"accepted_hits.bam"), samtools_path, os.path.join(cwd, genome, bname+"/"+"accepted_hits.reheader.bam") ),
                            "mv {} {}".format( os.path.join(cwd, genome, bname+"/"+"accepted_hits.reheader.bam"), os.path.join(cwd, genome, bname+"/"+"accepted_hits.bam") ),
                            "rm {}".format(os.path.join(cwd, genome, bname+"/"+"header.sam"), ),
                            ]

                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            ## SE
            else:
                for infile in infiles:
                    bname = re.sub(".fastq.gz$", "", os.path.basename(infile))

                    if not os.path.isdir( os.path.join(cwd, genome, bname) ):
                        os.mkdir( os.path.join(cwd, genome, bname) )

                    cmdl = "{}hisat {} -p {} -x {} {} -U {} --novel-splicesite-outfile {} --un-gz {} --al-gz {} --met-file {} 2> {} | {}samtools view -Sb - | {}samtools sort -@ {} -m {}G - {}"\
                                .format(hisat_path, hisat_opts, args.threads, available_genomes[genome]['hisat_index'], library_type, infile,
                                        os.path.join(cwd, genome, bname+"/"+"splice_sites.txt"),
                                        os.path.join(cwd, genome, bname+"/"+"un.fastq.gz"),         # --un
                                        os.path.join(cwd, genome, bname+"/"+"al.fastq.gz"),         # --al
                                        os.path.join(cwd, genome, bname+"/"+"metrics.txt"),
                                        os.path.join(cwd, genome, bname+"/"+"align_summary.txt"),
                                        samtools_path,
                                        samtools_path, samtools_threads, samtools_mem,
                                        os.path.join(cwd, genome, bname+"/"+"accepted_hits"),
                                        )

                    jobs = [cmdl,
                            # add command line to bam header
                            "cat <({}samtools view -H {}) <(echo '@PG\tCL:{}') > {}"\
                                .format(samtools_path, os.path.join(cwd, genome, bname+"/"+"accepted_hits.bam"), cmdl, os.path.join(cwd, genome, bname+"/"+"header.sam")),
                            "{}samtools reheader {} {} | {}samtools view -bS - > {}"\
                                .format(samtools_path, os.path.join(cwd, genome, bname+"/"+"header.sam"), os.path.join(cwd, genome, bname+"/"+"accepted_hits.bam"), samtools_path, os.path.join(cwd, genome, bname+"/"+"accepted_hits2.bam") ),
                            "mv {} {}".format( os.path.join(cwd, genome, bname+"/"+"accepted_hits2.bam"), os.path.join(cwd, genome, bname+"/"+"accepted_hits.bam") ),
                            "rm {}".format(os.path.join(cwd, genome, bname+"/"+"header.sam")),
                            ]

                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            print
            q.join()
            if is_error:
                exit(is_error)

            ## Generate links in main HISAT output folder and index files
            for infile in infiles:
                if args.paired:
                    bname = re.sub("_R*[1|2].fastq.gz$", "", os.path.basename(infile[0]))
                else:
                    bname = re.sub(".fastq.gz$", "", os.path.basename(infile))
                tophat_file = os.path.join(cwd, genome, bname, "accepted_hits.bam") #"{}/accepted_hits.bam".format(bname)
                os.symlink(tophat_file, "{}___{}.bam".format(genome, bname) )
                subprocess.call("{}samtools index {}___{}.bam".format(samtools_path, genome, bname), shell=True)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)



def run_project_report(args, q):
    """
    Compile a report on the Project.
    """
    analysis_name = "project_report"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "CWD:", cwd

        if os.path.isdir(os.path.join(cwd,"figure")):
            shutil.rmtree(os.path.join(cwd,"figure"))  ## knit exits with error if this folder is present!!!

        jobs = ["[ -f {} ] || ( {} cat {}htsqc/Report_table.R | {}R --vanilla --quiet --args {} {} {} {})".format(os.path.join(cwd,"Report.tsv"), R_libraries_export, script_path, R_path, args.main_indir, args.main_outdir, args.fastq_downsample, args.paired),
                # """( echo 'library(knitr); tsvfile="{tsv}"; knit("{rmd}",output="{md}",envir=parent.frame())' | {R} --vanilla --quiet ) \
                # && {pandoc} {md} --to html --output {html} --css {css} --template {tmpl} --table-of-contents""".format(
                #                         rmd=os.path.join(script_path,"htsqc","htsqc_project.RMD"),
                #                         md=os.path.join(cwd,"Report.md"),
                #                         R=os.path.join(R_path, "R"),
                #                         tsv=os.path.join(cwd,"Report.tsv"),
                #                         figdir=os.path.join(cwd,"figure"),
                #                         pandoc=os.path.join(pandoc_path,"pandoc"),
                #                         html=os.path.join(cwd,"Report.html"),
                #                         css=os.path.join(script_path,"htsqc","TOC_forRmd.css"),
                #                         tmpl=os.path.join(script_path,"htsqc","default.html"),
                # ),

                # {pandoc} {md} --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures --output {html} --smart --email-obfuscation none --self-contained --standalone --section-divs --table-of-contents --toc-depth 3 --template {tmpl} --css {css} --variable 'theme:cerulean'".format(pandoc=os.path.join(pandoc_path,"pandoc"),
                #                         md=os.path.join(cwd,"Report.md"),
                #                         html=os.path.join(cwd,"Report.html"),
                #                         tmpl=os.path.join(script_path,"htsqc","default.html"),
                #                         css=os.path.join(script_path,"htsqc","TOC_forRmd.css"),
                #                         figdir=os.path.join(cwd,"figure")
                # ),

                "{rlib} cat {rmd2html} | {R} --vanilla --quiet --args {rmd} {infile} {html}".format( html=os.path.join(cwd, "Report.hmtl"),
                                                                                                     rlib=R_libraries_export,
                                                                                                     rmd2html=os.path.join(script_path, "htsqc", "RMD2html.R"),
                                                                                                     R=os.path.join(R_path, "R"),
                                                                                                     rmd=os.path.join(script_path, "htsqc", "htsqc_project.RMD"),
                                                                                                     infile=os.path.join(cwd, "Report.tsv"),
                                                                                                    ),
                "rm -R {} {}".format(os.path.join(cwd, "figure"), os.path.join(cwd, "htsqc_project.md"))
                ]
        print "\n".join(jobs)
        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True))
        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)



#### MAIN PROGRAM ######################################################################################################

def main():
    """
    Main program.
    """
    args = parse_args()
    ##print "Args:", args

    print "\n{} htsqc start".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'))
    if args.verbose:
        print "Temp dir:", temp_dir
        # print "Genome:", args.genome
        print "Input dir:", args.indir
        print "Host:", socket.gethostname()
        print "Threads per process:", args.threads
        print "Files in parallel:", args.parallel
        print "Seed (random):", args.seed
        print "Mapping to:", ",".join(args.genomes)
        print
        try:
            print "PATH:", os.environ["PATH"]
        except:
            print ""

    ### Output dir
    ##print "Outdir:", args.outdir
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    os.chdir(args.outdir)

    args.analysis_counter = 0   # Counter for analyzes conducted
    indir = args.indir

    check_for_paired_infiles(args, indir, ".fastq.gz", verbose=args.verbose)      # sets the args.paired to True if sequences are paired end

    q = Queue()
    for i in range(args.parallel):
        worker = Thread(target=queue_worker, args=(q, args.verbose))
        worker.setDaemon(True)
        worker.start()

    available_genomes = get_genomes(args.genomes)

    ## FASTQ downsampling
    t1 = datetime.datetime.now()
    if args.fastq_downsample:
        indir = run_fastq_downsampling(args, q, args.indir)
    ### print "Output folder:", indir
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    # ## Run FastQC
    # t1 = datetime.datetime.now()
    # run_fastqc(args, q, indir)
    # t2 = datetime.datetime.now()
    # print "Duration:", t2-t1

    ## RUN HISAT
    t1 = datetime.datetime.now()
    bam_dir = run_hisat(args, q, indir, available_genomes)
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    ## Run project_report
    t1 = datetime.datetime.now()
    run_project_report(args, q)
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    return args.outdir


if __name__ == "__main__":
    start = datetime.datetime.now()
    #print "os.environ['PATH']:    ", os.environ['PATH']
    #print "os.system('echo $PATH): ", os.system("echo $PATH")
    #print "Args:", sys.argv
    outdir = main()
    print "\n{} htsqc finished (runtime: {})".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), datetime.datetime.now() - start)
    print "Output stored in: {}\n".format(outdir)
