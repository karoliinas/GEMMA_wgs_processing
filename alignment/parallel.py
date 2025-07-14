#!/usr/bin/python3

"""

A tool for running parallel jobs on locally. Reads a list of items
from standard input and runs the specified command on each item. Variable $x
contains the item name in each command invocation. The command must be given
inside single tick marks. The output of each job is stored by default
in a log file under '~/.jobs/'. Variable ${x} can be used in given commands,
all different (white-space separated) values for ${x} are read from stdin.

Usage:
  parallel2 [options] [<command>]
  parallel2 [options] [-a=<FILE>]... [-C=<FILE>]... [-n=N] [<command>]

Examples:
  echo *.bam | parallel2 'samtools sort $x ${x/.bam/} && samtools index $x'
  echo *_1.fq.gz | parallel2 -n8  \\
    'tophat2 tophat-indexes/hg19 $x ${x/_1.fq/_2.fq} -o ${x/_1.fq.gz/}'

Options:
  -h --help             Display this help.
  -n --workers=N        Number of workers to run in parallel.
  -c --cpus=N           How many CPUs to allocate for each job [default: 1].
                        Does not affect locally running jobs.
  -m --memory=I         Show memory usage during the run of job with I index
                        or log directory
  -t --time=N           Time limit in hours for the full analysis [default: 72].
  -j --job-name=NAME    Job name shown by SLURM and --show [default: job].
  -J --max-jobs=N       Maximum jobs to process before quitting. By default
                        maximum is unlimited [default: 0].
  -P --partition=NAME   SLURM partition on which to run jobs [default: local].
  -f --file             <command> is a file to read commands from (one per line)
  -l --logdir=DIR       Set directory for log files [default:~/.jobs/]
  -L --linked=DIR       Create the job log dir as a symbolic link that can
                        be shared by multiple VMs (DIR needs to be on nfs
                        filesystem that supports exclusive opening of files).
                        When linking to a running job, if -n is not specified
                        the same number of workers is used as the creator of
                        the linked job folder is using.
  -R --reset            Delete all previously existing log files
  -M --maintenance      Delete all logs older than 2 months
  -i --increase=I       Increase the number of workers by one for job with I
                        index or log directory. It is safe use this flag
                        multiple times.
  -s --show=N           Show paths of previous N job log directories in the
                        default log directory (~/.jobs/).
  -S --show-jobs        Show paths of previous 10 job log directories in the
                        default log directory (~/.jobs/).
  -v --view=I           Show tail of finished on ongoing tasks log files for  
                        job with index or log dir I (see --show).
  -p --path=I           Output the path of job with index I.                        
  -V --view-job         Show tail of finished and ongoing tasks log files for
                        job with index 1 (see --show).
  -k --shutdown=I       Stops processing after current tasks finish for
                        job with I index or log directory.
  -r --resume=I         Cancels issued --shutdown command for
                        job with I index or log directory.
  -W --wait-all         Wait for all current parallel runs to finish
                        before starting this job.
  -w --wait=I           Wait for job with index or log dir I to finish.
  -a --wait-file=FILE   If FILE exists, wait for it to be removed before
                        starting the job. If FILE does not exists it is
                        created instantly, and removed once the job has
                        finished (status done or interrupted) can be used
                        multiple times. The file acts as a semafore for
                        processes, only the process that is able to write
                        this file will be allowed to run.
  -C --file-comp=FILE   Wait for file to be completed (file size is
                        non-zero and does not change for an hour) before
                        starting job.
  -d --disarrange       Randomize the order that jobs are processed in.
  -I --id               Print out machine id and quit.
  -T --runtimes=I      Display runtimes of finished job

  
"""

from __future__ import print_function
import docopt
import subprocess, sys, os, datetime, time, re
from pypette import shell, open_exclusive, natural_sorted
import atexit
import uuid
import random
import shutil
from glob import glob

# Set version number to worker thread call when updating
VERSION = "2_95"

sbatch_template = '''#!/bin/bash -l
#SBATCH -p %s
#SBATCH -J %s
#SBATCH -n 1
#SBATCH -c %d
#SBATCH --mem-per-cpu=%d
#SBATCH -t %d
#SBATCH -o %s/worker.err
#SBATCH -e %s/worker.err
#SBATCH --open-mode=append
parallel worker %s
'''

UNEXPECTED_EXIT = True
JOB_DIR = ""
SYMLINK_TARGET_DIR = "" # Remote job dir
WORKER_THREAD = False
WORKER_FILE = None
RANDOMIZE_JOBS = False
READ_WAIT_FILES=[]
WRITE_WAIT_FILES=[]
WRITE_WAIT_HANDLES=[]
WATCHED_FILES=[]
RUNNING_DOCKERS=[]
LINKED_JOB=False

def info(message):
    print("INFO: " + message, file=sys.stderr)
    
def error(message):
    print("ERROR: " + message, file=sys.stderr)
    sys.exit(-1)

def warning(message):
    print("WARNING: " + message, file=sys.stderr)

def sanitize_path(path):
    path = path.replace('../', '')
    if path.startswith('./'): path = path[2:]
    return path.replace('/', '_')

def parallel(commands, job_name, max_workers, cpus, partition,
             time_limit, log_dir, wait_to_run=None, linked_dir=None, max_jobs=0):

    global JOB_DIR, SYMLINK_TARGET_DIR, UNEXPECTED_EXIT, READ_WAIT_FILES, WRITE_WAIT_FILES, \
           WRITE_WAIT_HANDLES, LINKED_JOB, RUNNING_DOCKERS, VERSION

    n_exports = 0

    if len( commands) == 0:
        error('No commands specified.')
        sys.exit(-2)

    # Lines (commands) ending with ";"" are merged with the previous line
    # and each line will be run in parallel with others
    edited_commands = []
    merge_with_previous = False


    # Commands
    for c in range( len( commands)):
        elem = str( commands[ c]).replace('\n', ' ').strip()
        #print("ELEM:",elem)
        #sys.exit( 0)
        if len( elem) == 0: continue
        if elem.startswith('#'): continue

        dockers = re.findall(r"docker exec ([^\s;\\]+)", elem)
        RUNNING_DOCKERS += dockers

        multiline = False
        if elem.endswith("\\"):
            elem = elem[:-1]
            multiline = True

        if merge_with_previous: edited_commands[-1] += str(" " + elem)
        else: edited_commands.append( elem)

        merge_with_previous = False
        if elem.endswith(";"): merge_with_previous = True
        elif elem.endswith("&"): merge_with_previous = True # either && or &
        elif elem.endswith("|"): merge_with_previous = True 
        elif multiline: merge_with_previous = True


    commands = edited_commands[:]
    if len( RUNNING_DOCKERS):
        info( "Found %i docker exec commands. Issuing 'docker kill container' if process is interrupted." % len(RUNNING_DOCKERS))

    if sys.stdin.isatty():
        # If the user did not provide any input, just run the command once.
        # The command must not contain $x.
        for command in commands:
            if '$x' in command or '${x' in command:
                error('Command contains $x but no targets provided.')
        targets = []
    else:
        # Parse whitespace-delimited target items from standard input.
        targets = []
        for line in sys.stdin:
            targets += line.split() #The str.split() method without an argument splits on whitespace
        targets = [t.replace('\n', '') for t in targets]

        # Check is targets are really required in any command
        targets_required = False
        for command in commands:
            if '$x' in command or '${x' in command:
                targets_required = True
                break

        if not targets and targets_required: error('Command requires targets but none provided.')

    if len(set(targets)) < len(targets):
        error('Target list contains multiple instances of the following targets:\n' + '\n'.join(s for s in set(targets)
            if targets.count(s) > 1))

    n_jobs = (len(commands)-n_exports) * (1 if len(targets) == 0 else len(targets))

    if n_jobs <= 0:
        error('No jobs to run found.')
        sys.exit( -1)

    if max_workers > n_jobs: max_workers = n_jobs
    if max_jobs > 0 and max_workers > max_jobs: max_workers = max_jobs


    if log_dir == None or len( log_dir) < 1 or log_dir.find(".jobs/") >= 0:
        log_dir = os.path.expanduser('%s%s_%s' % (log_dir, "job", datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')))

    original_thread = True

    # Create or join job at linked dir
    if linked_dir != None:
        if not linked_dir.endswith("/"): linked_dir += "/"

        if os.path.isfile( linked_dir + "tasks"):
            # Joining linked job
            info("Joining running linked job at '%s'. (%s)" % (linked_dir, count_progress( linked_dir)))
            original_thread = False
            LINKED_JOB = True
            load_basedir( linked_dir) # Switch to original working dir
            update_max_jobs_count( linked_dir, max_jobs)
        else:
            # Creating linked job (shared folder)
            try:
                os.makedirs( linked_dir)
                info("Created directory for linked jobs at '%s'." % linked_dir)
                save_basedir( linked_dir)
                save_owner_info( linked_dir)
            except Exception as ex:
                error( "Could not create linked job directory at '%s'. (%s)" % (linked_dir, str(ex)))        

    local_dir = log_dir

    while True:

        try:
            if linked_dir == None:
                # Non-linked job
                os.makedirs( log_dir)
            else:
                # Linked job
                # Log dir is local dir, linked dir is network dir
                os.symlink( linked_dir, log_dir) # os.symlink(src, dst)
                #info("Linked dir:'%s' log dir: '%s'" % ( linked_dir, log_dir))                
                # Use created link as path for job
                SYMLINK_TARGET_DIR = linked_dir
                # workers from original thread are started directly to 
                # SYMLINK_TARGET_DIR so they can know they are from the 
                # original thread
                if original_thread: log_dir = linked_dir

        except FileExistsError:
            # Add numbering to folder name, if it already exists
            m = re.search('_(\d+)$', log_dir)
            if m:
                try_num = int( m.group( 1)) + 1
                if try_num > 99 :
                    error( "Could not create log directory '%s'." % log_dir)
                    UNEXPECTED_EXIT = False
                    sys.exit( -2)
                log_dir = re.sub( '_\d+$', "_%i" % try_num, log_dir)
            else:
                log_dir = log_dir + "_2"
            continue # Try creating dir again
        except Exception as ex:
                sys.stderr.write( "ERROR: Could not create %slog directory '%s'.\n" % (("LINKED " if LINKED_JOB else ""), log_dir))
                sys.stderr.write( "       An exception of type {0} occured.\n       Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
                UNEXPECTED_EXIT = False 
                sys.exit( -2) 
        # Folder created succesfully
        break

    JOB_DIR = log_dir    

    if partition != 'local':
        info('Distributing %d %s named "%s" on %s partition '
            '(with %d %s and X GB of memory per job).' % (
            n_jobs, 'jobs' if n_jobs != 1 else 'job',
            job_name, partition, cpus, 'CPUs' if cpus != 1 else 'CPU'))
    else:
        info('Starting %d %s named "%s" on local machine. (log at:\'%s\')' % (
            n_jobs, 'jobs' if n_jobs != 1 else 'job', job_name, log_dir))
        if job_name != None and len( job_name) > 0 and original_thread:
            try:
                with open( '%s/name' % log_dir, 'w') as f:
                    f.write( "%s\n" % job_name)
            except:
                sys.stderr.write("Warning: Could not write name file for job '%s'.\n" % log_dir)

    # Write "tasks" file
    if original_thread:        
        try:
            tasks_filepath = '%s/tasks' % log_dir            
            with open( tasks_filepath, 'w') as f:
                for cmd in commands:
                    if not cmd or len( cmd) == 0 or cmd.startswith("#"): continue
                    #if cmd.startswith("export"): 
                    #    f.write('EXP:%s\n' % cmd)
                    #else: 
                    f.write('CMD:%s\n' % cmd)
                    for target in targets: f.write('TRG:%s\n' % target)

            # Write recommended number of workers so that
            # VMs linking to job do not need specify it
            workers_filepath = '%s/n_workers' % log_dir
            with open( workers_filepath, 'w') as f:
                f.write('%i\n' % max_workers)         

            update_max_jobs_count( log_dir, max_jobs)   

        except Exception as ex:
                error( "Could not create task file in job directory at '%s'. (%s)" % (log_dir, str(ex)))


    if partition == 'local':

        # Jobs to wait for
        n_jobs_running = 0

        if wait_to_run != None:

            if wait_to_run == "ALL": 
                special_dir = local_dir
                exclude_dir = True # exlude only self
            else:
                special_dir = wait_to_run
                exclude_dir = False # Include chosen dir only

            # Jobs
            #info("Checking if other parallel tasks are running... (--wait)"

            n_jobs_running, running_job_names = get_n_jobs_running(special_dir=special_dir, exclude=exclude_dir)
            if n_jobs_running > 0:
                info("Waiting for %i currently running job%s to finish. (%s)" % (n_jobs_running, ("s" if n_jobs_running > 1 else ""), (",".join(running_job_names))))

            #info( "special:" + special_dir)
            #info( "exclude:" + str( exclude_dir))
            #info( "running:" + str(running_job_names))

        n_jobs_to_wait = n_jobs_running

        # Files
        try:
            waited_stop = 0
            waited_go = 0
            n_watched = len( WATCHED_FILES)
            watched_file_sizes = n_watched * [0]
            while waited_go < 60:
                new_file_sizes = [ (int(os.stat(f).st_size) if os.path.isfile(f) else 0) for f in WATCHED_FILES ]
                changed_or_zero_size = [ (new_file_sizes[ i] == 0 or new_file_sizes[ i] != watched_file_sizes[ i]) for i in range( n_watched) ]                
                if any( changed_or_zero_size):
                    waited_go = 0; 
                    if waited_stop >= 60*24*7: #Week
                      error("Max file completion wait time reached.")
                    time.sleep( 60)
                    waited_stop += 1;
                else:
                    waited_go += 1

                watched_file_sizes = new_file_sizes

        except Exception as ex:
            error( "Exception occurred while waiting for files to complete. (%s)" % str( ex))

        waiting_filepath = '%s/waiting' % log_dir
        remaing_files = process_wait_files()
        n_files = len( remaing_files)
        if n_files > 0: 
            info("Waiting for %i file%s to be removed before starting job." % (n_files, ("" if n_files == 1 else "s")))
        wait_files_used = n_files > 0

        waited = 0
        if n_jobs_running != 0 or n_files != 0:
            try: open( waiting_filepath, 'a').close()
            except: pass

        # Waiting loop
        while n_jobs_running != 0 or n_files != 0:
            time.sleep( 60) # wait 1 minute
            waited += 1
            # Jobs
            if n_jobs_to_wait > 0:                
                n_jobs_running, running_job_names = get_n_jobs_running(special_dir=special_dir, exclude=exclude_dir)
                #info( "special:" + special_dir)
                #info( "exclude:" + str(exclude))
                #info( "running:" + str(running_job_names))
            # Wait files
            if wait_files_used:
                remaing_files = process_wait_files()
                n_files = len( remaing_files)

            if n_jobs_to_wait > 0 and n_jobs_running != 0 and waited % 60 == 0:
                info("Waiting for %i current job(s) to finish: %s  (wait time %.1fh)." % (n_jobs_running, (",".join(running_job_names)), waited/60.0))
            if wait_files_used and n_files != 0 and waited % 60 == 0:
                info("Waiting for %i file(s) to be removed: %s  (wait time %.1fh)." % (n_files, (",".join(remaing_files)), waited/60.0))                    
            if waited > 60*24*7: # week
                error("Max wait time reached. Quitting.")
                try: os.remove( waiting_filepath)
                except: pass
                sys.exit( 0)

        if waited > 0:
            if n_jobs_to_wait > 0: info("All previous jobs have finished.")
            if wait_files_used: info("All wait-files have been removed.")
            info(  "Starting jobs...")
        elif n_jobs_to_wait > 0:
            info("No parallel jobs are currently running (--wait). Starting new job.")
        elif wait_files_used:
            info("No wait-files to wait for found. Starting new job.")

        try: os.remove( waiting_filepath)
        except: pass

        try: 
            parallel_executable = os.path.realpath(__file__)
        except:
            parallel_executable = "/data/scripts/parallel%s.py" % VERSION

        if not os.path.isfile( parallel_executable):
            parallel_executable = "parallel2"

        worker_cmd = [ parallel_executable, 'worker', log_dir]
        workers = [subprocess.Popen(worker_cmd) for w in range( max_workers)]
        all_done = False
        moar_errors = 0
        done_filepath = '%s/done' % log_dir
        # Machine specific files for linked folders
        moar_filepath = '%s/moar_%s' % (log_dir, get_machine_id())
        memory_filepath = '%s/memory_usage_%s' % (log_dir, get_machine_id())

        moar_retries = 2
        n_new_workers = 0
        n_started = 0
        worker_wait = 0
        max_mem_usage = 0.0
        log_mem_usage( memory_filepath, len( workers))

        while not all_done:

            # Log memory usage every 10 mins
            if worker_wait % 60 == 0:
                log_mem_usage( memory_filepath, len( workers))

            try:
                # Wait for all workers to complete
                # Check for requests to start more workers
                # every 10 seconds
                for w in workers: w.wait( 10)
                all_done = True
                
            except subprocess.TimeoutExpired:                
                pass    
            except KeyboardInterrupt:         
                return False
         
            worker_wait += 1

            if LINKED_JOB and is_shutting_down( JOB_DIR): 
                error( "Job owner has interrupted processing.")
                return False

            # Check for existence of special file for requesting more workers in job dir
            if not all_done and moar_errors < moar_retries and os.path.isfile( moar_filepath):
                try:
                    n_new_workers = 0

                    with open( moar_filepath, "r") as mfile:
                        for line in mfile:
                            line = line.strip()
                            if line.startswith("#"): continue
                            elif len( line): n_new_workers = int( line)

                    # Remove file before starting new workers
                    # to avoid starting infinite number of threads
                    # in case removal of moar file fails
                    os.remove( moar_filepath)

                    # Start moar workers
                    for n in range( n_new_workers):
                        workers.append( subprocess.Popen(worker_cmd))
                        n_started += 1

                    sys.stderr.write("INFO: Started %i new worker%s.\n" % (n_new_workers, ("s" if n_new_workers > 1 else "")))
                    moar_errors = 0

                except Exception as ex:
                    sys.stderr.write( "ERROR: Could not process file '%s'.\n" % moar_filepath)
                    template = "ERROR: An exception of type {0} occured. Arguments: {1!r}\n"                    
                    err_message = template.format( type(ex).__name__, ex.args)
                    sys.stderr.write( err_message)
                    if n_started > 0: sys.stderr.write( "INFO: New workers started: %i.\n" % n_started)
                    moar_errors += 1
                    if moar_errors >= moar_retries: sys.stderr.write( "ERROR: No more retries.\n")
                    else: sys.stderr.write( "ERROR: Retrying in 10 seconds...\n")


    else:
        # Run the job steps on a SLURM cluster using sbatch.
        # Required memory is given in GB per job step. Convert to MB per CPU.
        mem_per_cpu = round(float(memory) / cpus * 1000)
        sbatch_script = sbatch_template % (partition, job_name, cpus,
            mem_per_cpu, 60 * time_limit, log_dir, log_dir, log_dir)
        workers = [subprocess.Popen(['sbatch', '-Q'], stdin=subprocess.PIPE) for p in range(max_workers)]

        for w in workers:
            w.stdin.write(sbatch_script.encode('utf-8'))
            w.stdin.close()
        for w in workers: w.wait()
        all_done = True

    UNEXPECTED_EXIT = False

    
    # Check if linked directory workers are still running
    # by counting ".working" files in the log dir
    n_running = 0
    if linked_dir != None: 
        try:
            progress_msg = count_progress( log_dir)
            info("Progress: '%s'" % progress_msg)
            search = re.search('n(\d+)$', progress_msg)
            if search != None: n_running = int( search.group(1))
        except Exception as ex: warning( "Could not determine number of running linked workers. (%s)" % str(ex))


    # Create an empty file to signal jobs have finished
    try: 
        if all_done and n_running == 0: 
            open( done_filepath, 'a').close()
            if linked_dir != None: info("All linked jobs have finished.")
        elif all_done and n_running > 0: info("Job finished, but %i linked job%s still running." % (n_running, ("" if n_running != 1 else "s")))
    except Exception: 
        pass

    clean_up_wait_files()


def log_mem_usage( memory_filepath, n_workers):
    try:
        if not hasattr( log_mem_usage, "max_mem_usage"): log_mem_usage.max_mem_usage=0.0 # Init once
        used_mem_percentage = float( subprocess.check_output( "free | grep Mem | awk '{printf \"%.1f\", $3/$2 * 100.0}'", shell=True))
        log_mem_usage.max_mem_usage = max( log_mem_usage.max_mem_usage, used_mem_percentage)
        with open( memory_filepath, "a") as mf:
            mf.write("%s\tn%i\t%.1f%%\tmax:%.1f%%\n" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), n_workers, used_mem_percentage, log_mem_usage.max_mem_usage))
    except: pass

def get_machine_id():
    # Python 3.6 and Python 3.8 return different machine IDs!
    return str(uuid.UUID(int=uuid.getnode()))[-12:]


# progress: 0 to read, -1 to update, >0 initialize
# return -2 for error
# return -1 for unlimited
# return >= 0 remaining job count (0 means no more jobs should be run)
def update_max_jobs_count( folder, progress):

    if not folder.endswith("/"): folder += "/"
    maxjobs_filepath = "%s%s_%s" % ( folder, "maxjobs", get_machine_id())    

    file_exists = os.path.isfile( maxjobs_filepath)

    # max jobs file not used
    if progress <= 0 and not file_exists: return -1

    # Create new file
    if progress > 0 and not file_exists:
        try:
            with open( maxjobs_filepath, 'w') as f:
                f.write('%i\n' % progress)
        except Exception as ex:
            error("Could not write max job count to file '%s'. (%s)" % (maxjobs_filepath, str(ex)))
            return -2

        return progress


    # Read and update count
    retries = 4
    success = False
    prev_progress = -1

    while not success and retries > 0:            

        try:
            fd = os.open(maxjobs_filepath, os.O_RDWR | os.O_EXCL)
            # Read existing job count            
            prev_progress = int( os.read(fd, 100))            

            # No need to update
            if prev_progress == 0 and progress <= 0:
                os.close(fd)
                return 0

            os.lseek(fd, 0, 0) # rewind to beginning

            new_progress = max( 0, progress + prev_progress)

            # Update existing job count
            fo = os.fdopen(fd, "w+")
            fo.write( "%i\n" % new_progress)
            fo.close()
            success = True

        except ValueError as ex:
            warning( "invalid integer value in file '%s'. (%s)" % (maxjobs_filepath, str(ex)))
            retries = 0
        except:
            retries -= 1
            if retries > 0: time.sleep(random.randint(0, 2))

    if success: return prev_progress
    return -2




def save_owner_info( folder):

    if not folder.endswith("/"): folder += "/"
    owner_filepath = "%s%s_%s" % ( folder, "owner", get_machine_id())

    # If working directory has been deleted or renamed
    # getcwd() can throw an exception
    try: cwd = os.getcwd()
    except FileNotFoundError: cwd = os.path.expanduser( "~/")    

    try:
        with open( owner_filepath, "w") as f:
            f.write( cwd)
    except Exception as ex:
        error( "Could not write owner info file '%s'. (%s)" % (owner_filepath, str(ex)))

def is_owned( folder):

    if not folder.endswith("/"): folder += "/"
    owner_filepath = "%s%s_%s" % ( folder, "owner", get_machine_id())
    #info("Owner path: '%s'" % owner_filepath)
    return os.path.isfile( owner_filepath)

def save_basedir( folder):

    if not folder.endswith("/"): folder += "/"
    basedir_filepath = "%s%s" % ( folder, "basedir")

    # If working directory has been deleted or renamed
    # getcwd() can throw an exception
    try: cwd = os.getcwd()
    except FileNotFoundError: cwd = os.path.expanduser( "~/")    

    try:
        with open( basedir_filepath, "w") as f:
            f.write( cwd)
    except Exception as ex:
        error( "Could not write basedir file '%s'. (%s)" % (basedir_filepath, str(ex)))

def load_basedir( folder):

    if not folder.endswith("/"): folder += "/"
    basedir_filepath = "%s%s" % ( folder, "basedir")

    try:
        with open( basedir_filepath, "r") as f:
            basedir = f.read()
    except Exception as ex:
        error( "Could not read basedir file '%s'. (%s)" % (basedir_filepath, str(ex)))

    if not os.path.isdir( basedir): error( "Basedir '%s' not found.")
    os.chdir(basedir)



def parallel_worker(log_dir):

    global WORKER_THREAD, WORKER_FILE, RANDOMIZE_JOBS, SYMLINK_TARGET_DIR

    WORKER_THREAD = True
    WORKER_FILE = None
    commands = []
    exports = []

    if os.path.islink( log_dir):
        LINKED_JOB = True
        SYMLINK_TARGET_DIR = os.readlink( log_dir)
        if SYMLINK_TARGET_DIR.endswith("/"): SYMLINK_TARGET_DIR = SYMLINK_TARGET_DIR[:-1]         

    with open('%s/tasks' % log_dir, 'r') as f:

        line_num = 0

        for line in f:
            line_num += 1
            line = line.strip()
            if len( line) < 5: continue

            line_type = line[:4]
            line_str = line[4:]

            if line_type == "CMD:":
                commands.append( [line_str, list(), list( exports[:])])
            elif line_type == "TRG:":
                if len( commands) == 0:
                    sys.stderr.write("ERROR: No command specified for target '%s'.\n" % line[4:])
                    return
                commands[ -1][ 1].append( line_str)
            elif line_type == "EXP:":
                if line_str[ -1] == ";": line_str = line_str[:-1] #remove semicolon if there           
                export_var = line_str.partition("=")[ 0]
                # Find previous exports matching until the "=" sign
                matching = [s for s in range( len( exports)) if export_var in exports[ s]]
                if len( matching) == 0:
                    # Add new export
                    exports.append( line_str)
                else:
                    # Replace existing export(s)
                    for m in matching: exports[ m] = line_str
            elif line_type[ 0] == "#":
                #Comment line
                continue
            else:
                sys.stderr.write("Warning: Unknown line type '%s' on line %i.\n" % (line_type, line_num))
                continue
        #command = next(f).strip()[4:]
        #targets = [target.strip()[4:] for target in f]
        
    if len( commands) < 1: 
        sys.stderr.write("ERROR: No commands specified.\n")
        return

    cmd_index = 0    
    cmd_index_str = ""
    n_commands = len( commands) 
    remaining_jobs = -1 # -1 for unlimited jobs

    if RANDOMIZE_JOBS:        
        random.shuffle( commands)   

    for cmd_arr in commands:
    
        cmd_index += 1
        if n_commands > 1: cmd_index_str = "%s_" % cmd_index
        command = cmd_arr[ 0]
        targets = cmd_arr[ 1]
        exports = "".join( [exp+"; " for exp in cmd_arr[ 2]])
        # Create single empty target if none specified
        if len( targets) == 0: targets = ['']

        for target in targets:

            if len( target) == 0: target_in_fn = "cmd"
            else: target_in_fn = sanitize_path( target)

            # Shutdown requested
            #info("Is '%s' shutting down? %i" % (log_dir, is_shutting_down( log_dir)))
            #info("Is symlinked dir '%s' shutting down? %i" % (SYMLINK_TARGET_DIR, is_shutting_down( SYMLINK_TARGET_DIR)))
            if is_shutting_down( log_dir) or is_shutting_down( SYMLINK_TARGET_DIR): break

            remaining_jobs = update_max_jobs_count( log_dir, -1) # Worker ready to start new job substracts one            
            if remaining_jobs == 0: break # do not start new jobs

            filebase = '%s/%s%s' % (log_dir, cmd_index_str, target_in_fn)
            out_file = filebase + ".out"
            work_file = filebase + ".working"
            done_file = filebase + ".done"
            error_file = filebase + ".error"

            if os.path.isfile( done_file) and not os.path.isfile( out_file):
                warning( "File '%s' already exists! (skipping)" % done_file)
                continue

            if os.path.isfile( error_file) and not os.path.isfile( out_file): 
                warning( "File '%s' already exists! (skipping)" % error_file)
                continue
              
            out = open_exclusive( out_file)
            if not out: continue
            cmd_with_target = '''export x='%s';%s%s''' % (target, exports, command)
            out.write('%s\n%s\n' % (cmd_with_target, '-'*80))
            out.flush()

            WORKER_FILE = out_file

            try: open( work_file, 'a').close()            
            except Exception as ex: sys.stderr.write("ex setting working:%s\n" % str( ex))            

            # Create a file if running docker exec
            docker_containers = re.findall(r"docker exec ([^\s;\\]+)", cmd_with_target)
            docker_files = []
            for cont in docker_containers:
                try: 
                    docker_filename = "docker_%s_%s.container" % (get_machine_id(), cont)
                    docker_filename = ("%s/%s" % (log_dir, docker_filename)) 
                    open( docker_filename, 'a').close()            
                    docker_files.append( docker_filename)
                except Exception as ex: out.write("ERROR: Could not create docker file '%s'. (%s)\n" % (docker_filename, str( ex)))


            retval = False # Error if stopped during processing
            was_interrupted = True

            try:
                start_time = datetime.datetime.now()
                # RUN TASK
                retval = shell(cmd_with_target, stdout=out, stderr=out) 

                was_interrupted = False

                end_time = datetime.datetime.now()
                duration = end_time - start_time
                #out.write('%s\nJOB FINISHED. ELAPSED TIME WAS %s.\n' % ('-'*80, duration.strftime("%Hh%Mm%Ss")))
                # Duration str -> 00:00:00.0000 HH:MM:SS.decimals 

                if remaining_jobs >= 0: out.write( "remaining jobs: %i\n" % remaining_jobs)

                out.write('%s\nJOB FINISHED. ELAPSED TIME WAS %s\n' % ('-'*80, str( duration).split(".")[0 ]))
                out.close()
                
            except:
                pass
            finally:                
                WORKER_FILE = None

            # Done file for each individual task
            try: open( done_file, 'a').close()
            except Exception as ex: sys.stderr.write("exception setting done: %s\n" % str( ex))

            try: os.remove( work_file)
            except: pass

            
            if not was_interrupted:

                #info("removing docker files: %s" % str(docker_files ))
                for df in docker_files:
                    try: os.remove( df)
                    except: pass                

            if retval == False or was_interrupted:
                try: open( error_file, 'a').close()
                except Exception as ex: sys.stderr.write("ex setting err:%s\n" % str( ex))

                #info( "ISLINKED: %i" % LINKED_JOB)
                if len(SYMLINK_TARGET_DIR) and was_interrupted and LINKED_JOB:
                    if os.path.isfile( out_file ):
                        new_job_log = out_file.replace(".out", ".linked_out")
                        # Remove job status files so that it can be reprocessed
                        try: os.remove( new_job_log)                            
                        except: pass                        
                        try: os.remove( work_file)                            
                        except: pass
                        try: os.remove( done_file)                            
                        except: pass
                        try: os.remove( error_file)                            
                        except: pass

                        try: 
                            os.rename( out_file, new_job_log)
                            info( "Moved file job log to '%s'" % ( new_job_log))
                        except: pass

                        #info_str = "Remove these files to reprocess interrupted linked job: '%s'." % out_file.replace(".out",".*")                   
                        #info( info_str)


        if is_shutting_down( log_dir) or is_shutting_down( SYMLINK_TARGET_DIR) or remaining_jobs == 0: break




def reset_logs( log_dir, older_than_a_month_only=False):

    n_removed = 0
    n_skipped = 0
    n_failed = 0
    p = re.compile('.*?_\d\d\d\d-\d\d-\d\d_\d\d.\d\d.\d\d')
    today_str = time.strftime("job_%Y-%m-%d")

    for jobdir in glob( log_dir + "*/"):

        if not p.match( jobdir): continue

        if jobdir.find( today_str) == 0:
            # Do not remove jobs from the same day
            n_skipped += 1
            continue

        # As an extra safety feature, require jobdirs to 
        # have a "tasks" file
        if not os.path.isfile( jobdir + "tasks"): 
            
            # Linked jobs without a tasks file could have been removed from
            # the shared folder, clean up (unlink) these
            if os.path.islink( jobdir[:-1]):
                # Linked jobs
                rp = os.path.realpath( jobdir[:-1])
                try: 
                    os.unlink( jobdir[:-1])
                    n_removed += 1
                    info( "Unlinked '%s' -> '%s'" % (jobdir[:-1], rp))
                except: 
                    warning( "Could not unlink '%s' -> '%s'" % (jobdir[:-1], rp))
                    n_skipped += 1
            else: 
                # Clean up empty directories on local server
                # Non-linked jobs    
                try: 
                    # Throws OSError if directory is not empty
                    os.rmdir( jobdir[:-1])
                    n_removed += 1
                except: 
                    warning( "Could not remove dir '%s'" % jobdir[:-1])
                    n_skipped += 1
            continue   

        if older_than_a_month_only:
             tasks_file_mod_time = os.stat( jobdir + "tasks").st_mtime                
             last_time = (time.time() - tasks_file_mod_time) / float(60*60*24*30) # n_months
             #info("last time: '%s'" % str(last_time))
             if last_time < 2.0: # two months
                #info("Keeping: '%s'" % jobdir)
                n_skipped += 1
                continue

        try: 
            shutil.rmtree( jobdir)
            n_removed += 1
        except: 
            n_failed += 1


    if sum( [n_removed, n_skipped, n_failed]) > 0:
        sys.stderr.write(
            "INFO: Logs reset: %i removed (%i skipped, %i failed).\n" % 
            (n_removed,n_skipped,n_failed))

def get_job_status( aJobDir):

    if aJobDir == None or len( aJobDir) == 0:
        sys.stderr.write( "ERROR: No job log directory specified.\n")
        return "invalid"

    if aJobDir[ -1] != '/': aJobDir += '/'
    if not os.path.isdir( aJobDir): aJobDir = os.path.expanduser( "~/.jobs") + "/" + aJobDir

    if not os.path.isdir( aJobDir): return "invalid"            
    if not os.path.isfile( aJobDir + "tasks"): return "invalid" 

    if os.path.isfile( aJobDir + "done"): 
        files = os.listdir( aJobDir )
        n_err_files = len([ file for file in files if file.endswith(".error")])
        n_out_files = len([ file for file in files if file.endswith(".out")])
        #n_out_files = len( os.listdir( aJobDir + '*.out'))
        return ("done, e%i/%i" % (n_err_files, n_out_files))
    if os.path.isfile( aJobDir + "interrupted"): return "interrupted" 
    if os.path.isfile( aJobDir + ("interrupted_%s" % get_machine_id())): return "interrupted" 
    if os.path.isfile( aJobDir + "shutdown"): return "shutting down"
    if os.path.isfile( aJobDir + "waiting"): return ("waiting %s" % count_progress( aJobDir))

    return ("running %s" % count_progress( aJobDir))


def count_progress( aJobdir): 

    if aJobdir[ -1] != '/': aJobdir += '/'
    #if not os.path.isdir( aJobDir): aJobDir = os.path.expanduser( "~/.jobs") + "/" + aJobDir

    try:

        if not os.path.isdir( aJobdir): return "?/? (dir not found) n0"
        if not os.path.isfile( aJobdir + "tasks"): return "?/? (tasks file not found) n0" 
    
        n_tasks = 0
        n_cmds = 0
    
        with open(aJobdir + "tasks", "r") as f:
            for line in f:
                if line.startswith( "TRG:"): n_tasks += 1
                elif line.startswith( "CMD:"): n_cmds += 1

        if n_tasks == 0: n_tasks = n_cmds
        if n_tasks == 0: return "0/0 (no tasks) n0"
    
        n_done = len( [name for name in os.listdir(aJobdir) if name.endswith(".done")])
        n_error = len( [name for name in os.listdir(aJobdir) if name.endswith(".error")])
        n_working = len( [name for name in os.listdir(aJobdir) if name.endswith(".working")])
        error_str = "" if n_error < 1 else ('\x1b[2;31;40me%i\x1b[0m ' % n_error)

        return "%i/%i %s(%.1f%%) n%i" % (n_done,n_tasks, error_str, (n_done/float(n_tasks)*100.0),n_working)

    except Exception as ex: 
        return "?/? (%s) n0" % str( ex)
        


def check_job_status( aJobDir):

    if aJobDir == None or len( aJobDir) == 0:
        sys.stderr.write( "ERROR: No job log directory specified.\n")
        sys.exit( -2)

    if aJobDir[ -1] != '/': aJobDir += '/'

    if not os.path.isdir( aJobDir):

        at_default_path = os.path.expanduser( "~/.jobs") + "/" + aJobDir
        if not os.path.isdir( at_default_path):        
            sys.stderr.write( "ERROR: Path '%s' is not valid directory.\n" % aJobDir)
            sys.exit( -1) 
        else:
            aJobDir = at_default_path
        
    if not os.path.isfile( aJobDir + "tasks"):
        sys.stderr.write( "ERROR: Directory '%s' does not appear to be a valid log folder.\n" % aJobDir)
        sys.exit( -3) 

    if os.path.isfile( aJobDir + "done"):
        sys.stderr.write( "INFO: This job appears to have finished already.\n")
        sys.exit(0) 

    if os.path.isfile( aJobDir + "interrupted"):
        sys.stderr.write( "INFO: This job appears to have been interrupted.\n")
        sys.exit(0) 

    return aJobDir

def cancel_shutdown( aJobDir):

    aJobDir = check_job_status( aJobDir)
    shutdown_file = "%sshutdown" % aJobDir 

    if os.path.isfile( shutdown_file ):
        try:
            os.remove( shutdown_file)
            info( "Job shutdown cancelled.")
        except:
            sys.stderr.write( "ERROR: Could cancel shutdown for job '%s'.\n" % aJobDir)
    else:      
        info( "Job '%s' is not shutting down." % aJobDir)  

def request_shutdown( aJobDir):

    if len( aJobDir) == 0: return
    aJobDir = check_job_status( aJobDir)
    shutdown_file = "%sshutdown" % aJobDir 

    if os.path.isfile( shutdown_file ): 
        info( "Job already shutting down after current tasks finish.")
        return

    try:
        with open( shutdown_file, "w") as f:
            f.write("Job shutting down after current tasks finish.\n")
        sys.stdout.write("INFO: Job shutting down after current tasks finish.\n")
    except:
        sys.stderr.write( "ERROR: Could not write to directory '%s'.\n" % aJobDir)        


def is_interrupted( aJobDir) :

    if not aJobDir or len( aJobDir) == 0: return False

    try:
        if aJobDir[ -1] != '/': aJobDir += '/'    
        interrupted_file = "%sinterrupted" % aJobDir
        interrupted_mid_file = "%sinterrupted_%s" % (aJobDir, get_machine_id())
        return os.path.isfile( interrupted_file) or \
               os.path.isfile( interrupted_mid_file)
    except: pass

    return False


def is_shutting_down( aJobDir) :

    if not aJobDir or len( aJobDir) == 0: return False

    try:
        if aJobDir[ -1] != '/': aJobDir += '/'    
        shutdown_file = "%sshutdown" % aJobDir
        interrupted_file = "%sinterrupted" % aJobDir
        interrupted_mid_file = "%sinterrupted_%s" % (aJobDir, get_machine_id())
        done_file = "%sdone" % aJobDir
        return (os.path.isfile( shutdown_file) or \
               os.path.isfile( interrupted_file) or \
               os.path.isfile( interrupted_mid_file) or \
               os.path.isfile( done_file))
    except: pass

    return False

    
def request_more_workers( aDir):

    aDir = check_job_status( aDir)

    if is_shutting_down( aDir):
        sys.stderr.write("ERROR: This job has been marked for shutdown.\n")
        return

    n_new_workers = 1
    n_retries = 2
    comments = []
    success = False
    moar_filepath = ("%smoar_%s" % (aDir, get_machine_id()))

    while not success and n_retries > 0:

        # Read existing
        try:
            with open( moar_filepath, "r") as mfile:
                for line in mfile:
                    line = line.strip()
                    if line.startswith("#"): comments.append( line)
                    elif len( line): n_new_workers = int( line)+1
        except Exception:
            pass

        if n_new_workers < 1: n_new_workers = 1
        if n_new_workers > 100: n_new_workers = 100

        # Write new
        try:            
            with open( moar_filepath, "w") as mfile:
                if len( comments):
                    mfile.write( "\n".join( comments))
                    mfile.write( "\n")
                mfile.write( "%i\n" % n_new_workers)
            success = True

        except Exception:
            n_retries -= 1
            time.sleep( 1) # wait 1 sec


    if not success:
        error("Failed to request more workers for job '%s'." % aDir)
    else:
        info("INFO: Request submitted for %i new worker%s." % (n_new_workers, ("" if n_new_workers == 1 else "s")))


def sort_job_dirs( dirs):

    times = []

    for path in dirs:
        
        dt = None

        # Parse creation time from folder name
        # Note: Does not account for multiple jobs created on the same second! 
        #try:            
        #    folder = os.path.basename(os.path.normpath( path))
        #    #datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')
        #    dt = datetime.datetime.strptime(folder, 'job_%Y-%m-%d_%H.%M.%S')
        #    #info("path: '%s'" % path)
        #    #info("folder: '%s'" % folder)
        #except ValueError:
        #    # Non-standard folder name
        #    pass
        #if dt == None:

        try:                
            # Task file modification date should be job creation time
            dt = datetime.datetime.fromtimestamp( os.path.getmtime( path + "/tasks"))
        except: 
            dt = datetime.datetime.fromtimestamp( os.path.getmtime( path))                  
            

        times.append( dt)


    sorted_dirs = [x for _,x in sorted(zip(times,dirs), reverse=True)]
    return sorted_dirs


def log_dir_for_index( aJobIndex):    

    try:
        index = int( aJobIndex) 
        if index < 1: error( "Index '%i' is not a valid job index." % index)
    except Exception:
        # Accept directory instead of index
        if os.path.isdir( aJobIndex): 
            return aJobIndex
        else:
            error( "Path '%s' is not valid directory." % aJobIndex)
            return None

    default_dir = os.path.expanduser( "~/.jobs")
    if not os.path.isdir( default_dir):
        error( "Path '%s' is not valid directory.\n" % default_dir)
        sys.exit( -1) 

    try: saved_dir = os.getcwd()
    except FileNotFoundError: saved_dir = os.path.expanduser( "~/")   

    os.chdir( default_dir)  
    files = filter( os.path.isdir, os.listdir( default_dir))
    files = [os.path.join( default_dir, f) for f in files] # add path to each file
    #print( [f + str(os.path.isdir(f)) for f in os.listdir( default_dir)])
    if index > len( files):
        error( "No log directory found for index %i" % index)

    if len( files) == 0:
        error( "No log directories found in '%s'.\n" % default_dir)
        sys.exit( -1)

    # NOTE: getmtime of the folder gets updated when files are updated, so order of 
    #       jobs in listing may change between consecutive calls
    #files.sort( key=lambda x: os.path.getmtime(x), reverse=True)
    files = sort_job_dirs( files)

    #i = 0
    #for f in files:
    #    i += 1
    #    info( "%i: %s" % (i, f))
    #sys.exit( 0)

    view_dir = files[index-1] # listing has 1-based numbering


    os.chdir( saved_dir) 
    return view_dir

def show_runtimes( runtimes ):

    view_dir = log_dir_for_index( runtimes)
    outfiles = [name for name in os.listdir( view_dir) if name.endswith(".out")]

    #outfiles.sort( key=lambda x: os.path.getmtime(view_dir + "/" + x), reverse=True)
    outfiles.sort()

    for of in outfiles:
        cmd = "tail -n 10 \"%s/%s\" | egrep -o \"ELAPSED TIME WAS.*\"" % (view_dir, of)
        #info( "CMD: " + cmd)
        #result = subprocess.run( cmd, capture_output=True, shell=True, stderr=subprocess.DEVNULL, text=True)
        result = None
        try:
            result = subprocess.check_output( cmd, shell=True, stderr=subprocess.DEVNULL)
        except Exception: pass
        
        if result and len( result):
            sys.stdout.write("%s: %s\n" % (of, result.decode("utf-8").rstrip()))
        else:
            sys.stdout.write("%s: Timing data missing\n" % (of))
        #os.system( cmd)


def view_job_logs( aJobIndex):

    view_dir = log_dir_for_index( aJobIndex)
    info("Viewing job 1 output files at '%s'" % view_dir)

    outfiles = [name for name in os.listdir( view_dir) if name.endswith(".out")]
    # Most recent 
    outfiles.sort( key=lambda x: os.path.getmtime(view_dir + "/" + x), reverse=False)
    #if len( outfiles) > 1000: outfiles = outfiles[-1000:]

    cmds_and_targets = []
    if "1_cmd.out" in outfiles:
        try:
            with open( view_dir + "/tasks", "r") as tf:
                for line in tf:
                    if line.startswith( "CMD:") or line.startswith("TRG:"):
                        cmds_and_targets.append( line[:70].rstrip())
                        # Add ellipsis if cut
                        if len( line) > 70: cmds_and_targets[ -1] = cmds_and_targets[ -1] + "..."
        except:
            error("Could not read task file for job '%s'." % view_dir)


    #print tails
    n_printed = 0
    for outfile in outfiles:
        if n_printed == 1000: 
            warning( "Only the first 1000/%i log files were tailed." % len( outfiles))
            break        
        task_str = ""
        if len( cmds_and_targets) > n_printed:
            try:
                digits = re.findall(r'(\d+)_cmd\.out', outfile)
                if len( digits) > 0: 
                    digit = int( digits[ 0])
                    task_str = "(task %i: %s)" % (digit,cmds_and_targets[digit-1])
            except: pass        
            

        sys.stderr.write( "###\n")
        sys.stderr.write( "### LOGFILE: '%s' %s\n" % (outfile, task_str))
        sys.stderr.write( "###        : '%s/%s'\n" % (view_dir, outfile))
        sys.stderr.write( "###\n" )
        try:
            shell("tail %s/%s" % (view_dir, outfile))#, stdout=stdout=subprocess.PIPE, stderr=out)
            shell("if [[ $(grep -m 1 -c -i error %s/%s) == 1 ]] ; then echo 'INFO: Errors found:'; grep -i error %s/%s; else echo 'INFO: No errors found'; fi" % (view_dir, outfile, view_dir, outfile))            
        except Exception as ex:
            error( "tail or grep cmd failed for file '%s': '%s'" % (outfile,str( ex)))
            break
        #info( "####END LOGFILE: '%s'" % outfile)
        n_printed += 1


    # print progress
    status = "?"
    try: status = get_job_status( view_dir)
    except: pass

    name = ""
    try: name = open("%s/name" % view_dir, 'r').readline().strip()
    except: pass
    if len(name) == 0: name = view_dir

    sys.stdout.write( "\n### Job '%s' status: [%s]\n" % (name, status))  
    

def get_n_jobs_running( special_dir=None, exclude=True):

    default_dir = os.path.expanduser( "~/.jobs")
    if not os.path.isdir( default_dir):
        sys.stderr.write( "ERROR: Path '%s' is not valid directory.\n" % default_dir)
        sys.exit( -1) 

    try: saved_dir = os.getcwd()
    except FileNotFoundError: saved_dir = os.path.expanduser( "~/")   
    os.chdir( default_dir)    
    files = filter( os.path.isdir, os.listdir( default_dir))
    files = [os.path.join( default_dir, f) for f in files] # add path to each file
    #print( [f + str(os.path.isdir(f)) for f in os.listdir( default_dir)])

    if len( files) == 0: return 0

    n_running = 0
    job_names = []

    for jobdir in files:

        if special_dir != None:
            if exclude: 
                if jobdir == special_dir: continue
            else:
                if jobdir != special_dir: continue

        if os.path.isfile( jobdir + "/done") or \
           os.path.isfile( jobdir + "/interrupted") or \
           os.path.isfile( jobdir + "/shutdown"): continue
        elif not os.path.isfile( jobdir + "/tasks"): continue # Requires tasks file
        
        n_running += 1

        name = "unnamed"
        try: name = open("%s/name" % jobdir, 'r').readline().strip()
        except: pass
        job_names.append( name)

    os.chdir( saved_dir)
    return (n_running, job_names)
        


def show_previous_log_dirs( aN=1):

    try:
        n_dirs = int( aN)
    except Exception:
        sys.stderr.write( "ERROR: '--show' expected an integer argument, got '%s'.\n" % aN)
        sys.exit( -2)
    
    default_dir = os.path.expanduser( "~/.jobs")
    if not os.path.isdir( default_dir):
        sys.stderr.write( "ERROR: Path '%s' is not valid directory.\n" % default_dir)
        sys.exit( -1) 

    # If working directory has been deleted or renamed
    # getcwd() can throw an exception
    try: saved_dir = os.getcwd()
    except FileNotFoundError: saved_dir = os.path.expanduser( "~/")

    os.chdir( default_dir)  
    files = filter( os.path.isdir, os.listdir( default_dir))
    files = [os.path.join( default_dir, f) for f in files] # add path to each file
    #print( [f + str(os.path.isdir(f)) for f in os.listdir( default_dir)])

    if len( files) == 0:
        warning( "No log directories found in '%s'.\n" % default_dir)
        os.chdir( saved_dir)
        return # Nothing to show

    #Getting file creation dates, on the other hand, is fiddly and platform-dependent, differing even between the three big OSes:
    #On Linux, this is currently impossible, at least without writing a C extension for Python
    #files.sort( key=lambda x: os.path.getmtime(x), reverse=True)
    files = sort_job_dirs( files)

    dirs = files[:n_dirs]
    i = 0
    max_len = 5
    buf = [] 

    for d in dirs:
        i += 1
        name = "unnamed"        
        try: name = open("%s/name" % d, 'r').readline().strip()
        except: pass

        link_indicator = " "   
        if os.path.islink( d): link_indicator = "@" if is_owned( d) else ">"

        status = "?"
        try: status = get_job_status( d)
        except: pass

        # Add some colors
        if status.startswith("done"):
            if status.find("e0") > 0:
                # No errors
                status = status.replace( "done", '\x1b[2;32;40m' + "done" + '\x1b[0m')
            else:
                # Some errors
                status = status.replace( "done", '\x1b[2;31;40m' + "done" + '\x1b[0m')
        else:
            status = status.replace( "running", '\x1b[1;32;40m' + "running" + '\x1b[0m')
        #status = status.replace( "interrupted", '\x1b[2;36;40m' + "interrupted" + '\x1b[0m')
        
        #i_str = '\x1b[1;35;40m' + ("%2i" % i) + "." + '\x1b[0m'
        #d     = '\x1b[2;35;40m' + d + '\x1b[0m'
        link_indicator = '\x1b[1;37;40m' + link_indicator + '\x1b[0m'
        #name = '\x1b[2;35;40m' + "('" + '\x1b[1;35;40m' + name + '\x1b[2;35;40m' + "')" + '\x1b[0m'
        name = "('" + name + "')" 
        max_len = max( max_len, len( name))

        #sys.stdout.write( "%s %s  %s  ('%s') [%s]\n" % (i_str,d, link_indicator, name, status))
        #sys.stdout.write( "%2i. %s  %s  %-55s [%s]\n" % (i,d, link_indicator, name, status))
        buf.append( (i,d, link_indicator, name, status))

    max_len = min( 55, max_len)
    template = "%%2i. %%s  %%s  %%-%is [%%s]" % max_len # Set required space for job names

    sys.stdout.write( "\n".join( [template % t for t in buf ]) + "\n")
    os.chdir( saved_dir)

    #sys.stdout.write( "\n".join( files[:n_dirs]) + "\n")


def SystemSpecificFolderSeparators( aPath):

    if os.name == 'nt': return aPath.replace( "/", "\\")
    return aPath.replace( "\\", "/")

def InsertFolderSeparator( folder):

    if len( folder) == 0: return folder

    folder_sep = "/"
    if folder.find( "\\") >= 0: folder_sep = "\\"

    #Make sure folder ends with a folder separator character
    if folder[ -1] != folder_sep:
        folder += folder_sep

    return SystemSpecificFolderSeparators( folder)

def ChangeToFolder( filepath, folder):

    path, filename = os.path.split( filepath)
    return InsertFolderSeparator( folder) + filename


def KillDockers():

    if WORKER_THREAD: return

    # Signals do not propage into containers
    # https://github.com/moby/moby/issues/9098    
    n_containers = len( RUNNING_DOCKERS)        
    if n_containers == 0: return

    try:
        
        mid = get_machine_id()
        n_killed = 0
        # Kill only executing containers on this VM
        info("Stopping active docker containers, please wait...")
        for docker_file in glob( JOB_DIR+"/docker_"+mid+"_*.container"):
            try:
                containers = re.findall(r"docker_"+mid+"_(.+)\.container", docker_file)
                if len( containers) > 0:
                    #shell( "docker kill %s" % containers[ 0])
                    shell( "docker stop %s" % containers[ 0])                
            except Exception as ex: warning( "Could not kill docker container '%s'. (%s)" % (docker_file, str( ex)))            
            n_killed += 1
            try: os.remove( docker_file)
            except: pass

        info("Stopped %i containers." % n_killed)
        if n_killed == 0: info("All running containers can be stopped with 'docker stop $(docker ps -aq)'.")

        #time.sleep( 2)

    except Exception as ex: warning( "Could not kill all docker containers. (%s)" % str( ex))      


@atexit.register
def exit_handler():

    global JOB_DIR, UNEXPECTED_EXIT, WORKER_THREAD, WRITE_WAIT_FILES, WRITE_WAIT_HANDLES, \
           WORKER_FILE, LINKED_JOB, SYMLINK_TARGET_DIR, RUNNING_DOCKERS

    #if not WORKER_THREAD:
    #    sys.stderr.write( "Exitting...\n")
    #    sys.stderr.write( "unexp: %s, jobdir %s\n" % (str(UNEXPECTED_EXIT), str(JOB_DIR)))

    if not WORKER_THREAD and UNEXPECTED_EXIT and len( JOB_DIR) > 0:

        sys.stderr.write("\nERROR: Processing interrupted by user.\n")

        try: os.remove( "%s/waiting" % JOB_DIR)
        except: pass

        retries = 1

        if not LINKED_JOB:
            while retries > 0:
                try:             
                    open( ("%s/interrupted" % JOB_DIR), 'a').close()            
                    retries = 0
                except Exception as ex:
                    sys.stderr.write("Could not write 'interrupted' file. Ex:%s\n" % str( ex))
                    retries -= 1
                    time.sleep( 1)
        else:
            open( ("%s/interrupted_%s" % (JOB_DIR, get_machine_id())), 'a').close() 
            info("Stopping this linked job does not stop processing on other machines.")

    KillDockers()

    if WORKER_THREAD and WORKER_FILE != None:
        info( "Worker '%s' quitting." % WORKER_FILE)
        worker_file = WORKER_FILE.replace(".out", ".working")        
        try: os.remove( worker_file)
        except Exception as ex: warning("Could not remove worker file: '%s' (%s)" % (worker_file, str( ex)))

        if LINKED_JOB:
            try: # Rename interrupted linked work files to ".link"
                move_dest_base = WORKER_FILE.replace(".out", ("_"+get_machine_id()+".link"))
                move_dest = move_dest_base[:]
                retries = 0
                while os.path.isfile( move_dest) and retries < 10:
                    move_dest = move_dest_base.replace(".link", ("_%i.link" % retries))
                    retries += 1
                if LINKED_JOB and UNEXPECTED_EXIT: os.rename( WORKER_FILE, move_dest)
            except Exception as ex: warning("Could not move worker file: '%s' (%s)" % (WORKER_FILE, str( ex)))
        WORKER_FILE = None

    if not WORKER_THREAD: clean_up_wait_files()

# Not used
def CreateLocalFolderForLinkedJob():

    # Clean up symlink to linked dir if processing interrupted
    if LINKED_JOB and not WORKER_THREAD and UNEXPECTED_EXIT and len( SYMLINK_TARGET_DIR) > 0:

        try:
            if os.path.islink( JOB_DIR): 
                # JOB_DIR == real_job_dir                
                real_job_dir = os.readlink( JOB_DIR)
                if real_job_dir.endswith("/"): real_job_dir = real_job_dir[:-1]               
                # Change symlinked dir into real local dir
                info( "Creating a local copy of linked dir '%s' to '%s'." % (real_job_dir, JOB_DIR))
                os.unlink( JOB_DIR)
                os.mkdir( JOB_DIR)
               
                try: open( ("%s/interrupted" % JOB_DIR), 'a').close()
                except: warning("Could not set linked dir 'interrupted'.")
                try: shutil.copyfile( real_job_dir+"tasks", JOB_DIR+"/tasks" )
                except: warning("Could not copy 'tasks' from linked job.")
                try: shutil.copyfile( real_job_dir+"name", JOB_DIR+"/name" )
                except: warning("Could not copy 'name' from linked job.")
                # Copy output from this VM
                mid = get_machine_id()
                n_copied = 0
                info("Copy str '%s'" % (real_job_dir + "*"+mid+"*.link")) 
                for linked_worker_file in glob( real_job_dir + "*"+mid+"*.link"):                
                    # strip machine_id, .link -> .out
                    info("Copying: '%s' -> '%s'" % (linked_worker_file, dest_file)) 
                    dest_file = ChangeToFolder( linked_worker_file, JOB_DIR).replace(".link", ".out").replace( "_"+mid, "")
                    try: 
                        shutil.copyfile( linked_worker_file, dest_file)
                        n_copied += 1
                    except: warning("Copy failed '%s' -> '%s'" % (linked_worker_file, dest_file))                

                info( "%i output log files (.out) copied to local folder '%s'." % (n_copied, JOB_DIR))

            else: warning("Linked dir is not a symlink: '%s'." % (JOB_DIR))
        except Exception as ex: 
            warning("Could not unlink linked dir: '%s' (%s)." % (JOB_DIR, str( ex)))    



def CheckWaitFiles( aWaitFiles):

    global READ_WAIT_FILES, WRITE_WAIT_FILES, WRITE_WAIT_HANDLES

    wait_file_retries = 3
    wait_status = None

    if aWaitFiles and len( aWaitFiles) != 0:
        for wf in aWaitFiles:

            wait_type = None
            file_handle = None
            
            while wait_file_retries > 0:
                
                if os.path.isfile( wf):                    
                    info("Waiting for file '%s' to be removed before starting tasks." % wf)
                    wait_type = "R"
                    break

                file_handle = open_exclusive( wf)

                if not file_handle: 
                    time.sleep( 1)
                    wait_file_retries -= 1
                    continue # Retry
                else:
                    wait_type = "W"
                    break
                            
            if not wait_type:
                error("Could not create wait file '%s'." % wf)
            elif wait_type == "W":
                WRITE_WAIT_FILES.append( wf)
                WRITE_WAIT_HANDLES.append( file_handle)
                info("Removing file '%s' after job finishes." % wf)
            elif wait_type == "R":
                READ_WAIT_FILES.append( wf)
            else:
                error("Unknown wait type")





def read_tasks_file_as_input( tasks_file):

    #import StringIO # Python 2
    from io import StringIO

    ret_arr = []
    input_arr = []
    with open( tasks_file, "r") as tf:
        for line in tf:
            if line.startswith("CMD:"):
                ret_arr.append( line[4:].strip())
            elif line.startswith("TRG:"):
                input_arr.append( line[4:].strip())

    if len( input_arr):
        #info("inputting to stdin")
        sys.stdin = StringIO( " ".join(input_arr))

    return ret_arr


def set_up_wait_files( wait_files):

    global WRITE_WAIT_FILES, WRITE_WAIT_HANDLES, READ_WAIT_FILES
   
    # Nothing to process
    if len( wait_files) == 0: return

    # Make sure no duplicates
    wait_files = list( set( wait_files))

    wait_file_retries = 3

    if wait_files and len( wait_files) != 0:

        for wf in wait_files:

            wait_type = None
            file_handle = None
            
            while wait_file_retries > 0:
                
                if os.path.isfile( wf):                    
                    info("Waiting for file '%s' to be removed before starting tasks." % wf)
                    wait_type = "R"
                    break

                file_handle = open_exclusive( wf)

                if not file_handle: 
                    time.sleep( 1)
                    wait_file_retries -= 1
                    continue # Retry
                else:
                    wait_type = "W"
                    break
                            
            if not wait_type:
                error("Could not create wait file '%s'." % wf)
            elif wait_type == "W":
                WRITE_WAIT_FILES.append( wf)
                WRITE_WAIT_HANDLES.append( file_handle)
                info("Removing file '%s' after job finishes." % wf)
            elif wait_type == "R":
                READ_WAIT_FILES.append( wf)
            else:
                error("Unknown wait type")




            
def process_wait_files():

    global READ_WAIT_FILES

    if len( READ_WAIT_FILES) == 0: return []

    isfiles = [os.path.isfile( rwf) for rwf in READ_WAIT_FILES]
    #remaining_files = [x for x, y in zip(READ_WAIT_FILES, isfiles) if y == True]
    removed_files = [x for x, y in zip(READ_WAIT_FILES, isfiles) if y == False]
    READ_WAIT_FILES = [x for x in READ_WAIT_FILES if x not in removed_files]

    # Change READ files to WRITE files
    # Also re-add to READ_WAIT_FILES if another process using the same file was started
    set_up_wait_files( removed_files)


    return READ_WAIT_FILES




def clean_up_wait_files():

    global WRITE_WAIT_FILES, WRITE_WAIT_HANDLES, READ_WAIT_FILES

    for wi in range( len(WRITE_WAIT_FILES)):

        try: WRITE_WAIT_HANDLES[ wi].close()
        except Exception as ex: warning("Could not close wait file: '%s'. (%s)" % (WRITE_WAIT_FILES[ wi], str( ex)))            

        try: os.remove( WRITE_WAIT_FILES[ wi])            
        except Exception as ex: warning("Could not remove wait file: '%s'. (%s)" % (WRITE_WAIT_FILES[ wi], str( ex)))

    WRITE_WAIT_FILES=[]
    WRITE_WAIT_HANDLES=[]
    READ_WAIT_FILES=[]



#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':

    # Users do not need to know about the "parallel worker" invocation. 
    if len(sys.argv) >= 3 and sys.argv[1] == 'worker':
        parallel_worker(log_dir=sys.argv[2])
        exit()

    args = docopt.docopt(__doc__)

    # DEBUG
    #info( "args" + str(args))
    #sys.exit( 0)

    partition = args['--partition']
    if partition == 'env': 
        try:
            partition = os.environ['PARALLEL_PARTITION']
        except KeyError:
            sys.stderr.write( "WARNING: Environment variable 'PARALLEL_PARTITION' not set, using partition: 'local'.\n")
            partition = local

    log_dir = args['--logdir']
    if log_dir != None and len( log_dir) > 0:
         if log_dir[ -1] != '/': log_dir += '/'
    else:
        log_dir = "~/.jobs/"
    log_dir = os.path.expanduser( log_dir)
    #print( "LOGDIR: '%s'" % log_dir)


    # Job control commands
    inc_workers_in_dir = args['--increase']
    if inc_workers_in_dir and len( inc_workers_in_dir):
        request_more_workers( log_dir_for_index( inc_workers_in_dir))
        UNEXPECTED_EXIT = False
        sys.exit( 0)        

    show_dirs = args['--show']
    if show_dirs and len( show_dirs):
        show_previous_log_dirs( show_dirs)
        UNEXPECTED_EXIT = False
        sys.exit( 0)   

    show_dirs = args['--show-jobs']
    if show_dirs == True:
        show_previous_log_dirs( "10")
        UNEXPECTED_EXIT = False
        sys.exit( 0) 

    try:
        max_jobs = int(args['--max-jobs'])
        if max_jobs < 0: raise ValueError
        if max_jobs > 0: info("Processing a maximum of %i jobs." % max_jobs)
    except:
        error( "The --max-jobs argument should be an positive integer value, 0 for unlimited.")


    show_id = args['--id']
    if show_id == True:
        sys.stdout.write("%s\n" % get_machine_id())
        UNEXPECTED_EXIT = False
        sys.exit( 0) 

    if bool( args['--disarrange']):
        RANDOMIZE_JOBS = True

    view_logs = args['--view']
    if view_logs and len( view_logs):
        view_job_logs( view_logs)
        UNEXPECTED_EXIT = False
        sys.exit( 0)  

    get_path = args['--path']
    if get_path and len( get_path):
        sys.stdout.write( log_dir_for_index( get_path))
        UNEXPECTED_EXIT = False
        sys.exit( 0)     

    view_logs = args['--view-job']
    if view_logs == True:
        view_job_logs( "1")
        UNEXPECTED_EXIT = False
        sys.exit( 0) 

    shutdown_dir = args['--shutdown']
    if shutdown_dir and len( shutdown_dir):
        request_shutdown( log_dir_for_index( shutdown_dir))
        UNEXPECTED_EXIT = False
        sys.exit( 0)  

    shutdown_dir = args['--resume']
    if shutdown_dir and len( shutdown_dir):
        cancel_shutdown( log_dir_for_index( shutdown_dir))
        UNEXPECTED_EXIT = False
        sys.exit( 0)  

    mem_usage_dir = args['--memory']
    if mem_usage_dir and len( mem_usage_dir):   
        mu_file = "memory_usage_%s" % get_machine_id()
        mu_file_path = "%s/%s" % (log_dir_for_index( mem_usage_dir), mu_file)
        try:            
            shell("tail %s" % mu_file_path)
            info( "File: '%s'" % mu_file_path)
        except Exception as ex:
            error( "tail cmd failed for file '%s': '%s'" % (mu_file, str( ex)))
        UNEXPECTED_EXIT = False
        sys.exit( 0)

    runtimes = args['--runtimes']
    if runtimes and len( runtimes):
        show_runtimes( runtimes )
        sys.exit( 0)

    if bool( args['--reset']): 
        reset_logs( log_dir, older_than_a_month_only=False)    
        UNEXPECTED_EXIT = False
        sys.exit( 0)


    if bool( args['--maintenance']): 
        reset_logs( log_dir, older_than_a_month_only=True)    
        UNEXPECTED_EXIT = False
        sys.exit( 0)


    view_logs = args['--view']
    if view_logs and len( view_logs):
        view_job_logs( view_logs)
        UNEXPECTED_EXIT = False
        sys.exit( 0)

    wait_to_run = args['--wait']
    # Change index number to a directory
    # Do it here, so this new job-to-be-started
    # does not change indexing
    if wait_to_run and len( wait_to_run):
        wait_to_run = log_dir_for_index( wait_to_run)


    if bool( args['--wait-all']):
        wait_to_run = "ALL"

    # --wait-file can be specified multiple times
    wait_files = args['--wait-file']
    set_up_wait_files( wait_files)
                

    wait_comps = args['--file-comp']

    if wait_comps and len( wait_comps) != 0:
        for wf in wait_comps:
            WATCHED_FILES.append( wf)
            info("Waiting for file '%s' to be completed before starting tasks." % wf)

    n_workers = args['--workers']

    # Linking directories on network drives to be shared by multiple
    # instance of parallel
    linked_dir = args['--linked'] 
    attached_to_linked = False    
    if linked_dir and len( linked_dir):  

        if not linked_dir.endswith("/"): linked_dir += "/"    

        if os.path.isdir( linked_dir):  
            try: os.remove( "%sinterrupted_%s" % (linked_dir, get_machine_id()))
            except: pass

            linked_job_status = get_job_status( linked_dir)
            if not linked_job_status.startswith("running") or linked_job_status.startswith("waiting"):
                error( "Linked job '%s' is not running. (%s)" % (linked_dir, linked_job_status))                     
        existing_tasks_file = linked_dir + "tasks"
        # Read tasks from linked folder
        if os.path.isfile( existing_tasks_file):
            if (args['<command>'] != None and len( args['<command>']) != 0) or bool( args['--file']):
                error("Linked folder '%s' already has tasks." % linked_dir)
            if not sys.stdin.isatty():
                error("Cmd modifiers will be read from an existing task file only, not stdin.")

            args['<command>'] = read_tasks_file_as_input( existing_tasks_file)

            if args['<command>'] == None or len( args['<command>']) == 0:
                error("No tasks found in linked folder '%s'." % linked_dir)
            attached_to_linked = True

            n_workers_file = linked_dir + "n_workers"
            if n_workers == None: 
                try:
                    with open(n_workers_file, "r") as f: 
                        n_workers = f.readline().rstrip()
                    info("Using %s workers as recommended by linked job." % n_workers)
                except: pass


    try:
        screen = os.environ['STY']
        if screen == None or len( screen) < 1: raise ValueError        
        if args['--job-name'] == "job": args['--job-name'] = screen
    except:
        if linked_dir and len( linked_dir): error("Linked jobs should be run inside screen only.")

    # Debug
    #info( "WF:")
    #info( wait_file)
    #info( "ARGS: '%s'" % str(args['<command>']))
    #isatty = sys.stdin.isatty()
    #info( "Isatty: %i" % isatty)
    #if not isatty:
    #    si = []
    #    for line in sys.stdin: si.append( line)
    #    info( "stdin: %s" % " ".join( si))
    #sys.exit( 0)


    # Starting a new job
    if bool( args['--file']): 
        if attached_to_linked: error("Flag '--file' is not compatible with '-L'.")
        filename = args['<command>']
        if not filename or len( filename) == 0:
            sys.stderr.write( "ERROR: File not specified.\n")
            sys.exit(-1)
        elif not os.path.isfile( filename):
            sys.stderr.write( "ERROR: File '%s' not found.\n" % filename)
            sys.exit(-1)

        try:
            with open( filename) as f:
                commands = f.readlines()
        except Exception as ex:
            sys.stderr.write( "ERROR: Could not read command file '%s'.\n" % filename)
            sys.stderr.write( "       An exception of type {0} occured.\n".format( type( ex).__name__))
            sys.stderr.write( "       Arguments: {1!r}\n".format( ex.args))
            UNEXPECTED_EXIT = False
            sys.exit( -2)
    elif attached_to_linked:
        # args['<command>'] have been read from a tasks file
        commands = args['<command>']
    else:
        # commands specified on command line
        if args['<command>'] == None or len( args['<command>']) == 0:
            #sys.stderr.write( "No command given.\n")
            show_previous_log_dirs( "10")
            UNEXPECTED_EXIT = False
            sys.exit( -1)
        commands = [args['<command>']]

    #atexit.register( exit_handler)

    if n_workers == None:
        error("Please set the number of parallel workers to use (with --workers or -n).")

    try:
        n_workers = int( n_workers)
        if n_workers < 0: raise ValueError
    except:
        error( "The --workers argument should be an positive, non-zero integer value.")


    parallel( commands, job_name=args['--job-name'],
        max_workers=n_workers, cpus=int(args['--cpus']), # memory=int(args['--memory']),
        partition=partition, time_limit=int(args['--time']),log_dir=log_dir, 
        wait_to_run=wait_to_run, linked_dir=linked_dir, max_jobs=max_jobs)

