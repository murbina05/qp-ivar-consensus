# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
# from operator import length_hint
# from sys import stdlib_module_names
import pandas as pd

from os import environ
from os.path import basename, join
from glob import glob
from itertools import zip_longest

from qiita_client import ArtifactInfo

MEMORY = '16g'
WALLTIME = '30:00:00'
FINISH_MEMORY = '10g'
FINISH_WALLTIME = '10:00:00'
MAX_RUNNING = 8
#  -i    (Required) Sorted bam file,
# with aligned reads, to consensus primers and quality

#  -b    (Required) BED file with
# primer sequences and positions

#  -m    Minimum length of read to
#  retain after consensusming (Default: 30)

#  -q    Minimum quality threshold
#  for sliding window to pass (Default: 20)

#  -s    Width of sliding window (Default: 4)

#  -e    Include reads with no primers.
# By default, reads with no primers are excluded

QC_REFERENCE = environ["QC_REFERENCE"]
IVAR_consensus_BASE = 'ivar consensus -x 5 -e -i %s -b %s -p %s -m %s -q %s -s %s'

IVAR_consensus_CMD = ' '.join([IVAR_consensus_BASE, ' -o {out_dir}/%s -O {out_dir}/%s'])

# i dont think i need this part


def get_dbs_list():

    folder = QC_REFERENCE
    # skip human database
    list = [basename(f) for f in glob(f'{folder}/*.bed')]
    return list


# might need to add envrionment var, passes database
# however not need due to not using minimap2 :/

def _generate_commands(BAM_file, prefix, out_dir,
                       min_length=30, min_quality=20, slideing_window_width=4):
    """Helper function to generate commands and facilite testing"""
    files = BAM_file
    cmd = IVAR_consensus_CMD
#        if database is not None:
#            cmd = COMBINED_CMD
#    else:
#        cmd =_CMD_SINGLE
#        if database is not None:
#            cmd = COMBINED_CMD_SINGLE
    command = cmd.format(BAM_file=BAM_file, primer=QC_REFERENCE, prefix=prefix,
                         min_length=min_length, min_quality=min_quality,
                         slideing_window_width=slideing_window_width,
                         out_dir=out_dir)

    out_files = []
    commands = []
    for BAM_file in files:
        fname = basename(BAM_file)
        out_files.append((f'{out_dir}/{fname}',
                         'consensusmed'))  # might be consensusmed
#        if rev_fp:
#            rname = basename(rev_fp)
#            out_files.append((f'{out_dir}/{rname}',
#            'raw_reverse_seqs'))
#            cmd = command % (fwd_fp, rev_fp, fname, rname)
#        else:
        cmd = command % (BAM_file, fname)
        commands.append(cmd)

    return commands, out_files


def ivar_consensus(qclient, job_id, parameters, out_dir):
    """Run ivar consensus with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """

    qclient.update_job_step(
        job_id, "Step 3 of 4: Finishing Ivar consensus")

    ainfo = []
    # Generates 2 artifacts: one for the ribosomal
    # reads and other for the non-ribosomal reads
    out_files = []
    with open(f'{out_dir}/{job_id}.out_files.tsv') as f:
        for line in f.readlines():
            fp, ft = line.split()
            out_files.append((fp, ft))

    # Step 4 generating artifacts
    msg = "Step 4 of 4: Generating new artifact"
    qclient.update_job_step(job_id, msg)
    ainfo = [ArtifactInfo('consensusmed files', 'per_sample_FASTQ', out_files)]
#   ^^^^ looks like the  part might need to change ^^^^
#       might need to change to fastq for multisample pipeline
    return True, ainfo, ""


def ivar_consensus_to_array(files, out_dir, params, prep_info, url, job_id):
    """Creates qsub files for submission of bam

    Parameters
    ----------
    files : dict
        The dictionary of files to process, raw_forward_seqs/raw_reverse_seqs
    out_dir : str
        The output directory
    params : dict
        The parameter values to run Ivar consensus
    prep_info : str
        The path to prep_info
    url : str
        The url to send info to
    job_id : str
        The job id

    Returns
    -------
    str, str, str
        The paths of the main_qsub_fp, finish_qsub_fp, out_files_fp
    """
    database = None
    if params['primer'] != 'None':
        database = [join(QC_REFERENCE, f'{db}')
                    for db in get_dbs_list()
                    if params['primer'] in db][0]

    bam_reads = sorted(files['unconsensusmed_sorted_bam'])
#    if 'raw_reverse_seqs' in files:
#        rev_seqs = sorted(files['raw_reverse_seqs'])
#    else:
#        rev_seqs = []
    df = pd.read_csv(prep_info, sep='\t', dtype='str',
                     na_values=[], keep_default_na=True)
    df.set_index('sample_name', inplace=True)
    if 'run_prefix' not in df.columns:
        raise ValueError('Missing run_prefix column in your preparation')

    # Note that for processing we don't actually need the run_prefix so
    # we are not going to use it and simply loop over the ordered
    # fwd_seqs/rev_seqs
    commands, out_files = _generate_commands(
        bam_reads, database, params['threads'], out_dir)

    # writing the job array details
    details_name = join(out_dir, 'ivar_consensus.array-details')
    with open(details_name, 'w') as details:
        details.write('\n'.join(commands))
    n_jobs = len(commands)

    # all the setup pieces
    PPN = params['threads']
    lines = ['#!/bin/bash',
             '#PBS -M qiita.help@gmail.com',
             f'#PBS -N {job_id}',
             f'#PBS -l nodes=1:ppn={PPN}',
             f'#PBS -l walltime={WALLTIME}',
             f'#PBS -l mem={MEMORY}',
             f'#PBS -o {out_dir}/{job_id}' + '_${PBS_ARRAYID}.log',
             f'#PBS -e {out_dir}/{job_id}' + '_${PBS_ARRAYID}.err',
             f'#PBS -t 1-{n_jobs}%{MAX_RUNNING}',
             '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh',
             'set -e',
             f'cd {out_dir}',
             f'{params["environment"]}',
             'date',  # start time
             'hostname',  # executing system
             'echo ${PBS_JOBID} ${PBS_ARRAYID}',
             'offset=${PBS_ARRAYID}',
             'step=$(( $offset - 0 ))',
             f'cmd=$(head -n $step {details_name} | tail -n 1)',
             'eval $cmd',
             'set +e',
             'date']
    main_qsub_fp = join(out_dir, f'{job_id}.qsub')
    with open(main_qsub_fp, 'w') as job:
        job.write('\n'.join(lines))
        job.write('\n')

    # finish job
    lines = ['#!/bin/bash',
             '#PBS -M qiita.help@gmail.com',
             f'#PBS -N finish-{job_id}',
             '#PBS -l nodes=1:ppn=1',
             f'#PBS -l walltime={FINISH_WALLTIME}',
             f'#PBS -l mem={FINISH_MEMORY}',
             f'#PBS -o {out_dir}/finish-{job_id}.log',
             f'#PBS -e {out_dir}/finish-{job_id}.err',
             '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh',
             'set -e',
             f'cd {out_dir}',
             f'{params["environment"]}',
             'date',  # start time
             'hostname',  # executing system
             'echo $PBS_JOBID',
             f'finish_qp_ivar_consensus {url} {job_id} {out_dir}\n'
             "date"]
    finish_qsub_fp = join(out_dir, f'{job_id}.finish.qsub')
    with open(finish_qsub_fp, 'w') as out:
        out.write('\n'.join(lines))
        out.write('\n')

    out_files_fp = join(out_dir, f'{job_id}.out_files.tsv')
    with open(out_files_fp, 'w') as out:
        out.write('\n'.join([f'{fp}\t{ft}'for fp, ft in out_files]))

    return main_qsub_fp, finish_qsub_fp, out_files_fp
