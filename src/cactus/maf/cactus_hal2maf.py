#!/usr/bin/env python3

"""
This is a spiritual successor to hal2mafMP.py (from hal).  putting it here in cactus as cactus already has the setup tools and docker command stuff
needed toil autoscale, and it'll be easier to keep all dependencies managed in submodules.  
"""

import os, sys
from argparse import ArgumentParser
import copy
import timeit, time
import math

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.version import cactus_commit
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("halFile", help = "HAL file to convert to MAF")
    parser.add_argument("outputMAF", help = "Output MAF (will be gzipped if ends in .gz)")
    parser.add_argument("--batchSize", type=int, help = "Number of chunks for each hal2maf batch", default=None)
    parser.add_argument("--batchCount", type=int, help = "Number of hal2maf batches [default 1 unless --batchSize set]", default=None)
    parser.add_argument("--batchCores", type=int, help = "Number of cores for each hal2maf batch.")
    parser.add_argument("--chunkSize", type=int, help = "Size of chunks to operate on.", required=True)
    parser.add_argument("--raw", action="store_true", help = "Do not run taf-based normalization on the MAF")

    # pass through a subset of hal2maf options
    parser.add_argument("--refGenome", required=True,
                        help="name of reference genome (root if empty)",
                        default=None)
    parser.add_argument("--rootGenome",
                        help="name of root genome (none if empty)",
                        default=None)
    parser.add_argument("--targetGenomes",
                        help="comma-separated (no spaces) list of target "
                        "genomes (others are excluded) (vist all if empty)",
                        default=None)
    parser.add_argument("--maxRefGap",
                        help="maximum gap length in reference", type=int,
                        default=None)
    parser.add_argument("--noDupes",
                        help="ignore paralogy edges",
                        action="store_true",
                        default=False)
    parser.add_argument("--onlyOrthologs",
                        help="make only orthologs to the reference appear in the MAF blocls",
                        action="store_true",
                        default=False)
    parser.add_argument("--noAncestors",
                        help="don't write ancestral sequences. IMPORTANT: "
                        "Must be used in conjunction with --refGenome"
                        " to set a non-ancestral genome as the reference"
                        " because the default reference is the root.",
                        action="store_true",
                        default=False)
    
    # pass through taf_add_gap_bases options
    parser.add_argument("--gapFill",
                        help="use TAF tools to fill in reference gaps up to this length (currently more reliable than --maxRefGap) [default=50]",
                        type=int,
                        default=50)

    # pass through taf_norm options
    parser.add_argument("--maximumBlockLengthToMerge",
                        help="Only merge together any two adjacent blocks if one or both is less than this many bases long, [default=200]",
                        type=int,
                        default=200)
    parser.add_argument("--maximumGapLength",
                         help="Only merge together two adjacent blocks if the total number of unaligned bases between the blocks is less than this many bases, [default=3]",
                         type=int,
                         default=3)
    parser.add_argument("--fractionSharedRows",
                        help="The fraction of rows between two blocks that need to be shared for a merge, [default=0.6]",
                        type=float,
                        default=0.6)                            
    
    #Progressive Cactus Options
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    if options.batchSize and options.batchCount:
        raise RuntimeError('Only one of --batchSize and --batchCount can be specified')

    # apply cpu override                
    if options.batchCores is None:
        if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
            options.batchCores = cpu_count()
            logger.info('Setting batchCores to {}'.format(options.batchCores))
        else:
            raise RuntimeError('--batchCores must be specified for batch systems other than singleMachine')
    
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    if not options.batchCount and not options.batchSize:
        logger.info('Using default batch count of 1')
        options.batchCount = 1
                    
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            maf_id = toil.restart()
        else:
            logger.info("Importing {}".format(options.halFile))
            hal_id = toil.importFile(options.halFile)
            maf_id = toil.start(Job.wrapJobFn(hal2maf_workflow, hal_id, options))

        #export the maf
        toil.exportFile(maf_id, makeURL(options.outputMAF))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-hal2maf has finished after {} seconds".format(run_time))


def hal2maf_workflow(job, hal_id, options):

    hal2maf_ranges_job = job.addChildJobFn(hal2maf_ranges, hal_id, options, cores=1, disk=hal_id.size)    
    hal2maf_all_job = hal2maf_ranges_job.addFollowOnJobFn(hal2maf_all, hal_id, hal2maf_ranges_job.rv(), options)
    hal2maf_merge_job = hal2maf_all_job.addFollowOnJobFn(hal2maf_merge, hal2maf_all_job.rv(), options, disk=hal_id.size)
    return hal2maf_merge_job.rv()

def hal2maf_ranges(job, hal_id, options):
    """ get the ranges (in reference *sequence* coordinates) for each hal2maf job """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    RealtimeLogger.info("Computing range information")

    res = cactus_call(parameters=['halStats', hal_path, '--sequenceStats', options.refGenome], check_output=True)
    ref_sequence_stats = []
    for line in res.strip().split('\n')[1:]:
        tokens = line.strip().split(",")
        if len(tokens) == 4:
            ref_sequence_stats.append((tokens[0], int(tokens[1]), int(tokens[2]), int(tokens[3])))

    chunks = []
    for ref_stats in ref_sequence_stats:
        ref_name = ref_stats[0]
        ref_len = ref_stats[1]
        start = 0
        while start < ref_len:
            end = min(start + options.chunkSize, ref_len)
            chunks.append((ref_name, start, end))
            start = end

    assert chunks
    return chunks

def hal2maf_all(job, hal_id, chunks, options):
    """ make a job for each batch of chunks """
    num_batches = options.batchCount
    if not num_batches:        
        num_batches = int(len(chunks) / options.batchCores)
        remainder = len(chunks) % options.batchCores
        if remainder:
            num_batches += 1
        RealtimeLogger.info('Setting batchCount to {}'.format(num_batches))
            
    batch_size = options.batchSize
    if not batch_size:
        batch_size = math.ceil(len(chunks) / num_batches)
        RealtimeLogger.info('Setting batchSize to {}'.format(batch_size))
        
    chunks_left = len(chunks)
    batch_results = []        
    for i in range(num_batches):
        cur_chunk = i * batch_size
        cur_batch_size = min(chunks_left, batch_size)
        batch_results.append(job.addChildJobFn(hal2maf_batch, hal_id, chunks[cur_chunk:cur_chunk+cur_batch_size], options,
                                               disk=math.ceil((1 + 1 / num_batches)*hal_id.size), cores=options.batchCores).rv())
        chunks_left -= cur_batch_size
    assert chunks_left == 0
    
    return batch_results

def hal2maf_cmd(hal_path, chunk, chunk_num, options):
    """ make a hal2maf command for a chunk """

    # we are going to run this relative to work_dir
    hal_path = os.path.basename(hal_path)
    
    cmd = 'set -eo pipefail && hal2maf {} stdout --refGenome {} --refSequence {} --start {} --length {}'.format(hal_path, options.refGenome, chunk[0], chunk[1], chunk[2]-chunk[1])
    if options.rootGenome:
        cmd += ' --rootGenome {}'.format(options.rootGenome)
    if options.targetGenomes:
        cmd += ' --rootGenome {}'.format(options.targetGenomes)
    if options.maxRefGap:
        cmd += ' --maxRefGap {}'.format(options.maxRefGap)
    if options.noDupes:
        cmd += ' --noDupes'
    if options.onlyOrthologs:
        cmd += ' --onlyOrthologs'
    if options.noAncestors:
        cmd += ' --noAncestors'
    if not options.raw:
        # todo: we can parameterize these guys in the cactus config
        cmd += ' | maf_to_taf | taf_add_gap_bases -a {} -m {} | taf_norm -k -m {} -n {} -q {}'.format(hal_path, options.gapFill,
                                                                                                      options.maximumBlockLengthToMerge,
                                                                                                      options.maximumGapLength,
                                                                                                      options.fractionSharedRows)
    if chunk[1] != 0:
        cmd += ' | grep -v ^#'
    if options.outputMAF.endswith('.gz'):
        cmd += ' | bgzip'        
    cmd += ' > {}.maf'.format(chunk_num)
    return cmd

def hal2maf_batch(job, hal_id, batch_chunks, options):
    """ run hal2maf on a batch of chunks in parallel """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))
    job.fileStore.readGlobalFile(hal_id, hal_path)

    cmds = [hal2maf_cmd(hal_path, chunk, i, options) for i, chunk in enumerate(batch_chunks)]
    
    # do it with parallel
    cmd_path = os.path.join(work_dir, 'hal2maf_cmds.txt')
    with open(cmd_path, 'w') as cmd_file:
        for cmd in cmds:
            cmd_file.write(cmd + '\n')
    parallel_cmd = [['cat', cmd_path],
                    ['parallel', '-j', str(job.cores), '{}']]
    RealtimeLogger.info('First of {} commands in parallel batch: {}'.format(len(cmds), cmds[0]))
    cactus_call(parameters=parallel_cmd, work_dir=work_dir)

    # merge up the results
    maf_path = os.path.join(work_dir, 'out.maf')
    catFiles([os.path.join(work_dir, '{}.maf'.format(i)) for i in range(len(batch_chunks))], maf_path)

    return job.fileStore.writeGlobalFile(maf_path)
            

def hal2maf_merge(job, maf_ids, options):
    """ just cat the results """
    work_dir = job.fileStore.getLocalTempDir()
    merged_path = os.path.join(work_dir, 'merged.maf')
    in_paths = []
    for i, maf_id in enumerate(maf_ids):
        in_paths.append(job.fileStore.readGlobalFile(maf_id))
    catFiles(in_paths, merged_path)
    return job.fileStore.writeGlobalFile(merged_path)
