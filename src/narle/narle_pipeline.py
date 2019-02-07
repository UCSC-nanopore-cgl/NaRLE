#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import os
import multiprocessing
import subprocess
import sys
import textwrap
import tarfile
from urlparse import urlparse
import math
import shutil
import glob
import time
import datetime
import logging
from contextlib import closing
from docker.errors import ContainerError
import pysam

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil import physicalMemory, physicalDisk
from toil.job import Job
from toil.lib.docker import apiDockerCall, dockerCall
from toil_lib.jobs import map_job
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.urls import download_url, s3am_upload


# input / output schemes
SCHEMES = ('http', 'https', 'file', 's3', 'ftp')

# filenames
DEFAULT_CONFIG_NAME = 'config-toil-nanopore.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-nanopore.tsv'

# docker images
DOCKER_SAMTOOLS_IMG = "quay.io/ucsc_cgl/samtools"
DOCKER_SAMTOOLS_TAG = "1.8--cba1ddbca3e1ab94813b58e40e56ab87a59f1997"

DOCKER_SAMTOOLS_SORT_IMG = "tpesout/samtools_sort"
DOCKER_SAMTOOLS_SORT_TAG = "latest"
DOCKER_SAMTOOLS_SORT_OUT = "samtools_sort.bam"
DOCKER_SAMTOOLS_SORT_LOG = "samtools_sort.log"

DOCKER_SAMTOOLS_VIEW_IMG = "tpesout/samtools_view"
DOCKER_SAMTOOLS_VIEW_TAG = "latest"
DOCKER_SAMTOOLS_VIEW_OUT = "samtools_view.out"
DOCKER_SAMTOOLS_VIEW_LOG = "samtools_view.log"

DOCKER_MINIMAP2_IMG = "tpesout/minmap2"
DOCKER_MINIMAP2_TAG = "latest"
DOCKER_MINIMAP2_OUT = "minimap2.sam"
DOCKER_MINIMAP2_LOG = "minimap2.log"

DOCKER_MARGIN_POLISH_IMG = "tpesout/margin_polish"
DOCKER_MARGIN_POLISH_TAG = "latest"
DOCKER_MARGIN_POLISH_LOG = "marginPolish.log"

DOCKER_CUSTOM_STATS_IMG = "tpesout/custom_assembly_stats"
DOCKER_CUSTOM_STATS_TAG = "latest"
DOCKER_CUSTOM_STATS_LOG = "custom_assembly_stats.log"

DOCKER_QUAST_IMG = "tpesout/quast"
DOCKER_QUAST_TAG = "latest"
DOCKER_QUAST_LOG = "quast.log"

DOCKER_RTG_IMG = "tpesout/rtg_tools"
DOCKER_RTG_TAG = "latest"
DOCKER_RTG_LOG = "rtg.log"

# resource
# DEFAULT_CPU = 8

# assembly
ASM_READS_MEM_FACTOR = 1
ASM_READS_DSK_FACTOR_kb2b = 1024

# alignment
ALN_READS_MEM_FACTOR = 1
ALN_READS_DSK_FACTOR_kb2b = 1024 * 2

# polishing
POLISH_READS_MEM_FACTOR = 1
POLISH_READS_DSK_FACTOR_kb2b = 1024

# resource estimation functions
estimate_aln_dsk_usage = lambda reads, config: int(min(config.maxDisk, int(reads) * ALN_READS_DSK_FACTOR_kb2b))
estimate_aln_mem_usage = lambda reads, config: int(min(config.maxDisk, int(reads) * ALN_READS_MEM_FACTOR))
estimate_polish_dsk_usage = lambda reads, config: int(min(config.maxDisk, int(reads) * POLISH_READS_DSK_FACTOR_kb2b))
estimate_polish_mem_usage = lambda reads, config: int(min(config.maxMemory, int(reads) * POLISH_READS_MEM_FACTOR))
estimate_asm_dsk_usage = lambda reads, config: int(min(config.maxDisk, int(reads) * ASM_READS_DSK_FACTOR_kb2b))
estimate_asm_mem_usage = lambda reads, config: int(min(config.maxMemory, int(reads) * ASM_READS_MEM_FACTOR))

# job data
JD_READS_FA = 'reads_fa'
JD_READS_FQ = 'reads_fq'
JD_TRUE_REF = 'true_ref'
JD_ASSEMBLY = 'assembly'
JD_ASSEMBLY_STATS = 'assembly_stats'
JD_ASSEMBLY_ALIGN = 'assembly_align'
JD_MARGIN_POLISH_ASSEMBLY = 'margin_polish_assembly'
JD_MARGIN_POLISH_ASSEMBLY_STATS = 'margin_polish_assembly_stats'
JD_MARGIN_POLISH_PARAMS = 'margin_polish_params'

# file types
FT_FASTA = ['fa', 'fasta']
FT_FASTQ = ['fq', 'fastq']
FT_GZ = ['gz']
is_type = lambda input, type_suffixes: any(map(
    lambda t_suffix: input.endswith(t_suffix) or
                     any(map(lambda tz_suffix: input.endswith("{}.{}".format(t_suffix, tz_suffix)), FT_GZ)),
    type_suffixes))
is_fasta = lambda input: is_type(input, FT_FASTA)
is_fastq = lambda input: is_type(input, FT_FASTQ)
is_zipped = lambda input: is_type(input, FT_GZ)


# for debugging
DEBUG = True
DOCKER_LOGGING = True



###########################################################
#                    WORKFLOW STEPS                       #
###########################################################


def prepare_input(job, sample, config):
    # time tracking
    start = time.time()

    # job prep
    config = argparse.Namespace(**vars(config))
    uuid, reads_url, assembly_url = sample
    config.uuid = uuid
    config.read_filetype = FT_FASTA[0] if is_fasta(reads_url) else (FT_FASTQ[0] if is_fastq(reads_url) else None)
    if is_zipped(reads_url):
        config.read_filetype = "{}.{}".format(config.read_filetype, FT_GZ[0])
    work_dir = job.tempDir

    # set up output dir
    config.output_dir = os.path.join(config.output_dir, uuid)
    if urlparse(config.output_dir).scheme != "s3":
        mkdir_p(config.output_dir)

    # set up intermediate dir
    if config.intermediate_file_location is not None:
        config.intermediate_file_location = os.path.join(config.intermediate_file_location, uuid)
        mkdir_p(config.intermediate_file_location)

    # loggit
    log(job, "START:{}".format(datetime.datetime.now()), uuid)
    log(job, "Preparing input with reads_url:{}, assembly_url:{}".format(reads_url, assembly_url),
        uuid, "prepare_input")

    # resource estimation
    config.maxCores = min(config.maxCores, multiprocessing.cpu_count())
    config.defaultCores = config.maxCores #min(DEFAULT_CPU, config.maxCores)
    config.maxMemory = min(config.maxMemory, int(physicalMemory() * .95))
    config.maxDisk = min(config.maxDisk, int(physicalDisk(None, toilWorkflowDir=work_dir) * .95))

    ### download files, save to jobstore ###

    # reads
    download_url(reads_url, work_dir=work_dir)
    reads_basename = os.path.basename(reads_url)
    reads_location = os.path.join(work_dir, "{}.original.{}".format(uuid, reads_basename))
    os.rename(os.path.join(work_dir, reads_basename), reads_location)
    reads_fileid = job.fileStore.writeGlobalFile(reads_location)
    reads_size = os.stat(reads_location).st_size
    reads_fa_fileid = None
    reads_fq_fileid = None
    if is_fasta(reads_basename):
        reads_fa_fileid = reads_fileid
        if config.intermediate_file_location is not None:
            copy_files(file_paths=[reads_location], output_dir=config.intermediate_file_location)
        raise UserError("FASTA filetype is not supported for input reads (yet)")
    else:
        reads_fq_fileid = reads_fileid
        #TODO convert to fa

    # true reference
    true_reference_fileid = None
    true_reference_size = 0
    if config.true_reference is not None:
        download_url(config.true_reference, work_dir=work_dir)
        true_reference_basename = os.path.basename(config.true_reference)
        true_reference_location = os.path.join(work_dir, "{}.true_reference.{}".format(uuid, true_reference_basename))
        os.rename(os.path.join(work_dir, true_reference_basename), true_reference_location)
        true_reference_fileid = job.fileStore.writeGlobalFile(true_reference_location)
        true_reference_size = os.stat(true_reference_location).st_size
        if config.intermediate_file_location is not None:
            copy_files(file_paths=[true_reference_location], output_dir=config.intermediate_file_location)

    # assembly
    assembly_fileid = None
    if assembly_url is not None:
        download_url(assembly_url, work_dir=work_dir)
        assembly_basename = os.path.basename(assembly_url)
        assembly_location = os.path.join(work_dir, "{}.assembly.{}".format(uuid, assembly_basename))
        os.rename(os.path.join(work_dir, assembly_basename), assembly_location)
        assembly_fileid = job.fileStore.writeGlobalFile(assembly_location)
        if config.intermediate_file_location is not None:
            copy_files(file_paths=[assembly_location], output_dir=config.intermediate_file_location)

    # params
    download_url(config.margin_params, work_dir=work_dir)
    params_basename = os.path.basename(config.margin_params)
    params_location = os.path.join(work_dir, params_basename)
    params_fileid = job.fileStore.writeGlobalFile(params_location)
    if config.intermediate_file_location is not None:
        copy_files(file_paths=[params_location], output_dir=config.intermediate_file_location)

    # job data
    job_data = {
        JD_READS_FA: reads_fa_fileid,
        JD_READS_FQ: reads_fq_fileid,
        JD_TRUE_REF: true_reference_fileid,
        JD_ASSEMBLY: assembly_fileid,
        JD_ASSEMBLY_ALIGN: None,
        JD_MARGIN_POLISH_ASSEMBLY: None,
        JD_MARGIN_POLISH_PARAMS: params_fileid,
    }

    # for final consolidation
    return_values = list()

    # do we need to enqueue assembly job?
    assembly_job = job
    assembly_job_data = job_data
    if assembly_url is None:
        # alignment resource usage
        asm_mem = estimate_asm_mem_usage(reads_size, config)
        asm_dsk = estimate_asm_dsk_usage(reads_size, config)

        # add jobs
        assembly_job = job.addChildJobFn(generate_shasta_assembly, config, job_data,
                                         cores=config.defaultCores, memory=asm_mem, disk=asm_dsk)
        assembly_job_data = assembly_job.rv()

    # enqueue alignment job
    aln_mem = estimate_aln_mem_usage(reads_size, config)
    aln_dsk = estimate_aln_dsk_usage(reads_size, config)
    assembly_align_job = assembly_job.addChildJobFn(align, config, assembly_job_data, JD_ASSEMBLY, JD_ASSEMBLY_ALIGN,
                                       cores=config.defaultCores, memory=aln_mem, disk=aln_dsk)

    # enqueue polish job
    pol_mem = estimate_aln_mem_usage(reads_size, config)
    pol_dsk = estimate_polish_mem_usage(reads_size, config)
    margin_polish_job = assembly_align_job.addChildJobFn(run_margin_polish, config, assembly_align_job.rv(),
                                           cores=config.defaultCores, memory=pol_mem, disk=pol_dsk)
    return_values.append(margin_polish_job.rv())

    # assembly qc
    if config.true_reference is not None:
        # add for assembly
        assembly_stats_job = assembly_job.addChildJobFn(gather_assembly_stats, config, assembly_job_data,
                                                        JD_ASSEMBLY, JD_ASSEMBLY_STATS, cores=1)
        return_values.append(assembly_stats_job.rv())
        # add for polishing
        margin_polish_stats_job = margin_polish_job.addChildJobFn(gather_assembly_stats, config, margin_polish_job.rv(),
                                                  JD_MARGIN_POLISH_ASSEMBLY, JD_MARGIN_POLISH_ASSEMBLY_STATS, cores=1)
        return_values.append(margin_polish_stats_job.rv())

    # submit final consolidation job
    consolidate_output_job = job.addFollowOnJobFn(consolidate_output, config, return_values)

    # log
    log_generic_job_debug(job, config.uuid, 'prepare_input', work_dir=work_dir)
    log_time(job, "prepare_input", start, config.uuid)

    return job_data


def align(job, config, job_data, ref_in_descriptor, align_out_descriptor):
    # prep
    start = time.time()
    uuid = config.uuid
    job_data = consolidate_job_data(job, config, job_data, 'align')
    read_filetype = FT_FASTQ[0] if job_data[JD_READS_FQ] is not None else FT_FASTA[0]
    read_fileid = job_data[JD_READS_FQ] if job_data[JD_READS_FQ] is not None else job_data[JD_READS_FA]
    ref_fileid = job_data[ref_in_descriptor]
    work_dir = job.tempDir
    log(job, "job_data: {}, {}".format(ref_in_descriptor, job_data), uuid, "align")

    # download files
    input_reads_filename = "{}.{}.{}".format(uuid, ref_in_descriptor, read_filetype)
    input_reads_location = os.path.join(work_dir, input_reads_filename)
    job.fileStore.readGlobalFile(read_fileid, input_reads_location)
    input_ref_filename = "{}.{}.reference.fa".format(uuid, ref_in_descriptor)
    input_ref_location = os.path.join(work_dir, input_ref_filename)
    job.fileStore.readGlobalFile(ref_fileid, input_ref_location)

    # do alignment
    alignment_filename = "{}.{}.unsorted.sam".format(uuid, ref_in_descriptor)
    alignment_location = os.path.join(work_dir, alignment_filename)
    docker_minimap2(job, config, work_dir, input_reads_filename, input_ref_filename, output_filename=alignment_filename)
    if not os.path.isfile(alignment_location):
        raise UserError("Alignment file {} not created for {}".format(alignment_filename, uuid))
    if config.intermediate_file_location is not None:
        log_location = os.path.join(work_dir, "minimap2.{}.log".format(input_reads_filename))
        os.rename(os.path.join(work_dir, DOCKER_MINIMAP2_LOG), log_location)
        copy_files(file_paths=[alignment_location, log_location], output_dir=config.intermediate_file_location)

    # sort it
    sorted_aln_filename = "{}.{}.bam".format(uuid, ref_in_descriptor)
    sorted_aln_location = os.path.join(work_dir, sorted_aln_filename)
    docker_samtools_sort(job, config, work_dir, alignment_filename, output_filename=sorted_aln_filename)
    if not os.path.isfile(sorted_aln_location):
        raise UserError("Sorted alignment file {} not created for {}".format(sorted_aln_filename, uuid))

    # save it
    if config.intermediate_file_location is not None:
        log_location = os.path.join(work_dir, "samtools_sort.{}.log".format(alignment_filename))
        os.rename(os.path.join(work_dir, DOCKER_SAMTOOLS_SORT_LOG), log_location)
        copy_files(file_paths=[sorted_aln_location, log_location], output_dir=config.intermediate_file_location)
    move_or_upload(job, config, [sorted_aln_location])

    # save to jobstore
    sorted_aln_fileid = job.fileStore.writeGlobalFile(sorted_aln_location)
    job_data[align_out_descriptor] = sorted_aln_fileid

    # log
    log_generic_job_debug(job, uuid, 'align', work_dir=work_dir)
    log_time(job, "align", start, uuid)
    return job_data


def gather_assembly_stats(job, config, job_data, assembly_in_descriptor, stats_out_descriptor):
    # prep
    start = time.time()
    uuid = config.uuid
    job_data = consolidate_job_data(job, config, job_data, 'gather_assembly_stats')
    work_dir = job.tempDir

    # get files
    assembly_filename = "{}.assembly.fa".format(uuid)
    assembly_location = os.path.join(work_dir, assembly_filename)
    job.fileStore.readGlobalFile(job_data[assembly_in_descriptor], assembly_location)
    true_ref_filename = "{}.true_ref.fa".format(uuid)
    true_ref_location = os.path.join(work_dir, true_ref_filename)
    job.fileStore.readGlobalFile(job_data[JD_TRUE_REF], true_ref_location)

    # get stats
    log(job, "In gather_assembly_stats: {}->{}, {}".format(assembly_in_descriptor, stats_out_descriptor, job_data),
        config.uuid, 'gather_assembly_stats')
    custom_stats_files = docker_custom_assembly_stats(job, config, work_dir, assembly_filename, true_ref_filename,
                                                      "{}.{}.custom_assembly_stats".format(uuid, assembly_in_descriptor))
    quast_files = docker_quast(job, config, work_dir, assembly_filename, true_ref_filename,
                              "{}.{}.quast".format(uuid, assembly_in_descriptor),
                              None if config.quast_params is None else config.quast_params.split())

    # tarball result
    files_to_tar = []
    files_to_tar.extend(custom_stats_files)
    files_to_tar.extend(quast_files)
    output_tarball = "{}.assembly_stats.{}.tar.gz".format(uuid, assembly_in_descriptor)
    tarball_files(output_tarball, list(map(lambda x: os.path.abspath(os.path.join(work_dir, x)), files_to_tar)), work_dir)
    if not os.path.isfile(os.path.join(work_dir, output_tarball)):
        err = "Failed to create tarball of assembly_stats {}".format(output_tarball)
        log(job, "{}.  Desired files: {}".format(err, files_to_tar), config.uuid, 'gather_assembly_stats')
        raise UserError(err)

    # handle file
    stats_tarball = os.path.join(work_dir, output_tarball)
    if config.intermediate_file_location is not None:
        cas_log_location = os.path.join(work_dir, "custom_assembly_stats.{}.log".format(assembly_in_descriptor))
        os.rename(os.path.join(work_dir, DOCKER_CUSTOM_STATS_LOG), cas_log_location)
        quast_log_location = os.path.join(work_dir, "quast.{}.log".format(assembly_in_descriptor))
        os.rename(os.path.join(work_dir, DOCKER_QUAST_LOG), quast_log_location)
        copy_files(file_paths=[stats_tarball, cas_log_location, quast_log_location],
                   output_dir=config.intermediate_file_location)
    move_or_upload(job, config, [stats_tarball])
    stats_tarball_fileid = job.fileStore.writeGlobalFile(stats_tarball)
    job_data[stats_out_descriptor] = stats_tarball_fileid

    # log
    log_generic_job_debug(job, uuid, 'gather_assembly_stats', work_dir=work_dir)
    log_time(job, "gather_assembly_stats", start, uuid)
    return job_data


def run_margin_polish(job, config, job_data):
    # prep
    start = time.time()
    uuid = config.uuid
    job_data = consolidate_job_data(job, config, job_data, 'run_margin_polish')
    work_dir = job.tempDir

    # get files
    assembly_filename = "{}.marginPolish.fa".format(uuid)
    assembly_location = os.path.join(work_dir, assembly_filename)
    job.fileStore.readGlobalFile(job_data[JD_ASSEMBLY], assembly_location)
    alignment_filename = "{}.marginPolish.bam".format(uuid)
    alignment_location = os.path.join(work_dir, alignment_filename)
    job.fileStore.readGlobalFile(job_data[JD_ASSEMBLY_ALIGN], alignment_location)
    params_filename = "{}.marginPolish.params.json".format(uuid)
    params_location = os.path.join(work_dir, params_filename)
    job.fileStore.readGlobalFile(job_data[JD_MARGIN_POLISH_PARAMS], params_location)

    # run margin polish
    log(job, "Running run_margin_polish: {}".format(job_data), config.uuid, 'run_margin_polish')
    polished_assembly_filename = docker_marginPolish(job, config, work_dir,
                                                     alignment_filename, assembly_filename, params_filename,
                                                     output_base="{}.marginPolish.out".format(uuid))
    polished_assembly_location = os.path.join(work_dir, polished_assembly_filename)
    if not os.path.isfile(polished_assembly_location):
        raise UserError("Polished assembly file {} not created for {}".format(polished_assembly_location, uuid))

    # handle file
    if config.intermediate_file_location is not None:
        log_location = os.path.join(work_dir, "marginPolish.log".format(alignment_filename))
        os.rename(os.path.join(work_dir, DOCKER_MARGIN_POLISH_LOG), log_location)
        copy_files(file_paths=[polished_assembly_location, log_location], output_dir=config.intermediate_file_location)
    move_or_upload(job, config, [polished_assembly_location])

    # save to jobstore
    polished_assembly_fileid = job.fileStore.writeGlobalFile(polished_assembly_location)
    job_data[JD_MARGIN_POLISH_ASSEMBLY] = polished_assembly_fileid

    # log
    log_generic_job_debug(job, uuid, 'run_margin_polish', work_dir=work_dir)
    log_time(job, "run_margin_polish", start, uuid)
    return job_data


def generate_shasta_assembly(job, config, job_data):
    # prep
    start = time.time()
    uuid = config.uuid
    job_data = consolidate_job_data(job, config, job_data, 'generate_shasta_assembly')
    work_dir = job.tempDir

    log(job, "In generate_shasta_assembly: {}".format(job_data), uuid, 'generate_shasta_assembly')

    # log
    log_generic_job_debug(job, uuid, 'generate_shasta_assembly', work_dir=work_dir)
    log_time(job, "generate_shasta_assembly", start, uuid)

    #TODO
    raise UserError("generate_shasta_assembly is not implemented!")
    return job_data


def consolidate_output(job, config, job_data_list):
    # prep
    start = time.time()
    uuid = config.uuid
    final_job_data = consolidate_job_data(job, config, job_data_list, "consolidate_output")

    log(job, "In consolidate_output: {}".format(final_job_data), uuid, 'consolidate_output')

    return final_job_data


###########################################################
#                 JOB UTILITY FUNCTIONS                   #
###########################################################


def consolidate_job_data(job, config, job_data_list, function_name):
    # prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.tempDir

    # consolodate all data
    final_job_data = dict()
    def update_final_job_data(job_data, lvl_string):
        job_data_idx = 0
        log(job, "JOB_DATA_{}.{}: {}".format(lvl_string, job_data_idx, job_data), uuid, function_name)
        if type(job_data) == list:
            for jd in job_data:
                update_final_job_data(jd, "{}.{}".format(lvl_string, job_data_idx))
                job_data_idx += 1
        elif type(job_data) == dict:
            for jdk in job_data.keys():
                jdv = job_data[jdk]
                if jdv is None: continue
                if jdk in final_job_data and final_job_data[jdk] != jdv:
                    log(job, "Found different non-null values during data consolidation: {} -> [{},{}]".format(
                        jdk, jdv, final_job_data[jdk]), uuid, function_name)
                final_job_data[jdk] = jdv
    update_final_job_data(job_data_list, "0")

    # loggit
    log(job, "Consolidated job data: {}".format(final_job_data), uuid, function_name)

    return final_job_data


def docker_samtools_sort(job, config, work_dir, input_filename, output_filename=None, extra_args=None):
    # prep
    data_location = os.path.join("/data", input_filename)
    args = [data_location, "-@", str(job.cores)]
    if extra_args is not None:
        duplicated_args = list(filter(lambda x: x in args, extra_args))
        if len(duplicated_args) > 0:
            log(job, "Duplicated args in call to samtools_sort: {}".format(duplicated_args), config.uuid,
                'samtools_sort')
        args.extend(extra_args)

    # call
    docker_call(job, config, work_dir, args, DOCKER_SAMTOOLS_SORT_IMG, DOCKER_SAMTOOLS_SORT_TAG)

    # loggit
    log_file = os.path.join(work_dir, DOCKER_SAMTOOLS_SORT_LOG)
    log_debug_from_docker(job, log_file, config.uuid, 'samtools_sort',
                          input_file_locations=[os.path.join(work_dir, input_filename)])

    # sanity check
    require_docker_file_output(job, config, work_dir, [DOCKER_SAMTOOLS_SORT_OUT], 'samtools_sort',
                               DOCKER_SAMTOOLS_SORT_LOG)

    # rename file to appropriate output
    if output_filename is not None:
        os.rename(os.path.join(work_dir, DOCKER_SAMTOOLS_SORT_OUT), os.path.join(work_dir, output_filename))
        return output_filename
    else:
        return DOCKER_SAMTOOLS_SORT_OUT


def docker_samtools_view(job, config, work_dir, input_filename, selection_str=None, output_filename=None):
    # prep
    data_location = os.path.join("/data", input_filename)
    args = ['-hb', data_location]
    if selection_str is not None:
        args.append(selection_str)

    # call
    docker_call(job, config, work_dir, args, DOCKER_SAMTOOLS_VIEW_IMG, DOCKER_SAMTOOLS_VIEW_TAG)

    # sanity check
    require_docker_file_output(job, config, work_dir, [DOCKER_SAMTOOLS_VIEW_OUT], 'samtools_view',
                               DOCKER_SAMTOOLS_VIEW_LOG)

    # rename file to appropriate output
    if output_filename is not None:
        os.rename(os.path.join(work_dir, DOCKER_SAMTOOLS_VIEW_OUT), os.path.join(work_dir, output_filename))
        return output_filename
    else:
        return DOCKER_SAMTOOLS_VIEW_OUT


def docker_samtools_merge(job, config, work_dir, sam_input_filenames, sam_output_filename):
    # prep
    data_input_names = [os.path.join("/data", sam) for sam in sam_input_filenames]
    data_output_name = os.path.join("/data", sam_output_filename)
    args = ['merge', '--output-fmt', 'SAM', '-@', str(job.cores), data_output_name]
    args.extend(data_input_names)

    # call
    docker_call(job, config, work_dir, args, DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG)

    # sanity check
    require_docker_file_output(job, config, work_dir, [sam_output_filename], 'samtools_merge')


def docker_minimap2(job, config, work_dir, input_filename, ref_filename, output_filename=None, extra_args=None):
    # prep
    data_location = os.path.join("/data", input_filename)
    ref_location = os.path.join("/data", ref_filename)
    args = ["-t", str(job.cores), '-a']
    if extra_args is not None:
        duplicated_args = list(filter(lambda x: x in args, extra_args))
        if len(duplicated_args) > 0:
            log(job, "Duplicated args in call to minimap: {}".format(duplicated_args), config.uuid, 'minimap2')
        args.extend(extra_args)
    args.append(ref_location)
    args.append(data_location)

    # call
    docker_call(job, config, work_dir, args, DOCKER_MINIMAP2_IMG, DOCKER_MINIMAP2_TAG)

    # loggit
    log_debug_from_docker(job, os.path.join(work_dir, DOCKER_MINIMAP2_LOG), config.uuid, 'minimap2',
                          input_file_locations=[os.path.join(work_dir, input_filename),
                                                os.path.join(work_dir, ref_filename)])

    # sanity check
    minimap_out = os.path.join(work_dir, DOCKER_MINIMAP2_OUT)
    require_docker_file_output(job, config, work_dir, [DOCKER_MINIMAP2_OUT], 'minimap2', DOCKER_MINIMAP2_LOG)
    with open(minimap_out) as hopefully_sam_file:
        top_lines = []
        for line in hopefully_sam_file:
            if line.startswith("@"): break
            top_lines.append(line)
            if len(top_lines) > 3: break
        if len(top_lines) != 0:
            raise UserError("Output file {} for minimap2 does not appear to be in SAM format: {}".format(
                DOCKER_MINIMAP2_OUT, list(map(lambda x: x.strip(), top_lines))))

    # rename file to appropriate output
    if output_filename is not None:
        os.rename(minimap_out, os.path.join(work_dir, output_filename))
        return output_filename
    else:
        return DOCKER_MINIMAP2_OUT


def docker_index_bam(job, config, work_dir, bam_filename):
    # prep
    data_bam_location = os.path.join("/data", bam_filename)
    docker_params = ["index", data_bam_location]

    # call
    docker_call(job, config, work_dir, docker_params, DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG)

    # sanity check
    if not os.path.isfile(os.path.join(work_dir, bam_filename + ".bai")):
        raise UserError("BAM index file not created for {}".format(bam_filename))


def docker_get_contig_names(job, config, work_dir, bam_filename):
    # prep
    data_bam_location = os.path.join("/data", bam_filename)
    samtools_params = ["view", '-H', data_bam_location]

    # call
    bam_header = docker_call(job, config, work_dir, samtools_params, DOCKER_SAMTOOLS_IMG, DOCKER_SAMTOOLS_TAG,
                             stdout=True)
    bam_header_lines = bam_header.split("\n")
    contigs = list(set(map(lambda x: x.strip().split("\t")[1].replace("SN:",''),
                       filter(lambda x: x.startswith("@SQ"), bam_header_lines))))

    # sanity check
    if len(contigs) == 0:
        present_tags = set(map(lambda x: x.split("\t")[0], bam_header_lines))
        log(job, "Could not find contigs in {}, with {} header lines and present tags {}".format(
            bam_filename, len(bam_header_lines), present_tags), config.uuid, 'docker_get_contig_names')
        raise UserError("Could not find contigs in {}".format(bam_filename))

    return contigs


def docker_merge_vcfs(job, config, work_dir, vcf_input_names, vcf_output_name):
    # zip and tabix
    for in_vcf in vcf_input_names:
        docker_call(job, config, work_dir, ['bgzip', os.path.join("/data", in_vcf)], DOCKER_RTG_IMG, DOCKER_RTG_TAG)
        docker_call(job, config, work_dir, ['index', os.path.join("/data", "{}.gz".format(in_vcf))], DOCKER_RTG_IMG, DOCKER_RTG_TAG)

    # prep
    data_input_names = [os.path.join("/data", "{}.gz".format(vcf)) for vcf in vcf_input_names]
    data_output_name = os.path.join("/data", vcf_output_name)
    docker_params = ["vcfmerge", '-o', data_output_name]
    docker_params.extend(data_input_names)

    # call
    docker_call(job, config, work_dir, docker_params, DOCKER_RTG_IMG, DOCKER_RTG_TAG)

    # loggit
    log_debug_from_docker(job, os.path.join(work_dir, DOCKER_RTG_LOG), config.uuid, 'rtg_merge_vcfs',
                          input_file_locations=[os.path.join(work_dir, "{}.gz".format(vcf)) for vcf in vcf_input_names])

    # sanity check
    require_docker_file_output(job, config, work_dir, [vcf_output_name], 'rtg_merge_vcfs', DOCKER_RTG_LOG)


def docker_marginPolish(job, config, work_dir, input_bam_filename, input_ref_filename, input_params_filename,
                        output_base=None, extra_args=None):
    # prep
    bam_location = os.path.join("/data", input_bam_filename)
    ref_location = os.path.join("/data", input_ref_filename)
    params_location = os.path.join("/data", input_params_filename)
    output_base_filename = "out" if output_base is None else output_base
    output_base_location = os.path.join("/data", output_base_filename)
    args = [bam_location, ref_location, params_location, "-t", str(job.cores), '-o', output_base_location]
    if extra_args is not None:
        duplicated_args = list(filter(lambda x: x in args, extra_args))
        if len(duplicated_args) > 0:
            log(job, "Duplicated args in call to marginPolish: {}".format(duplicated_args), config.uuid, 'marginPolish')
        args.extend(extra_args)

    # verify index
    if not os.path.isfile(os.path.join(work_dir, "{}.bai".format(input_bam_filename))):
        docker_index_bam(job, config, work_dir, input_bam_filename)

    # call
    docker_call(job, config, work_dir, args, DOCKER_MARGIN_POLISH_IMG, DOCKER_MARGIN_POLISH_TAG)

    # loggit
    output_filename = "{}.fa".format(output_base_filename)
    log_debug_from_docker(job, os.path.join(work_dir, DOCKER_MARGIN_POLISH_LOG), config.uuid, 'marginPolish',
                          input_file_locations=[os.path.join(work_dir, input_bam_filename),
                                                os.path.join(work_dir, input_ref_filename)],
                          output_file_locations=[os.path.join(work_dir, output_filename)])

    # sanity check
    require_docker_file_output(job, config, work_dir, [output_filename], 'marginPolish', DOCKER_MARGIN_POLISH_LOG)

    # return
    return output_filename


def docker_custom_assembly_stats(job, config, work_dir, assembly_filename, true_ref_filename,
                        output_base=None, extra_args=None):
    # prep
    assembly_location = os.path.join("/data", assembly_filename)
    true_ref_location = os.path.join("/data", true_ref_filename)
    output_base_filename = "out" if output_base is None else output_base
    output_base_location = os.path.join("/data", output_base_filename)

    args = ['--sequences', assembly_location, '--ref', true_ref_location, '--output_dir', output_base_location]
    if extra_args is not None:
        duplicated_args = list(filter(lambda x: x in args, extra_args))
        if len(duplicated_args) > 0:
            log(job, "Duplicated args in call to custom_assembly_stats: {}".format(duplicated_args), config.uuid,
                'custom_assembly_stats')
        args.extend(extra_args)

    # call
    docker_call(job, config, work_dir, args, DOCKER_CUSTOM_STATS_IMG, DOCKER_CUSTOM_STATS_TAG)

    # loggit
    filename_no_ext = lambda x: '.'.join(os.path.basename(x).split('.')[:-1])
    output_file_base = "{}_VS_{}".format(filename_no_ext(assembly_location).replace(".", "_"),
                                            filename_no_ext(true_ref_filename).replace(".", "_"),)
    expected_file_suffixes = ['sam', 'sorted.sam', 'sorted.bam', 'sorted.bam.bai', 'sorted.png']
    output_files = list(map(lambda x: "{}.{}".format(output_file_base, x), expected_file_suffixes))
    output_files.append("summary_{}_VS_{}.sorted.csv".format(filename_no_ext(assembly_location).replace(".", "_"),
                                                             filename_no_ext(true_ref_filename).replace(".", "_")))
    output_files = list(map(lambda x: os.path.join(output_base_filename, x), output_files))
    log_debug_from_docker(job, os.path.join(work_dir, DOCKER_CUSTOM_STATS_LOG), config.uuid, 'custom_assembly_stats',
                          input_file_locations=[os.path.join(work_dir, assembly_filename),
                                                os.path.join(work_dir, true_ref_filename)],
                          output_file_locations=list(map(lambda x: os.path.join(work_dir, x), output_files)))

    # sanity check
    require_docker_file_output(job, config, work_dir, output_files,
                               'custom_assembly_stats', DOCKER_CUSTOM_STATS_LOG)

    # remove sam files
    for file in filter(lambda x: x.endswith("sam"), output_files):
        os.remove(os.path.join(work_dir,file))

    # return the base output folder
    return [output_base_filename]


def docker_quast(job, config, work_dir, assembly_filename, true_ref_filename,
                        output_base=None, extra_args=None):
    # prep
    assembly_location = os.path.join("/data", assembly_filename)
    true_ref_location = os.path.join("/data", true_ref_filename)
    output_base_filename = "out" if output_base is None else output_base
    output_base_location = os.path.join("/data", output_base_filename)

    args = ["-o", output_base_location, '-r', true_ref_location]
    if extra_args is not None:
        duplicated_args = list(filter(lambda x: x in args, extra_args))
        if len(duplicated_args) > 0:
            log(job, "Duplicated args in call to quast: {}".format(duplicated_args), config.uuid,
                'quast')
        args.extend(extra_args)
    args.append(assembly_location)

    # call
    docker_call(job, config, work_dir, args, DOCKER_QUAST_IMG, DOCKER_QUAST_TAG)

    # loggit
    expected_files = ['report.pdf', 'report.html'] #todo add more?
    output_files = list(map(lambda x: os.path.join(output_base, x), expected_files))
    log_debug_from_docker(job, os.path.join(work_dir, DOCKER_QUAST_LOG), config.uuid, 'quast',
                          input_file_locations=[os.path.join(work_dir, assembly_filename),
                                                os.path.join(work_dir, true_ref_filename)],
                          output_file_locations=list(map(lambda x: os.path.join(work_dir, x), output_files)))

    # sanity check
    require_docker_file_output(job, config, work_dir, output_files,
                               'quast', DOCKER_QUAST_LOG)

    # return
    return [output_base]


###########################################################
#                 TOIL UTILITY FUNCTIONS                  #
###########################################################


def move_or_upload(job, config, files, uuid=None, function="move_or_upload"):
    log(job, "Attempting to save {} to {}".format(files, config.output_dir, uuid, function))
    if urlparse(config.output_dir).scheme == 's3':
        for f in files:
            s3am_upload(fpath=f, s3_dir=config.output_dir)
    elif urlparse(config.output_dir).scheme != 's3':
        copy_files(file_paths=files, output_dir=config.output_dir)


def docker_call(job, config, work_dir, params, image, tag, detach=False, stdout=None, stderr=None):
    tagged_image = "{}:{}".format(image, tag)
    if DOCKER_LOGGING:
        log(job, "Running '{}' with parameters: {}".format(tagged_image, params), config.uuid, 'docker')
    return apiDockerCall(job, tagged_image, working_dir=work_dir, parameters=params, user="root", detach=detach,
                         stdout=stdout, stderr=stderr)


def require_docker_file_output(job, config, work_dir, output_filenames, function_id, log_filename=None,
                               max_directory_contents=None, max_log_lines=None):
    missing_filenames = list(filter(lambda x: not os.path.isfile(os.path.join(work_dir, x))
                                               or os.stat(os.path.join(work_dir, x)).st_size == 0, output_filenames))
    if len(missing_filenames) > 0:
        # document missing
        log(job, "Missing files after docker call: ", config.uuid, function_id)
        for missing in missing_filenames:
            log(job, "\t{}{}".format(missing, " (empty)" if os.path.isfile(os.path.join(work_dir, missing)) else ""),
                config.uuid, function_id)

        # document contents
        directory_contents = os.listdir(work_dir)
        if max_directory_contents is not None and len(directory_contents) > max_directory_contents:
            directory_contents = directory_contents[0:max_directory_contents]
            directory_contents.append("[{} items total]".format(len(directory_contents)))
        log(job, "Current files in work_dir: {}".format(work_dir), config.uuid, function_id)
        for existing in directory_contents:
            log(job, "\t{}".format(existing), config.uuid, function_id)

        # document log
        if log_filename is not None:
            log_location = os.path.join(work_dir, log_filename)
            if os.path.isfile(log_location):
                log(job, "Log file contents: {}".format(log_filename), config.uuid, function_id)
                log_lines = 0
                with open(log_location) as log_stream:
                    for ll in log_stream:
                        if max_log_lines is None or log_lines < max_log_lines:
                            log(job, "\t{}".format(ll.rstrip()), config.uuid, function_id)
                        log_lines += 1
                if max_log_lines is not None and log_lines <= max_log_lines:
                    log(job, "\t[{} lines total]".format(log_lines), config.uuid, function_id)
        else:
            log(job, "Log file {} was not found".format(log_filename), config.uuid, function_id)

        # die
        raise UserError("Missing or empty files after running {} on {}: {}".format(function_id, config.uuid,
                                                                                   missing_filenames))


def log_debug_from_docker(job, log_file_location, identifier, function,
                          input_file_locations=None, output_file_locations=None):
    # sanity check
    if not os.path.isfile(log_file_location):
        log(job, "Logfile missing!", identifier, function)
        return

    # file size logging
    if input_file_locations is not None:
        if type(input_file_locations) == str: input_file_locations = [input_file_locations]
        for input_file_location in input_file_locations:
            log(job, "DEBUG_INPUT_FILESIZE:{}:{}".format(
                os.stat(input_file_location).st_size, os.path.basename(input_file_location)),
                identifier, function)
    if output_file_locations is not None:
        if type(output_file_locations) == str: output_file_locations = [output_file_locations]
        for output_file_location in output_file_locations:
            if os.path.isfile(output_file_location):
                log(job, "DEBUG_OUTPUT_FILESIZE:{}:{}".format(
                    os.stat(output_file_location).st_size, os.path.basename(output_file_location)),
                    identifier, function)
            else:
                log(job, "DEBUG_OUTPUT_FILESIZE:-1:{}".format(os.path.basename(output_file_location)),
                    identifier, function)

    # any logging from logfile
    with open(log_file_location) as log_file:
        for line in log_file:
            line = line.strip()
            if line.startswith("DEBUG"):
                log(job, line, identifier, function)


def log_generic_job_debug(job, identifier, function, work_dir=None):

    # job resource logging
    log(job, "DEBUG_JOB_CORES:{}".format(job.cores), identifier, function)
    log(job, "DEBUG_JOB_MEM:{}".format(job.memory), identifier, function)
    log(job, "DEBUG_JOB_DISK:{}".format(job.disk), identifier, function)

    # workdir logging
    if work_dir is not None:
        workdir_abspath = os.path.abspath(work_dir)
        try:
            du_line = subprocess.check_output(['du', '-s', workdir_abspath])
            directory_filesize = du_line.split()[0]
            log(job, "DEBUG_WORKDIR_FILESIZE:{}:{}".format(directory_filesize, workdir_abspath),
                identifier, function)
        except Exception, e:
            log(job, "Exception ({}) finding size of directory {}: {}".format(type(e), workdir_abspath, e),
                identifier, function)


def log(job, message, identifier=None, function=None):
    if function is not None:
        message = "{}: {}".format(function, message)
    if identifier is not None:
        message = "{}:{}{}".format(identifier, "" if message.startswith(str(function)) else " ", message)
    job.log(message)


def log_time(job, function_name, start_time, sample_identifier=''):
    log(job, "TIME:{}".format(int(time.time() - start_time)), sample_identifier, function_name)



###########################################################
#             WORKFLOW MANAGEMENT FUNCTIONS               #
###########################################################


def generate_config():
    return textwrap.dedent("""
        # nanopore Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-marginphase run"
        #
        # URLs can take the form: {scheme}
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
        ##############################################################################################################

        # Required: Output location of sample. Can be full path to a directory (with or without file:// protocol) or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        # Notes: files will be grouped by their UUID
        output-dir: file:///tmp/output

        # Required: URL {scheme} for marginPolish parameter file
        margin-params: https://raw.githubusercontent.com/benedictpaten/marginPhase/polisher/params/allParams.np.json

        # Optional: URL {scheme} for true FASTA reference. Setting this value will trigger QC data generation
        # Example:
        #   true-reference: http://server.location/of/reference.fa
        true-reference: 

        # Optional: Extra parameters for QUAST execution
        # Example (recommended for human):
        #   quast-params: --large
        quast-params: 

        # Optional: intermediate files will be saved to this directory with a job timestamp (only accepts file:// scheme)
        # Example:
        #   save-intermediate-files: file:///tmp/intermediate
        save-intermediate-files: 

    """.format(scheme=", ".join([x + '://' for x in SCHEMES])))


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #
        #   There are 3 tab-separated columns: UUID, READS_URL, ASSEMBLY_URL
        #
        #   UUID            Required    A unique identifier for the sample to be processed.
        #   READS_URL       Required    A URL [{scheme}] pointing to the sample fastq
        #   ASSEMBLY_URL    Optional    A URL [{scheme}] pointing to assembly (will skip assembly steps)
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1\tfile:///path/to/reads.fastq
        #   UUID_2\ts3://path/to/reads.fastq\ts3://path/to/assembly.fa
        #   UUID_3\tfile://path/to/reads.fastq\thttp://path/to/assembly.fa
        #
        #   Place your samples below, one per line.
        """.format(scheme="'" + ("', '".join([x + '://' for x in SCHEMES])) + "'"))


def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message

    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def parse_samples(config, path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """

    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')

            # validate structure
            if len(sample) < 2:
                raise UserError('Bad manifest format! Required at least 2 tab-separated columns, got: {}'.format(sample))
            if len(sample) > 3:
                raise UserError('Bad manifest format! Required at most 3 tab-separated columns, got: {}'.format(sample))

            # extract sample parts
            uuid = sample[0]
            reads_url = sample[1]
            assembly_url = None
            if len(sample) > 2: assembly_url = sample[2]

            # sanity checks
            if urlparse(reads_url).scheme not in SCHEMES:
                raise UserError("Reads URL {} is missing scheme. Expected: {}".format(assembly_url, [x + "://" for x in SCHEMES]))
            if not is_fasta(reads_url) and not is_fastq(reads_url):
                raise UserError("Input type (fasta or fastq) of sample {} could not be determined: {}".format(
                    uuid, reads_url))
            if assembly_url is not None and urlparse(assembly_url).scheme not in SCHEMES:
                raise UserError("Assembly URL {} is missing scheme. Expected: {}".format(assembly_url, [x + "://" for x in SCHEMES]))

            sample = [uuid, reads_url, assembly_url]
            samples.append(sample)
    return samples


def main():
    """
    NaRLE: Nanopore Run Length Encoding Analysis Pipeline

    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    """

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the NaRLE pipeline')
    group = parser_run.add_mutually_exclusive_group()
    parser_run.add_argument('--config', default=DEFAULT_CONFIG_NAME, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default=DEFAULT_MANIFEST_NAME, type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()

    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_NAME), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_MANIFEST_NAME), generate_manifest)

    # Pipeline execution
    elif args.command == 'run':
        # sanity check
        require(os.path.exists(args.config), '{} not found. Please run '
                                             '"toil-narle generate-config"'.format(args.config))
        require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                               'run "toil-narle generate-manifest"'.format(args.manifest))

        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxsize
        config.defaultCores = config.maxCores #int(min(DEFAULT_CPU, config.maxCores))
        config.maxDisk = int(args.maxDisk) if args.maxDisk else sys.maxint
        config.maxMemory = sys.maxint
        # fix parsing of GB to int
        if args.maxMemory:
            args.maxMemory = args.maxMemory.upper()
            if args.maxMemory.endswith('B'):
                args.maxMemory = args.maxMemory.rstrip('B')
            # actual parsing
            if args.maxMemory.endswith('G'):
                config.maxMemory = int(args.maxMemory.rstrip('G')) * 1024 * 1024 * 1024
            elif args.maxMemory.endswith('M'):
                config.maxMemory = int(args.maxMemory.rstrip('M')) * 1024 * 1024
            elif args.maxMemory.endswith('K'):
                config.maxMemory = int(args.maxMemory.rstrip('K')) * 1024
            else:
                config.maxMemory = int(args.maxMemory)

        # Config sanity checks
        required_config_values = ['output_dir', 'margin_params']
        missing_config_values = list(filter(lambda x: x not in config or parsed_config[x] is None or len(parsed_config[x]) == 0,
                                            required_config_values))
        if len(missing_config_values) != 0:
            raise UserError("Required config values are missing or empty: {}".format(missing_config_values))

        # management of output
        if urlparse(config.output_dir).scheme not in [None, "", "file", "s3"]:
            raise UserError("Config parameter 'output-dir' must be file:// or s3:// protocol")
        if urlparse(config.output_dir).scheme != "s3":
            config.output_dir = config.output_dir.replace("file://", "", 1)
            mkdir_p(config.output_dir)
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'

        # quast params
        if 'quast_params' not in config or not config.quast_params:
            config.quast_params = None

        # intermediate files
        if 'save_intermediate_files' not in config or not config.save_intermediate_files:
            config.intermediate_file_location = None
        elif urlparse(config.save_intermediate_files).scheme != 'file':
            raise UserError("Config parameter 'save_intermediate_files' must be used with local (file://) output directory")
        else:
            intermediate_location = os.path.join(config.save_intermediate_files.replace("file://", "", 1),
                                                 datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            if not intermediate_location.endswith('/'): config.output_dir += '/'
            mkdir_p(intermediate_location)
            config.intermediate_file_location = intermediate_location

        # margin params configuration
        if urlparse(config.margin_params).scheme not in SCHEMES:
            raise UserError("Margin Params URL {} is missing scheme. Expected: {}".format(config.margin_params, [x + "://" for x in SCHEMES]))
        if config.true_reference is not None and urlparse(config.true_reference).scheme not in SCHEMES:
            raise UserError("True Reference URL {} is missing scheme. Expected: {}".format(config.margin_params, [x + "://" for x in SCHEMES]))

        # get samples
        samples = parse_samples(config, args.manifest)

        # Program checks
        for program in ['docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        Job.Runner.startToil(Job.wrapJobFn(map_job, prepare_input, samples, config), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
