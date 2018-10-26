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

from toil_marginphase.marginphase_pipeline import prepare_input as mp_prepare_input
from toil_marginphase.marginphase_pipeline import run_margin_phase as mp_run_margin_phase
from toil_marginphase.marginphase_pipeline import merge_chunks as mp_merge_chunks
from toil_marginphase.marginphase_core import DOCKER_CPECAN_IMG_DEFAULT, DOCKER_CPECAN_TAG_DEFAULT, \
    DOCKER_MARGIN_PHASE_IMG_DEFAULT, DOCKER_MARGIN_PHASE_TAG_DEFAULT, TOIL_JOBSTORE_PROTOCOL, \
    CI_CHUNK_INDEX, CI_OUTPUT_FILE_ID, ID_MERGED, CI_UUID, TAG_HAPLOTYPE, percent




# input / output schemes
SCHEMES = ('http', 'file', 's3', 'ftp')

# identifiers for RLE timing
RLE_POST_ALIGN='post_align'
RLE_PRE_ALIGN='pre_align'
RLE_NONE='none'
RLE_ALIGN_TYPES=[RLE_POST_ALIGN, RLE_PRE_ALIGN, RLE_NONE]

# nucleotide types (for RLE)
NUCL_FASTA='fasta'
NUCL_FASTQ='fastq'

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

DOCKER_NAPPER_IMG = "tpesout/napper"
DOCKER_NAPPER_TAG = "latest"
DOCKER_NAPPER_LOG = "napper.log"

DOCKER_ALIGNQC_IMG = "tpesout/alignqc"
DOCKER_ALIGNQC_TAG = "latest"
DOCKER_ALIGNQC_LOG = "alignqc.log"

DOCKER_RTG_IMG = "tpesout/rtg_tools"
DOCKER_RTG_TAG = "latest"
DOCKER_RTG_LOG = "rtg.log"

# resource
DEFAULT_CPU = 8

ALN_FASTQ_MEM_FACTOR = 1
ALN_FASTQ_DSK_FACTOR = 1
ALN_REF_MEM_FACTOR = 1
ALN_REF_DSK_FACTOR = 1

AQC_BAM_MEM_FACTOR = 1
AQC_BAM_DSK_FACTOR = 1
AQC_REF_MEM_FACTOR = 1
AQC_REF_DSK_FACTOR = 1

# job data
JD_RAW_FASTQ = 'raw_fastq'
JD_RLE_FASTQ = 'rle_fastq'
JD_RLE_FASTQ_RL = 'rle_fastq_rl'
JD_RAW_REF = 'raw_ref'
JD_RLE_REF = 'rle_ref'
JD_RLE_REF_RL = 'rle_ref_rl'
JD_RAW_BAM = 'raw_bam'
JD_RLE_BAM = 'rle_bam'
JD_ALIGNQC = 'alignqc'
JD_MP_VCF = 'mp_vcf'
JD_MP_SAM = 'mp_sam'

# haplotyped sam
HS_CHROM = 'chrom'
HS_PHASE_BLOCK = 'phase_block'
HS_HAP1_FILEID = 'h1_fid'
HS_HAP2_FILEID = 'h2_fid'

# file types
FT_FASTA = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
FT_FASTQ = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
is_type = lambda input, type_suffixes: any(map(lambda t_suffix: input.endswith(t_suffix), type_suffixes))
is_fasta = lambda input: is_type(input, FT_FASTA)
is_fastq = lambda input: is_type(input, FT_FASTQ)


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
    uuid, data_url, reference_url, rle_type = sample
    config.uuid = uuid
    config.reference_url = reference_url
    config.rle_type = rle_type
    config.read_filetype = NUCL_FASTA if is_fasta(data_url) else (NUCL_FASTQ if is_fastq(data_url) else None)
    if config.intermediate_file_location is not None:
        config.intermediate_file_location = os.path.join(config.intermediate_file_location, uuid)
        mkdir_p(config.intermediate_file_location)
    work_dir = job.tempDir
    log(job, "START:{}".format(datetime.datetime.now()), uuid)
    log(job, "Preparing input with URL:{}, reference_url:{}, rle_type:{}".format(data_url, reference_url, rle_type),
        uuid, "prepare_input")

    # resource estimation
    config.maxCores = min(config.maxCores, multiprocessing.cpu_count())
    config.defaultCores = min(DEFAULT_CPU, config.maxCores)
    config.maxMemory = min(config.maxMemory, int(physicalMemory() * .95))
    config.maxDisk = min(config.maxDisk, int(physicalDisk(None, toilWorkflowDir=work_dir) * .95))

    ### download files, save to jobstore ###
    #ref fasta
    download_url(reference_url, work_dir=work_dir)
    ref_genome_basename = os.path.basename(reference_url)
    ref_genome_location = os.path.join(work_dir, "{}.original.{}".format(uuid, ref_genome_basename))
    os.rename(os.path.join(work_dir, ref_genome_basename), ref_genome_location)
    ref_genome_fileid = job.fileStore.writeGlobalFile(ref_genome_location)
    ref_genome_size = os.stat(ref_genome_location).st_size

    # download data
    download_url(data_url, work_dir=work_dir)
    fastq_basename = os.path.basename(data_url)
    fastq_location = os.path.join(work_dir, "{}.original.{}".format(uuid, fastq_basename))
    os.rename(os.path.join(work_dir, fastq_basename), fastq_location)
    fastq_fileid = job.fileStore.writeGlobalFile(fastq_location)
    fastq_size = os.stat(fastq_location).st_size

    # save files
    if config.intermediate_file_location is not None:
        copy_files(file_paths=[ref_genome_location, fastq_location], output_dir=config.intermediate_file_location)

    # do RLE encoding (if appropriate)
    rle_reference_fileid = None
    rle_reference_rl_fileid = None
    rle_fastq_fileid = None
    rle_fastq_rl_fileid = None
    if config.rle_type != RLE_NONE:
        # do the encoding
        ref_rle_basename, ref_rle_rl_basename = docker_napper_rle(job, config, work_dir, os.path.basename(ref_genome_location),
                                                           "{}.rle.{}".format(uuid, ref_genome_basename))
        # get locations
        ref_rle_location = os.path.join(work_dir, ref_rle_basename)
        ref_rle_rl_location = os.path.join(work_dir, ref_rle_rl_basename)
        # save to jobstore
        rle_reference_fileid = job.fileStore.writeGlobalFile(ref_rle_location)
        rle_reference_rl_fileid = job.fileStore.writeGlobalFile(ref_rle_rl_location)
        ref_genome_size = os.stat(ref_rle_location).st_size
        # copy files if appropriate
        if config.intermediate_file_location is not None:
            log_location = os.path.join(work_dir, "napper.{}.log".format(ref_rle_basename))
            os.rename(os.path.join(work_dir, DOCKER_NAPPER_LOG), log_location)
            copy_files(file_paths=[ref_rle_location, ref_rle_rl_location, log_location],
                       output_dir=config.intermediate_file_location)
    if config.rle_type == RLE_PRE_ALIGN:
        # do the encoding
        fastq_rle_basename, fastq_rle_rl_basename = docker_napper_rle(job, config, work_dir, os.path.basename(fastq_location),
                                                           "{}.rle.{}".format(uuid, fastq_basename))
        # get locations
        fastq_rle_location = os.path.join(work_dir, fastq_rle_basename)
        fastq_rle_rl_location = os.path.join(work_dir, fastq_rle_rl_basename)
        # save to jobstore
        rle_fastq_fileid = job.fileStore.writeGlobalFile(fastq_rle_location)
        rle_fastq_rl_fileid = job.fileStore.writeGlobalFile(fastq_rle_rl_location)
        fastq_size = os.stat(fastq_rle_location).st_size
        # copy files if appropriate
        if config.intermediate_file_location is not None:
            log_location = os.path.join(work_dir, "napper.{}.log".format(fastq_rle_basename))
            os.rename(os.path.join(work_dir, DOCKER_NAPPER_LOG), log_location)
            copy_files(file_paths=[fastq_rle_location, fastq_rle_rl_location, log_location],
                       output_dir=config.intermediate_file_location)

    # job data
    job_data = {
        JD_RAW_FASTQ: fastq_fileid,
        JD_RLE_FASTQ: rle_fastq_fileid,
        JD_RLE_FASTQ_RL: rle_fastq_rl_fileid,
        JD_RAW_REF: ref_genome_fileid,
        JD_RLE_REF: rle_reference_fileid,
        JD_RLE_REF_RL: rle_reference_rl_fileid,
        JD_RAW_BAM: None,
        JD_RLE_BAM: None,
    }

    # alignment resource usage
    aln_cpu = config.defaultCores
    aln_mem = int(min(int(fastq_size * ALN_FASTQ_MEM_FACTOR + ref_genome_size * ALN_REF_MEM_FACTOR),
                      config.maxMemory))
    aln_dsk = int(min(int(fastq_size * ALN_FASTQ_DSK_FACTOR + ref_genome_size * ALN_REF_DSK_FACTOR),
                      config.maxDisk))

    # submit align job
    align_job = job.addChildJobFn(align, config, job_data, memory=aln_mem, cores=aln_cpu, disk=aln_dsk)

    # submit final consolidation job
    consolidate_output_job = align_job.addFollowOnJobFn(consolidate_output, config, align_job.rv())

    # log
    log_generic_job_debug(job, config.uuid, 'prepare_input', work_dir=work_dir)
    log_time(job, "prepare_input", start, config.uuid)


def align(job, config, job_data):
    # prep
    start = time.time()
    uuid = config.uuid
    read_filetype = config.read_filetype
    work_dir = job.tempDir
    rle_input = config.rle_type in [RLE_PRE_ALIGN]
    rle_identifier = "rle" if rle_input else "raw"

    # download files
    input_reads_filename = "{}.{}.{}".format(uuid, rle_identifier, read_filetype)
    input_reads_location = os.path.join(work_dir, input_reads_filename)
    job.fileStore.readGlobalFile(job_data[JD_RLE_FASTQ if rle_input else JD_RAW_FASTQ], input_reads_location)
    input_ref_filename = "{}.{}.reference.fa".format(uuid, rle_identifier)
    input_ref_location = os.path.join(work_dir, input_ref_filename)
    job.fileStore.readGlobalFile(job_data[JD_RLE_REF if rle_input else JD_RAW_REF], input_ref_location)

    # do alignment
    alignment_filename = "{}.{}.unsorted.sam".format(uuid, rle_identifier)
    alignment_location = os.path.join(work_dir, alignment_filename)
    docker_minimap2(job, config, work_dir, input_reads_filename, input_ref_filename, output_filename=alignment_filename,
                    kmer_size=18 if rle_input else 15)
    if not os.path.isfile(alignment_location):
        raise UserError("Alignment file {} not created for {}".format(alignment_filename, uuid))
    if config.intermediate_file_location is not None:
        log_location = os.path.join(work_dir, "minimap2.{}.log".format(input_reads_filename))
        os.rename(os.path.join(work_dir, DOCKER_MINIMAP2_LOG), log_location)
        copy_files(file_paths=[alignment_location, log_location], output_dir=config.intermediate_file_location)

    # sort it
    sorted_aln_filename = "{}.{}.bam".format(uuid, rle_identifier)
    sorted_aln_location = os.path.join(work_dir, sorted_aln_filename)
    docker_samtools_sort(job, config, work_dir, alignment_filename, output_filename=sorted_aln_filename)
    if not os.path.isfile(sorted_aln_location):
        raise UserError("Sorted alignment file {} not created for {}".format(sorted_aln_filename, uuid))
    if config.intermediate_file_location is not None:
        log_location = os.path.join(work_dir, "samtools_sort.{}.log".format(alignment_filename))
        os.rename(os.path.join(work_dir, DOCKER_SAMTOOLS_SORT_LOG), log_location)
        copy_files(file_paths=[sorted_aln_location, log_location], output_dir=config.intermediate_file_location)

    # save to jobstore
    sorted_aln_fileid = job.fileStore.writeGlobalFile(sorted_aln_location)
    job_data[JD_RLE_BAM if rle_input else JD_RAW_BAM] = sorted_aln_fileid

    # next steps
    return_values = list()
    return_values.append(job_data)

    # do alignment quality control (if appropriate)
    do_align_qc = False # todo move this to config
    if do_align_qc:
        # resource estimation
        sorted_aln_size = os.stat(sorted_aln_location).st_size
        ref_genome_size = os.stat(input_ref_location).st_size
        aqc_cpu = 1
        aqc_mem = int(min(int(sorted_aln_size * AQC_BAM_MEM_FACTOR + ref_genome_size * AQC_REF_MEM_FACTOR),
                          config.maxMemory))
        aqc_dsk = int(min(int(sorted_aln_size * AQC_BAM_DSK_FACTOR + ref_genome_size * AQC_REF_DSK_FACTOR),
                          config.maxDisk))
        # alignqc job
        alignqc_job = job.addChildJobFn(alignqc, config, job_data, memory=aqc_mem, cores=aqc_cpu, disk=aqc_dsk)
        return_values.append(alignqc_job.rv())

    # todo add marginphase
    margin_phase_job = job.addChildJobFn(enqueue_margin_phase, config, job_data)
    return_values.append(margin_phase_job.rv())

    # log
    log_generic_job_debug(job, uuid, 'align', work_dir=work_dir)
    log_time(job, "align", start, uuid)
    return return_values


def alignqc(job, config, job_data):
    # prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.tempDir
    rle_input = config.rle_type in [RLE_PRE_ALIGN]
    rle_identifier = "rle" if rle_input else "raw"

    # download files
    input_bam_filename = "{}.{}.bam".format(uuid, rle_identifier)
    input_bam_location = os.path.join(work_dir, input_bam_filename)
    job.fileStore.readGlobalFile(job_data[JD_RLE_BAM if rle_input else JD_RAW_BAM], input_bam_location)
    input_ref_filename = "{}.{}.ref.fa".format(uuid, rle_identifier)
    input_ref_location = os.path.join(work_dir, input_ref_filename)
    job.fileStore.readGlobalFile(job_data[JD_RLE_REF if rle_input else JD_RAW_REF], input_ref_location)

    # do alignqc
    alignqc_filename = "{}.{}.alignqc.xhtml".format(uuid, rle_identifier)
    alignqc_location = os.path.join(work_dir, alignqc_filename)
    success = docker_alignqc(job, config, work_dir, input_bam_filename, input_ref_filename, alignqc_filename)
    log_generic_job_debug(job, uuid, 'alignqc', work_dir=work_dir)

    # sanity check (failure is ok here)
    if not success:
        log(job, "AlignQC Failure!", uuid, 'alignqc')
        log_time(job, "alignqc", start, uuid)
        return job_data

    # save output
    if config.intermediate_file_location is not None:
        log_location = os.path.join(work_dir, "alignqc.{}.log".format(input_bam_filename))
        os.rename(os.path.join(work_dir, DOCKER_ALIGNQC_LOG), log_location)
        copy_files(file_paths=[alignqc_location, log_location], output_dir=config.intermediate_file_location)

    # save to jobstore
    alignqc_fileid = job.fileStore.writeGlobalFile(alignqc_location)
    job_data[JD_ALIGNQC] = alignqc_fileid

    # log
    log_generic_job_debug(job, uuid, 'alignqc', work_dir=work_dir)
    log_time(job, "alignqc", start, uuid)
    return job_data


def enqueue_margin_phase(job, config, job_data):
    # prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.tempDir
    rle_input = config.rle_type in [RLE_PRE_ALIGN]
    rle_identifier = "rle" if rle_input else "raw"

    # download files
    input_bam_filename = "{}.{}.bam".format(uuid, rle_identifier)
    input_bam_location = os.path.join(work_dir, input_bam_filename)
    job.fileStore.readGlobalFile(job_data[JD_RLE_BAM if rle_input else JD_RAW_BAM], input_bam_location)
    input_ref_filename = "{}.{}.ref.fa".format(uuid, rle_identifier)
    input_ref_location = os.path.join(work_dir, input_ref_filename)
    job.fileStore.readGlobalFile(job_data[JD_RLE_REF if rle_input else JD_RAW_REF], input_ref_location)

    # get contigs in bam
    docker_index_bam(job, config, work_dir, input_bam_filename)
    bam_contigs = docker_get_contig_names(job, config, work_dir, input_bam_filename)
    log(job, "Contigs found in {}: {}".format(input_bam_filename, bam_contigs))

    # separate bam into contigs
    margin_phase_sample_map = dict()
    for contig in bam_contigs:
        bam_contig_filename = "{}.{}.{}.bam".format(uuid, rle_identifier, contig)
        docker_samtools_view(job, config, work_dir, input_bam_filename, selection_str=contig,
                             output_filename=bam_contig_filename)

        # save files
        bam_contig_location = os.path.join(work_dir, bam_contig_filename)
        if config.intermediate_file_location is not None:
            copy_files(file_paths=[bam_contig_location], output_dir=config.intermediate_file_location)
        bam_contig_fileid = job.fileStore.writeGlobalFile(bam_contig_location)
        margin_phase_sample_map[contig] = \
            ["{}.{}".format(uuid, contig), "{}{}".format(TOIL_JOBSTORE_PROTOCOL, bam_contig_fileid), contig]

    # get contigs in reference
    all_ref_contigs = list()
    with open(input_ref_location) as ref_in:
        # how to handle contig when finished with it
        def save_ref_contig(contig, location, out):
            if out is None: return
            out.close()
            if config.intermediate_file_location is not None:
                copy_files(file_paths=[location], output_dir=config.intermediate_file_location)
            ref_contig_fileid = job.fileStore.writeGlobalFile(location)
            margin_phase_sample_map[contig].append("{}{}".format(TOIL_JOBSTORE_PROTOCOL, ref_contig_fileid))

        # data we need to track
        contig_name = None
        ref_contig_location = None
        ref_contig_out = None

        # handle all lines in ref
        try:
            for line in ref_in:
                # new contig
                if line.startswith(">"):
                    # save old contig (if exists)
                    if ref_contig_out is not None:
                        save_ref_contig(contig_name, ref_contig_location, ref_contig_out)
                    # get current contig
                    contig_name = line.lstrip(">").strip().split()[0]
                    all_ref_contigs.append(contig_name)
                    # clear old data
                    ref_contig_location = None
                    ref_contig_out = None
                    # intialize new data (if appropriate)
                    if contig_name in margin_phase_sample_map:
                        ref_contig_location = os.path.join(work_dir, "{}.{}.{}.ref.fa".format(
                            uuid, rle_identifier, contig_name))
                        ref_contig_out = open(ref_contig_location, 'w')

                # write header or line if we care about contig
                if ref_contig_out is not None:
                    ref_contig_out.write(line)
                pass
        # we don't expect problems here
        except Exception, e:
            log(job, "Got {} exception handling contig {} in {}: {}".format(type(e), contig_name, input_ref_filename, e),
                uuid, 'enqueue_margin_phase')
            raise e
        # be sure to close file if unexpected problem
        finally:
            if ref_contig_out is not None: ref_contig_out.close()

        # save the last ref_contig
        if ref_contig_out is not None:
            save_ref_contig(contig_name, ref_contig_location, ref_contig_out)

    # get all contigs we have bam and ref for (should be all)
    margin_phase_contigs = list(filter(lambda x: len(margin_phase_sample_map[x]) == 4, margin_phase_sample_map.keys()))
    margin_phase_contigs.sort()
    if len(margin_phase_contigs) != len(margin_phase_sample_map):
        missing_margin_phase_contigs = list(filter(lambda x: x not in margin_phase_contigs, margin_phase_sample_map.keys()))
        missing_margin_phase_contigs.sort()
        all_ref_contigs.sort()
        log(job, "Missing contigs in ref which were found in aligned bam!", uuid, 'enqueue_margin_phase')
        log(job, "\tFound contigs: {}".format(margin_phase_contigs), uuid, 'enqueue_margin_phase')
        log(job, "\tMissing contigs: {}.".format(missing_margin_phase_contigs), uuid, 'enqueue_margin_phase')
        log(job, "\tAll ref contigs: {}. ".format(all_ref_contigs), uuid, 'enqueue_margin_phase')

    # prep for marginPhase enqueue
    margin_phase_samples = [margin_phase_sample_map[c] for c in margin_phase_contigs]

    # add params to everything
    mp_params_location = os.path.join(work_dir, "{}.mp_params.json".format(uuid))
    write_margin_phase_params(mp_params_location)
    if config.intermediate_file_location is not None:
        copy_files(file_paths=[mp_params_location], output_dir=config.intermediate_file_location)
    mp_params_fileid = job.fileStore.writeGlobalFile(mp_params_location)
    for mps in margin_phase_samples:
        mps.append("{}{}".format(TOIL_JOBSTORE_PROTOCOL, mp_params_fileid))

    # get job configuration
    margin_phase_config = argparse.Namespace(**{
        'output_dir': None, #should not be used
        'partition_size': 2000000,
        'partition_margin': 50000,
        'intermediate_file_location': config.intermediate_file_location,
        'margin_phase_image': DOCKER_MARGIN_PHASE_IMG_DEFAULT,
        'margin_phase_tag': DOCKER_MARGIN_PHASE_TAG_DEFAULT,
        'cpecan_image': DOCKER_CPECAN_IMG_DEFAULT,
        'cpecan_tag': DOCKER_CPECAN_TAG_DEFAULT,
        'minimal_output': True,
        'minimal_cpecan_output': True,
        'cpecan_probabilities': True,
        'maxCores': config.maxCores,
        'maxDisk': config.maxDisk,
        'maxMemory': config.maxMemory,
    })

    # prep for return values
    margin_phase_rvs = list()

    # enqueue jobs
    for mps in margin_phase_samples:
        log(job, "Enqueing MP job: {}".format(mps), uuid, 'enqueue_margin_phase')
        mp_job = job.addChildJobFn(mp_prepare_input, mps, margin_phase_config, False)
        margin_phase_rvs.append(mp_job.rv())

    # enqueue consolidation job
    cmp_job = job.addFollowOnJobFn(consolidate_margin_phase_jobs, config, job_data, margin_phase_rvs)

    # log
    log_generic_job_debug(job, uuid, 'enqueue_margin_phase', work_dir=work_dir)
    log_time(job, "enqueue_margin_phase", start, uuid)
    return cmp_job.rv()


def write_margin_phase_params(file_location):
    with open(file_location, 'w') as out:
        out.write('{"alphabet":["Aa","Cc","Gg","Tt","-"],"wildcard":"Nn","haplotypeSubstitutionModel":[0.9993,0.0001,'
                  '0.0004,0.0001,0.0001,0.0001,0.9993,0.0001,0.0004,0.0001,0.0004,0.0001,0.9993,0.0001,0.0001,0.0001,'
                  '0.0004,0.0001,0.9993,0.0001,0.0001,0.0001,0.0001,0.0001,0.9996],"readErrorModel":[0.95956038,'
                  '0.00768823,0.02313562,0.00951578,0.00009999,0.01216431,0.94530676,0.00747140,0.03495754,0.00009999,'
                  '0.03353173,0.00747909,0.94699127,0.01189791,0.00009999,0.00950825,0.02377039,0.00805501,0.95856636,'
                  '0.00009999,0.00961502,0.00961502,0.00961502,0.00961502,0.96153992],"maxNotSumTransitions":true,'
                  '"minPartitionsInAColumn":100,"maxPartitionsInAColumn":100,"minPosteriorProbabilityForPartition":0.0,'
                  '"maxCoverageDepth":64,"minReadCoverageToSupportPhasingBetweenHeterozygousSites":2,'
                  '"onDiagonalReadErrorPseudoCount":1000,"offDiagonalReadErrorPseudoCount":10,"trainingIterations":0,'
                  '"estimateReadErrorProbsEmpirically":false,"filterBadReads":false,"filterMatchThreshold":0.8,'
                  '"useReferencePrior":true,"includeInvertedPartitions":true,"filterLikelyHomozygousSites":true,'
                  '"minSecondMostFrequentBaseFilter":0,"minSecondMostFrequentBaseLogProbFilter":20,'
                  '"gapCharactersForDeletions":false,"filterAReadWithAnyOneOfTheseSamFlagsSet":3840,"compareVCFs":false,'
                  '"writeGVCF":false,"roundsOfIterativeRefinement":10,"verbose":0}\n')


def consolidate_margin_phase_jobs(job, config, job_data, margin_phase_rvs):
    # prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.tempDir
    rle_input = config.rle_type in [RLE_PRE_ALIGN]
    rle_identifier = "rle" if rle_input else "raw"

    # organize mp contigs
    mp_contig_tarballs = dict()
    for mprv in margin_phase_rvs:
        contig_uuid = "[missing]" if len(mprv) == 0 else mprv[0][CI_UUID]
        merged_files = list(filter(lambda x: x[CI_CHUNK_INDEX] == ID_MERGED, mprv))
        if len(merged_files) != 1:
            log(job, "Found {} merged files from {} marginPhase run".format(len(merged_files), contig_uuid))
            continue
        contig = contig_uuid.replace("{}.".format(uuid), '', 1)
        mp_contig_tarballs[contig] = merged_files[0][CI_OUTPUT_FILE_ID]
    log(job, "Got {}/{} merged tarballs from MarginPhase".format(len(mp_contig_tarballs), len(margin_phase_rvs)),
        uuid, 'consolidate_margin_phase_jobs')

    # prep for untarring
    untar_work_dir = os.path.join(work_dir, "untar")
    vcf_work_dir = os.path.join(work_dir, "vcf")
    sam_work_dir = os.path.join(work_dir, "sam")
    [os.mkdir(x) for x in [untar_work_dir, vcf_work_dir, sam_work_dir]]
    vcf_filenames = list()
    sam_filenames = list()

    # get untarred files
    for contig in mp_contig_tarballs.keys():
        contig_work_dir = os.path.join(untar_work_dir, contig)
        os.mkdir(contig_work_dir)
        tarball_filename = "{}.{}.tar.gz".format(uuid, contig)
        tarball_location = os.path.join(contig_work_dir, tarball_filename)
        job.fileStore.readGlobalFile(mp_contig_tarballs[contig], tarball_location)
        with tarfile.open(tarball_location) as tar:
            tar.extractall(path=contig_work_dir)
        for file in glob.glob(os.path.join(contig_work_dir, '*')):
            file_basename = os.path.basename(file)
            if file.endswith(".vcf"):
                os.rename(file, os.path.join(vcf_work_dir, file_basename))
                vcf_filenames.append(file_basename)
            elif file.endswith(".sam"):
                os.rename(file, os.path.join(sam_work_dir, file_basename))
                sam_filenames.append(file_basename)
            elif not file.endswith(".tar.gz"):
                log(job, "Unexpected file in {}: {}".format(tarball_filename, file_basename),
                    uuid, 'consolidate_margin_phase_jobs')
        shutil.rmtree(contig_work_dir)

    # merge vcf
    merged_vcf_filename = "{}.{}.marginPhase.vcf.gz".format(uuid, rle_identifier)
    merged_vcf_location = os.path.join(vcf_work_dir, merged_vcf_filename)
    docker_merge_vcfs(job, config, vcf_work_dir, vcf_filenames, merged_vcf_filename)
    if config.intermediate_file_location is not None:
        log_location = os.path.join(vcf_work_dir, "rtg_vcf_merge.{}.log".format(merged_vcf_filename))
        os.rename(os.path.join(vcf_work_dir, DOCKER_RTG_LOG), log_location)
        copy_files(file_paths=[merged_vcf_location, log_location], output_dir=config.intermediate_file_location)
    merged_vcf_fileid = job.fileStore.writeGlobalFile(merged_vcf_location)
    job_data[JD_MP_VCF] = merged_vcf_fileid
    log(job, "Merged {} VCF files into {}".format(len(vcf_filenames), merged_vcf_filename),
        uuid, 'consolidate_margin_phase_jobs')

    # merge sam
    merged_sam_filename = "{}.{}.marginPhase.sam".format(uuid, rle_identifier)
    merged_sam_location = os.path.join(sam_work_dir, merged_sam_filename)
    docker_samtools_merge(job, config, sam_work_dir, sam_filenames, merged_sam_filename)
    if config.intermediate_file_location is not None:
        copy_files(file_paths=[merged_sam_location], output_dir=config.intermediate_file_location)
    merged_sam_fileid = job.fileStore.writeGlobalFile(merged_sam_location)
    job_data[JD_MP_SAM] = merged_sam_fileid
    log(job, "Merged {} SAM files into {}".format(len(sam_filenames), merged_sam_filename),
        uuid, 'consolidate_margin_phase_jobs')

    # log
    log_generic_job_debug(job, uuid, 'consolidate_margin_phase_jobs', work_dir=work_dir)
    log_time(job, "consolidate_margin_phase_jobs", start, uuid)
    return job_data


def split_haplotyped_sam(job, config, job_data):
    HS_CHROM = 'chrom'
    HS_PHASE_BLOCK = 'phase_block'
    HS_HAP1_FILEID = 'h1_fid'
    HS_HAP2_FILEID = 'h2_fid'

    # prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.tempDir
    sam_work_dir = os.path.join(work_dir, 'split_sams')
    os.mkdir(sam_work_dir)
    rle_input = config.rle_type in [RLE_PRE_ALIGN]
    rle_identifier = "rle" if rle_input else "raw"

    # get marginPhase'd sam files
    haplotyped_sam_filename = "{}.{}.marginPhase.sam".format(uuid, rle_identifier)
    haplotyped_sam_location = os.path.join(work_dir, haplotyped_sam_filename)
    job.fileStore.readGlobalFile(job_data[JD_MP_SAM], haplotyped_sam_location)

    # get sam header
    haplotyped_header_filename = "{}.{}.marginPhase.header.txt".format(uuid, rle_identifier)
    haplotyped_header_location = os.path.join(work_dir, haplotyped_header_filename)
    with open(haplotyped_header_location, 'w') as output, open(haplotyped_sam_location, 'r') as input:
        for line in input:
            if line.startswith("@"):
                output.write(line)
            else:
                break

    # tracking for current phase blocks
    phase_block_files = dict()
    current_chrom = set()
    HAP1 = '1'
    HAP2 = '2'
    H1_OUT = '1.o'
    H2_OUT = '2.o'

    ### helper functions ###
    # this creates new files and file handles (closed when new chrom is encountered)
    def init_new_phase_block_files(chrom, phase_block):
        hap1_location = os.path.join(sam_work_dir, "{}.{}.{}.{}.1.sam".format(uuid, rle_identifier, chrom, phase_block))
        hap2_location = os.path.join(sam_work_dir, "{}.{}.{}.{}.2.sam".format(uuid, rle_identifier, chrom, phase_block))
        shutil.copy(haplotyped_header_location, hap1_location)
        shutil.copy(haplotyped_header_location, hap2_location)
        phase_block_files[chrom][phase_block] = {
            HAP1: hap1_location,
            HAP2: hap2_location,
            H1_OUT: open(hap1_location, 'w'),
            H2_OUT: open(hap2_location, 'w'),
        }
    # closes all open file handles for chromosome
    def close_phase_block_chrom(chrom):
        if chrom is None: return
        for pb in phase_block_files[chrom].values():
            if pb[H1_OUT] is not None:
                pb[H1_OUT].close()
                pb[H1_OUT] = None
            if pb[H2_OUT] is not None:
                pb[H2_OUT].close()
                pb[H2_OUT] = None
    # writes read to appropriate file handle (and creates if necessary
    def handle_phase_block_read(chrom, phase_block, is_hap1, line):
        # ensure chrom is accounted for
        if chrom not in current_chrom:
            [close_phase_block_chrom(c) for c in current_chrom]
            if chrom in phase_block_files:
                error = "Found read for inactive chrom {}! Handled: {}, current: {}".format(
                    chrom, phase_block_files.keys(), current_chrom)
                log(job, error, uuid, 'split_haplotyped_sam')
                raise UserError(error)
            phase_block_files[chrom] = dict()
            current_chrom.clear()
            current_chrom.add(chrom)
        # ensure phaseblock is accounted for
        if phase_block not in phase_block_files[chrom]:
            init_new_phase_block_files(chrom, phase_block)
        # write
        phase_block_files[chrom][phase_block][H1_OUT if is_hap1 else H2_OUT].write(line)
    # parse read and get phase block info
    def get_phase_blocks_from_read_line(line):
        # prep
        phase_blocks = list()
        line = line.rstrip().split("\t")
        read_id = line[0]
        chrom = line[2]
        # look for haplotype tag in all tags
        for tag in line[11:]:
            if not tag.startswith(TAG_HAPLOTYPE):
                continue
            # get the haplotype encodings
            tag = tag.split(":")
            encoded_pb = tag[2]
            for epb_part in encoded_pb.split(";"):
                epb_part = epb_part.split(",")
                # haplotype is 0 to indicate not happed, otherwise 1 or 2
                hap = int(epb_part[0].replace('h', ''))
                if hap not in [0,1,2]:
                    error = "Got unexpected haplotype tag in read {}: {}".format(read_id, encoded_pb)
                    log(job, error, uuid, 'split_haplotyped_sam')
                    raise UserError(error)
                if hap == 0:
                    continue
                is_hap1 = hap == 1
                # get phase block
                pb = int(epb_part[1].replace('p',''))
                # save relevant information
                phase_blocks.append([chrom, pb, is_hap1])
        # this may be empty
        return phase_blocks

    # having defined all the hard work above, now we just organize all reads in the sam file
    total_read_cnt = 0
    unphased_read_cnt = 0
    success = False
    try:
        with open(haplotyped_sam_location, 'r') as input:
            for line in input:
                # skip header lines
                if line.startswith("@"): continue
                # get phase blocks
                phase_blocks = get_phase_blocks_from_read_line(line)
                if len(phase_blocks) == 0:
                    unphased_read_cnt += 1
                for phase_block in phase_blocks:
                    chrom, phase_block, is_hap1 = phase_block
                    handle_phase_block_read(chrom, phase_block, is_hap1, line)
                total_read_cnt += 1
        success = True
    finally:
        # close the last set of files
        for chrom in current_chrom:
            close_phase_block_chrom(chrom)
        # lgogit
        log(job, "{} split haplotyped sam into {} phase blocks over {} contigs.".format(
            "Successfully" if success else "Failed", len(phase_block_files), sum(map(len, phase_block_files.values()))),
            uuid, 'split_haplotyped_sam')
        log(job, "Of {} total reads, {} were unphased".format(total_read_cnt, unphased_read_cnt,
                                                              percent(unphased_read_cnt, total_read_cnt)),
            uuid, 'split_haplotyped_sam')

    # set up intermediate savings
    if config.intermediate_file_location is not None:
        new_intermediate_location = os.path.join(config.intermediate_file_location, 'haplotyped_sams')
        os.mkdir(new_intermediate_location)
        def save_intermediate_file(file):
            copy_files(file_paths=[file], output_dir=new_intermediate_location)
    else:
        def save_intermediate_file(file):
            pass

    # save new files
    all_haplotyped_sams = list()
    for chrom in phase_block_files.keys():
        for phase_block in phase_block_files[chrom].keys():
            pb_sam_data = {
                HS_CHROM: chrom,
                HS_PHASE_BLOCK: phase_block,
            }
            save_intermediate_file(phase_block_files[chrom][phase_block][HAP1])
            save_intermediate_file(phase_block_files[chrom][phase_block][HAP2])
            pb_sam_data[HS_HAP1_FILEID] = job.filestore.writeGlobalFile(phase_block_files[chrom][phase_block][HAP1])
            pb_sam_data[HS_HAP2_FILEID] = job.filestore.writeGlobalFile(phase_block_files[chrom][phase_block][HAP2])
            all_haplotyped_sams.append(pb_sam_data)

    # log
    log_generic_job_debug(job, config.uuid, 'split_haplotyped_sam', work_dir=work_dir)
    log_time(job, "split_haplotyped_sam", start, uuid)
    return job_data


def consolidate_output(job, config, job_data_list):
    # prep
    start = time.time()
    uuid = config.uuid
    work_dir = job.tempDir

    # consolodate all data
    final_job_data = dict()
    for i, job_data in enumerate(job_data_list):
        log(job, "JOB_DATA_{}: {}".format(i, job_data), uuid, 'consolidate_output')
        for jdk in job_data.keys():
            jdv = job_data[jdk]
            if jdv is None: continue
            if jdk in final_job_data and final_job_data[jdk] != jdv:
                log(job, "Found different non-null values during data consolidation: {} -> [{},{}]".format(
                    jdk, jdv, final_job_data[jdk]), uuid, 'consolidate_output')
            final_job_data[jdk] = jdv

    # loggit
    log(job, "Consolidated job data: {}".format(final_job_data), uuid, 'consolidate_output')

    # log
    log_generic_job_debug(job, config.uuid, 'consolidate_output', work_dir=work_dir)
    log_time(job, "consolidate_output", start, uuid)
    return final_job_data


###########################################################
#                 JOB UTILITY FUNCTIONS                   #
###########################################################

def docker_alignqc(job, config, work_dir, bam_filename, ref_filename, output_filename, extra_args=None):
    # prep
    bam_location = os.path.join("/data", bam_filename)
    ref_location = os.path.join("/data", ref_filename)
    out_location = os.path.join("/data", output_filename)
    args = ['analyze', bam_location, '-g', ref_location, '--no_transcriptome', '-o', out_location]
    if extra_args is not None:
        duplicated_args = list(filter(lambda x: x in args, extra_args))
        if len(duplicated_args) > 0:
            log(job, "Duplicated args in call to alignqc: {}".format(duplicated_args), config.uuid, 'alignqc')
        args.extend(extra_args)

    try:
        # call
        docker_call(job, config, work_dir, args, DOCKER_ALIGNQC_IMG, DOCKER_ALIGNQC_TAG)

        # loggit
        log_file = os.path.join(work_dir, DOCKER_ALIGNQC_LOG)
        log_debug_from_docker(job, log_file, config.uuid, 'alignqc', input_file_locations=[
            os.path.join(work_dir, bam_filename), os.path.join(work_dir, ref_filename)])

        # sanity check
        require_docker_file_output(job, config, work_dir, [output_filename], 'alignqc', DOCKER_ALIGNQC_LOG)
    except Exception, e:
        log(job, "Error ({}) in AlignQC output: {}".format(type(e), e), config.uuid, 'alignqc')
        return None

    # return output location
    return output_filename


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


def docker_minimap2(job, config, work_dir, input_filename, ref_filename, output_filename=None, kmer_size=15, extra_args=None):
    # prep
    data_location = os.path.join("/data", input_filename)
    ref_location = os.path.join("/data", ref_filename)
    args = ["-k", str(kmer_size), "-t", str(job.cores), '-a']
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


def docker_napper_rle(job, config, work_dir, filename, output_filename=None, type=None):
    # try to infer type (if necessary)
    if type is None:
        if is_fasta(filename):
            type = NUCL_FASTA
        elif is_fastq(filename):
            type = NUCL_FASTQ
        else:
            raise UserError("Could not infer filetype of {} from: {}".format([NUCL_FASTA, NUCL_FASTQ], filename))

    # prep
    if output_filename is None:
        output_filename = "RLE.{}".format(filename)
    output_rl_filename = "{}.rl".format(output_filename)
    docker_params = [type, '--input', os.path.join("/data", filename), '--output', os.path.join("/data", output_filename)]

    # call
    docker_call(job, config, work_dir, docker_params, DOCKER_NAPPER_IMG, DOCKER_NAPPER_TAG)

    # loggit
    log_debug_from_docker(job, os.path.join(work_dir, DOCKER_NAPPER_LOG), config.uuid, 'napper_rle',
                          input_file_locations=[os.path.join(work_dir, filename)])

    # sanity check
    workdir_rl_location = os.path.join(work_dir, output_rl_filename)
    if not os.path.isfile(workdir_rl_location):
        raise UserError("RL file not created for {}".format(filename))

    # return filenames
    return output_filename, output_rl_filename


###########################################################
#                 TOIL UTILITY FUNCTIONS                  #
###########################################################


def docker_call(job, config, work_dir, params, image, tag, detach=False, stdout=None, stderr=None):
    tagged_image = "{}:{}".format(image, tag)
    if DOCKER_LOGGING:
        log(job, "Running '{}' with parameters: {}".format(tagged_image, params), config.uuid, 'docker')
    return apiDockerCall(job, tagged_image, working_dir=work_dir, parameters=params, user="root", detach=detach,
                         stdout=stdout, stderr=stderr)


def require_docker_file_output(job, config, work_dir, output_filenames, function_id, log_filename=None,
                               max_directory_contents=None, max_log_lines=None):
    missing_filenames = list(filter(lambda x: not os.path.isfile(os.path.join(work_dir, x)), output_filenames))
    if len(missing_filenames) > 0:
        # document missing
        log(job, "Missing files after docker call: ", config.uuid, function_id)
        for missing in missing_filenames:
            log(job, "\t{}".format(missing), config.uuid, function_id)

        # document contents
        directory_contents = os.listdir(work_dir)
        if max_directory_contents is not None and len(directory_contents) > max_directory_contents:
            directory_contents = directory_contents[0:max_directory_contents]
            directory_contents.append("[{} items total]".format(len(directory_contents)))
        log(job, "Current files in work_dir: {}".format(work_dir), config.uuid, function_id)
        for missing in directory_contents:
            log(job, "\t{}".format(missing), config.uuid, function_id)

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
        raise UserError("Missing files after running {} on {}: {}".format(function_id, config.uuid, missing_filenames))


def log_debug_from_docker(job, log_file_location, identifier, function, input_file_locations=None):
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

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        output-dir: file:///tmp

        # Optional: URL {scheme} for default FASTA reference
        default-reference: file://path/to/reference.fa

        # Optional: Default value for whether to perform RLE encoding on reference and fastq
        default-rle: none

        # Optional: for debugging, this will save intermediate files to the output directory (only works for file:// scheme)
        save-intermediate-files: False

    """.format(scheme=", ".join([x + '://' for x in SCHEMES])))


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #
        #   There are 6 tab-separated columns: UUID, URL, contig name, reference fasta URL, parameters URL, reference vcf URL
        #
        #   UUID            Required    A unique identifier for the sample to be processed.
        #   URL             Required    A URL [{scheme}] pointing to the sample fastq
        #   REFERENCE_URL   Optional    A URL [{scheme}] pointing to reference fasta file
        #   PERFORM_RLE     Optional    Identifier [{rle}] for when to perform run-length-encoding
        #
        #   For the two optional values, there must be a value specified either in this manifest, or in the configuration
        #   file.  Any value specified in the manifest overrides what is specified in the config file
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1\tfile:///path/to/file.fastq
        #   UUID_2\ts3://path/to/file.fastq\ts3://path/to/reference.fa
        #   UUID_3\ts3://path/to/file.fastq\ts3://path/to/reference.fa\tpost_align
        #   UUID_4\ts3://path/to/file.fastq\t\tpre_align
        #
        #   Place your samples below, one per line.
        """.format(scheme="'" + ("', '".join([x + '://' for x in SCHEMES])) + "'",
                   rle="'" + ("', '".join(RLE_ALIGN_TYPES)) + "'"))


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
            if len(sample) > 4:
                raise UserError('Bad manifest format! Required at most 6 tab-separated columns, got: {}'.format(sample))

            # extract sample parts
            uuid = sample[0]
            data_url = sample[1]
            reference_url, rle_type = "", ""
            if len(sample) > 2: reference_url = sample[2]
            if len(sample) > 3: rle_type = sample[3]

            # fill defaults
            if len(reference_url) == 0:
                reference_url = config.default_reference
            if len(rle_type) == 0:
                rle_type = config.default_rle
                if rle_type is None or len(rle_type) == 0:
                    rle_type = RLE_NONE

            # sanity checks
            if not is_fasta(data_url) and not is_fastq(data_url):
                raise UserError("Input type (fasta or fastq) of sample {} could not be determined: {}".format(
                    uuid, data_url))
            if reference_url is None or len(reference_url) == 0:
                raise UserError("Sample {} missing reference URL, no default specified".format(uuid))
            if rle_type not in RLE_ALIGN_TYPES:
                raise UserError("Sample {} has unknown RLE specification {}.  Expected one of {}"
                                .format(uuid, rle_type, RLE_ALIGN_TYPES))
            if urlparse(reference_url).scheme not in SCHEMES:
                raise UserError("Reference URL {} is missing scheme. Expected: {}".format(reference_url, [x + "://" for x in SCHEMES]))
            if urlparse(data_url).scheme not in SCHEMES:
                raise UserError("Data URL {} is missing scheme. Expected: {}".format(reference_url, [x + "://" for x in SCHEMES]))

            sample = [uuid, data_url, reference_url, rle_type]
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
        config.defaultCores = int(min(DEFAULT_CPU, config.maxCores))
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
        required_config_values = ['output_dir']
        missing_config_values = list(filter(lambda x: x not in config, required_config_values))
        if len(missing_config_values) != 0:
            raise UserError("Missing required config values: {}".format(missing_config_values))
        if urlparse(config.output_dir).scheme != "s3":
            config.output_dir = config.output_dir.replace("file://", "", 1)
            mkdir_p(config.output_dir)
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        if 'save_intermediate_files' not in config or not config.save_intermediate_files:
            config.intermediate_file_location = None
        elif urlparse(config.output_dir).scheme == "s3":
            raise UserError("Config parameter 'save_intermediate_files' cannot be used with s3 output directory")
        else:
            intermediate_location = os.path.join(config.output_dir, "intermediate", datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            mkdir_p(intermediate_location)
            config.intermediate_file_location = intermediate_location

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
