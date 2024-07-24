#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

 Author: Junjun Zhang <junjun.zhang@oicr.on.ca>
         Linda Xiang <linda.xiang@oicr.on.ca>
 """

import os
import sys
import json
import sys
import argparse
import subprocess
import json
import re
import hashlib
import uuid
import tarfile
from datetime import date
import copy
from glob import glob
import yaml
import io
import shutil
import csv
from math import log10, isnan

workflow_full_name = {
    'dna-seq-alignment': 'DNA Seq Alignment'
}

def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()


def get_aligned_seq_basename(files_to_upload):
    # get aligned bam/cram basename from '*.wgs.grch38.(cram|bam).qc_metrics.tgz'
    # or '*.wgs.grch38.(cram|bam).oxog_metrics.tgz'
    for f in files_to_upload:
        m = re.match(r'(.+?\.(cram|bam))\.(qc_metrics|oxog_metrics|duplicates_metrics)\.tgz$', f)
        if m: return(m.group(1))

    #sys.exit('Error: missing DNA alignment QC metrics or oxog metrics file with patten: *.{bam,cram}.{qc_metrics, oxog_metrics}.tgz')


def insert_filename_friendly_rg_id(metadata):
    filename_friendly_rg_ids = set()

    # let's loop it two times, first for the rg id actually doesn't need to convert
    for rg in metadata['read_groups']:
        submitter_read_group_id = rg['submitter_read_group_id']
        filename_friendly_rg_id = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in submitter_read_group_id ])

        if filename_friendly_rg_id == submitter_read_group_id:  # no change, safe to use
            rg['filename_friendly_rg_id'] = filename_friendly_rg_id
            filename_friendly_rg_ids.add(filename_friendly_rg_id)

    for rg in metadata['read_groups']:
        submitter_read_group_id = rg['submitter_read_group_id']
        filename_friendly_rg_id = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in submitter_read_group_id ])

        if filename_friendly_rg_id == submitter_read_group_id:  # no change, already covered
            continue

        if filename_friendly_rg_id in filename_friendly_rg_ids:  # the converted new friendly ID conflicts with existing one
            for i in range(len(metadata['read_groups'])):
                if not '%s_%s' % (filename_friendly_rg_id, i+1) in filename_friendly_rg_ids:
                    filename_friendly_rg_id += '_%s' % str(i+1)
                    break

        rg['filename_friendly_rg_id'] = filename_friendly_rg_id
        filename_friendly_rg_ids.add(filename_friendly_rg_id)


def get_rg_id_from_ubam_qc(tar, metadata):
    tar_basename = os.path.basename(tar)  # TEST-PR.DO250122.SA610149.D0RE2_1_.6cae87bf9f05cdfaa4a26f2da625f3b2.lane.bam.ubam_qc_metrics.tgz
    md5sum_from_filename = tar_basename.split('.')[-5]
    if not re.match(r'^[a-f0-9]{32}$', md5sum_from_filename):
        sys.exit('Error: ubam naming not expected %s' % tar_basename)

    for rg in metadata.get("read_groups"):
        rg_id_in_bam = rg.get("read_group_id_in_bam") if rg.get("read_group_id_in_bam") else rg.get("submitter_read_group_id")
        seq_file_name = rg.get("file_r1")
        bam_name = seq_file_name if seq_file_name.endswith('.bam') else ''
        md5sum_from_metadata = hashlib.md5(("%s %s" % (bam_name, rg_id_in_bam)).encode('utf-8')).hexdigest()
        if md5sum_from_metadata == md5sum_from_filename:
            return rg.get("filename_friendly_rg_id")

    # up to this point no match found, then something wrong
    sys.exit('Error: unable to match ubam qc metric tar "%s" to read group id' % tar_basename)


def get_dupmetrics(file_to_upload):
    library = {}
    library['libraries']=[]
    with tarfile.open(file_to_upload, 'r') as tar:
        for member in tar.getmembers():
            if member.name.endswith('.txt'):
                f = tar.extractfile(member)
                cols_name = []
                for r in f:
                    row = r.decode('utf-8')                
                    if row.startswith('LIBRARY'): 
                        cols_name = row.strip().split('\t')
                        continue
                    if cols_name:
                        if not row.strip(): break
                        metric = {}
                        cols = row.strip().split('\t')
                        for n, c in zip(cols_name, cols):
                            if n == "LIBRARY": metric.update({n: c})
                            elif '.' in c or 'e' in c: metric.update({n: float(c)}) 
                            else: metric.update({n: int(c)})
                        library['libraries'].append(metric)      
    return library

def get_alnmetrics(file_to_upload):
    metric = {}
    collected_sum_fields = {
        'raw total sequences': 'total_reads',
        'reads mapped': 'mapped_reads',
        'reads paired': 'paired_reads',
        'reads properly paired': 'properly_paired_reads',
        'pairs on different chromosomes': 'pairs_on_different_chromosomes',
        'total length': 'total_bases',
        'bases mapped (cigar)': 'mapped_bases_cigar',
        'mismatches': 'mismatch_bases',
        'error rate': 'error_rate',
        'bases duplicated': 'duplicated_bases',
        'insert size average': 'average_insert_size',
        'average length': 'average_length'
    }
    with tarfile.open(file_to_upload, 'r') as tar:
        for member in tar.getmembers():
            if member.name.endswith('stats'):
                f = tar.extractfile(member)
                cols_name = []
                for r in f:
                    row = r.decode('utf-8')          
                    if row.startswith('SN'): 
                        n = row.strip().split('\t')[1].replace(":","")
                        c = row.strip().split('\t')[2].split('\t')[0]
                        if n not in collected_sum_fields: continue
                        else:
                            if '.' in c or 'e' in c: metric[collected_sum_fields[n]]=float(c)
                            else: metric[collected_sum_fields[n]]=int(c)
    print(metric)
    return metric

def get_oxometrics(file_to_upload):
    library = []
    with tarfile.open(file_to_upload, 'r') as tar:
        for member in tar.getmembers():
            if member.name.endswith('.txt'):
                f = tar.extractfile(member)
                num = 0
                NTOT = 0
                NALTOXO = 0
                NALTNON = 0
                oxoQ = float('NaN')
                oxoReader = csv.DictReader(filter(lambda row: not re.match(r'^#|\s*$', row), io.TextIOWrapper(f, encoding='utf-8')), delimiter='\t')

                for line in oxoReader:
                    CTXT = line['CONTEXT']
                    if re.match('CCG', CTXT):
                        num = num + 1
                        NTOT = NTOT + int(line['TOTAL_BASES'])
                        NALTOXO = NALTOXO + int(line['ALT_OXO_BASES'])
                        NALTNON = NALTNON + int(line['ALT_NONOXO_BASES'])
                        oxoQ = float(line['OXIDATION_Q']) if '?' not in line['OXIDATION_Q'] else float('NaN')

                if num > 1:
                    if NTOT > 0:
                        er = float(max(NALTOXO - NALTNON, 1.0001)) / float(NTOT)
                        oxoQ = -10.0 * log10(er)
                    else:
                        oxoQ = float('NaN')
    library.append({"context": "CCG","oxoQ_score": float('%.4f' % oxoQ) if not isnan(oxoQ) else None})
    return library

def get_readgroupinfo(file,rg_id,size):
    cols_name = []
    for r in file:
        row = r.decode('utf-8')                
        if row.startswith('TOTAL_READS'): 
            cols_name = row.strip().split('\t')
            continue
        if cols_name:
            if not row.strip(): break
            metric = {}
            cols = row.strip().split('\t')
            for n, c in zip(cols_name, cols):
                if n == "LIBRARY": metric.update({n: c})
                elif '.' in c or 'e' in c: metric.update({n.lower(): float(c)}) 
                else: metric.update({n.lower(): int(c)})
    metric['file_size']=size
    metric['read_group_id']=rg_id  
    return {x:metric[x] for x in ['total_reads','pf_reads','read_length','file_size','read_group_id']}

def get_files_info(file_to_upload, seq_experiment_analysis_dict):
    print(file_to_upload)
    file_info = {
        'fileName': os.path.basename(file_to_upload),
        'fileType': file_to_upload.split(".")[-1].upper(),
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'info': {
            'data_category': 'Quality Control Metrics',
            'data_subtypes': None,
            'files_in_tgz': []
        }
    }

    if re.match(r'.+?\.ubam_qc_metrics\.tgz$', file_to_upload):
        file_info.update({'dataType': 'Sequencing QC'})
        file_info['info']['data_subtypes'] = ['Read Group Metrics']
        file_info['info'].update({'description': 'Read group level QC metrics generated by Picard CollectQualityYieldMetrics.'})
        file_info['info'].update({'analysis_tools': ['Picard:CollectQualityYieldMetrics']})
    elif re.match(r'.+?\.(cram|bam)\.qc_metrics\.tgz$', file_to_upload):
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Alignment Metrics']
        file_info['info'].update({'description': 'Alignment QC metrics generated by Samtools stats.'})
        file_info['info'].update({'analysis_tools': ['Samtools:stats']})
    elif re.match(r'.+?\.duplicates_metrics\.tgz$', file_to_upload):
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Duplicates Metrics']
        file_info['info'].update({'description': 'Duplicates metrics generated by biobambam2 bammarkduplicates2.'})
        file_info['info'].update({'analysis_tools': ['biobambam2:bammarkduplicates2']})
    elif re.match(r'.+?\.oxog_metrics\.tgz$', file_to_upload):
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['OxoG Metrics']
        file_info['info'].update({'description': 'OxoG metrics generated by GATK CollectOxoGMetrics.'})
        file_info['info'].update({'analysis_tools': ['GATK:CollectOxoGMetrics']})
    else:
        sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)

    with tarfile.open(file_to_upload, 'r') as tar:
        for member in tar.getmembers():
            file_info['info']['files_in_tgz'].append(os.path.basename(member.name))

    # retrieve duplicates metrics from the file
        if file_info['info']['data_subtypes'][0] == 'Alignment Metrics':
            extra_info = get_alnmetrics(file_to_upload)
        if file_info['info']['data_subtypes'][0] == 'OxoG Metrics':
            extra_info = get_oxometrics(file_to_upload)
        if file_info['info']['data_subtypes'][0] == 'Duplicates Metrics':
            extra_info = get_dupmetrics(file_to_upload)
        if file_info['info']['data_subtypes'][0] == 'Read Group Metrics':
            f = tar.extractfile(member)
            extra_info = get_readgroupinfo(f,".".join(member.name.split(".")[3:-1]),member.size)
    if extra_info:
         file_info['info'].update({'metrics': extra_info})

    return file_info


def get_sample_info(sample_list):
    samples = copy.deepcopy(sample_list)
    for sample in samples:
        for item in ['info', 'sampleId', 'specimenId', 'donorId', 'studyId']:
            sample.pop(item, None)
            sample['specimen'].pop(item, None)
            sample['donor'].pop(item, None)

    return samples


def main(args):
    with open(args.seq_experiment_analysis, 'r') as f:
        seq_experiment_analysis_dict = json.load(f)

    payload = {
        'analysisType': {
            'name': 'qc_metrics'
        },
        'studyId': seq_experiment_analysis_dict.get('studyId'),
        'info': {},
        'workflow': {
            'workflow_name': workflow_full_name.get(args.wf_name, args.wf_name),
            'workflow_version': args.wf_version,
            'genome_build': args.genome_build,
            'run_id': args.wf_run,
            'session_id': args.wf_session,
            'inputs': [
                {
                    'analysis_type': 'sequencing_experiment',
                    'input_analysis_id': seq_experiment_analysis_dict.get('analysisId')
                }
            ]
        },
        'files': [],
        'experiment': seq_experiment_analysis_dict.get('experiment'),
        'samples': get_sample_info(seq_experiment_analysis_dict.get('samples'))
    }

    # pass `info` dict from seq_experiment payload to new payload
    if 'info' in seq_experiment_analysis_dict and isinstance(seq_experiment_analysis_dict['info'], dict):
        payload['info'] = seq_experiment_analysis_dict['info']
    else:
        payload.pop('info')

    if 'library_strategy' in payload['experiment']:
        experimental_strategy = payload['experiment'].pop('library_strategy')
        payload['experiment']['experimental_strategy'] = experimental_strategy

    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    insert_filename_friendly_rg_id(seq_experiment_analysis_dict)

    aligned_seq_basename = get_aligned_seq_basename(args.files_to_upload)

    # get file of the payload
    for f in sorted(args.files_to_upload):
        # renmame duplicates_metrics file to have the same base name as the aligned seq
        if re.match(r'.+\.duplicates_metrics\.tgz$', f):
            new_name = '%s.duplicates_metrics.tgz' % aligned_seq_basename
            dst = os.path.join("out", new_name)
            os.symlink(os.path.abspath(f), dst)
            f = new_name
        # renmame ubam_qc_metrics file to have the same base name as the aligned seq
        elif re.match(r'.+?\.lane\.bam\.ubam_qc_metrics\.tgz$', f):
            rg_id = get_rg_id_from_ubam_qc(f, seq_experiment_analysis_dict)
            new_name = '%s.%s.ubam_qc_metrics.tgz' % (re.sub(r'\.aln\.(cram|bam)$', '', aligned_seq_basename), rg_id)
            dst = os.path.join("out", new_name)
            os.symlink(os.path.abspath(f), dst)
            f = new_name
        else:
            dst = os.path.join(new_dir, f)
            os.symlink(os.path.abspath(f), dst)

        payload['files'].append(get_files_info(f, seq_experiment_analysis_dict))

    with open("%s.dna_seq_qc.payload.json" % str(uuid.uuid4()), 'w') as f:
        f.write(json.dumps(payload, indent=2))




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool: payload-dnaseq-alignment')
    parser.add_argument("-f", "--files_to_upload", dest="files_to_upload", type=str, required=True,
                        nargs="+", help="Aligned reads files to upload")
    parser.add_argument("-a", "--seq_experiment_analysis", dest="seq_experiment_analysis", required=True,
                        help="Input analysis for sequencing experiment", type=str)
    parser.add_argument("-u", "--read_group_ubam_analysis", dest="read_group_ubam_analysis", default=[],
                        help="Input payloads for the analysis", type=str, nargs='+')
    parser.add_argument("-w", "--wf_name", dest="wf_name", required=True, help="Workflow name")
    parser.add_argument("-v", "--wf_version", dest="wf_version", required=True, help="Workflow version")
    parser.add_argument("-r", "--wf_run", dest="wf_run", required=True, help="workflow run ID")
    parser.add_argument("-s", "--wf_session", dest="wf_session", required=True, help="workflow session ID")
    parser.add_argument("-b", "--genome_build", dest="genome_build", default="", help="Genome build")
    parser.add_argument("-p", "--pipeline_yml", dest="pipeline_yml", required=False, help="Pipeline info in yaml")
    args = parser.parse_args()

    main(args)