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
Authors:
  Linda Xiang <linda.xiang@oicr.on.ca>
  Junjun Zhang <junjun.zhang@oicr.on.ca>
"""

import os
import subprocess
import sys
import json
import re
import hashlib
from argparse import ArgumentParser
from multiprocessing import cpu_count
import glob


def check_readgroups(seq_experiment_analysis):
    new_dict=[]
    # we assume information in seq_experiment_analysis has gone through
    # all validation checks in sequencing experiment submission
    # note: since advanced SONG validation is not ready, here we still validate uniqueness of
    #       read_group_id_in_bam and submitter_read_group_id

    fastq_files=[]
    read_group_id_in_bam_list=[]
    submitter_read_group_id_list=[]
    for rg in seq_experiment_analysis.get('read_groups'):
        new_dict.append({"read_groups":rg})
        new_dict[-1]['read_groups']['experiment']=seq_experiment_analysis.get('experiment')
        if rg.get('file_r2'):
            if rg.get('file_r1').split('.')[-1]!=rg.get('file_r2').split('.')[-1]:
                sys.exit('Error found: File1:%s and File2:%s are different file types' % (rg.get('file_r1'),rg.get('file_r2')))
        
        if rg.get("file_r1").endswith('.bam'):
            new_dict[-1]['format']='BAM'

            # make sure no duplicate of read_group_id_in_bam (when populated) within the same bam
            
            if rg.get('read_group_id_in_bam') in read_group_id_in_bam_list:
                sys.exit('Error found: read_group_id_in_bam duplicated %s' % (rg.get('read_group_id_in_bam')))
            else:
                read_group_id_in_bam_list.append(rg.get('read_group_id_in_bam'))
        else:
            new_dict[-1]['format']='FASTQ'

            if rg.get("file_r1") in fastq_files:
                sys.exit('Error found: same FASTQ %s must not be associated to more than one read group\n' % (rg.get("file_r1")))
            else:
                fastq_files.append(rg.get("file_r1"))
            if rg.get("file_r2") in fastq_files:
                sys.exit('Error found: same FASTQ %s must not be associated to more than one read group\n' % (rg.get("file_r2")))
            else:
                fastq_files.append(rg.get("file_r2"))

        # make sure no duplicate of submitter_read_group_id
        if rg.get('submitter_read_group_id') in submitter_read_group_id_list:
            sys.exit('Error found: submitter_read_group_id duplicated: %s' % rg['submitter_read_group_id'])
        else:
            submitter_read_group_id_list.append(rg.get('submitter_read_group_id'))

    return new_dict



def main(args):
    with open(args.metadata_json, 'r') as f:
        seq_experiment_analysis = json.load(f)

    study_id = seq_experiment_analysis['studyId']
    donor_id = seq_experiment_analysis['samples'][0]['donor']['donorId']
    specimen_id = seq_experiment_analysis['samples'][0]['specimenId']
    sample_id = seq_experiment_analysis['samples'][0]['sampleId']
    specimenType = seq_experiment_analysis['samples'][0]['specimen']['specimenType']
    tumourNormalDesignation = seq_experiment_analysis['samples'][0]['specimen']['tumourNormalDesignation']

    readgroups = check_readgroups(seq_experiment_analysis)

    ### Make 
    os.mkdir("out")
    for rg,content in enumerate(readgroups):
        print(readgroups[rg])

        
        readgroups[rg]['info']={}
        readgroups[rg]['info']['studyId']=study_id
        readgroups[rg]['info']['donorId']=donor_id
        readgroups[rg]['info']['specimenId']=specimen_id
        readgroups[rg]['info']['sampleId']=sample_id
        readgroups[rg]['info']['specimenType']=specimenType
        readgroups[rg]['info']['tumourNormalDesignation']=tumourNormalDesignation
        submitter_read_group_id=readgroups[rg]['read_groups']['submitter_read_group_id']
        os.mkdir("out/%s" % (submitter_read_group_id))

        if readgroups[rg]['read_groups']['is_paired_end']:
            dict_keys=["file_r1","file_r2"]
        else:
            dict_keys=["file_r1"]

        if readgroups[rg]['format']=='FASTQ':
            for key in dict_keys:
                file=readgroups[rg]['read_groups'][key]
                if file.endswith('.bz2'):
                    cmd = 'bunzip2 -k -c %s | gzip -c > out/%s/%s.%s.%s.%s.%s.fastq.gz' % (file,submitter_read_group_id,study_id,donor_id,sample_id,submitter_read_group_id,"R1" if key=='file_r1' else "R2")
                    subprocess.run(cmd, check=True)
                else:
                    os.symlink(
                        os.path.abspath(file),
                        "out/%s/%s.%s.%s.%s.%s" % (submitter_read_group_id,study_id,donor_id,sample_id,submitter_read_group_id,file.replace(".bz2",".gz")))                    

            with open("out/%s/%s.%s.%s.%s.json" % (submitter_read_group_id,study_id,donor_id,sample_id,submitter_read_group_id), "w") as outfile:
                json.dump(readgroups[rg], outfile,indent=2)

        else:
            with open("out/%s/%s.%s.%s.%s.json" % (submitter_read_group_id,study_id,donor_id,sample_id,submitter_read_group_id), "w") as outfile:
                json.dump(readgroups[rg], outfile,indent=2)

            file=readgroups[rg].get('read_groups').get("file_r1")
            #md5= [f['fileMd5sum'] for f in seq_experiment_analysis["files"] if f['fileName']==file][0]
            os.symlink(
                os.path.abspath(file),
                "out/%s/%s.%s.%s.%s" % (submitter_read_group_id,study_id,donor_id,sample_id,file))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--seq-files", dest="seq_files", required=True,
                        help="Seq files to process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata-json", dest="metadata_json", required=True,
                        help="JSON file containing sequencing_experiment analysis")
    args = parser.parse_args()

    main(args)