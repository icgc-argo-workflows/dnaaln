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

def main(args):

    os.mkdir("out")

    for metadata_json in args.metadata_json:
        print(metadata_json)
        with open(metadata_json, 'r') as f:
            seq_experiment_analysis = json.load(f)

            study_id = metadata_json.split(".")[0]
            donor_id = metadata_json.split(".")[1]
            sample_id = metadata_json.split(".")[2]

            submitter_read_group_id=seq_experiment_analysis['read_groups']['submitter_read_group_id']
            read_group_id_in_bam=seq_experiment_analysis['read_groups']["read_group_id_in_bam"]

        return_file=None
        for ind,file in enumerate(args.seq_files):
            if ("%s.bam" % read_group_id_in_bam) == file.split("/")[-1]:
                return_file=file
                break
            if ("%s.bam" % submitter_read_group_id) == file.split("/")[-1]:
                return_file=file
                break

        if return_file==None:
            sys.exit("Read group %s could not be mapped to RG level Bams" % (submitter_read_group_id))

        os.mkdir("out/%s" % (submitter_read_group_id))
        os.symlink(
            os.path.abspath(return_file),
            "out/%s/%s.%s.%s.%s" % (submitter_read_group_id,study_id,donor_id,sample_id,file)
        )
        os.symlink(
            os.path.abspath(metadata_json),
            "out/%s/%s.%s.%s.%s.json" % (submitter_read_group_id,study_id,donor_id,sample_id,submitter_read_group_id)
        )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--seq-files", dest="seq_files", required=True,
                        help="Seq files to process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata-json", dest="metadata_json", required=True,
                        help="JSON file containing sequencing_experiment analysis", type=str, nargs='+')
    args = parser.parse_args()

    main(args)