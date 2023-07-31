# nfcore/dnaseqaln: Documentation

The nfcore/dnaseqaln documentation is split into the following pages:

- [Usage](usage.md)
Alignment starting from local sequencing payload
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --local_sequencing_json test/data/local_sequencing.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta  -resume --tools aln,oxo_qc,aln_qc,rg_qc,up_qc,up_aln,cleanup
```
Uploading QC starting from local aligment payload
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --local_alignment_json test/data/local_alignment.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools up_aln
```
Uploading QC starting from local QC payload
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN--local_qc_json test/data/local_qc.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools up_qc
```
Alignment starting from SONG/SCORE
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --tools aln,oxo_qc,aln_qc,rg_qc,up_qc,up_aln,cleanup --reference_fasta test/reference/GRCh38_Verily_v1.fasta
```
Alignment starting from SONG/SCORE w/ Indexing
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --tools index,aln,oxo_qc,aln_qc,rg_qc,up_qc,up_aln,cleanup --reference_fasta test/reference/GRCh38_Verily_v1.fasta
```
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
