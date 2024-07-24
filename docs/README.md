# nfcore/dnaseqaln: Documentation

The nfcore/dnaseqaln documentation is split into the following pages:

- [Usage](usage.md)
Alignment starting from local sequencing payload via BWA
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --local_sequencing_json test/data/local_sequencing.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem_aln,markdup,cleanup
```
Alignment starting from local sequencing payload via BWA2
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --local_sequencing_json test/data/local_sequencing.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem2_aln,markdup,cleanup
```
Alignment and upload starting from local sequencing payload via both BWA and BWA2
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --local_sequencing_json test/data/local_sequencing.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem2_aln,bwamem2_aln,markdup,up_qc,up_aln,cleanup
```
Alignment and upload starting from local sequencing payload via both BWA and BWA2 without markdup
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --local_sequencing_json test/data/local_sequencing.json --local_data_directory test/data --local true --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem2_aln,bwamem2_aln,up_aln,cleanup
```
Alignment starting from SONG/SCORE
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --tools bwamem2_aln,bwamem2_aln,up_aln,cleanup --reference_fasta test/reference/GRCh38_Verily_v1.fasta
```
Alignment starting from SONG/SCORE w/ Indexing
```
nextflow run main.nf -profile debug_qa,docker --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --tools index,bwamem2_aln,bwamem2_aln,up_aln,cleanup --reference_fasta test/reference/GRCh38_Verily_v1.fasta
```
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
