# nfcore/dnaseqaln: Documentation

The nfcore/dnaseqaln documentation is split into the following pages:

- [Usage](usage.md)
Alignment starting from local sequencing via BWA
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem_aln,markdup,cleanup --max_cpus 2 --max_memory 2GB --samplesheet test/data/local_sequencing.tsv
```
Alignment starting from local sequencing via BWA2
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem2_aln,markdup,cleanup --max_cpus 2 --max_memory 2GB --samplesheet test/data/local_sequencing.tsv
```
Alignment and upload starting from local sequencing via both BWA and BWA2
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem_aln,bwamem2_aln,markdup,cleanup --max_cpus 2 --max_memory 2GB --samplesheet test/data/local_sequencing.tsv
```
Alignment and upload starting from local sequencing via both BWA and BWA2 without markdup
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem_aln,bwamem2_aln,cleanup --max_cpus 2 --max_memory 2GB --samplesheet test/data/local_sequencing.tsv
```
Alignment starting using data from SONG/SCORE
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem_aln,markdup,cleanup --max_cpus 2 --max_memory 2GB --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN
```
Alignment starting using data from SONG/SCORE w/ no upload
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools bwamem_aln,markdup,cleanup --max_cpus 2 --max_memory 2GB --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN
```
Alignment starting from SONG/SCORE w/ Indexing
```
nextflow run main.nf -profile docker,rdpc_qa --reference_fasta test/reference/GRCh38_Verily_v1.fasta --tools index,bwamem_aln,markdup,cleanup --max_cpus 2 --max_memory 2GB --analysis_id 026e7dbd-8a7b-4ee1-ae7d-bd8a7b0ee120 --study_id TEST-QA --api_token NNNNNNNN-NNNN-NNNN-NNNN-NNNNNNNNNNNN
```
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
