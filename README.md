# biopython_pipeline
### This pipeline utilizes biopython and other bioinformatics software tools to execute multiple tasks with one command
#### These tasks include:
- Downloading SRR files from the NCBI database
- Isolating coding sequence (CDS) features from a RefSeq file to make a CDS-only fasta file
- Quantifying the transcripts per million (TPM) of CDSs using Kallisto (1)
- Finding differentially expressed genes between two time points with Sleuth (2)
- Mapping reads to a reference genome using Bowtie2 (3)
- Generating assemblies from multiple reads with SPAdes (4)
- Aligning an assembly to a nucleotide database with blast+ (5)

#### Pipeline execution steps and code
1. This pipeline works with both pre-existing fastq files and SRR IDs that can be downloaded from the NCBI database
- Downloading SRR fastq files from NCBI in the pipeline execution is the default setting
- To use fastq files that are already downloaded: use `--no-wget` flag in the command
  - If `--no-wget` is used a path to the fastq files must be specified with `--reads-dir`
  - all fastq files must be named sampleid_1.fastq.gz and sampleid_2.fastq.gz where 'sampleid' is specified in the input file
2. A working directory where a folder containing all ouputs will be added must be specified with `--working-dir`
3. A path to the comma-delimited input file containing donor_id, sample_timepoint, sampleid/SRR_accession number must be specified with `--input`
- If any category is not applicable, replace with NA
- See example_input for format example
4. An accession number of the RefSeq organism must be specified with `--accession`
- ex. NC_006273.2
5. K-mer size to use in SPAdes assembly must be specified with `--kmer-size`
6. The organism to align assembly with blast+ must be specified with `--organism`

The following must be installed prior to pipeline running: Biopython, Kallisto, Sleuth, Bowtie2, SPAdes, Blast, and Argparse

*Unix Code with NCBI search:*

`bash {/path/to/00wrapper.sh} --working-dir {/path/to/output} --input {/path/to/input_file} --accession {ncbi_accession#} --kmer-size {#} --organism {organism_name}`

Example: `bash ~/biopython_pipeline/scripts/00wrapper.sh --working-dir ~/ --input ~/input_file --accession NC_006273.2 --kmer-size 77 --organism Betaherpesvirinae`

*Unix Code with pre-existing fastq files:*

`bash {/path/to/00wrapper.sh} --working-dir {/path/output} --input {/path/input_file} --accession {ncbi_accession#} --kmer-size {#} --organism {organism_name} --no-wget --reads-dir {path}`

-> An output file - `/path/to/output/pipelineproject_claudia_krueger/PipelineProject.log` - will contain output information from the pipeline

#### Sample test data
Sample data fastq files are included within the sample_data folder with an output log file to test the pipeline and compare output. To run:

`bash {/path/to/00wrapper.sh} --working-dir {/path/output} --input {path/to/biopython_pipeline/sample_data/input_file_mini} --accession NC_006273.2 --kmer-size 77 --organism Betaherpesvirinae --no-wget --reads-dir {path/to/biopython_pipeline/sample_data}`

Example: `bash ~/biopython_pipeline/scripts/00wrapper.sh --working-dir ~/ --input ~/biopython_pipeline/sample_data/input_file_mini --accession NC_006273.2 --kmer-size 77 --organism Betaherpesvirinae --no-wget --reads-dir ~/biopython_pipeline/sample_data`

#### Citations
Test and sample data: Cheng, S., Caviness, K., Buehler, J., Smithey, M., Nikolich-Žugich, J., & Goodrum, F. (2017). Transcriptome-wide characterization of human cytomegalovirus in natural infection and experimental latency. Proceedings of the National Academy of Sciences of the United States of America, 114(49), E10586–E10595. https://doi.org/10.1073/pnas.1710522114

(1) Bray, N., Pimentel, H., Melsted, P. et al. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525–527 (2016). https://doi.org/10.1038/nbt.3519

(2) Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. Nature methods, 14(7), 687–690. https://doi.org/10.1038/nmeth.4324

(3) Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923

(4) Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology : a journal of computational molecular cell biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

(5) Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and applications. BMC Bioinformatics 10, 421 (2009). https://doi.org/10.1186/1471-2105-10-421
