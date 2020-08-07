# low_pass_WGS
Karyotyping using low coverage whole genome sequencing.

This pipeline utilises BWA (http://bio-bwa.sourceforge.net/) to map short paired-end sequencing data to a genome sequence to determine a sample's karyotpye. After postprocessing with samtools (http://samtools.sourceforge.net/), the data is further analysed using the QDNAseq R package (https://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html).

The input data should be trimmed, e.g. using `TrimGalore` or any other comparable quality trimming tool. Unpaired reads resulting from this should be retained and concatenated into one additional read file, e.g. `<sample_name>_unpaired.fastq.gz`.
