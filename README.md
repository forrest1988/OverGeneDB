# OverGeneDB - Overlapping Genes Detection

Scripts located in OverlappingGenesDetection.py file were used to identify overlapping gene pairs in human and mouse genomes. Detection is based on the alternative transcription start sites (TSSs) coordinates from DBTSS and reference genome annotations from RefSeq. Scripts were written in Python 2.7.

Run example:
An example of how to identify overlapping genes in human testis with the minimal expression level set to 5 ppm was implemented in OverlappingGenesDetection.py. It may be run using the following command:
python OverlappingGenesDetection.py

For that purpose three input files, located in "InputFiles" folder, are required:
- DBTSS_Human_AdultTestis.tab - TSS details file downloaded from the DBTSS website. File is in Tab Separated Value format (TSV).
- Human_GeneNames.tsv - List of all human transcripts and their associated gene names from RefSeq. Example file was provided for human hg38 genome reference. File is in Tab Separated Value format (TSV).
- Human_RefSeqTranscripts.tsv - Human RefSeq transcripts in BED12 format, where the columns are as follows:
		1. chromosome
		2. transcript start
		3. transcript end
		4. transcript accession number
		5. score
		6. strand
		7. CDS start
		8. CDS end
		9. color
		10. number of exons
		11. exon sizes
		12. exon starts (in relation to transcript start)

Example output:
OverlappingGenePairs_Human_AdultTestis_5ppmCUTOFF.tsv is the example output file for human testis.
		
More TSS data:
Additional Raw TSS data may be downloaded from DBTSS FTP server at ftp://ftp.hgc.jp/pub/hgc/db/dbtss/. Data for six mouse tissues - brain, heart, kidney, liver, spleen and thymus - were stored at http://rhesus.amu.edu.pl/ovrGenes/dbtss. 
