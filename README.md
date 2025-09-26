Filtration Strategy:
Mutect, Strelka and Varscan were used to call SNVs.
Two or more callers called the SNV and VAF>0.1 and Read Depth>20
gnomAD_genome_ALL<0.1 & ExAC_ALL<0.1

Average VAF and Read Depth is calculated across all the callers that the mutation is called in and stored in the column name ‘vaf’ and ‘depth’.

There is a binary flag column added that indicates if the mutation was called in each caller (e.g., mutect_flag)
