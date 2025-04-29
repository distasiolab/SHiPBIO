# Preparing GWAS data for single cell analysis with [SCDRS](https://github.com/martinjzhang/scDRS)


- Download GWAS data from GWAS Catalog FTP: (e.g. http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003219/Fritsche-26691988.txt.gz)
- Unzip tsv file, rename to **.pval

- To make a n_total column, summing the
awk -F'\t' 'BEGIN {OFS = FS} NR == 1 {print $0, "n_total"} NR > 1 {print $0, $25+$26}' **.pval > out.pval
or
awk -F'\t' 'BEGIN {OFS = FS} {RS=ORS="\r\n"} NR == 1 {print $0, "n_total"} NR > 1 {print $0, $c1 +$c2}' c1="n_cas" c2="n_con" **.pval > out.pval

- For Chromosome position to rsid:
      - awk 'BEGIN {FS="\t"; OFS=""} {print "select chrom,chromStart,name from snp138 where chrom = \"chr",$1,"\" and chromStart + 1 = ",$2,";"}' Alzheimers_GCST013197_PMID_34493870.pval > snps.sql
      - mysql -h genome-mysql.cse.ucsc.edu -u genome -A -D hg19 --skip-column-names < snps.sql > rsids.txt

- To replace empty values with NA:
     - awk -i inplace -F'\t' -v OFS='\t' '{for(i=1; i<=NF; i++) if($i=="") $i="NA"; print}' **.tsv

- Use MAGMA to get zscores (RunMAGMA_001.sh), resulting in **.genes.out

magma --bfile [DATA] --gene-annot [ANNOT].genes.annot --pval [PVAL_FILE] N=[N]  
magma --bfile [DATA] --gene-annot [ANNOT].genes.annot --pval [PVAL_FILE] ncol=[N_COL]


- Delete all columns except GeneI (i.e. Entrez accession number) and Z-score (which you rename to the disease name): Use awk to select 'GENE' and 'ZSTAT' columns:
	   awk 'BEGIN { FS = OFS } { print $1, $8 }' **.genes.out > * genes.zscore.tsv # Use only columns 1 and 8 (indexing starts with 1)
	   awk -i inplace -v OFS="\t" '$1=$1' *.genes.zscore.tsv 
	   
Use https://biit.cs.ut.ee/gprofiler/convert to convert Entrez accession IDs to Gene Names
       - Target namespace is ENTREZGENE; Numeric IDs treated as ENTREZGENE_ACC
       	        - Download csv and replace commas with tabs to create 'gProfiler_hsapiens_*.tsv'
		awk -F ',' 'BEGIN {OFS="\t"} {$1=$1}1' gProfiler_hsapiens_*.csv > gProfiler_hsapiens_*.tsv

		- Remove quotes:
		cat gProfiler_hsapiens_**.tsv | tr -d '"' > gProfiler_hsapiens_**.tsv2
		dos2unix gProfiler_hsapiens_**.tsv2
		
       	 	- Do inner join to get gene names for each gene ID (accession #) in *.genes.stats.tsv
		~/Programs/tsv-utils/bin/tsv-join -H -f gProfiler_hsapiens_**.tsv2 -a name -k initial_alias -z -d GENE Glaucoma_GlobalBiobank_36777996_GCST90399726.genes.stats.tsv > COMBINED.tsv
		
		- Remove duplicate gene names in COMBINED:
		awk -i inplace '!seen[$3]++' COMBINED.tsv
		
		  	 
- Select and rename columns, resulting in **.genes.zscore.tsv (which just has column names 'GENE' and disease name (e.g. 'AMD', 'MS', 'Alzheimers')
	   awk -F'\t' 'BEGIN {OFS = FS} NR == 1 {print "GENE", "AMD"} NR > 1 {print $3, $2}' COMBINED.tsv > **.genes.zscore.tsv

- Run SCDRS Munge (RunSCDRS_Munge_001.sh) to get **.gs

      	   ~/miniconda3/envs/cellcharter-env/bin/scdrs munge-gs --out-file **.gs --zscore-file **.genes.zscore.tsv --weight zscore --n-max 1000

