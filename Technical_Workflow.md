# ✅ Variant Calling Preperation
This script calls the germline variants in a human WGS paired end reads 2X100bp

## 1. Download data 
```
wget -P /home/dmatute/Projects/Variant_Calling/data ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P /home/dmatute/Projects/Variant_Calling/data ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
```

## 2. Preprare Files 

### 2.1. Download reference genome 
```
wget - P /home/dmatute/Projects/Variant_Calling/data https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/Desktop/demo/supporting_files/hg38/hg38.fa.gz
```

### 2.2. Index reference genome 
The reference is indexed to increase the speed and efficiency of aligning fragmented genomic sequence.
```
samtools faidx /home/dmatute/Projects/Variant_Calling/data/hg38.fa 
```

### 2.3. Create a sequence dictionary for a reference sequence
This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, <mark>  which is required by many processing and analysis tools</mark>. The output file contains a header but no SAMRecords, and the header contains only sequence records.

[gatk CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360036712531-CreateSequenceDictionary-Picard)
```
gatk CreateSequenceDictionary R=/home/dmatute/Projects/Variant_Calling/data/hg38.fa O=/home/dmatute/Projects/Variant_Calling/data/hg38.dict
```

### 2.4. Download known sites files for BQSR from GATK resource bundle
Known variant from GATK’s hg38 bundle (dbSNP build 138) will be used for base recalibration and variant annotation. 
```
wget -P . https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P . https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

## 4. Read Quality Control
```
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/
```
No trimming required, quality looks okay.

# ✅ Variant Calling

Constant Variables
|variable|path|
|---|---|
|ref|"/home/dmatute/Projects/Variant_Calling/data/hg38.fa"|
|known_sites|"/home/dmatute/Projects/Variant_Calling/data/known_varinats/Homo_sapiens_assembly38.dbsnp138.vcf"|
|aligned_reads|"/home/dmatute/Projects/Variant_Calling/aligned_reads"|
|reads|"/home/dmatute/Projects/Variant_Calling/reads"|
|results|"/home/dmatute/Projects/Variant_Calling/results"|
|data |"/home/dmatute/Projects/Variant_Calling/data"|

## 5. Map to reference using BWA-MEM
[BWA Documentation](https://bio-bwa.sourceforge.net/bwa.shtml)
### 5.1. Index Reference with BWI
- Index the reference genome for fast BWA alignment.
```
bwa index ${ref}
```
- The following files are produced:

    |File|Purpose|
    |---|---|
    |.bwt	|Burrows–Wheeler Transform — a compressed, searchable version of your genome.|
    |.pac	|The packed (binary) version of the genome sequence.|
    |.ann, .amb	|Annotation and ambiguity files that store sequence names and ambiguous bases (N’s).|
    |.sa	|Suffix array — helps locate matches quickly.|

### 5.2. BWA alignment
- The BWA mode selected is *mem*; it is the latest, fastest and most accrurare. 
- The algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).
- The BWA-MEM algorithm performs local alignment.
```
bwa mem -t 8 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam
```
<img src="https://github.com/DanyMatute/Variant_Calling/blob/main/Screenshot%202025-10-23%20192954.png" alt="BWA_MAM">

[GPU Acceleration of BWA-MEM DNA Sequence Alignment](https://doi.org/10.13140/RG.2.2.33045.19687)

## 6. SAM/BAM file manipulation
```
samtools view -bS aligned_reads/SRR062634.paired.sam > aligned_reads/SRR062634.bam
samtools sort -o aligned_reads/SRR062634_sorted.bam aligned_reads/SRR062634.bam
samtools index aligned_reads/SRR062634_sorted.bam
```

## 7. Mark Duplicates and Sort with GATK4 
[Markduplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360057439771-MarkDuplicates-Picard)
```
gatk MarkDuplicates   -I SRR062634_sorted.bam   -O SRR062634_dedup_reads.bam   -M SRR062634_metrics.txt
```
- This tool locates and tags **duplicate reads** in a BAM or SAM file, where duplicate reads are defined as originating from a **single fragment of DNA.**
- Duplicates can come from:
    1. sample preparation
    2. from simple amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as **optical duplicates.**
- ```MarkDuplicates``` tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file
- After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).
- The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. it does not identify the type of duplicate.

<details> 
<summary>Metrics.txt</summary>

|Field	|Description|
|---|---|
|LIBRARY	|The library on which the duplicate marking was performed.|
|UNPAIRED_READS_EXAMINED	|The number of mapped reads examined which did not have a mapped mate pair, either because the read is unpaired, or the read is paired to an unmapped mate.|
|READ_PAIRS_EXAMINED	|The number of mapped read pairs examined. (Primary, non-supplemental)|
|SECONDARY_OR_SUPPLEMENTARY_RDS|	The number of reads that were either secondary or supplementary|
|UNMAPPED_READS|	The total number of unmapped reads examined. (Primary, non-supplemental)|
|UNPAIRED_READ_DUPLICATES|	The number of fragments that were marked as duplicates.|
|READ_PAIR_DUPLICATES	|The number of read pairs that were marked as duplicates.|
|READ_PAIR_OPTICAL_DUPLICATES|	The number of read pairs duplicates that were caused by optical duplication. Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.|
|PERCENT_DUPLICATION	|The fraction of mapped sequence that is marked as duplicate.|
|ESTIMATED_LIBRARY_SIZE|	The estimated number of unique molecules in the library based on PE duplication.|

</details>

## 8. Base quality recalibration 
[gatk BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)
- Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. 
- This allows us to get more accurate base qualities overall, which in turn improves the accuracy of our variant calls. 
- we can't correct the base calls themselves, but we can at least tell the variant caller more accurately how much it can trust a base. 
- The base recalibration process involves two key steps: 
    - first the ```BaseRecalibrator``` tool builds a model of covariation based on the input data and a set of known variants, producing a recalibration file; 
  	- then the ```ApplyBQSR``` tool adjusts the base quality scores in the data based on the model, producing a new BAM file.
### 8.1. Build the Model
```
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
```
### 8.1. Apply the model to adjust the quality scores
```
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam
```

## 9. Collect Alignment & Insert Size Metrics
- **CollectAlignmentSummaryMetrics**
	- Produces a summary of alignment metrics from a SAM or BAM file.
 	- This tool takes a SAM/BAM file input and produces metrics detailing the *quality of the read alignments* as well as the proportion of the reads that *passed machine signal-to-noise threshold quality filters*.
- **CollectInsertSizeMetrics**
	- This tool provides useful metrics for *validating library construction* including the insert size distribution and read orientation of paired-end libraries.
 	- The CollectInsertSizeMetrics tool outputs the *percentages of read pairs* in each of the three orientations (FR, RF, and TANDEM) as a histogram.
  	- In addition, the insert size distribution is output as both a histogram (.insert_size_Histogram.pdf) and as a data table (.insert_size_metrics.txt).
```
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf
```

# ✅ Call Variants 
## 12. HaplotypeCaller 
- The HaplotypeCaller identify single nucleotide polymorphisms (SNPs) and insertions/deletions (Indels) from genome alignment.
```
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf
```
<img width="1002" height="641" alt="image" src="https://github.com/user-attachments/assets/5f668294-81ab-488a-be69-e16d5e96eea6" />

## 11. Extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf

# ✅ Variant Analysis
## 12. Quality Check the aligned raw_variants, snps & indels vcf files.
Parses VCF or BCF through ```bcftools stats``` and to obtain alignment metrics, these stats can then be graphed with ```plot-vcfstats```. 
```
bcftools stats raw_variants.vcf > qc/raw_variants_stats.txt
bcftools stats raw_indels.vcf > qc/raw_indels_stats.txt
bcftools stats raw_snps.vcf > qc/raw_snps_stats.txt
```

## 13. Filtering
### 13.1 Filter SNPs & Indels
[gatk VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)

```
# Filter SNPs
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"

# Filter INDELS
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"
```
### 13.2 Select Variants that PASS filters
```
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis-ready-snps.vcf

gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis-ready-indels.vcf
```

### 13.3 Exclude variants that failed genotype filters. 
```
cat analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > analysis-ready-snps-filteredGT.vcf
cat analysis-ready-indels.vcf| grep -v -E "DP_filter|GQ_filter" > analysis-ready-indels-filteredGT.vcf
```

## 14. Quality Check the filtered snps & indels vcf
Parses VCF or BCF and produces text file stats which is suitable for machine processing and can be plotted using plot-vcfstats. 

```
bcftools stats ../filtering/analysis-ready-snps-filteredGT.vcf > analysis-ready-snps-filteredGT_stats.txt

bcftools stats ../filtering/analysis-ready-indels-filteredGT.vcf > analysis-ready-indels-filteredGT_stats.txt

plot-vcfstats -p filtered_snps_visual analysis-ready-snps-filteredGT_stats.txt 

plot-vcfstats -p filtered_indels_visual analysis-ready-indels-filteredGT_stats.txt 
```

# ✅ Annotation
## 16.1 VEP
- The [Ensembl](https://www.ensembl.org/Homo_sapiens/Tools/VEP?tl=UtCOFFs56h7fKewR-11380108) can do free VEP, this was used, and bellow is the command that the website used. 
```
./vep --af --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --polyphen b --pubmed --regulatory --show_ref_allele --sift b --species homo_sapiens --symbol --transcript_version --tsl --uploaded_allele --cache --input_file [input_data] --output_file [output_file]

./vep --af --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --polyphen b --pubmed --regulatory --show_ref_allele --sift b --species homo_sapiens --symbol --transcript_version --tsl --uploaded_allele --cache --input_file [input_data] --output_file [output_file]
```

## 16.2 Annotation Annotation Summary 
```
results/visuaization$ grep -v '^#' ../annotation/VEP/annotation_VEP_indel.vcf_ | cut -f8 | sort | uniq -c | sort -nr | head > indel_vcf_summary.txt
results/visuaization$ grep -v '^#' ../annotation/VEP/annotation_VEP_snps.vcf_ | cut -f8 | sort | uniq -c | sort -nr | head > snp_vcf_summary.txt
```
<table border="1" cellspacing="0" cellpadding="4">
<caption>Example indel_vcf_summary.txt</caption>
  <thead>
    <tr>
      <th>AC</th>
      <th>AF</th>
      <th>AN</th>
      <th>DP</th>
      <th>ExcessHet</th>
      <th>FS</th>
      <th>MLEAC</th>
      <th>MLEAF</th>
      <th>MQ</th>
      <th>QD</th>
      <th>SOR</th>
      <th>CSQ</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2</td>
      <td>1.00</td>
      <td>2</td>
      <td>61</td>
      <td>0.0000</td>
      <td>0.000</td>
      <td>1</td>
      <td>0.500</td>
      <td>41.84</td>
      <td>29.55</td>
      <td>0.941</td>
      <td>TCTCTCTC | intergenic_variant | MODIFIER | ... | G/GTCTCTCTC</td>
    </tr>
    <tr>
      <td>2</td>
      <td>1.00</td>
      <td>2</td>
      <td>51</td>
      <td>0.0000</td>
      <td>0.000</td>
      <td>2</td>
      <td>1.00</td>
      <td>55.44</td>
      <td>34.15</td>
      <td>0.836</td>
      <td>- | intergenic_variant | MODIFIER | ... | TTC/T</td>
    </tr>
    <tr>
      <td>2</td>
      <td>1.00</td>
      <td>2</td>
      <td>50</td>
      <td>0.0000</td>
      <td>0.000</td>
      <td>1</td>
      <td>0.500</td>
      <td>55.54</td>
      <td>23.68</td>
      <td>0.741</td>
      <td>- | intergenic_variant | MODIFIER | ... | CT/C</td>
    </tr>
  </tbody>
</table>

<details> 
<summary>Field Key</summary>

| Field         | Meaning               | Interpretation                                                                                                                           |
| ------------- | --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| **AC**        | Allele Count          | Number of alternate alleles observed in called samples. Since your data shows `AC=2`, it likely means both alleles (homozygous variant). |
| **AF**        | Allele Frequency      | Frequency of the alternate allele. `AF=1.00` → variant is fixed (homozygous ALT).                                                        |
| **AN**        | Allele Number         | Total number of alleles in all called genotypes (usually 2 per diploid sample).                                                          |
| **DP**        | Depth                 | Total read depth at the locus. Higher = more confidence. Your variants have 30–350 reads → good coverage.                                |
| **MQ**        | Mapping Quality       | Average mapping quality of reads supporting the call. 40–60 = reliable alignments.                                                       |
| **QD**        | Quality by Depth      | Variant confidence normalized by depth. >2 is generally good.                                                                            |
| **SOR**       | Strand Odds Ratio     | Detects strand bias; near 1 = balanced; >3 = possible strand bias.                                                                       |
| **FS**        | Fisher Strand         | Another strand bias measure; `0.000` = no bias detected.                                                                                 |
| **ExcessHet** | Excess Heterozygosity | High values may indicate sequencing or alignment errors; your values are all `0.0000` → fine.                                            |
| **CSQ**       | Consequence           | Annotation from VEP or Funcotator: shows gene impact, type (intergenic, upstream, downstream, etc.).                                     |

</details>

- Most variants are **intergenic** meaning they  fall outside of any known gene, likely non functional or neutral. But there are some up/downstreem variants, influencing gene regulation. 
- the variants are around mitochondrial genome variants/
- if "Modifier" then the variation will cause no coding change 
- if you want to find functionally significant variant look for these words **“missense_variant”, “frameshift_variant”, “stop_gained”**
  
