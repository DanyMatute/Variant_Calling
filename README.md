# Project Overview
🎯 This project's aim is to identify high-confidence germline INDEL and SNPs variant calls, following GATK best practices for human whole-genome data.

- The dataset consists of human WGS paired-end Illumina reads (2×100 bp) from the 1000 Genomes Project (sample ID: HG00096), sequenced at ~30× coverage.
- Two paired-end FASTQ files were downloaded for the sample HG00096:
  - SRR062634_1.filt.fastq.gz
  - SRR062634_2.filt.fastq.gz
- The reads were generated on an Illumina HiSeq 2000 platform with a read length of 100 bp. The estimated mean coverage is 30× across the genome.

Reads were aligned to the GRCh38/hg38 reference genome using BWA-MEM, followed by sorting, duplicate marking, and known variant resources from GATK’s hg38 bundle (dbSNP build 138) were used for base recalibration and variant annotation. Variants were called using GATK HaplotypeCaller and filtered based on standard hard-filtering parameters. The resulting SNPs and indels were evaluated using bcftools statistics and annotated with Ensembl VEP to assess potential functional consequences.

~3.5 M variants were identified; most were intergenic, with a few near mitochondrial genes.

# Data & Methods Summary 
| Step            | Tool                                | Key Purpose                   |
| --------------- | ----------------------------------- | ----------------------------- |
| Read QC         | FastQC                              | Assess raw read quality       |
| Alignment       | BWA-MEM                             | Map reads to reference (hg38) |
| Post-processing | Samtools, GATK MarkDuplicates, BQSR | Improve mapping accuracy      |
| Variant Calling | GATK HaplotypeCaller                | Identify SNPs & Indels        |
| Filtering       | GATK VariantFiltration              | Remove low-quality variants   |
| Annotation      | VEP                     | Predict variant consequence   |

# Key Analytical Results
## Alignment Summary Metrics
| Metric                             | Result    |
| ---------------------------------- | --------- |
| TOTAL_READS                | 48,297,986 |
| % of reads that successfully aligned                | 99.65%|
| Fraction of properly paired reads               | 99.8% |
| Fraction of mismatched bases                 | 0.5% |
| Mean read length                | 100bp|
| Fraction of chimeric reads                | 0.17% |
| Base-level error rate in aligned reads                 | 0.45% |

- ~99.6% of reads aligned → your reference genome matches your sequencing sample well.
- Low mismatch & error rates → good sequencing quality.
- Proper pairing >99% → library prep and read pairing are consistent.
- Low chimera and adapter rates → minimal contamination or sequencing artifacts.
- 2.3% improper pairs → slightly above perfect but completely normal (typically under 5% is fine).
- **Overall, this alignment is clean, high-quality, and ready for variant calling.**

## DuplicationMetrics
| Metric                             | Result    |
| ---------------------------------- | --------- |
| PERCENT_DUPLICATION                | 0.9% |
| ESTIMATED_LIBRARY_SIZE (number of unique DNA fragments)       | 1,383,600,948 |
| Base-level error rate in aligned reads                 | 0.45% |

<p>
  <img src="https://github.com/DanyMatute/Variant_Calling/blob/main/Insert_size_histogram.png" alt="Alt Text" style="width:50%; height:auto; ">
  
  <em>The insert size distribution for SRR062634_sorted_dedup_bqsr_reads.bam shows a single peak around ~180 bp with a narrow standard deviation, indicating consistent library preparation and proper fragment length. No secondary peaks or abnormal broadening were observed.</em>
</p>



## Variant Calling Matrics
| Metric                             | Result    |
| ---------------------------------- | --------- |
| Total raw variants                 | 1,091,729 |
| Number of samples |1|
| Number of High-Quality SNPs         | 1,384     |
| number of multiallelic sites  | 27 |
| number of multiallelic SNP sites | 11 |
| Transition-to-Transversion ratio | 1.80 |
| Indels NUmber of high-quality indels  | 967  |

<details>
  
  - *multiallelic site A genomic position with more than one alternate allele, regardless of type (SNP, indel, etc.). Example REF = A; ALT = C,CG*
  - *Multiallelic SNP site - Specifically, a multiallelic site where all alternates are SNPs (single base substitutions). Example REF = A; ALT = C,G*
  - *Transitions = A↔G or C↔T (changes between similar base types: purine↔purine or pyrimidine↔pyrimidine). Transversions = A↔C, A↔T, G↔C, or G↔T (purine↔pyrimidine).*
  - *Transition-to-Transversion ratio(ts/tv) Ratio of transitions to transversions. Only Snps have ts/tv trations*

</details>

**The Ts/Tv ratio = 1.08, which is lower than expected for high-quality human data.**
  
    Typically:
    - ~2.0–3.0 for high-quality germline calls.
    - ~1.0–1.2 for low-quality or somatic calls, or poorly filtered data.
    
    So this might mean:
    - The dataset has some false positives.
      - A false positive means the variant caller called a position as a variant, but in reality, it’s not a true biological change. 
      - It’s usually due to: 
        - Sequencing errors (e.g., low-quality bases or reads), 
        - Misalignment (reads mapped to the wrong region), 
        - Poor filtering (keeping low-quality calls).
    - Or, it’s not human (some non-model organism).
    - Or, you’re still before final filtering, but we have done filtering
<p>
  <img src="https://github.com/DanyMatute/Variant_Calling/blob/main/snps_substitution_types.png" alt="Alt Text" style="width:30%; height:auto; ">
  <img src="https://github.com/DanyMatute/Variant_Calling/blob/main/indel_distribution.png" alt="Alt Text" style="width:35%; height:auto; ">
  <img src="https://github.com/DanyMatute/Variant_Calling/blob/main/snps_ts-tv_statified.png" alt="Alt Text" style="width:30%; height:auto; ">
  <em></em>
</p>

## Annotation Metrics

| Indel Metric                             | Result    |
| ---------------------------------- | --------- |
|Variants processed|	967|
|Variants filtered out|	0|
|Novel / existing variants|	379 (39.2) / 588 (60.8)|
|Overlapped genes	|144|
|Overlapped transcripts|	1042|
|Overlapped regulatory features|	26|

# Discussion/Conclusion

The majority of high-confidence variants were non-coding, suggesting no direct impact on protein structure. However, several upstream variants were detected in mitochondrial genes involved in the oxidative phosphorylation pathway (MT-ND1, MT-ND2). While likely benign, regulatory effects cannot be ruled out without functional validation.

This small project demonstrates the variant calling workflow, particularly around data preprocessing and variant quality filtration. Future improvements could include integrating automated annotation (Funcotator/SnpEff), building a Jupyter report for summary statistics, and containerizing the entire workflow with Docker or Nextflow.
