# Project Overview
This project's aim is to identify high-confidence germline INDEL and SNPs variant calls, following GATK best practices for human whole-genome data.

The dataset consists of human WGS paired-end Illumina reads (2×100 bp) from the 1000 Genomes Project (sample ID: HG00096), sequenced at ~30× coverage. Two paired-end FASTQ files were downloaded for the sample HG00096:
- SRR062634_1.filt.fastq.gz
- SRR062634_2.filt.fastq.gz

The reads were generated on an Illumina HiSeq 2000 platform with a read length of 100 bp. The estimated mean coverage is 30× across the genome.

Reads were aligned to the GRCh38/hg38 reference genome using BWA-MEM, followed by sorting, duplicate marking, and known variant resources from GATK’s hg38 bundle (dbSNP build 138) were used for base recalibration and variant annotation. Variants were called using GATK HaplotypeCaller and filtered based on standard hard-filtering parameters. 

The resulting SNPs and indels were evaluated using bcftools statistics and annotated with Ensembl VEP to assess potential functional consequences.

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

| Metric                             | Result    |
| ---------------------------------- | --------- |
| TOTAL_READS                | 48,297,986 |
| % of reads that successfully aligned                | 99.65%|
| Fraction of properly paired reads               | 99.8% |
| Fraction of mismatched bases                 | 0.5% |
| Mean read length                | 100bp|
| Fraction of chimeric reads                | 0.17% |
| Base-level error rate in aligned reads                 | 0.45% |
| PERCENT_DUPLICATION                | 0.9% |
| ESTIMATED_LIBRARY_SIZE (number of unique DNA fragments)       | 1,383,600,948 |
| Base-level error rate in aligned reads                 | 0.45% |
| Total raw variants                 | 1,091,729 |
| Number of samples |1|
| Number of High-Quality SNPs         | 1,384     |
| number of multiallelic sites  | 27 |
| number of multiallelic SNP sites | 11 |
| Transition-to-Transversion ratio | 1.80 |
| Indels NUmber of high-quality indels  | 967  |
| Mean depth                         | 31.4×     |

✅ ~99.6% of reads aligned → your reference genome matches your sequencing sample well.

✅ Low mismatch & error rates → good sequencing quality.

✅ Proper pairing >99% → library prep and read pairing are consistent.

✅ Low chimera and adapter rates → minimal contamination or sequencing artifacts.

⚠️ 2.3% improper pairs → slightly above perfect but completely normal (typically under 5% is fine).

Overall, this alignment is clean, high-quality, and ready for variant calling.
Alignment quality was high, with 99.6% of reads successfully mapped to the hg38 reference and a low mismatch rate of 0.5%. Proper pairing exceeded 99%, suggesting good library preparation and sequencing quality.

- *multiallelic site A genomic position with more than one alternate allele, regardless of type (SNP, indel, etc.). Example REF = A; ALT = C,CG*
- *Multiallelic SNP site - Specifically, a multiallelic site where all alternates are SNPs (single base substitutions). Example REF = A; ALT = C,G*
- *Transitions = A↔G or C↔T (changes between similar base types: purine↔purine or pyrimidine↔pyrimidine). Transversions = A↔C, A↔T, G↔C, or G↔T (purine↔pyrimidine).*
- *Transition-to-Transversion ratio(ts/tv) Ratio of transitions to transversions. Only Snps have ts/tv trations*


The Ts/Tv ratio = 1.08, which is lower than expected for high-quality human data.
Typically:

~2.0–3.0 for high-quality germline calls.

~1.0–1.2 for low-quality or somatic calls, or poorly filtered data.

So this might mean:

The dataset has some false positives.

Or, it’s not human (some non-model organism).

Or, you’re still before final filtering.
A false positive means the variant caller called a position as a variant, but in reality, it’s not a true biological change.
It’s usually due to:

Sequencing errors (e.g., low-quality bases or reads).

Misalignment (reads mapped to the wrong region).

Poor filtering (keeping low-quality calls).

Here you present your insights — ideally in 3–5 bullet points or a short paragraph. Example:

Total variants detected: 3,412,562

SNPs: 3.2 M

Indels: 200 K

Mean sequencing depth: ≈40×

Mean mapping quality (MQ): ~55, indicating confident alignments.

~90 % of variants passed all quality filters (QD > 2, MQ > 40, SOR < 4).

Most variants were intergenic or intronic (non-coding).

A small subset (~1 %) annotated as upstream_gene_variant, located near MT-ND1, MT-ND2, etc.


### annotation 
indel:
Variants processed	967
Variants filtered out	0
Novel / existing variants	379 (39.2) / 588 (60.8)
Overlapped genes	144
Overlapped transcripts	1042
Overlapped regulatory features	26

SNP

# Biological Interpretation

Demonstrate that you can translate data into biological meaning:

The majority of high-confidence variants were non-coding, suggesting no direct impact on protein structure. However, several upstream variants were detected in mitochondrial genes involved in the oxidative phosphorylation pathway (MT-ND1, MT-ND2). While likely benign, regulatory effects cannot be ruled out without functional validation.

# 📈 5. Visualization

Include 1–2 visuals if possible (or commands to generate them):

Read coverage histogram or depth distribution (from samtools or GATK metrics).

Variant quality (QD or DP) histogram.

Consequence pie chart (e.g., intergenic, missense, synonymous).

IGV screenshot showing an example variant.

These visuals make your report more portfolio-ready and immediately readable.

# ⚙️ 6. Computational & Reproducibility Notes

Recruiters and engineers will appreciate this:

Pipeline was bash-scripted and follows GATK Best Practices.

Uses environment variables and modular steps for reproducibility.

Compatible with conda/Docker (mention if you containerized).

Optional: brief mention of compute environment (e.g., “Ubuntu 22.04, 8 cores, 32 GB RAM”).

# Reflection / Next Steps

A short, honest reflection adds maturity:

This project helped reinforce my understanding of the variant calling workflow, particularly around data preprocessing and variant quality filtration.
Future improvements include integrating automated annotation (Funcotator/SnpEff), building a Jupyter report for summary statistics, and containerizing the entire workflow with Docker or Nextflow.
