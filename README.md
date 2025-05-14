# Practical Haplotype Graph v2.4+ Analysis Scripts

This repository contains scripts to assist with analyzing outputs from the Practical Haplotype Graph (PHG) version 2.4 and above.

The [PHG](https://github.com/maize-genetics/phg_v2), developed by the Buckler Lab, is a computational tool for creating, storing, and using a pangenome for genomic imputation. For more information and full documentation, please visit the [official PHG website](https://phg.maizegenetics.net/).

While this repository provides custom tools and workflows for PHG-based analysis, it is not officially affiliated with the PHG software or its authors. These scripts were developed for specific use cases but are shared here in the hope that they may be useful for others working with PHG data.



## Contents

- `scripts/imputed_parents_merged_vcf.py`: Generates a merged VCF of imputed genotypes using a merged parent VCF and imputed parents files output from `phg find-paths`.
- More scripts to come.



## Script: `imputed_merged_vcf.py`

This script merges a PHG-derived merged parent VCF with a directory of imputed parents files, producing a final genotype matrix where each sample has one genotype call per reference range.



### Step 1: Generate the Merged Parent VCF

Before running `imputed_parents_merged_vcf.py`, first use the PHG tool `merge-gvcfs` to combine all founder GVCFs into a single VCF.

**Note:** The reference genome can be used during imputation by adding its hVCF to your vcf_files directory, but the reference gVCF should be excluded from merging for this script. See script options below.

```
export JAVA_OPTS="-Xmx12g"

phg merge-gvcfs \
    --input-dir phg_v2/output/vcf_files \
    --output-file phg_v2/merged_parents/merged_parents.vcf
```



### Step 2 (Optional): Filter Variants with `bcftools`

Filtering can improve downstream analysis performance and focus results on biallelic, informative SNPs.

#### Keep only:
- **Biallelic SNPs**
- **Sites with ALT allele in ≥2 and ≤36 samples**
- **Sites with no missing data**

```
bcftools view \
    -m2 -M2 \
    --types snps \
    -c 2 -C 36 \
    -Ov -o phg_v2/merged_parents/filtered_merged_parents.vcf \
    phg_v2/merged_parents/merged_parents.vcf

bcftools view -g ^miss \
    phg_v2/merged_parents/filtered_merged_parents.vcf \
    > phg_v2/merged_parents/non_miss_filtered_merged_parents.vcf
```



### Step 3: Run the Script

```
python3 imputed_parents_merged_vcf.py \
    --ref_ranges_file phg_v2/output/ref_ranges.bed \
    --merged_parents_vcf_path phg_v2/merged_parents/merged_parents.vcf \
    --out_parents_dir phg_v2/output/read_mappings/vcf_files
```

Optional:
```
    --reference_sample_name B73 \
    --merged_imputed_vcf_path /path/to/output/merged_imputed.vcf
```



### Script Features

- Designed for **homozygous/haploid imputation** output from phg find-paths, using the --path-type haploid option and specifying --out-parents-dir.
- Marks missing haplotypes as "." and reference alleles as "0".
- If a reference sample was used during imputation, provide the reference sample name as it appears in the imputed parents files using the optional --reference_sample_name argument.



### Performance Notes

- Tested on a 1.7 GB merged parent VCF (37 parents, 16 million variants) with 2,000 imputed samples.
- Max RAM usage: ~70 GB.
- Runtime: ~6 hours on a high-memory server.



### Dependencies

- Python 3.8+
- Standard libraries: `gzip`, `os`, `argparse`, `glob`
- [bcftools](https://samtools.github.io/bcftools/) for optional pre-processing



## Future Scripts

### Step 4 (Upcoming): Calculate Founder Contribution

Planned features:
- Quantify founder haplotype contributions for each sample.
- Generate summary tables and genome-wide visualizations.
- Tools for specific founder panels (e.g., Zea mays).



## Acknowledgements

This project was developed by **Micah Kelleher** as part of work in the **Baxter Lab** at the **Donald Danforth Plant Science Center**. The code reflects workflows tailored to *Zea mays* analysis and PHG-based imputation.

Special thanks to:

- **Dr. Ivan Baxter & Lab** – For providing the environment, resources, and insights critical to this project.  
- **Dr. Sherry Flint-Garcia** – For her guidance and for providing the *Zea mays* population used in development and testing.  
- **Dr. Jeff Ross-Ibarra** – For his help with performance benchmarking and interpreting PHG outputs.  
- **The Buckler Lab** – For creating and maintaining the PHG, and for helpful guidance on usage.  
- **Tim Kosfeld** – For foundational discussions and early support for this work.
