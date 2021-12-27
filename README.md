# circos
R script to make Circos figure based on BioCircos package

## Usage
### create and activate circos in anaconda environment and install BioCircos via conda
1. conda create circos
2. conda activate circos
3. conda install -c r r-biocircos
4. conda install -c r r-htmlwidgets

### To run with sample data
```
Rscript script.R --args vcf_snp_file="test_data/snp_gatk_mutect2.vcf" bedpe_SV_file="test_data/sv_manta.bedpe" cnv_file="test_data/cnv_accucopy.tsv"
```
```
Rscript script.R --args vcf_snp_file="test_data/snp_gatk_mutect2.vcf" bedpe_SV_file="test_data/sv_manta.bedpe" cnv_file="test_data/cnv_accucopy.tsv" AF_label=0.2 DP_label1=30 DP_label2=50 tumor_vcf_column=10 
```

### Parameters
1. vcf_snp_file: variant file (e.g. from GATK pipeline)
2. bedpe_SV_file: structural variation file (e.g. from manta software)
3. cnv_file: copy number variation file (e.g. from accucopy software)
4. AF_labe: minimum allelic fraction threshold for highlighting snp (default: 0.2)
5. DP_label1: minimum read depth threshold for highlighting snp in red color (default: 30)
6. DP_label2: minimum read depth threshold for highlighting snp in green color (default: 50)
7. tumor_vcf_column: column in VCF file that contain tumor read depth information (default: 10) 
8. out: output file name
