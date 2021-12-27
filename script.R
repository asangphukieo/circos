#######################
# Developer: Apiwat Sangphukieo
# Email: apiwat.sang@cmu.ac.th
# 
# ref: https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html#multi-track-example
#######################

library(BioCircos)
library(htmlwidgets)


###### Usage  ######
### create and activate conda circos and install BioCircos via conda
#conda create circos
#conda activate circos
#conda install -c r r-biocircos
#conda install -c r r-htmlwidgets

#to run:
#Rscript script.R --args vcf_snp_file="test_data/snp_gatk_mutect2.vcf" bedpe_SV_file="test_data/sv_manta.bedpe" cnv_file="test_data/cnv_accucopy.tsv"
#Rscript script.R --args vcf_snp_file="test_data/snp_gatk_mutect2.vcf" bedpe_SV_file="test_data/sv_manta.bedpe" cnv_file="test_data/cnv_accucopy.tsv" AF_label=0.2 DP_label1=30 DP_label2=50 tumor_vcf_column=10 

##############################################################
##############################################################

#Note: use only variants in vcf file that contain "PASS" only
###1.1 vcf for SNP file (from GATK pipeline)
###1.2 bedpe file for structural variants (from manta)
###1.3 vcf file for copy number alteration (from accucopy])

args <- commandArgs(trailingOnly = TRUE)

vcf_snp_file="INPUT_circos/OS0130_fresh_pass.vcf"
bedpe_SV_file="INPUT_circos/somaticSV_fresh.bedpe"
cnv_file="INPUT_circos/cnv.accucopy_fresh.tsv" 
AF_label=0.2    #minimum allelic fraction threshold for highlighting snp
DP_label1=30    #minimum read depth threshold for highlighting snp (red)
DP_label2=50    #minimum read depth threshold for highlighting snp (green)
tumor_vcf_column=10 #column in VCF file that contain tumor read depth information of tumor 
out="plot"

for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

##########################
#  load file in  #
##########################

vcf_snp=read.table(vcf_snp_file, comment="#",header=F)
vcf_snp = vcf_snp[vcf_snp[,7]=="PASS",] #use only "PASS" variants

bedpe_SV=read.table(bedpe_SV_file, comment="#",header=F)
bedpe_SV = bedpe_SV[bedpe_SV[,12]=="PASS",] #use only "PASS" variants
print(table(bedpe_SV[,11])) #check SV type

#copy number from Accucopy
cnv=read.table(cnv_file, comment="#",header=T)


##########################
# optional:
# select candidate SNP with 30X and 50x read depth and 0.2 AF , in OS0130_Tissue_tumor in column 10
##########################
vcf_snp[,c("AF","DP")]=do.call(rbind,strsplit(as.character(vcf_snp[,as.numeric(tumor_vcf_column)]),split=":"))[,3:4]
snp_AF= vcf_snp[vcf_snp[,"AF"]>=AF_label,] #apply AF filter >=0.2
snp_30x=snp_AF[snp_AF[,"DP"]>=DP_label1,] #apply DP filter >=30
snp_50x=snp_AF[snp_AF[,"DP"]>=DP_label2,] #apply DP filter >=50


##########################
#  SNP track  #
##########################

SNP_chr=as.character(gsub("chr","",vcf_snp[,1]))
SNP_pos=round(vcf_snp[,2])
SNP_col=runif(length(SNP_chr), 2, 12) #random snp vertical position on the track

snp_minRadius=0.55
snp_maxRadius=0.70
tracks = BioCircosSNPTrack("testSNP2", SNP_chr, 
  SNP_pos, SNP_col, color = "#474747", 
  minRadius = snp_minRadius, maxRadius = snp_maxRadius, range = c(2,12), borderColors = "#808080",opacities = 0.6,size = 1)


# Overlap point of interest on previous track, fix range to use a similar scale
SNP_chr=as.character(gsub("chr","",snp_30x[,1]))
SNP_pos=round(snp_30x[,2])
SNP_col=runif(length(snp_30x), 2, 12) #random snp vertical position on the track

tracks = tracks + BioCircosSNPTrack("testSNP3", SNP_chr, 
  SNP_pos, SNP_col, minRadius = snp_minRadius, maxRadius = snp_maxRadius, size = 2,range = c(2,12), opacities = 0.8,color="#fa7f7f")

### 50x
SNP_chr=as.character(gsub("chr","",snp_50x[,1]))
SNP_pos=round(snp_50x[,2])
SNP_col=runif(length(snp_50x), 2, 12) #random snp vertical position on the track

tracks = tracks + BioCircosSNPTrack("testSNP3", SNP_chr, 
  SNP_pos, SNP_col, minRadius = snp_minRadius, maxRadius = snp_maxRadius, size = 2,range = c(2,12), opacities = 0.8, color="#2bab2f")


# Background and text tracks
tracks = tracks + BioCircosBackgroundTrack("testBGtrack1",
  borderColors = "#FFFFFF", borderSize = 0.6, minRadius = snp_minRadius, maxRadius = snp_maxRadius)    

tracks = tracks + BioCircosBackgroundTrack("testBGtrack2", borderColors = "#FFFFFF", 
  fillColor = "#FFEEEE", borderSize = 0.6, minRadius = snp_minRadius, maxRadius = snp_maxRadius)

##########################


##########################
#  SV track  #
##########################

sv_minRadius=0.70
sv_maxRadius=0.77

#Tandem dup, Deletion, Insertion use only start and end position in CI position
###### DUPLICATION ######
SV_list=bedpe_SV[bedpe_SV[,11]=="DUP",][,1:6] 

if(nrow(SV_list)>=1){
arcsEnds = round(SV_list[,6])
arcsLengths = round(SV_list[,6]-SV_list[,2])
arcsChr = as.character(gsub("chr","",SV_list[,1]))

#arcsEnds = round(runif(7, 50000001, 133851895))
#arcsLengths = round(runif(7, 1, 50000000))
#arcsChr=as.character(sample(1:12, 7, replace=T))

tracks = tracks + BioCircosArcTrack("fredTestArc", arcsChr, 
  starts = arcsEnds - arcsLengths, ends = arcsEnds, labels = 1:7, 
  minRadius =sv_minRadius, maxRadius = sv_maxRadius, opacities = 0.6)
}

###### DELETION ######
SV_list=bedpe_SV[bedpe_SV[,11]=="DEL",][,1:6] 

if(nrow(SV_list)>=1){
arcsEnds = round(SV_list[,6])
arcsLengths = round(SV_list[,6]-SV_list[,2])
arcsChr = as.character(gsub("chr","",SV_list[,1]))

tracks = tracks + BioCircosArcTrack("fredTestArc", arcsChr, 
  starts = arcsEnds - arcsLengths, ends = arcsEnds, labels = 1:7, 
  minRadius =sv_minRadius, maxRadius = sv_maxRadius,colors ="#a00000", opacities = 0.6)
}

###### INSERTION ######
SV_list=bedpe_SV[bedpe_SV[,11]=="INS",][,1:6] 

if(nrow(SV_list)>=1){
arcsEnds = round(SV_list[,6])
arcsLengths = round(SV_list[,6]-SV_list[,2])
arcsChr = as.character(gsub("chr","",SV_list[,1]))

tracks = tracks + BioCircosArcTrack("fredTestArc", arcsChr, 
  starts = arcsEnds - arcsLengths, ends = arcsEnds, labels = 1:7, 
  minRadius =sv_minRadius, maxRadius = sv_maxRadius,colors ="#30a000", opacities = 0.6)
}
print("DEL=Red, DUP=Blue, INS=Green")

##########################
#  SV translocation track  #
##########################

SV_trans_list=bedpe_SV[bedpe_SV[,11]=="BND",][,1:6] 


sv_trans_maxRadius=0.55

chr1=as.character(gsub("chr","",SV_trans_list[,1]))
chr2=as.character(gsub("chr","",SV_trans_list[,4]))
gene1Starts_sv=round(SV_trans_list[,2])
gene1Ends_sv=round(SV_trans_list[,3])
gene2Starts_sv=round(SV_trans_list[,5])
gene2Ends_sv=round(SV_trans_list[,6])

tracks = tracks + BioCircosLinkTrack("testLink", gene1Chromosomes = chr1, 
  gene1Starts = gene1Starts_sv, gene1Ends = gene1Ends_sv, gene2Chromosomes = chr2,
  color = "#FF6666", labels = paste(chr1, chr2, sep = "*"), displayLabel = F,
  gene2Starts = gene2Starts_sv, gene2Ends = gene2Ends_sv, maxRadius = sv_trans_maxRadius,borderSize = 0.1, width = 0.5)


##########################
#  CNA track  #
##########################

cna_minRadius=0.77
cna_maxRadius=0.98

snvChr = as.character(gsub("chr","",cnv[,"chr"]))
snvStart = cnv[,"start"]
snvEnd = cnv[,"end"]
snvValues = cnv[,"cp"] # for redian

# Create CNV track
tracks = tracks + BioCircosCNVTrack('cnv_track', as.character(snvChr), snvStart, snvEnd, snvValues, 
  color = "#CC0000", range = c(0,25))

# Add background
tracks = tracks + BioCircosBackgroundTrack("arcs_background", 
    color = "#FFFFFF", minRadius = cna_minRadius, 
    maxRadius = cna_maxRadius,borderColors = "#FFFFFF", 
    borderSize = 0.6)

##########################
# save to html file #
##########################
plot <- BioCircos(tracks, genomeFillColor = "Spectral", yChr = F, chrPad = 0, displayGenomeBorder = F, 
  genomeTicksLen = 1, genomeTicksTextSize = 0, genomeTicksScale = 50000000,
  genomeLabelTextSize = 15, genomeLabelDy = 0)
# Save it as an html file
saveWidget(plot, paste(out,".html",sep=""), selfcontained = FALSE)

# Output in screenshot.png
system(paste("/Applications/Google\\ Chrome.app/Contents/MacOS/Google\\ Chrome --headless --disable-gpu --screenshot ",out,".html",sep=""))
# Output in output.pdf
system(paste("/Applications/Google\\ Chrome.app/Contents/MacOS/Google\\ Chrome --headless --disable-gpu --print-to-pdf ",out,".html",sep=""))

print("Success!")
