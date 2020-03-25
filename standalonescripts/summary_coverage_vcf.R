###This script count the number of SNPs per individuals in a given vcf file (input at line 4)

###Parameters
vcf_file  =  #Path to your vcf i.e. "myvcf.vcf"

### Dependencies
library("VariantAnnotation") # install from https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html​

#Define the function
summary_coverage_vcf<-function(vcf_file){
  print ("Reading data...")
  data<-readVcf(vcf_file) # Read the data

  #​Define a function that counts the number of missing data
  numbermissing<-function(x){ 
    return(length(grep("\\./\\.",x)))
  }

  # Count the  missing and then infer the non-missing
  countsofmissing<-apply(geno(data)$GT,2,numbermissing)
  nonmissing<-dim( geno(data)$GT)[1]-countsofmissing
  print(paste("There are",dim(data)[1],"SNPs and ", dim(data)[2],"individuals"))
  print("The distribution of number of SNPs per individual is:")
  print(summary(nonmissing))
  print("This is the number of covered SNPs sample by sample:")
  sort(nonmissing) # Will output the number of SNPs for each sample from minimum to maxium
}

#Run it on your vcf
summary_coverage_vcf(vcf_file)
