library(data.table)
library(plyr)

setwd("/cluster/projects/kridelgroup/FLOMICS/DATA/TGL_BAM_RNASEQ")

count_files=list.files(pattern="ReadsPerGene.out.tab")

out="/cluster/projects/kridelgroup/FLOMICS/DATA/TGL_GENE_COUNTS/"

get_counts_clean = function(file_name){
    dat = fread(file_name)
    dat = dat[,1:2]
    dat=dat[5:nrow(dat),]
    dat$V1 = sapply(dat$V1, function(x){unlist(strsplit(x, "\\."))[1]})
    #add sample name
    dat$sample = unlist(strsplit(file_name, "Reads"))[[1]]
    return(dat)
}

all_counts=as.data.table(ldply(llply(count_files, get_counts_clean, .progress="text")))

#pull into one matrix
mat=(dcast(all_counts, V1 ~ sample, value.var="V2"))
colnames(mat)[1] = "gene"
write.table(mat, paste(out,"STAR_quantmode_counts_matrix_FL_136_patients.txt",sep=""), quote=F, row.names=F, sep="\t")
