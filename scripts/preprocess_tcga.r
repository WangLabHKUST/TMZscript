########################################################################################################################################################################################################
### TCGA pan-glioma dataset
tcgaexp = read.delim("./GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt.gz", 
                     stringsAsFactors = F, row.names = 1, comment.char = "#")
idx = substr(names(tcgaexp),14,15) == "01" #only use the primary tumor samples
tcgaexp = tcgaexp[,idx]
names(tcgaexp) = substr(names(tcgaexp),start = 1,stop = 12)
dim(tcgaexp)

tcgaau = read.delim("./TCGA_GBMLGG_withRNA_CNVs_20200217.txt", stringsAsFactors = F, comment.char = "#", row.names = 1)
rownames(tcgaau) = gsub(rownames(tcgaau),pattern = "-",replacement = ".") # "-" ==> "." 
dim(tcgaau)

sharedSamples = intersect(names(tcgaexp), rownames(tcgaau))
tcgaexp = tcgaexp[,sharedSamples]
tcgaau = tcgaau[sharedSamples,]
stopifnot(identical(names(tcgaexp), rownames(tcgaau)))
  
write.csv(tcgaau,'./tcga.cna')
write.csv(tcgaexp,'./tcga.exp')
  





              
