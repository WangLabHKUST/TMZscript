

#data preprocessing
dataNor <-function(dfInput){
  df = dfInput
  df[is.na(df)] <- 0
  df = log2(df+1)
  df = scale(t(df))
  return (df)
}

infile = "./tcga.exp.csv"
ipt = read.csv(infile,row.names=1)

dmn = dim(ipt)
dtexp = dataNor(ipt)
write.table(t(dtexp),'./tcga.exp.nor.csv')
  

infile = "./200503.76.RPKM.csv"
ipt = read.csv(infile,row.names=1)

dmn = dim(ipt)
dtexp = dataNor(ipt)
write.table(t(dtexp),'./200503.76.RPKM.nor.csv')





              
