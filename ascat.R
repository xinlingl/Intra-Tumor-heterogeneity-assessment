library(ASCAT)
ascat.bc = ascat.loadData("tumorLogR.txt","tumorBAF.txt","normalLogR.txt", "normalBAF.txt")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.bc$Tumor_BAF_segmented[[2]]=ascat.bc$Tumor_BAF_segmented[[1]]
ascat.output = ascat.runAscat(ascat.bc, gamma=1)
write.table(ascat.output$aberrantcellfraction, file='tumor_purity.tsv',quote=FALSE, sep='\t', col.names = NA)
write.table(ascat.output$segments, file='major_minor_copy_number.tsv', quote=FALSE, sep='\t', col.names = NA)



