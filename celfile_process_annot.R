#Preprocessing of CEL files

source("http://bioconductor.org/biocLite.R")
biocLite('huex10sttranscriptcluster.db')

library(oligo)
library(huex10sttranscriptcluster.db)
library(sqldf)

#reading the CEL files
exonCELs <- list.celfiles()
affyExonFS <- read.celfiles(exonCELs)

#normalization and getting the expression values
input <- exprs(rma(affyExonFS, target = "core"))

#METHOD 1: USING sqldf
#annotating the transcript cluster IDs

Annot <- data.frame(SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=","),
                    DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=","),
                    ENTREZID=sapply(contents(huex10sttranscriptclusterENTREZID), paste, collapse=","))

names(Annot)

#write it to a file
write.csv(Annot,"Annot.csv")

new <- read.csv("Annot.csv")

#reading in the transcript cluster IDs
c <- read.csv("clusterIDs.csv")
head(c)
#check if its a data frame
is.data.frame(c)

#we have Annot table and another table containing TranscriptCluster IDs
final <- sqldf("select c.SYMBOL, new.* from c LEFT JOIN new WHERE new.SYMBOL==c.SYMBOL")
head(final)
write.table(final, file = "annotated file.txt")


###########################

# using annotationDB

#Extract transcript clusterIDs
probes=row.names(new)
probes

Symbols = unlist(mget(probes, huex10sttranscriptclusterSYMBOL, ifnotfound=NA))
head(Symbols)
Entrez_IDs = unlist(mget(probes, huex10sttranscriptclusterENTREZID, ifnotfound=NA))
head(Entrez_IDs)

#Combine gene annotations with raw data
new=cbind(probes,Symbols,Entrez_IDs,new)
head(new)
#Write RMA-normalized, mapped data to file
write.table(new, file = "input_wgcna.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

################################################
