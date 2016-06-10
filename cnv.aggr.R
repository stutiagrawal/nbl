#This code is revised from Zhenyu Zhang's code here:http://github.com/ZhenyuZ/eqtl/blob/master/cnv/GetGeneLevelCNA.r

options(stringsAsFactors=F)
library("GenomicRanges")
library(DESeq)
library(optparse)

option_list <- list(
    make_option("--locfile", type="character", help="path to locations file"),
    make_option("--patients", type="character", help="path to uuids"),
    make_option("--segfiles", type="character", help="path to datasets"),
    make_option("--outfile", type="character", help="name of output file")
    )

parser <- OptionParser(usage="main.R [options] file", option_list=option_list)
args <- parse_args(parser)

#get gene locations
get_gene_loc <- function(geneloc){
    
    loc = GRanges(seqnames = geneloc$chr,
                  ranges = IRanges(start=as.numeric(geneloc$start), end=as.numeric(geneloc$end)),
                  strand = "+",
                  name = geneloc$gene_id)
    #gene.length = width(loc) #vector of lengths of each gene id.
    #n = length(loc) #no. of geneids
    return(loc)
}


get_seg_loc <- function(segfile){
    
    cnvseg <- data.frame(read.table(segfile, header=T, colClasses="character", sep="\t"))
    seg <- GRanges(seqnames = paste("chr", cnvseg$Chromosome, sep=""),
                   ranges = IRanges(start=as.numeric(cnvseg$Start), end=as.numeric(cnvseg$End)),
                   strand = "+",
                   Num_probes = as.numeric(cnvseg$Num_Probes),
                   Segment_Mean = as.numeric(cnvseg$Segment_Mean))
    return(seg)
    
}

get_copy_number<- function(loc, seg){
    
    segm <- values(seg)$Segment_Mean
    Hits <- findOverlaps(loc,seg)
    gene_to_seg_map <- data.frame(cbind(queryHits(Hits), subjectHits(Hits))) #contains the indices of the genes 
                                                                             #from the locfile mapped to
                                                                             #the indices of the segfile.
    colnames(gene_to_seg_map ) <- c("loc.index", "seg.index")
    cnv <- numeric(length(loc))
    
    for (gene_index in 1:nrow(gene_to_seg_map)){
        cnv[gene_index] <- 2 * (2^segm[gene_to_seg_map$seg.index[gene_index]])
    }

    return(cnv)
    
}
    
locfile = args$locfile
geneloc = read.table(locfile, h=F, sep="\t", colClasses="character")
colnames(geneloc) <- c("gene_id", "chr", "start", "end")
w <- which(!duplicated(geneloc$gene_id))
geneloc <- geneloc[w,]
loc <- get_gene_loc(geneloc)

patients <- read.table(args$patients, header=T, colClasses="character")

all_cnv = matrix(nrow=nrow(geneloc), ncol=0)
col_counter = 0

for (i in 1:nrow(patients)){
    segfile=paste(args$segfiles, patients[i,1], ".seg.txt", sep="")
    if (file.exists(segfile)){
        col_counter = col_counter + 1
        print(paste("Getting copy number for", patients[i,1]))
        seg <- get_seg_loc(segfile)
        cnv <- get_copy_number(loc, seg)
        all_cnv <- cbind(all_cnv, as.numeric(cnv))
        colnames(all_cnv)[col_counter] <- patients[i,1]
                
    }
}
rownames(all_cnv) <- geneloc$gene_id
write.table(all_cnv, file=args$outfile, quote=FALSE, sep='\t')




