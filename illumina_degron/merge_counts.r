#!/usr/local/bin/Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)

options(width=180, stringsAsFactors=F)

# Collect files with raw counts and frequencies
# ---------------------------------------------

lib = read.table("ref.seq")
colnames(lib) = c("name","dna")
stopifnot(length(lib$dna) == length(unique(lib$dna)))


raw = lib
raw$count = 1
file_list = list()
n_genotypes = c()
n_lib_genotypes = c()
n = 0
for (countfile in args) {
    stopifnot(file.exists(countfile))
    stopifnot(substring(countfile, nchar(countfile)-10, nchar(countfile)) == "_counts.txt")

    # File and directory names
    l = strsplit(countfile, "/")[[1]]
    cfile_bn = l[length(l)]
    if (length(l) > 1) {
        cfile_dn = paste(l[1:(length(l)-1)], collapse="/")
    } else {
        cfile_dn="."
    }
    cfile_id = substring(cfile_bn, 1, nchar(cfile_bn)-11)
    # cfile_id = substring(cfile_bn, 1, nchar(cfile_bn)-4)
    file_list[[cfile_id]] = countfile

    # Read count file
    # print(sprintf("Read %s",countfile))
    d = read.table(countfile)
    colnames(d) = c("dna","count","rpm")

    n = n+1
    thres = 1
    i = which(d$count > thres)
    n_lib_gt = length(intersect(d$dna,lib$dna))
    print(sprintf("File %3d/%-3d %-20s has %5d (%4d) genotypes with an overlap of %4d (%4d) to original library (counts > %d). Merging with %6d genotypes overlap %4d",
          n, length(args), cfile_id, nrow(d), length(i),
          n_lib_gt, length(intersect(d[i,"dna"],lib$dna)), thres, nrow(raw), length(intersect(d[i,"dna"],raw$dna))))
    n_genotypes = c(n_genotypes, nrow(d))
    n_lib_genotypes = c(n_lib_genotypes, n_lib_gt)

    # merge while keeping all dna sequences. raw$count gets suffix "" and d$count gets suffix ".file_id"
    raw = merge(raw, d[,c("dna","count")], by="dna", all=T, suffixes=c("",paste0(".",cfile_id)))
}

print(sprintf("In total %d genotypes",nrow(raw)))

print("Summary of genotypes per file")
print(summary(n_genotypes))

print(sprintf("Summary of genotype overlap with original library of %d",nrow(lib)))
print(summary(n_lib_genotypes))

# order raw according to dna length and secondarily sequence (poor mans alignment)
raw$len = nchar(raw$dna)
raw = raw[order(raw$len,raw$dna),]

# names of colunms with counts
cns = paste0("count.",names(file_list))

# replace NA with zero
for (cn in cns) { raw[ which( is.na(raw[,cn]) ),cn ] = 0 }

# new data.frame without the 'len' and merge-practical 'count' column
raw = raw[,c("name","dna",cns)]

# rename columns
colnames(raw) = c("name","dna",names(file_list))

# dump
save(raw, file="counts_zt.rda")

