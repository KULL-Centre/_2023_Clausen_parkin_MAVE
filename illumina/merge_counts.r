options(width=160, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    files = c("aspa_toxicity_call/17_S17_L002_counts.txt","aspa_toxicity_call/19_S19_L003_counts.txt")
    # files = c("test_counts.txt")
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript process.r  <counts1.txt>  [counts2.txt  ...]")
    quit(save="no")
} else {
    files = args
}

bc_file = "barcode_map.rda"
load(bc_file)
print(sprintf("Loaded barcode map for %s from %s",wt[["name"]],bc_file))

raw = data.frame(barcode = bc_map$barcode, var_aa = bc_map$var_aa, var_dna = bc_map$var_dna)

for (file in files) {
    d = read.table(file)
    fl = strsplit(file, "/")[[1]]
    stopifnot(substr(fl[length(fl)], nchar(fl[length(fl)])-10, nchar(fl[length(fl)])) == "_counts.txt")
    file_id = substr(fl[length(fl)], 1, nchar(fl[length(fl)])-11)
    raw[,file_id] = d[match(raw$barcode,d$V1),"V2"]
    na_mask = is.na(raw[,file_id])
    raw[na_mask,file_id] = 0
    # Report
    mapped_bc = sum(! na_mask)
    mapped_reads = sum(raw[which(! na_mask),file_id])
    print(sprintf("Mapped %.2f%% of barcodes and %.2f%% of reads from %s (%d barcodes and %d reads)",
                  mapped_bc/nrow(d)*100, mapped_reads/sum(d$V2)*100, file_id, mapped_bc, mapped_reads))
    
}

save(raw, wt, file="raw.rda")

# colnames(raw) = c("counts","barcode")
# raw$var = bc_map[match(raw$barcode,bc_map$barcode),"var"]
# raw$aa_subst = bc_map[match(raw$barcode,bc_map$barcode),"aa_subst"]

# n_na = sum(is.na(raw$var))
# print(sprintf("Barcodes not mapped: %d of %d, %.2f%%, total %d counts", n_na, nrow(raw), n_na/nrow(raw)*100.0, sum(raw[is.na(raw$var),"counts"])))

# n_wt = sum(raw$var=="", na.rm=T)
# print(sprintf("Barcodes mapped to WT: %d of %d, %.2f%%, total %d counts", n_wt, nrow(raw), n_wt/nrow(raw)*100.0, sum(raw[which(raw$var==""),"counts"])))

# n_1m = sum(raw$aa_subst==1, na.rm=T)
# print(sprintf("Barcodes mapped to single mutants: %d of %d, %.2f%%, total %d counts",
#   n_1m, nrow(raw), n_1m/nrow(raw)*100.0, sum(raw[which(raw$aa_subst==1),"counts"])))


