args = commandArgs(trailingOnly=TRUE)

# Interactive debugging
# args = c("counts_zt.rda")

options(width=200)

print("Loading samples.rda")
load("samples.rda")

filename = args[1]
stopifnot(substr(filename,nchar(filename)-3,nchar(filename)) == ".rda")
filename_noext = substr(filename,1,nchar(filename)-4)

print(sprintf("Loading file with raw counts %s",filename))
stopifnot(file.exists(filename))
load(filename)
# make sure something is loaded
stopifnot("raw" %in% ls())

# Lane merged data.frame (lm)
raw_lm = data.frame(name=raw$name, dna=raw$dna)
lane_rp = c()
lane_nreads = c()

snames = samples[which(! is.na(samples$lane)),"file"]
snames_nolane = unique( substr(snames, 1, nchar(snames)-5) )

# loop over all samples and sum lanes
for (sname in snames_nolane) {
    cn_lanes = paste0(sname,"_L00",seq(4))

    # Sum lane counts
    raw_lm[,sname] = apply(raw[,cn_lanes], MARGIN=1, sum)

    # Correlations
    rp = c(cor(raw[,cn_lanes[1]], raw[,cn_lanes[2]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[1]], raw[,cn_lanes[3]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[1]], raw[,cn_lanes[4]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[2]], raw[,cn_lanes[3]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[2]], raw[,cn_lanes[4]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[3]], raw[,cn_lanes[4]], method="pearson"))
    print(sprintf("Lane correlations of %s: %s", sname, paste(rp,collapse=", ")))
    lane_rp = c(lane_rp, rp)

    # Read counts per lane
    nreads = apply(raw[,cn_lanes], MARGIN=2, sum)
    lane_nreads = c(lane_nreads, nreads)
    print(sprintf("Number of reads for lanes: %s", paste(nreads, collapse=", ")))
}
print(sprintf("Lane correlations %.2f sd %.2f min %.2f max %.2f",mean(lane_rp),sd(lane_rp),min(lane_rp),max(lane_rp)))
print(sprintf("Average reads per lane %.0f sd %.0f min %d max %d",mean(lane_nreads),sd(lane_nreads),min(lane_nreads),max(lane_nreads)))

# add columns that are already lane merged
snames = samples[which( is.na(samples$lane)),"file"]
raw_lm[,snames] = raw[,snames]

raw = raw_lm
outfile = paste0(filename_noext,"_lanemerge.rda")
print(sprintf("Dump %s",outfile))
save(raw, file=outfile)
