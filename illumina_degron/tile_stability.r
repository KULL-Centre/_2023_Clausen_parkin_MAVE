#####################################################################
##     Load read counts, sample information and library information
#####################################################################
args = commandArgs(trailingOnly=TRUE)

filename = args[1]
print(sprintf("Loading file with raw counts %s",filename))
load(filename)

options(width=200)
settings = list()
settings[["version"]] = 1.2
settings[["input"]] = filename

# Load information on samples (FASTQ files)
print("Loading samples.rda")
load("samples.rda")

# remove lane info - counts should be merged
samples_lanes = samples
samples = samples_lanes[which(samples_lanes$lane==1 | is.na(samples_lanes$lane)),]
samples$name = samples$file
i = which(! is.na(samples$lane))
samples[i,"name"] = substr(samples[i,"file"], 1, nchar(samples[i,"file"])-5)

# Load library
print("Loading library.rda")
load("library.rda")

# assign a residue position per tile
library$resi = (library$tile_last + library$tile_first)/2
library$name2 = paste(library$protein, library$tile_first, library$tile_last, sep="_")
i = which(is.na(library$tile_first))
library[i,"name2"] = library[i,"name"]

# columns containing read counts
cns = colnames(raw)[3:ncol(raw)]
stopifnot(setequal(samples$name, cns))


#####################################################################
##     Clean
#####################################################################

# remove reads that are not in any library
ri_out = which(is.na(raw$name))
for (cn in cns) {
    n_out = sum(raw[ri_out,cn])
    n_all = sum(raw[,cn])
    print(sprintf("Removing %d of %d read counts (%.1f%%) from %s that are outside all libraries", n_out, n_all, n_out/n_all*100.0, cn))
}
ri_in = which(! is.na(raw$name))
print(sprintf("Only keeping %d of %d unique reads that match library",length(ri_in),nrow(raw)))
raw_clean = raw[ri_in,]

# zero read counts that does not belong to O, E or CT library
coverage = c()
for (cn in cns) {
    cn_lib = samples[match(cn,samples$name),"lib"]
    ir_outoflib = which(! (raw_clean$name %in% library[which(library$pool==cn_lib | library$pool=="all"),"name"]))
    n_out = sum(raw_clean[ir_outoflib,cn])
    n_all = sum(raw_clean[,cn])
    print(sprintf("Removing %d of %d read counts (%.1f%%) from %s that does not belong to library %s", n_out, n_all, n_out/n_all*100.0, cn, cn_lib ))
    raw_clean[ir_outoflib,cn] = 0
    
    # coverage based on cleaned counts
    coverage = c(coverage, sum(raw_clean[,cn])/sum(library$pool==cn_lib | library$pool=="all"))
}

# plot reads per tile
quartz(width=10, height=6)
lib_fac = as.factor(samples[match(cns,samples$name),"lib"])
barplot(coverage, ylab="Reads per tile", ylim=c(0,10000), col=as.integer(lib_fac)+1) #, legend.text=levels(lib_fac), args.legend=list(col=seq_along(levels(lib_fac))))
legend("topright", levels(lib_fac), col=seq_along(levels(lib_fac))+1, pch=15)
quartz.save("coverage.png", type="png")

# order rows as in lib
raw_clean = raw_clean[match(library$name,raw_clean$name),]
stopifnot( all( raw_clean$dna == library$dna ) )


#####################################################################
##     Merge replicates
#####################################################################
# Technical replica merged (tm)
trep_nreads = c()
trep_rp = c()
raw_tm = data.frame(name=raw_clean$name)
par(ask=T)
for (si in which(samples$tech_rep==1)) {
    # sample row of other tech_rep 
    si2 = which(samples$lib == samples[si,"lib"]  &  samples$gate == samples[si,"gate"]  &  samples$bio_rep == samples[si,"bio_rep"]  &  
               samples$facs_rep == samples[si,"facs_rep"]  &  samples$tech_rep == 2)
    stopifnot(length(si2) == 1)

    cn_rep = c(samples[si,"name"], samples[si2,"name"])
    cn = sprintf("%s%d_%s_fasc%d", samples[si,"lib"], samples[si,"gate"], samples[si,"bio_rep"], samples[si,"facs_rep"])
    
    raw_tm[,cn] = apply(raw_clean[,cn_rep], MARGIN=1, sum)

    # Correlations of relevant library members
    # ri = which(! is.na(raw_clean$name))
    ri = which(raw_clean$name %in% library[which(library$pool==samples[si,"lib"]),"name"])
    rp = cor(raw_clean[ri,cn_rep[1]], raw_clean[ri,cn_rep[2]], method="pearson")
    rs = cor(raw_clean[ri,cn_rep[1]], raw_clean[ri,cn_rep[2]], method="spearman")
    rp_log = cor(log(raw_clean[ri,cn_rep[1]]+1e-6), log(raw_clean[ri,cn_rep[2]]+1e-6), method="pearson")
    trep_rp = c(trep_rp, rp)
    print(sprintf("Pearson between technical replica %s and %s of %d lib %s members %.2f (log %.2f, spearman %.2f). Counts summaries:",
                  cn_rep[1],cn_rep[2],length(ri),samples[si,"lib"],rp,rp_log,rs))
    # print(summary(raw_clean[ri,cn_rep[1]]))
    # print(summary(raw_clean[ri,cn_rep[2]]))

    # plot(raw_clean[ri,cn_rep[1]], raw_clean[ri,cn_rep[2]], xlab=cn_rep[1], ylab=cn_rep[2])

    # Read counts per replica
    nreads = apply(raw_clean[,cn_rep], MARGIN=2, sum)
    trep_nreads = c(trep_nreads, nreads)
    print(sprintf("Number of reads for technical replica: %s", paste(nreads, collapse=", ")))
}
print("Tech replica correlations")
print(summary(trep_rp))
print("Reads per tech replica")
print(summary(trep_nreads))


#####################################################################
##     Calculate PSI
#####################################################################
# Function that normalize all columns to RPM
reads_per_million = function(df) {
    col_sums = apply(df, MARGIN=2, sum)
    df * matrix(10^6/col_sums, nrow=dim(df)[1], ncol=dim(df)[2], byrow=TRUE)
}

# Calculate count frequencies with psudo_counts
psudo_counts = 0
rpm = data.frame(name=raw_tm$name)
rpm = cbind(rpm, reads_per_million(raw_tm[,2:ncol(raw_tm)] + psudo_counts))


# PSI calculated per tile and replica
psi = data.frame(name=library$name2, aa=library$aa)
rpm_sum = data.frame(name=raw_clean$name)
for (si in which(samples$gate==1 & samples$tech_rep==1)) {
    cn_gates = sprintf("%s%d_%s_fasc%d", samples[si,"lib"], seq(4), samples[si,"bio_rep"], samples[si,"facs_rep"])
    cn = sprintf("%s_%s_facs%d", samples[si,"lib"], samples[si,"bio_rep"], samples[si,"facs_rep"])
    print(sprintf("Calculate PSI for %s",cn))

    psi[,cn] = apply(rpm[,cn_gates], MARGIN=1, function(v){vsum=sum(v); sum( v/vsum * c(1,2,3,4) )})
    rpm_sum[,cn] = apply(rpm[,cn_gates], MARGIN=1, sum)
}

# Average PSI over replicates for each sub-library
bn = unique(samples$bio_rep)
libs = c("O","E","CT")
all_rps = c()
for (lib in libs) {
    cns_rep = sprintf("%s_%s_facs%d", lib, bn, rep(seq(2), each=length(bn)))
    cn = paste0(lib,".mean")
    psi[,cn] = apply(psi[,cns_rep], MARGIN=1, mean, na.rm=T)
    cn = paste0(lib,".sd")
    psi[,cn] = apply(psi[,cns_rep], MARGIN=1, sd, na.rm=T)

    print(sprintf("Summary of SD between PSI values of biological replicates for %s library",lib))
    print(summary(psi[,cn]))

    # correlations - these are independt of normalization
    rp_mat = cor(psi[,cns_rep], method="pearson", use="pairwise.complete.obs")
    rps = rp_mat[upper.tri(rp_mat)]
    print(sprintf("Lib %s PSI Pearson of %d replica mean %.4f range %.4f - %.4f", lib, length(rps), mean(rps), min(rps), max(rps)))
    all_rps = c(all_rps, rps)

    # plot replica correlations    
    quartz(width=8, height=8)
    plot(psi[,cns_rep], main=sprintf("Library %s",lib), upper.panel=NULL)
    quartz.save(sprintf("replica_correlations_%s.png",lib), type="png")
}
print(sprintf("All lib PSI Pearson of %d replica mean %.4f range %.4f - %.4f", length(all_rps), mean(all_rps), min(all_rps), max(all_rps)))


# plot psi distribution per library
quartz(width=8, height=5)
breaks = seq(1,4,0.5)
h_odd = hist(psi[,"O.mean"], breaks=breaks, plot=F)
h_even = hist(psi[,"E.mean"], breaks=breaks, plot=F)
h_ct = hist(psi[,"CT.mean"], breaks=breaks, plot=F)
str_odd = sprintf("Odd, %d tiles",sum(library$pool %in% c("O","all")))
str_even = sprintf("Even, %d tiles",sum(library$pool %in% c("E","all")))
str_ct = sprintf("CT, %d tiles",sum(library$pool %in% c("CT","all")))
bp = barplot(matrix(c(h_odd$counts,h_even$counts,h_ct$counts), ncol=length(breaks)-1, byrow=T), beside=T,
             main="PSI distribution per library", xlab="PSI", ylab="Tiles",
             legend.text=c(str_odd,str_even,str_ct), args.legend=list(x="top", ncol=3))
axis(1, seq(bp[1,1]-1, bp[nrow(bp),ncol(bp)]+1, nrow(bp)+1), seq(1,4,0.5))
quartz.save("lib_psi.png", type="png")


#####################################################################
##     Normalize sub-libraries and merge PSI values
#####################################################################
# Normalize even and ct to odd library - in this case seems to have lowest uncertainty
# Also odd lib has lowest PSI for the stable control which means that other values are scaled down and not up, potentially above 4.
ip_stab = which(psi$name == "APPY_WT_1")
ip_intm = which(psi$name == "APPY_RAAA_1")
ip_dest = which(psi$name == "APPY_DAAA_1")
appy = c(ip_stab,ip_intm,ip_dest)

print("Per library replica average")
print(psi[appy,c("name",paste0(libs,rep(c(".mean",".sd"),each=length(libs))))])

# Normalize even
fit_even = lm( psi[c(ip_stab,ip_dest),"O.mean"] ~ psi[c(ip_stab,ip_dest),"E.mean"] )
print(sprintf("Normalization of even lib: even_norm = %.4f + even * %.4f", coef(fit_even)[1], coef(fit_even)[2]))
psi$E.mean_norm = psi[,"E.mean"] * coef(fit_even)[2] + coef(fit_even)[1]
psi$E.sd_norm = psi[,"E.sd"] * coef(fit_even)[2]

# Normalize CT
fit_ct = lm( psi[c(ip_stab,ip_dest),"O.mean"] ~ psi[c(ip_stab,ip_dest),"CT.mean"] )
print(sprintf("Normalization of CT lib: ct_norm = %.4f + ct * %.4f", coef(fit_ct)[1], coef(fit_ct)[2]))
psi$CT.mean_norm = psi[,"CT.mean"] * coef(fit_ct)[2] + coef(fit_ct)[1]
psi$CT.sd_norm = psi[,"CT.sd"] * coef(fit_ct)[2]

# merge into a common column
single_non_na = function(v) { i=which(! is.na(v)); if (length(i)==1) {v[i]} else {NA} }
psi$mean_merge = apply(psi[,c("O.mean","E.mean_norm","CT.mean_norm")], MARGIN=1, single_non_na)
psi$sd_merge = apply(psi[,c("O.sd","E.sd_norm","CT.sd_norm")], MARGIN=1, single_non_na)

# positive and negative control peptides are by normalization definition the value from the odd lib
psi[appy,"mean_merge"] = psi[appy,"O.mean"]
psi[appy,"sd_merge"] = psi[appy,"O.sd"]

# For intermediate control, average normalized values
psi[ip_intm,"mean_merge"] = mean(unlist(psi[ip_intm,c("O.mean","E.mean_norm","CT.mean_norm")]))
psi[ip_intm,"sd_merge"] = sd(unlist(psi[ip_intm,c("O.mean","E.mean_norm","CT.mean_norm")]))
stopifnot(all(! is.na(psi$mean_merge)))
stopifnot(all(! is.na(psi$sd_merge)))

# Report normalization
print("Normalized APPY mean PSI:")
print(psi[appy,c("name","O.mean","E.mean","CT.mean","E.mean_norm","CT.mean_norm")])

print("Summary of standard deviations between PSI values of biological replicates")
print(summary(psi$sd_merge))

# Plot normalization parameters
bn = unique(samples$bio_rep)
repn = sprintf("%s_facs%d", bn, rep(seq(2), each=length(bn)))
quartz(width=6, height=6)
plot(0,0,col=0, xlim=c(1,4), ylim=c(1,4), main="PSI normalization between sub-libraries", xlab="PSI", ylab="PSI, odd lib")
points(unlist(psi[appy,paste0("E_",repn)]), unlist(psi[appy,paste0("O_",repn)]), pch=rep(seq_along(appy),times=length(repn)), col=2, cex=.7)
# plot odd
abline(c(0,1))
# plot even
points(psi[appy,"E.mean"], psi[appy,"O.mean"], pch=seq_along(repn), col=2, cex=1.5, lwd=2)
abline(coef(fit_even), col=2)
points(unlist(psi[appy,paste0("CT_",repn)]), unlist(psi[appy,paste0("O_",repn)]), pch=rep(seq_along(appy),times=length(repn)), col=3, cex=.7)
# plot ct
points(psi[appy,"CT.mean"], psi[appy,"O.mean"], pch=seq(3), col=3, cex=1.5, lwd=2)
abline(coef(fit_ct), col=3)
# legend
str_even = sprintf("Even %.2fx+%.2f",coef(fit_even)[2],coef(fit_even)[1])
str_ct = sprintf("CT %.2fx+%.2f",coef(fit_ct)[2],coef(fit_ct)[1])
legend("topleft", c("Odd 1.00x+0.00",str_even,str_ct,"APPY RLLL","APPY RAAA","APPY DAAA"), pch=c(NA,NA,NA,1,2,3), lty=c(1,1,1,NA,NA,NA), col=c(1,2,3,1,1,1))
quartz.save("norm.png", type="png")

# Plot overall PSI distribution before and after normalization
quartz(width=8, height=5)
# psi can move outside [1,4] when normalized
breaks = seq(0.5,4.5,0.5)
h1 = hist(apply(psi[,c("O.mean","E.mean","CT.mean")], MARGIN=1, mean, na.rm=T), breaks=breaks, plot=F)
h2 = hist(apply(psi[,c("O.mean","E.mean_norm","CT.mean_norm")], MARGIN=1, mean, na.rm=T), breaks=breaks, plot=F)
h3 = hist(psi$mean_merge, breaks=breaks, plot=F)
bp = barplot(matrix(c(h1$counts,h2$counts), ncol=length(h1$breaks)-1, byrow=T), beside=T, xlab="PSI", ylab="Tiles",
             legend.text=c("Un-normalized","Normalized"), args.legend=list(x="top"))
axis(1, seq(bp[1,1]-1, bp[nrow(bp),ncol(bp)]+1, nrow(bp)+1), breaks)
quartz.save("norm_psi.png", type="png")

# Plot control PSI variation between odd, even and ct libraries
quartz(height=6, width=5)
bp = barplot(as.matrix(psi[appy,c("O.mean","E.mean","CT.mean")]), beside=T, ylab="PSI", ylim=c(1.0,4.0), xpd=F,
             names.arg=c("Odd","Even","CT"), legend.tex=psi[appy,"name"], args.legend=list(x="top"))
y0 = psi[appy,c("O.mean","E.mean","CT.mean")] + psi[appy,c("O.sd","E.sd","CT.sd")]
y1 = psi[appy,c("O.mean","E.mean","CT.mean")] - psi[appy,c("O.sd","E.sd","CT.sd")]
arrows(as.vector(bp), y0=as.vector(as.matrix(y0)), y1=as.vector(as.matrix(y1)), code=3, length=.1, angle=90)
quartz.save("appy_psi.png", type="png")

# Plot profiles per protein
for (pn in unique(library$protein)) {
    pil = which(library$protein==pn)
    if (length(pil) > 1) {
        quartz(height=4, width=10)
	ntic = 20
        tic_dist = ceiling((length(pil)*12+12) / ntic)
	# print(sprintf("%s with %d tiles gets %d tics distanced %.1f",pn,length(pil),ntic,tic_dist), cex.axis=.5)
        plot(library[pil,"resi"], psi[pil,"mean_merge"], type="b", ylim=c(1,4), main=pn, xlab="Residue number", ylab="PSI", xaxp=c(0,ntic*tic_dist,tic_dist))
	psi_nonorm = apply(psi[pil,paste0(c("O","E","CT"),".mean")], MARGIN=1, sum, na.rm=T)
        points(library[pil,"resi"], psi_nonorm, pch=20, col=2)
        quartz.save(paste0("profile_",pn,".png"), type="png")
    }
}


#####################################################################
##      Dump result
#####################################################################
save(psi, samples, library, proteins, settings, file="tile_stability.rda")

df = data.frame(name=psi$name, aa=psi$aa, psi=psi$mean_merge, psi_sd=psi$sd_merge)

write.table(df, sep=";", row.names=F, quote=F, file="tile_stability.csv")

