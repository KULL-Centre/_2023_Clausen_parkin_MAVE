load("../non_released_output/raw.rda")

# Samples are named with a sample id, a sample number and a lane number separated by underscore, 17-3-C3-1_S83_L002
# The sample id is composed by day-month-Xfacs-tech, where X is the facs bin A-D, facs is the facs replicate and tech
# is the technical replicate number. Gate A has lowest fluorescence
dates     = rep(c("07-4","17-3","21-4","21-5"), each=24)
gates     = rep(rep(c("A","B","C","D"), each=6), 4)
facs_rep  = rep(rep(seq(3), each=2), 16)
tech_rep  = rep(seq(2), 48)

sample = data.frame(date=dates, gate=gates, facs_rep=facs_rep, tech_rep=tech_rep)

# exception, this day had 3 technical replicates and number one failed so replace 1 by 3
sample[which(sample$date=="21-4" & sample$tech_rep==1),"tech_rep"] = 3
sample[which(sample$date=="21-5" & sample$tech_rep==2),"tech_rep"] = 3

# generate sample id's including the exceptions
sample$sample_id = sprintf("%s-%s%d-%d", sample$date, sample$gate, sample$facs_rep, sample$tech_rep)

# change exceptions back for automated replica identificantion 
sample[which(sample$date=="21-4" & sample$tech_rep==3),"tech_rep"] = 1
sample[which(sample$date=="21-5" & sample$tech_rep==3),"tech_rep"] = 2

# match sample number to sample id's
cn_list = strsplit(colnames(raw)[4:ncol(raw)], "_")
cn_samples = sapply(cn_list, '[[', 1)
sample$sample_n = sapply(cn_list, '[[', 2)[match(sample$sample_id,cn_samples)]

# Lane merged (lm)
raw_lm = data.frame(barcode=raw$barcode, var=raw$var_aa)
lane_rp = c()
lane_nreads = c()
# loop over all samples
for (i in seq(nrow(sample))) {
    cn_lanes = paste0(sample[i,"sample_id"],"_",sample[i,"sample_n"],"_L00",seq(4))
    raw_lm[,sample[i,"sample_id"]] = apply(raw[,cn_lanes], MARGIN=1, sum)

    # Correlations
    rp = c(cor(raw[,cn_lanes[1]], raw[,cn_lanes[2]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[1]], raw[,cn_lanes[3]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[1]], raw[,cn_lanes[4]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[2]], raw[,cn_lanes[3]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[2]], raw[,cn_lanes[4]], method="pearson"))
    rp = append(rp, cor(raw[,cn_lanes[3]], raw[,cn_lanes[4]], method="pearson"))
    print(sprintf("Lane correlations of %s: %s", sample[i,"sample_id"], paste(rp,collapse=", ")))
    lane_rp = c(lane_rp, rp)

    # Read counts per lane
    nreads = apply(raw[,cn_lanes], MARGIN=2, sum)
    lane_nreads = c(lane_nreads, nreads)
    print(sprintf("Number of reads for lanes: %s", paste(nreads, collapse=", ")))
}
print(sprintf("Lane correlations %.2f sd %.2f min %.2f max %.2f",mean(lane_rp),sd(lane_rp),min(lane_rp),max(lane_rp)))
print(sprintf("Average reads per lane %.0f sd %.0f min %d max %d",mean(lane_nreads),sd(lane_nreads),min(lane_nreads),max(lane_nreads)))

# Function that normalize all columns to RPM
reads_per_million = function(df) {
    col_sums = apply(df, MARGIN=2, sum)
    df * matrix(10^6/col_sums, nrow=dim(df)[1], ncol=dim(df)[2], byrow=TRUE)
}

# Calculate count frequencies with psudo_counts
psudo_counts = 0
rpm = data.frame(barcode=raw$barcode, var=raw$var_aa)
rpm = cbind(rpm, reads_per_million(raw_lm[,seq(ncol(raw_lm)-nrow(sample)+1,ncol(raw_lm))] + psudo_counts))


# Technical replica merged (rm)
rep_nreads = c()
rep_rp_raw = c()
rep_rp_rpm = c()
raw_rm = data.frame(barcode=raw$barcode, var=raw$var_aa)
rpm_rm = data.frame(barcode=raw$barcode, var=raw$var_aa)
for (i in which(sample$tech_rep==1)) {
    # sample row of other tech_rep 
    i2 = which(sample$date     == sample[i,"date"]      &  sample$gate     == sample[i,"gate"] &
               sample$facs_rep == sample[i,"facs_rep"]  &  sample$tech_rep == 2)
    stopifnot(length(i2) == 1)
    cn_rep = c(sample[i,"sample_id"], sample[i2,"sample_id"])
    cn = sprintf("%s-%s%d", sample[i,"date"], sample[i,"gate"], sample[i,"facs_rep"])
    raw_rm[,cn] = apply(raw_lm[,cn_rep], MARGIN=1, sum)
    rpm_rm[,cn] = apply(rpm[,cn_rep], MARGIN=1, mean)

    # Correlations
    rp_raw = cor(raw_lm[,cn_rep[1]], raw_lm[,cn_rep[2]], method="pearson")
    rep_rp_raw = c(rep_rp_raw, rp_raw)
    rp_rpm = cor(rpm[,cn_rep[1]], rpm[,cn_rep[2]], method="pearson")
    rep_rp_rpm = c(rep_rp_rpm, rp_rpm)    
    print(sprintf("Pearson between technical replica %s and %s raw %.2f and rpm %.2f",cn_rep[1],cn_rep[2],rp_raw,rp_rpm))
    
    # Read counts per replica
    nreads = apply(raw_lm[,cn_rep], MARGIN=2, sum)
    rep_nreads = c(rep_nreads, nreads)
    print(sprintf("Number of reads for technical replica: %s", paste(nreads, collapse=", ")))
}
print(sprintf("Replica raw correlations %.2f sd %.2f min %.2f max %.2f",mean(rep_rp_raw),sd(rep_rp_raw),min(rep_rp_raw),max(rep_rp_raw)))
print(sprintf("Replica rpm correlations %.2f sd %.2f min %.2f max %.2f",mean(rep_rp_rpm),sd(rep_rp_rpm),min(rep_rp_rpm),max(rep_rp_rpm)))
print(sprintf("Average reads per replica %.0f sd %.0f min %d max %d",mean(rep_nreads),sd(rep_nreads),min(rep_nreads),max(rep_nreads)))

# # Dump counts per barcode
# write.table(raw_rm, file="raw_bc.csv", sep=";", row.names=F, quote=F)
# write.table(rpm_rm, file="rpm_bc.csv", sep=";", row.names=F, quote=F)


# Function to calculate PSI based on all columns of a data frame
protein_stability_index = function(df, name) {
    # First gate (gate 1) is stable with index 4, and the late (gate 4) unstable with index 1
    psi = t(apply(df, MARGIN=1, function(v){vsum=sum(v); c(vsum, sum( v/vsum * c(4,3,2,1) ))}))
    # put in a data frame and give column name
    ret_df = data.frame(psi)
    colnames(ret_df) = paste0(name,c(".rpm_sum",".psi"))
    return(ret_df)
}

# PSI calculated per barcode
psi_bc = data.frame(barcode=raw$barcode, var=raw$var_aa, var_dna=raw$var_dna)
for (i in which(sample$gate=="A" & sample$tech_rep==1)) {
    cn_gates = sprintf("%s-%s%d", sample[i,"date"], c("A","B","C","D"), sample[i,"facs_rep"])
    cn = sprintf("run%s_facs%d",sample[i,"date"],sample[i,"facs_rep"])
    print(sprintf("Calculate PSI for %s",cn))
    psi_bc = cbind(psi_bc, protein_stability_index(rpm_rm[,cn_gates], cn))
}

# # Dump counts per barcode
# write.table(psi_bc, file="psi_bc.csv", sep=";", row.names=F, quote=F)

# Aggregate per variant
print("Aggregate barcodes per variant")
agg = aggregate(seq(nrow(rpm_rm)), list(rpm_rm$var), paste, collapse=":")
var_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(var_indices) = agg[,1]

psi = data.frame(var=names(var_indices), n_bc=sapply(var_indices,length))
subst_list = strsplit(psi$var, ":")
psi$n_subst = sapply(subst_list, length)
for (i in which(sample$gate=="A" & sample$tech_rep==1)) {
    cn = sprintf("run%s_facs%d",sample[i,"date"],sample[i,"facs_rep"])
    print(sprintf("Aggregate %s",cn))
    psi[,cn] = sapply(var_indices, function(v) { mean(psi_bc[v,paste0(cn,".psi")], na.rm=T) })
    psi[,paste0("sd_",cn)] = sapply(var_indices, function(v) { sd(psi_bc[v,paste0(cn,".psi")], na.rm=T) })
}

# order variants
psi$resi = sapply(subst_list, function(v){ if (length(v)>0) as.numeric(substr(v[1],2,nchar(v[1])-1)) else 0 })
psi$mut_aa = sapply(subst_list, function(v){ if (length(v)>0) substr(v[1],nchar(v[1]),nchar(v[1])) else "" })
psi = psi[order(psi$n_subst, psi$resi, psi$mut_aa, decreasing=F),]
# var_rpm = var_rpm[order(var_rpm$n_subst, var_rpm$resi, var_rpm$mut_aa, decreasing=F),]

# Calculate a mean of replica
cn = unique( sprintf("run%s_facs%d",dates,facs_rep) )
psi$mean.psi = apply(psi[,cn], MARGIN=1, mean, na.rm=T)
psi$sd.psi = apply(psi[,cn], MARGIN=1, sd, na.rm=T)

# Min-max normalization
psi_nonsense = median(psi[which(psi$mut_aa == "*"),"mean.psi"], na.rm=T)
scale = 1/(psi[1,"mean.psi"] - psi_nonsense)
psi$mean.psi.norm = (psi$mean.psi - psi_nonsense) * scale
psi$sd.psi.norm = psi$sd.psi * scale
print(sprintf("Normalized PSI to WT = %.3f and nonsense = %.3f",psi[1,"mean.psi"],psi_nonsense))

date_rp = c()
date_mae = c()
for (date in unique(dates)) {
    cn = sprintf("run%s_facs%d",date,unique(sample[which(sample$date==date),"facs_rep"]))
    
    rp = cor(psi[,cn[1]], psi[,cn[2]], method="pearson", use="complete.obs")
    mae = mean( abs(psi[,cn[1]]-psi[,cn[2]]), na.rm=T)
    if (length(cn) > 2) {
        rp = c(rp, cor(psi[,cn[1]], psi[,cn[3]], method="pearson", use="complete.obs"))
        rp = c(rp, cor(psi[,cn[2]], psi[,cn[3]], method="pearson", use="complete.obs"))
        mae = c(mae, mean( abs(psi[,cn[1]]-psi[,cn[3]]), na.rm=T))
        mae = c(mae, mean( abs(psi[,cn[2]]-psi[,cn[3]]), na.rm=T))
    }
    print(sprintf("Pearson correlations of variant PSI from %s: %s",sample[i,"date"],paste(rp,collapse=", ")))
    date_rp = c(date_rp, rp)
    
    print(sprintf("Mean absolute error of variant PSI from %s: %s",sample[i,"date"],paste(mae,collapse=", ")))
    date_mae = c(date_mae, mae)
}
print(sprintf("Replica correlations %.2f sd %.2f min %.2f max %.2f",mean(date_rp),sd(date_rp),min(date_rp),max(date_rp)))
print(sprintf("Replica MAE %.2f sd %.2f min %.2f max %.2f",mean(date_mae),sd(date_mae),min(date_mae),max(date_mae)))

facs_rp = c()
facs_mae = c()
for (facs in unique(facs_rep)) {
    cn = sprintf("run%s_facs%d",unique(sample[which(sample$facs_rep==facs),"date"]),facs)
    
    rp = cor(psi[,cn[1]], psi[,cn[2]], method="pearson", use="complete.obs")
    mae = mean( abs(psi[,cn[1]]-psi[,cn[2]]), na.rm=T)    
    if (length(cn) > 2) {
        rp = c(rp, cor(psi[,cn[1]], psi[,cn[3]], method="pearson", use="complete.obs"))
        rp = c(rp, cor(psi[,cn[2]], psi[,cn[3]], method="pearson", use="complete.obs"))
        mae = c(mae, mean( abs(psi[,cn[1]]-psi[,cn[3]]), na.rm=T))
        mae = c(mae, mean( abs(psi[,cn[2]]-psi[,cn[3]]), na.rm=T))
    }
    print(sprintf("Pearson correlations of variant PSI for FACS rep %d: %s",sample[i,"facs_rep"],paste(rp,collapse=", ")))
    facs_rp = c(facs_rp, rp)

    print(sprintf("Mean absolute error of variant PSI from %s: %s",sample[i,"date"],paste(mae,collapse=", ")))
    facs_mae = c(date_mae, mae)
}
print(sprintf("Replica correlations %.2f sd %.2f min %.2f max %.2f",mean(facs_rp),sd(facs_rp),min(facs_rp),max(facs_rp)))
print(sprintf("Replica MAE %.2f sd %.2f min %.2f max %.2f",mean(facs_mae),sd(facs_mae),min(facs_mae),max(facs_mae)))

# coverage, which single substitutions are present
i1 = which(psi$n_subst==1)   # & psi$mut_aa != "*")
agg = aggregate(psi[i1,"mut_aa"], list(psi[i1,"resi"]), paste0, collapse="")

residue = data.frame(resi=seq(nchar(wt[["aa"]])), wt=strsplit(wt[["aa"]],"")[[1]])
residue$mut_aa = ""
residue[agg[,1],"mut_aa"] = agg[,2]
residue$coverage = sapply(residue$mut_aa, nchar)
no_nonsense = psi$mut_aa != "*"
residue$median_psi = sapply(residue$resi, function(resi) { median(psi[which(psi$resi==resi & no_nonsense),"mean.psi"], na.rm=T) })
residue$median_vamp = sapply(residue$resi, function(resi) { median(psi[which(psi$resi==resi & no_nonsense),"mean.psi.norm"], na.rm=T) })


# Synonymous WT
# aggragate DNA variants of WT
wt_i = which(psi_bc$var=="")
agg = aggregate(wt_i, list(psi_bc[wt_i,"var_dna"]), paste, collapse=":")
var_indices = lapply( strsplit(agg[,2], ":"), as.numeric)
names(var_indices) = agg[,1]

# new data frame with WT DNA variants and PSI of all replicates
psi_wt_syn = data.frame(var=names(var_indices), n_bc=sapply(var_indices,length))
subst_list = strsplit(psi_wt_syn$var, ":")
psi_wt_syn$n_subst = sapply(subst_list, length)
for (i in which(sample$gate=="A" & sample$tech_rep==1)) {
    cn = sprintf("run%s_facs%d",sample[i,"date"],sample[i,"facs_rep"])
    print(sprintf("Aggregate %s",cn))
    psi_wt_syn[,cn] = sapply(var_indices, function(v) { mean(psi_bc[v,paste0(cn,".psi")], na.rm=T) })
    psi_wt_syn[,paste0("sd_",cn)] = sapply(var_indices, function(v) { sd(psi_bc[v,paste0(cn,".psi")], na.rm=T) })
}

# Average over replicates
cn = unique( sprintf("run%s_facs%d",dates,facs_rep) )
psi_wt_syn$mean.psi = apply(psi_wt_syn[,cn], MARGIN=1, mean, na.rm=T)
psi_wt_syn$sd.psi = apply(psi_wt_syn[,cn], MARGIN=1, sd, na.rm=T)

# Normalize synonymous WT
psi_wt_syn$mean.psi.norm = (psi_wt_syn$mean.psi - psi_nonsense) * scale
psi_wt_syn$sd.psi.norm = psi_wt_syn$sd.psi * scale


# save everything in R data file
save(psi, psi_bc, rpm_rm, raw_rm, sample, residue, wt, psi_nonsense, scale, psi_wt_syn, file="abundance.rda")

# dump vamp residue frame
write.table(residue, file="vamp_residues.csv", sep=";", row.names=F, quote=F)

# dump synonymous WT
psi_wt_syn[1,1] = "WT"
write.table(data.frame(var=psi_wt_syn$var, barcodes=psi_wt_syn$n_bc, vamp_score=psi_wt_syn$mean.psi.norm, vamp_std=psi_wt_syn$sd.psi.norm),
            file="wt_syn.csv", sep=";", row.names=F, quote=F)


# dump vamp scores only - biologist friendly excel-csv-format
i = which(psi$n_subst < 2)
df = data.frame(var=psi[i,"var"], barcodes=psi[i,"n_bc"], vamp_score=psi[i,"mean.psi.norm"], vamp_std=psi[i,"sd.psi.norm"])
stopifnot(df[1,1]=="")
df[1,1] = "WT"
write.table(df, file="abundance.csv", sep=";", row.names=F, quote=F)

# dump prism file
f = file("prism_mave_160_PRKN_vamp.txt", "wt")
write("# --------------------", f)
write("# version: 1", f)
write("# protein:", f)
write("#     name: PRKN", f)
write(sprintf("#     sequence: %s",wt[["aa"]]), f)
write("#     organism: Homo sapiens (Human)", f)
write("#     uniprot: O60260", f)
write("# mave:", f)
write("#     organism: Homo sapiens (Human)", f)
write("#     cloning: Landing pad", f)
write("#     expression: Overexpression", f)
write("#     technology: VAMPseq", f)
write("#     doi: unpublished", f)
write("#     year: 2022", f)
write("# variants:", f)
write(sprintf("#     number: %d",nrow(df)), f)
write(sprintf("#     coverage: %.2f",mean(residue$coverage > 0)), f)
write(sprintf("#     depth: %.2f",mean(nchar(gsub("*", "", residue[which(residue$coverage > 0),"mut_aa"], fixed=T)))), f)
write("#     width: single mutants", f)
write("# columns:", f)
write("#     barcodes: Number of different barcodes observed for variant", f)
write("#     vamp_score: VAMPseq score is protein stability index (PSI) normalized to wild-type and nonsense variants", f)
write("#     vamp_std: VAMPseq score standard deviation between 12 replicates (4 biological)", f)
write("# --------------------", f)
write("#", f)
write("# Unpublished version of March 2022 - kristoffer.johansson@bio.ku.dk", f)
write("# ", f)
write(sprintf("%5s  %8s  %10s  %8s", "var", "barcodes", "vamp_score", "vamp_std"), f)
write(sprintf("%5s  %8d  %10.4f  %8.4f",df$var, df$barcodes, df$vamp_score, df$vamp_std), f)
close(f)


# plot
quartz(width=8, height=6)
breaks = seq(-0.3, 1.6, 0.05)
h_wt = hist(psi_wt_syn[,"mean.psi.norm"], breaks=breaks, plot=F)
h_nons = hist(psi[which(psi$mut_aa=="*"),"mean.psi.norm"], breaks=breaks, plot=F)
h_lib = hist(psi[which(psi$mut_aa!="*" & psi$n_subst==1),"mean.psi.norm"], breaks=breaks, plot=F)
plot(0,0,col=0, xlim=c(-.1,1.3), ylim=c(0,10.0), xlab="VAMP score", ylab="Density")
lines(h_lib$mids, h_lib$density, col=1, lwd=2)
lines(h_nons$mids, h_nons$density, col=2, lwd=2)
lines(h_wt$mids, h_wt$density, col=3, lwd=2)
legend("top", c("WT syn.","Nonsense","Library"), lty=1, lwd=2, col=c(3,2,1))
quartz.save("vamp_distributions.png", type="png")
