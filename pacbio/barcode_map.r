options(width=160)

wt = list()
wt[["name"]] = "PARK2"
wt[["dna"]] = "ATGATCGTGTTTGTGCGCTTCAACAGCTCCCACGGCTTTCCAGTGGAGGTGGACTCTGATACCAGCATCTTCCAGCTGAAGGAGGTGGTGGCAAAGAGGCAGGGAGTGCCAGCCGACCAGCTGAGAGTGATCTTCGCCGGCAAGGAGCTGAGAAACGATTGGACAGTGCAGAATTGCGACCTGGATCAGCAGTCCATCGTGCACATCGTGCAGCGGCCCTGGAGAAAGGGACAGGAGATGAACGCAACCGGCGGCGACGATCCAAGGAATGCAGCAGGAGGATGTGAGAGGGAGCCTCAGTCTCTGACCAGAGTGGACCTGTCTAGCTCCGTGCTGCCTGGCGATAGCGTGGGCCTGGCCGTGATCCTGCACACCGACTCCAGGAAGGATTCTCCACCTGCAGGCAGCCCAGCAGGCCGGAGCATCTATAACTCCTTTTACGTGTATTGCAAGGGACCATGTCAGCGGGTGCAGCCTGGCAAGCTGAGAGTGCAGTGCAGCACCTGTAGGCAGGCCACACTGACCCTGACACAGGGCCCTTCCTGCTGGGACGATGTGCTGATCCCAAATAGAATGTCTGGCGAGTGCCAGAGCCCTCACTGTCCAGGCACATCCGCCGAGTTCTTTTTCAAGTGTGGCGCCCACCCCACCTCCGACAAGGAGACATCTGTGGCCCTGCACCTGATCGCCACCAACTCTCGGAATATCACATGCATCACCTGTACAGACGTGAGATCTCCTGTGCTGGTGTTTCAGTGCAACAGCAGGCACGTGATCTGCCTGGATTGTTTCCACCTGTATTGCGTGACCAGGCTGAATGACCGCCAGTTTGTGCACGATCCACAGCTGGGATACTCCCTGCCATGCGTGGCAGGATGTCCCAACTCTCTGATCAAGGAGCTGCACCACTTCAGGATCCTGGGCGAGGAGCAGTACAATCGCTATCAGCAGTACGGAGCAGAGGAGTGCGTGCTGCAGATGGGAGGCGTGCTGTGCCCACGCCCCGGCTGTGGAGCAGGCCTGCTGCCTGAGCCAGACCAGAGGAAGGTGACCTGCGAGGGAGGAAACGGCCTGGGATGTGGATTTGCCTTCTGCCGGGAGTGTAAGGAGGCCTATCACGAGGGCGAGTGCAGCGCCGTGTTTGAGGCATCCGGAACCACAACCCAGGCATACAGGGTGGATGAGAGAGCAGCAGAGCAGGCAAGATGGGAGGCAGCATCCAAGGAGACAATCAAGAAGACAACCAAGCCCTGCCCTAGGTGTCACGTGCCTGTGGAGAAGAATGGCGGCTGCATGCACATGAAGTGTCCACAGCCCCAGTGCAGACTGGAGTGGTGTTGGAACTGCGGCTGTGAGTGGAATAGGGTGTGCATGGGCGACCACTGGTTCGATGTG"
wt[["aa"]] = "MIVFVRFNSSHGFPVEVDSDTSIFQLKEVVAKRQGVPADQLRVIFAGKELRNDWTVQNCDLDQQSIVHIVQRPWRKGQEMNATGGDDPRNAAGGCEREPQSLTRVDLSSSVLPGDSVGLAVILHTDSRKDSPPAGSPAGRSIYNSFYVYCKGPCQRVQPGKLRVQCSTCRQATLTLTQGPSCWDDVLIPNRMSGECQSPHCPGTSAEFFFKCGAHPTSDKETSVALHLIATNSRNITCITCTDVRSPVLVFQCNSRHVICLDCFHLYCVTRLNDRQFVHDPQLGYSLPCVAGCPNSLIKELHHFRILGEEQYNRYQQYGAEECVLQMGGVLCPRPGCGAGLLPEPDQRKVTCEGGNGLGCGFAFCRECKEAYHEGECSAVFEASGTTTQAYRVDERAAEQARWEAASKETIKKTTKPCPRCHVPVEKNGGCMHMKCPQPQCRLEWCWNCGCEWNRVCMGDHWFDV"

infile = "bar_var_unq.txt"


##
## Helper functions
##

dna2aa = list()
# T          T                      C                    A                     G
dna2aa[["TTT"]] = "F"; dna2aa[["TCT"]] = "S"; dna2aa[["TAT"]] = "Y"; dna2aa[["TGT"]] = "C" # T
dna2aa[["TTC"]] = "F"; dna2aa[["TCC"]] = "S"; dna2aa[["TAC"]] = "Y"; dna2aa[["TGC"]] = "C" # C
dna2aa[["TTA"]] = "L"; dna2aa[["TCA"]] = "S"; dna2aa[["TAA"]] = "*"; dna2aa[["TGA"]] = "*" # A
dna2aa[["TTG"]] = "L"; dna2aa[["TCG"]] = "S"; dna2aa[["TAG"]] = "*"; dna2aa[["TGG"]] = "W" # G
# C
dna2aa[["CTT"]] = "L"; dna2aa[["CCT"]] = "P"; dna2aa[["CAT"]] = "H"; dna2aa[["CGT"]] = "R" # T
dna2aa[["CTC"]] = "L"; dna2aa[["CCC"]] = "P"; dna2aa[["CAC"]] = "H"; dna2aa[["CGC"]] = "R" # C
dna2aa[["CTA"]] = "L"; dna2aa[["CCA"]] = "P"; dna2aa[["CAA"]] = "Q"; dna2aa[["CGA"]] = "R" # A
dna2aa[["CTG"]] = "L"; dna2aa[["CCG"]] = "P"; dna2aa[["CAG"]] = "Q"; dna2aa[["CGG"]] = "R" # G
# A
dna2aa[["ATT"]] = "I"; dna2aa[["ACT"]] = "T"; dna2aa[["AAT"]] = "N"; dna2aa[["AGT"]] = "S" # A
dna2aa[["ATC"]] = "I"; dna2aa[["ACC"]] = "T"; dna2aa[["AAC"]] = "N"; dna2aa[["AGC"]] = "S" # C
dna2aa[["ATA"]] = "I"; dna2aa[["ACA"]] = "T"; dna2aa[["AAA"]] = "K"; dna2aa[["AGA"]] = "R" # A
dna2aa[["ATG"]] = "M"; dna2aa[["ACG"]] = "T"; dna2aa[["AAG"]] = "K"; dna2aa[["AGG"]] = "R" # G
# G
dna2aa[["GTT"]] = "V"; dna2aa[["GCT"]] = "A"; dna2aa[["GAT"]] = "D"; dna2aa[["GGT"]] = "G" # A
dna2aa[["GTC"]] = "V"; dna2aa[["GCC"]] = "A"; dna2aa[["GAC"]] = "D"; dna2aa[["GGC"]] = "G" # C
dna2aa[["GTA"]] = "V"; dna2aa[["GCA"]] = "A"; dna2aa[["GAA"]] = "E"; dna2aa[["GGA"]] = "G" # A
dna2aa[["GTG"]] = "V"; dna2aa[["GCG"]] = "A"; dna2aa[["GAG"]] = "E"; dna2aa[["GGG"]] = "G" # G

translate = function(dna) {
    n = nchar(dna)
    codons = substring(dna, seq(1, n-2, by=3), seq(3, n, by=3))
    if (length(codons)*3 != n) return(NA)
    paste0(dna2aa[codons], collapse="")
}

##
## Translate and compare to WT AA and DNA sequences
##

print(sprintf("Reading pacbio map %s for protein %s length %d amino acids",infile,wt[["name"]],nchar(wt[["aa"]])))
stopifnot(nchar(wt[["dna"]]) == 3*nchar(wt[["aa"]]))

# get unique pairs of barcodes and DNA variants, i.e. both barcode and variants lists are not unique by themself 
bc = read.table(infile)
colnames(bc) = c("pacbio","barcode","prot_dna")

print(sprintf("Translate %d variant DNA sequences to protein", length(bc$prot_dna)))
bc$prot_aa = sapply(bc$prot_dna, translate)

wt_dna_vec = strsplit(wt[["dna"]], "")[[1]]
wt_ndna = length(wt_dna_vec)

wt_aa_vec = strsplit(wt[["aa"]], "")[[1]]
wt_naa = length(wt_aa_vec)

get_subst = function(var, wt_vec) {
    var_vec = strsplit(var, "")[[1]]
    if (length(var_vec) != length(wt_vec)) return(NA)
    mut_resi = which(wt_vec != var_vec)
    paste0(wt_vec[mut_resi], mut_resi, var_vec[mut_resi])
}

# This is surprisingly fast!
subst_aa = lapply(bc$prot_aa, get_subst, wt_vec=wt_aa_vec)
subst_dna = lapply(bc$prot_dna, get_subst, wt_vec=wt_dna_vec)

stopifnot(all(! is.na(subst_aa)))
stopifnot(all(! is.na(subst_dna)))

bc$n_subst_aa = sapply(subst_aa, length)
bc$n_subst_dna = sapply(subst_dna, length)

# which barcodes correspondes to each variant
bc$var_aa = sapply(subst_aa, paste0, collapse=":")
bc$var_dna = sapply(subst_dna, paste0, collapse=":")

save(bc, file="bar_var_unq.rda")
# load("bar_var_unq.rda")

##
## Filtering on DNA substitutions
##

# Seems that 10 or more DNA substitutions is a resonable noise level (and anyway not interesting)
quartz(width=10, height=6)
par(mfrow=c(2,2), mar=c(4,4,2,2)+.1)
hist(bc$n_subst_dna, breaks=seq(0,1100)-.5)
hist(bc$n_subst_dna, breaks=seq(0,1100)-.5, ylim=c(0,200))
hist(bc$n_subst_dna, breaks=seq(0,1100)-.5, xlim=c(0,200), ylim=c(0,200))
hist(bc$n_subst_dna, breaks=seq(0,1100)-.5, xlim=c(0,20), ylim=c(0,200))
quartz.save("dna_substitutions.png", type="png")

# Remove barcode-variant pairs with this amount or more DNA substitutions in the variant
dna_subst_cut = 10

i_use = which(bc$n_subst_dna < dna_subst_cut)
i_rm = which(bc$n_subst_dna >= dna_subst_cut)
print(sprintf("Removing %d variants with %d or more DNA substitutions",length(i_rm),dna_subst_cut))
n = length(unique(bc$barcode))

n_filtered = length(unique(bc[i_use,"barcode"]))
print(sprintf("Unique barcodes reduced by %d (%.2f%%, %d -> %d)", n-n_filtered, (n-n_filtered)/n*100, n, n_filtered))

print("Distribution of amino acid substitutions among removed variants")
t = table(bc[i_rm,"n_subst_aa"])
print(t[which(as.numeric(names(t))<10)])

print("Distribution of number of variants per barcode BEFORE filtering on DNA substitutions")
print(table(table(bc$barcode)))
print("Distribution of number of variants per barcode AFTER filtering on DNA substitutions")
print(table(table(bc[i_use,"barcode"])))
# For PARK2, this filter heavily reduces the number of barcodes with 2 variants to 9.5% (48,089 -> 4,554)
# Now I only have 4,608 degenerata barcodes (8% of the previous 56,986)

# Distribution of pacbio reads
print("Distribution of pacbio reads BEFORE filtering on DNA substitutions")
t = table(bc[,"pacbio"])
print(t[order(as.numeric(names(t)))])
print("Distribution of pacbio reads AFTER filtering on DNA substitutions")
t = table(bc[i_use,"pacbio"])
print(t[order(as.numeric(names(t)))])
# There are still very many barcodes (t[1]/sum(t)*100 = 11%, before filter 29%) with only one read so dont filter on this

##
## Map of unique barcodes
##

# apply filter
bc_use = bc[i_use,]

# indices of bc dataframe for unique barcodes
bc_unq = aggregate(seq_along(bc_use$barcode), list(bc_use$barcode), paste0, collapse=":")
colnames(bc_unq) = c("barcode","bi")

# data.frame of unique barcodes (name to be filled later)
bc_map = data.frame(name = "bcx", barcode = bc_unq$barcode,
                    n = sapply(strsplit(bc_unq$bi,":"), length), bc_index = bc_unq$bi)

# maximum pacbio counts for each unique barcode
agg_max = aggregate(bc_use$pacbio, list(bc_use$barcode), max)
bc_map$pacbio_max = agg_max[match(bc_map$barcode,agg_max[,1]),2]

# sum of pacbio counts for each unique barcode
agg_sum = aggregate(bc_use$pacbio, list(bc_use$barcode), sum)
bc_map$pacbio_sum = agg_sum[match(bc_map$barcode,agg_sum[,1]),2]

# bc_index for variant with maximum pacbio counts for each unique barcode
bc_map$bi_max = bc_map$bc_index
i_degen = which(bc_map$n > 1)
bc_max = sapply(strsplit(bc_map[i_degen,"bc_index"], ":"),  function(v) { i=which.max(bc_use[as.numeric(v),"pacbio"]); v[i] })
bc_map[i_degen,"bi_max"] = bc_max

i = as.numeric(bc_map$bi_max)
bc_map$var_aa      = bc_use[i,"var_aa"]
bc_map$var_dna     = bc_use[i,"var_dna"]
bc_map$n_subst_aa  = bc_use[i,"n_subst_aa"]
bc_map$n_subst_dna = bc_use[i,"n_subst_dna"]
print(sprintf("Distribution of amino acid substitutions among %d barcodes in final map",nrow(bc_map)))
print(table(bc_map$n_subst_aa))

print("Distribution of pacbio dna variants per barcode:")
print(table(bc_map$n))
print(sprintf("Using variant with most reads for %d degenerate barcodes",length(i_degen)))

# reorder bc_map and name barcodes
resi = sapply(strsplit(bc_map$var_aa,":"), function(l) {if (length(l) > 0) as.numeric(substr(l[1],2,nchar(l[1])-1)) else 0})
mut_aa = sapply(strsplit(bc_map$var_aa,":"), function(l) {if (length(l) > 0) {nc=nchar(l[1]); substr(l[1],nc,nc)} else "wt"})
i = order(bc_map$n_subst_aa, resi, mut_aa, decreasing = FALSE)
bc_map = bc_map[i,]
bc_map$name = sprintf("bc%07d",seq(nrow(bc_unq))-1)

# barcodes of same substitution
print("Distribution of barcodes per unique amino acid variant")
t_var = table(bc_map$var_aa)
print(table(t))

print(sprintf("Made map of %d unique barcodes mapping to %d unique amino acid variants",nrow(bc_map),length(t_var)))

# dump data frame
save(bc_map, wt, dna_subst_cut, file="barcode_map.rda")

# dump text file
dump_df = data.frame(name=bc_map$name, barcode=bc_map$barcode, counts_pacbio=bc_map$pacbio_max, n_pacbio_var=bc_map$n,
                     subst=bc_map$var_aa, n_subst_aa=bc_map$n_subst_aa, n_subst_dna=bc_map$n_subst_dna)
write.table(dump_df, file="barcode_map.csv", sep=";", row.names=F, quote=F)

# dump FASTA file for BWA reference
fh = file("barcode_map.fasta", "w")
write(sprintf(">%s %s\n%s", bc_map$name, bc_map$var_aa, bc_map$barcode), fh)
close(fh)

