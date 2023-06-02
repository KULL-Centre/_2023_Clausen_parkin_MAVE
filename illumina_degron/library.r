options(width=180, stringsAsFactors=F)

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
## Assemble library
##

library = read.table("ref.seq")

colnames(library) = c("name","dna")
stopifnot(length(library$dna) == length(unique(library$dna)))

# Translate to protein sequence
print(sprintf("Translate %d variant DNA sequences to protein", length(library$dna)))
library$aa = sapply(library$dna, translate)

# assign tiles to odd, even or ct sublibraries
library$pool = "E" 
library[which(seq(nrow(library)) %% 2 == 1),"pool"] = "O"
ct_names = c("CYP2C9.0040","PTEN.0033","TPMT.0020","PARK2.0038","ASPA.0026","DHFR.0015","NUDT15.0013") 
library[match(ct_names,library$name),"pool"] = "CT"
library[which(grepl("APPY",library$name)),"pool"] = "all"

split_list = strsplit(library$name, ".", fixed=T)

library$protein = sapply(split_list, "[[", 1)
library$tile = sapply(split_list, function(v){ if (length(v)>1) {as.integer(v[2])} else {NA} })

##
## Assemble tiles into proteins
##

assemble_tiles = function(tiles) {
    # assume tiles are ordered and only have a single exact match of at least 4 amino acids
    full_seq = tiles[1]
    ret_df = data.frame(tile=tiles[1], nres=nchar(tiles[1]), tile_first=1, tile_last=nchar(tiles[1]), shift=0)
    for (it in seq(2,length(tiles))) {
        tile = tiles[it]
	shift = NA
        for (ir in seq(4,nchar(tile))) {
	    full_ct = substr(full_seq, nchar(full_seq)-ir+1, nchar(full_seq))
	    tile_nt = substr(tile, 1, ir)
	    if (tile_nt == full_ct) {
	        if (is.na(shift)) {
		    shift = ir
		} else {
		    print(sprintf("ERROR: Tile %d has two matches at shift %d and %d",it,shift,ir))
		    print(sprintf("CT of full seq: %s", full_ct))
		    print(sprintf("NT of tile:   : %s", tile_nt))
		    stopifnot(F)
		}
	    }
	}
	if (is.na(shift)) {
	    print(sprintf("ERROR: No match found for tile %d: %s",it,tile))
	    print(sprintf("Full seq: %s",full_seq))
	    stopifnot(F)
	} else {
	    nc = nchar(tile)
	    full_seq = paste0(full_seq, substr(tile, shift+1, nc))
	    ncf = nchar(full_seq)
	    ret_df = rbind(ret_df, data.frame(tile=tile, nres=nc, tile_first=ncf-nc+1, tile_last=ncf, shift=shift))
	}
    }
    
    # minor checks
    stopifnot(ret_df[nrow(ret_df),"tile_last"] == nchar(full_seq))
    nres_tab = table(ret_df$nres) 
    if (length(nres_tab) > 1) {
        print(sprintf("WARNING: Not all tiles have different length:"))
	print(nres_tab)
    }

    # return list
    ret = list()
    ret$df = ret_df
    ret$full_seq = full_seq
    return( ret )
}

prot_names = c()
prot_seqs = c()
library$tile_first = NA
library$tile_last = NA
for (protein in unique(library$protein)) {
    it = which(library$protein == protein)
    if (length(it) > 1) {
        at = assemble_tiles( library[it,"aa"] )
        prot_names = c(prot_names, protein)
        prot_seqs = c(prot_seqs, at$full_seq)
        stopifnot( all( at$df$tile == library[it,"aa"] ) )
        library[it,"tile_first"] = at$df$tile_first
        library[it,"tile_last"] = at$df$tile_last
    }
}
proteins = data.frame(name=prot_names, sequence=prot_seqs)

save(library, proteins, file="library.rda")

# Check
for (i in seq(nrow(library))) {
    ip = which(proteins$name == library[i,"protein"])
    if (length(ip) == 1) {
        tile = substr(proteins[ip,"sequence"], library[i,"tile_first"], library[i,"tile_last"])
	if (tile != library[i,"aa"]) {
	    print(sprintf("ERROR: Bad tile at row index %4d protein %6s: %s vs %s",i,library[i,"protein"],tile,library[i,"aa"]))
	}
    } else {
        print(sprintf("Cannot check protein not in proteins data frame: %s",library[i,"protein"]))
    }
}
