# Script that collects and analyse counts for different calls
options(width=200)

# load raw counts 
load("counts_zt.rda")
samples = data.frame(file=colnames(raw)[3:ncol(raw)])

# long complicated function to parse file names
parse_name = function(s) {
    ns = nchar(s)
    ret = list()
    ret$bio = 0
    ret$lib = "-"
    ret$facs = 0
    ret$sample = NA
    ret$lane = NA
    s_bio = substr(s,1,4)
    if (s_bio == "POI_") {
        ret$bio = "26-11-21"
	s_lib = substr(s,5,5)
	if (s_lib == "C") {
	    ret$lib="CT"
	    s_rest = substr(s,7,ns)
	} else {
	    ret$lib = s_lib
	    s_rest = substr(s,6,ns)
	}
	# print(sprintf("s_rest is %s",s_rest))
	s_facs = substr(s_rest,1,1)
	if (s_facs == "s") {
	    ret$facs = 2
	    s_rest = substr(s_rest,3,nchar(s_rest))
	} else {
	    ret$facs = 1
	    s_rest = substr(s_rest,2,nchar(s_rest))
	}
	# print(sprintf("s_rest is %s",s_rest))
	ret$gate = as.integer(substr(s_rest,1,1))
	s_tech = substr(s_rest,2,2)
        if (s_tech == "b") {
            ret$tech = 2
        } else {
            ret$tech = 1
        }
	ret$sample = strsplit(s_rest,"_")[[1]][2]
	
    } else if (s_bio == "14-6" | s_bio == "29-3") {
        if (s_bio == "14-6") { ret$bio = "14-6-22" } else { ret$bio = "29-3-22" }
	s_lib = substr(s,6,6)
	if (s_lib == "C") {
	    ret$lib="CT"
	    s_rest = substr(s,8,ns)
	} else {
	    ret$lib = s_lib
	    s_rest = substr(s,7,ns)
	}
	ret$gate = as.integer(substr(s_rest,1,1))
	s_facs = substr(s_rest,2,2)
	if (s_facs == "a") {
	    ret$facs = 1
	} else {
	    stopifnot(s_facs == "b")
	    ret$facs = 2
	}
        s_tech = substr(s_rest,3,3)
        if (s_tech == "r") {
            ret$tech = 2
        } else {
            ret$tech = 1
        }
        ret$sample = strsplit(s_rest,"_")[[1]][2]
        ret$lane = as.integer(substr(s,ns,ns))
    } else { return(NA) }
    
    return(ret)
}

# parse all names
l = lapply(samples$file, parse_name)

# assign name fields to all samples
samples$lib = sapply(l,"[[","lib")
samples$gate = sapply(l,"[[","gate")
samples$bio_rep = sapply(l,"[[","bio")
samples$facs_rep = sapply(l,"[[","facs")
samples$tech_rep = sapply(l,"[[","tech")
samples$lane = sapply(l,"[[","lane")
samples$sample = sapply(l,"[[","sample")

samples$reads = apply(raw[,samples$file], MARGIN=2, sum)

save(samples, file="samples.rda")

