# Reads in MEDIPS objects corresponding to refernce tissue identifies DMRs, which are used in subsequent steps

# Based on script by Keegan Korthauer

# this is derived from DMRs.R

suppressMessages(library(MEDIPS))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(annotatr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(ggdendro))
suppressMessages(library(matrixStats))
suppressMessages(library(stringr))
suppressMessages(library(circlize))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
source("scripts/limmamedips.R")

# TODO: see which of the above libraries can be removed

# TODO: remove variables that are no longer used (eg top)

# filelist = string with comma-separated list of medip object files
# lab1 = label for cases
# lab2 = label for controls
# sig.level = qvalue threshold for defining DMRs (optional, not used when top is specified)
# mdip.opt = logical whether to adjust for gc content bias in counts
# chrs = vector of chromosome names to include
# training = proportion of cases and tests to use as training set, if no validation sets provided  -- this is no longer used
# regulatory = logical whether to restrict to regulatory regions
# blacklist = bed file specifying regions to remove from consideration 
# top = threshold for number of DMRs (take top `top` DMRs ranked by qval) - will be appended to results files (heatmaps, results tables)
# colnames = logical whether to include sample names in heatmaps
# out.dir = directory location to save results
# iteration = a number that will be used to label subdirectory if running multiple iterations
# restrict_up = DMRs in the training data must overlap with this bed file and be UP in cases compared to controls
# restrict_down = DMRs in the training data must overlap with this bed file and be DOWN in cases compared to controls. Only used if restrict_up is set.
# depth_file = file with depth for each sample, will be included as a covariate in limma model if provided

# differential coverage
compute.diff <- function(filelist=files,
			metasheet=metasheet,
			lab1 = NULL, lab2 = NULL,
			sig.level = NULL,
			mdip.opt = FALSE,
			chrs =  paste0("chr", c(1:22, "X", "Y", "M")),
			regulatory = FALSE,
			blacklist = NULL,
			top = 300,
			colnames = TRUE,
			out.dir = "out", 
			sample_name = NULL,
			depth_file = NULL){

	set.seed(1)

	dir.create(out.dir, recursive=T)

        #read metasheet	
	met=read.table(metasheet, header=T, sep=",", as.is=T)
	if(!all(c("SampleName", "Class", "Type") %in% colnames(met))) stop("check metasheet file for required columns")
	met = subset(met, Class %in% c("reference_case", "reference_control"))
	if(nrow(met)==0) stop("no refernces_case or reference_control samples specified")

	#will batch info be included?
	if("Batch" %in% colnames(met)) show.batch=T

	#match files to samples in metasheet
        file=unlist(strsplit(filelist," "))
	SampleName = file %>% str_replace(".*/","") %>% str_replace(".medip.rds","") 
	#print(SampleName)
	print(met$SampleName)
	if(!(all(met$SampleName %in% SampleName))) stop("metasheet file does not match sample names")                     
                                                                                                                         
	message("excluding the following samples that are in the config file but not the metasheet file (if any):")       
	message(SampleName[!SampleName %in% met$SampleName])  
        
	met = merge(met, data.frame(cbind(SampleName, file)))  
        met$file = as.character(met$file)
 
        # load medips objects:                                                                                           
        message("reading in medips files") 
        obj1 <- lapply(met$file[met$Class=="reference_case"], readRDS) 
        obj2 <- lapply(met$file[met$Class=="reference_control"], readRDS) 
	print(obj2)
	# the coupling set contains counts of CG for each window, can be used for normalization
        if(file.exists("analysis/couplingset/couplingset.rds")) {
                CS = readRDS("analysis/couplingset/couplingset.rds")
        } else {
        	CS = MEDIPS.couplingVector(pattern = "CG", refObj = obj1[[1]]) 
        }

	#TODO: get read of anything here that isn't used:
        # set graph labels and file names
	diff.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".diff.rds"))
	bed.file <- file.path(out.dir, "reference_DMRs.bed")
	plot.q.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".q.lfc.pdf"))
	dmrs.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".dmrs.rds"))
	heatmap.file <- file.path(out.dir,paste0(lab1, ".", lab2, 
	  ".diff.heatmap.top", top,".pdf"))
	diff.file <- file.path(out.dir, paste0(lab1, ".", lab2, 
	  ".diff.reference.rds"))

        # identify DMRs:
        message("===== Identifying DMRs in reference set =====")

        message("Using the following cases:")
        l1 = unlist(lapply(obj1, function(x) {x@sample_name})) %>% str_replace(".dedup.bam", "") 
        message(paste(l1,delim=" "))                                                                                                  
        message("Using the following controls:")                                      
        l2 = unlist(lapply(obj2, function(x) {x@sample_name})) %>% str_replace(".dedup.bam", "") 
        message(paste(l2,delim=" ")) 

	n1=length(l1)
	n2=length(l2)

	if (!is.null(sig.level)){
		message("Using ", sig.level, " FDR cutoff instead of taking top ",
		  top, " regions.")
	}else{
		message("Using top ", top , " regions instead of FDR cutoff")
	}
	# uses a custom version of MEDIP.meth employing limma rather than ttest/edgeR

	if (!is.null(depth_file)){
		message("Including depth as covariates in the model")
		depth <- read.table(depth_file, stringsAsFactors=F, header=T)	
		#put depth dataframe in the same order as the medips objects
		depth <- depth[match(c(l1, l2), depth$sample),]
	} else {
		depth = NULL
	}
	diff = MEDIPS.meth(MSet1 = obj1, MSet2 = obj2,
	  CSet = CS, diff.method = "limma", chr = chrs,
	  p.adj = "BH", MeDIP = mdip.opt, minRowSum = 0.2*(n1+n2),
	  depth = depth)
	saveRDS(diff, file = diff.file)
	# diff has counts, rpkm, and adj p value for each window

	#NOTE: the way diff is calculated by MEDIPS.meth above, positive logFC means higher in CONTROLS (ie, MSet2 relative to MSet1)
	# reverse sign of logFC
	diff$logFC = -1 * diff$logFC

	# restrict to regulatory regions
	incl = !is.na(diff$P.Value)
	message("Trimming search regions. Starting with ", sum(incl), " windows")

	#restrict to windows with coupling factor > 0
	incl = incl & (CS@genome_CF > 0)
	message (sum(incl), " windows remain after restricting to windows with coupling factor > 0")

	if (regulatory){
      		# get cpg islands, shelves, and shores, plus enhancers
		annots <- c("hg19_cpgs", "hg19_enhancers_fantom")
		annotations = build_annotations(genome = 'hg19', annotations = annots)
		annotations <- annotations[annotations$type != "hg19_cpg_inter",]
		incl = incl & (makeGRangesFromDataFrame(diff) %over% annotations)
		message(sum(incl), " windows remain after restricting to CpGs/enhancers")
	}

	if (!is.null(blacklist)){
		bl=import(blacklist, format="BED")
		incl = incl & !(makeGRangesFromDataFrame(diff) %over% bl)
		message(sum(incl), " windows remain after removing blacklisted regions")
	}

	diff = diff[incl,]

	# recalculate the adjusted p values for differential methylation (originally calculated for all windows, not restrictedd set)
	diff$limma.adj.p.value = p.adjust(diff$P.Value, "BH")

        # plot p val vs logFC
        message("plotting pval vs logFC")
        plot = diff[!is.na(diff$limma.adj.p.value),]
	plot$logq = -1*log(plot$limma.adj.p.value,10)
	plot_notsig = subset(plot, logq < 6) 
	plot_notsig = plot_notsig[sample(1:nrow(plot_notsig),20000, replace=F),]

        ggplot()  + geom_point(data=subset(plot, logq > 6), aes(y=logq, x=logFC), color = "darkblue", pch = ".") +
          geom_point(data=plot_notsig, aes(y=logq, x=logFC), color="darkgray", pch = ".") +
          ylab("-log10 q-value") + xlab("log2 fold change (pos = higher in cases)") + theme_classic()
        ggsave(plot.q.file, height = 3, width=3)

        # create a bed file with significant DMRs at 0.005
        message("saving bed file with DMRs")
        d=diff[diff$limma.adj.p.value < 0.005 & !is.na(diff$limma.adj.p.value),]
        bed=data.frame(chr=d$chr,start=d$start,
          end=d$stop,case=d$MSets1.rpkm.mean,
          control=d$MSets2.rpkm.mean, LFC=d$logFC,Padj=d$limma.adj.p.value)
        write.table(bed,file=bed.file,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
	
        #also create bed files with only up regions and only down regions 
        write.table(bed[bed$LFC>0,1:3], paste0(bed.file,".case.up"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE) 
        write.table(bed[bed$LFC<0,1:3], paste0(bed.file,".case.down"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE) 
        
	
