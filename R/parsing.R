#' Load a VCF file into R
#'
#' @param vcfFile The filename (with necessary path) of the 
#'				  VCF file you want to read into R.
#' @param readDepth A Boolean argument indicating whether or not 
#'					to calculate the mean read depth for each individual 
#'					from the specified VCF file.  Default is \code{FALSE}.
#'				  VCF file you want to read into R.
#' @param outPath The file path prepended to all output objects.
#' @return A genotype matrix of 0s, 1s, and 2s denoting
#' 		 	the number of the counted allele in each individual 
#'			at each locus.  Missing data are indicated with \code{NA}.
#' @export
vcf2R <- function(vcfFile,readDepth=FALSE,outPath=NULL){
	`%>%` <- magrittr::`%>%`
	if(readDepth & is.null(outPath)){
		stop("\nyou must specify a filepath to save the read depth information\n")
	}
	vcf <- vcfR::read.vcfR(vcfFile, verbose = FALSE)
	#extract genotypes
	gt <- vcfR::extract.gt(vcf)
	# transpose and make into data.frame
	#NA here means SNP not called
	gt <- as.data.frame(t(gt), stringsAsFactors = FALSE) %>% as.matrix()
	#check which syntax is used to denote genotypes
	genos <- if(any(grepl("/",gt))){
				c("0/0","0/1","1/0","1/1")
			 } else if(any(grepl("|",gt))){
			 	c("0\\|0","0\\|1","1\\|0","1\\|1")
			 }
	gt <- gsub(genos[1],0,gt)	
	gt <- gsub(genos[2],1,gt)
	gt <- gsub(genos[3],1,gt)
	gt <- gsub(genos[4],2,gt)
	class(gt) <- "numeric"
	if(readDepth){
		meanDepth <- getReadDepth(vcf)
		utils::write.table(meanDepth,file=paste0(outPath,"_readDepth.txt"),row.names=FALSE,col.names=c("sampid","meanReadDepth"))
	}
	return(gt)
}

getReadDepth <- function(vcf){
	sampid <- read_depth <- NULL
	`%>%` <- magrittr::`%>%`
	dp <- vcfR::extract.gt(vcf, "DP", as.numeric = T)
	dp <- as.data.frame(t(dp), stringsAsFactors = FALSE)
	dp$sampid <- row.names(dp)
	dp <- dp %>% dplyr::select(sampid, tidyr::everything())
	#change dataframe from wide matrix to long list
	dp.long <- tidyr::pivot_longer(data = dp, names_to = "SNPid", values_to = "read_depth", cols = 2:ncol(dp))
	dp.mean <- dp.long %>% stats::na.omit() %>% dplyr::group_by(sampid) %>% dplyr::summarise(mean.dp = mean(read_depth))
	return(dp.mean)
}

#' Generate genotyped base-pair statistics
#'
#' @param stacksFAfile The filename (with necessary path) of the 
#'				  .fa file output by a STACKS run.
#' @param outPath The filename (with necessary path) of the 
#'					output object you want to generate.
#' @return This function does not return a value. Instead, a 
#'			named list is saved as a .Robj file at the location 
#'			designated by the \code{outPath} argument. The elements 
#'			of the list are:
#'			\itemize{
#'				\item \code{lociDistn} vector of length \emph{N}, 
#'						where \emph{N} is total number of individuals 
#'						in the dataset, for with the \emph{i}th element 
#'						gives the number of base pairs genotyped in exactly
#'						\emph{i} individuals. E.g., the 3rd element of 
#'						\code{lociDistn} gives the number of base pairs
#'						 genotyped in exactly 3 individuals.
#'				\item \code{coGeno} a symmetric matrix (dimensions N x N) 
#'						for which the \emph{i},\emph{j}th element gives 
#'						the number of base pairs genotyped in \emph{both} 
#'						samples \emph{i} and \emph{j}.
#'			}
#' @export
getBPstats <- function(stacksFAfile,outPath,propScoredinfilter){
	. <- V1 <- allele <- b <- clocus <- info <- label <- locus <- n.bp <- n_basepairs_in_locus <- n_samps_genoed <- sample.internal <- w <- x <- y <- z <- df.wide <- NULL	
	`%>%` <- magrittr::`%>%`
	df <- utils::read.table(stacksFAfile, header = F, skip = 1, sep = "\n")
	#add two dummy columns so we can rearrange single column of alternating data into two side-by-side cols
	df <- df %>% dplyr::mutate(label = as.character(rep(1:2, nrow(.)/2))) %>% 
				 dplyr::mutate(dummy.id = rep(1:(nrow(.)/2), each=2))
	#build the 2 cols
	df <- df %>% 
			 	tidyr::pivot_wider(names_from = label, values_from = V1) %>% 
			 	dplyr::rename("info" = "1") %>% 
			 	dplyr::rename("sequence" = "2")
	#count number of basepairs in each sequence
	df <- df %>% dplyr::mutate(n.bp = stringr::str_length(sequence))
	#pull out sample and locus ID info from info col into sep. cols.
	df <- df %>% tidyfast::dt_separate(., info, into = c("x","sample"), sep = " ", remove = F) %>% 
				 dplyr::mutate(sample = gsub("\\[","",sample)) %>% 
				 dplyr::mutate(sample = gsub("\\]","",sample)) %>% 
				 tidyfast::dt_separate(., x, 
				  					into = c("y","clocus","z","sample.internal","w","locus","b","allele"), 
				  					sep = "_", remove = T) %>% 
				 dplyr::select(-y,-z,-sample.internal,-w,-b)
	#keep just 1 haplotype per indiv and cols we want
	df <- df %>% dplyr::filter(allele == 0) %>% 
	  dplyr::select(clocus,sample,n.bp)
	#add number of individuals genotyped per every locus
	df <- df %>% dplyr::add_count(clocus, name = "n_samps_genoed")
	#get total number of individuals in dataset
	Nindivs <- length(unique(df$sample))
	#filter to just loci scored in at least X proportion of indivs
	nindivsthreshold <- round(Nindivs*propScoredinfilter,0)
	df <- df %>% dplyr::filter(n_samps_genoed>=nindivsthreshold)
	#get number of loci in dataset
		#this number matches that reported in populations.log at least for test dataset, yay!
	Nloci = df %>% dplyr::distinct(clocus) %>% nrow()
	#get total number of unique base pairs in dataset
		#this number matches that reported in populations.log at least for test dataset, yay!
	Nbp = df %>%  dplyr::select(clocus,n.bp) %>% 
				dplyr::distinct() %>% 
				dplyr::summarise(n=sum(n.bp))
	Nbp = Nbp$n
	#get number of individuals genotyped per every locus (and number of base pairs in locus)
	n_indiv_per_locus <- df %>% 
	  dplyr::select(-sample) %>% dplyr::distinct()
	#get number of bp genotyped for each number of individuals from 1:N
	lociDistn <- sapply(1:Nindivs,
	                    function(i){
	                      sum(n_indiv_per_locus$n.bp[which(n_indiv_per_locus$n_samps_genoed==i)])
	                    })
	#get matrix of number of cogenotyped basepairs for each pairs of individuals
	#make long data into wide and complete matrix (fill in missing data with zeros)
	df <- df %>% dplyr::select(-n_samps_genoed)
	df.wide <- stats::xtabs(n.bp ~ ., df)
	#take the crossproduct aka multiply number of basepairs by 1 if locus is cogenotyped and 0 if it's not, sum across all loci
	coGeno <- crossprod(df.wide, df.wide>0)
	
	BPstats <- list("lociDistn" = lociDistn,
					"coGeno" = coGeno,
					"nLoci" = Nloci,
					"Nbp" = Nbp,
					"Nindivs" = Nindivs)
	save(BPstats,file=paste0(outPath,"_BPstats.Robj"))
	return(invisible("BP stats generated!"))
}