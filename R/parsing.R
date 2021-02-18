#' Load a VCF file into R
#'
#' @param vcfFile The filename (with necessary path) of the 
#'				  VCF file you want to read into R.
#' @return A genotype matrix of 0s, 1s, and 2s denoting
#' 		 	the number of the counted allele in each individual 
#'			at each locus.  Missing data are indicated with \code{NA}.
#' @export
vcf2R <- function(vcfFile){
	`%>%` <- magrittr::`%>%`
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
	return(gt)
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
getBPstats <- function(stacksFAfile,outPath){
	. <- V1 <- allele <- b <- clocus <- info <- label <- locus <- n.bp <- n_basepairs_in_locus <- n_samps_genoed <- sample.internal <- w <- x <- y <- z <- NULL	
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
	df <- df %>% tidyr::separate(., info, into = c("x","sample"), sep = " ", remove = F) %>% 
				 dplyr::mutate(sample = gsub("\\[","",sample)) %>% 
				 dplyr::mutate(sample = gsub("\\]","",sample)) %>% 
				 tidyr::separate(., x, 
				  					into = c("y","clocus","z","sample.internal","w","locus","b","allele"), 
				  					sep = "_", remove = T) %>% 
				 dplyr::select(-y,-z,-sample.internal,-w,-b)
	#how many loci are in the dataset?
		#this number matches that reported in populations.log at least for test dataset, yay!
	Nloci = df %>% dplyr::group_by(clocus) %>% dplyr::summarise(n=dplyr::n()) %>% nrow()
	#how many total unique base pairs are in the dataset?
	#this number matches that reported in populations.log at least for test dataset, yay!
	Nbp = df %>% dplyr::filter(allele == 0) %>% 
				dplyr::select(-sequence,-allele,-locus) %>% 
				dplyr::group_by(clocus,n.bp) %>% 
				dplyr::summarise(n=dplyr::n()) %>% 
				dplyr::select(n.bp,clocus) %>% 
				dplyr::ungroup(.) %>% dplyr::summarise(n=sum(n.bp))
	Nbp = Nbp$n
	#keep just 1 haplotype per indiv and cols we want
	df <- df %>% dplyr::filter(allele == 0) %>% 
				dplyr::select(clocus,sample,n.bp)
	#make long data into wide and complete matrix (fill in missing data with zeros)
	df.wide <- stats::xtabs(n.bp ~ ., df)
	coGeno <- crossprod(df.wide, df.wide>0)				
	df.wide <- df.wide %>% 
				as.data.frame() %>% 
				dplyr::rename("n_basepairs_in_locus" = "Freq")
	#get number of individuals genotyped per every locus (and number of base pairs in locus)
	n_indiv_per_locus <- df.wide %>% 
									dplyr::filter(n_basepairs_in_locus > 0) %>% 
									dplyr::group_by(clocus) %>% 
									dplyr::mutate(n_samps_genoed = 1:dplyr::n()) %>% 
									dplyr::mutate(n_samps_genoed = max(n_samps_genoed)) %>% 
									dplyr::select(-sample) %>% dplyr::distinct()
	#total number of genotyped individuals
	N <- length(unique(df$sample))
	#number of bp genotyped for each number of individuals from 1:N
	lociDistn <- sapply(1:N,
					function(i){
						sum(n_indiv_per_locus$n_basepairs_in_locus[which(n_indiv_per_locus$n_samps_genoed==i)])
					})
	BPstats <- list("lociDistn" = lociDistn,
					"coGeno" = coGeno)
	save(BPstats,file=paste0(outPath,"_BPstats.Robj"))
	return(invisible("BP stats generated!"))
}