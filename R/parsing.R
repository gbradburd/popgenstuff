
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
getBPstats <- function(stacksFAfile,outPath,minPropIndivsScoredin,checkforlowcovsamps=FALSE,sampstodrop=NULL){
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
  #pull out sample and locus ID info from info col into sep. cols.
  df <- df %>% tidyfast::dt_separate(., info, into = c("x","sample"), sep = " ", remove = F) %>% 
    dplyr::mutate(sample = gsub("\\[","",sample)) %>% 
    dplyr::mutate(sample = gsub("\\]","",sample)) %>% 
    tidyfast::dt_separate(., x, 
                          into = c("y","clocus","z","sample.internal","w","locus","b","allele"), 
                          sep = "_", remove = T) %>% 
    dplyr::select(-y,-z,-sample.internal,-w,-b)
  
  # REMOVE SAMPLES WITH FEW CO-GENOTYPED LOCI (or others from a list)
  if (length(sampstodrop)>0){
    df <- df %>% dplyr::filter(!sample %in% sampstodrop)
  } else {
    print("no samples were dropped for this dataset")
  }
  
  #add number of individuals genotyped per every locus (divide by 2 bc still a row for each of 2 alleles at this point)
  df <- df %>% dplyr::add_count(clocus, name = "n_samps_genoed") %>% mutate(n_samps_genoed = n_samps_genoed/2)
  #get total number of individuals in dataset
  Nindivs <- length(unique(df$sample))
  #filter to just loci scored in at least X proportion of indivs
  nindivsthreshold <- round(Nindivs*minPropIndivsScoredin,0)
  df <- df %>% dplyr::filter(n_samps_genoed >= nindivsthreshold)
  #get new total number of individuals in dataset after filtering
  Nindivs <- length(unique(df$sample))
  #get position of Ns in sequences
  positions <- df %>% dplyr::mutate(position = stringr::str_locate_all(string = sequence, pattern = "N")) %>% 
    dplyr::select(sample,locus,allele,position) %>% tidyr::unnest(.,position) %>% as.data.frame()
  positions <- positions %>% dplyr::mutate(position_of_N = positions$position[,1]) %>% dplyr::select(-position)
  positions <- positions %>% dplyr::select(-allele) %>% dplyr::distinct() %>% dplyr::mutate(N_ID = paste(locus,position_of_N,sep=":"))
  #print percent of total bp positions that are N in at least one indiv
  n_N_positions = positions %>% dplyr::distinct(locus,position_of_N) %>% nrow()
  n_bps = df %>% dplyr::mutate(n.positions = stringr::str_length(sequence)) %>% dplyr::distinct(clocus,n.positions) %>% dplyr::summarise(n=sum(n.positions))
  print(paste("Percent of base positions with N for at least one indiv (after removing lowcov indivs and loci in <",
              minPropIndivsScoredin*100,"% of indivs): ",
              round((n_N_positions/n_bps$n)*100,2)," %",
              sep=""))
  #keep just 1 haplotype per indiv
  df <- df %>% dplyr::filter(allele == 0) 
  
  # build matrix of number of unique positions that are N in one or both indivs
  positions <- positions %>% dplyr::select(N_ID,sample) %>% mutate(n = 1)
  #check for samples that didn't have any Ns
  sampswithnoNs = setdiff(df %>% dplyr::distinct(sample), positions %>% dplyr::distinct(sample))
  if(nrow(sampswithnoNs)!=0) {
    print("these samples have no Ns:")
    print(sampswithnoNs)
    sampswithnoNs <- sampswithnoNs %>% mutate(N_ID = "dummy", n = 0) %>% dplyr::select(N_ID, sample, n)
    positions <- rbind(positions,sampswithnoNs)
  }
  positions.wide <- stats::xtabs(n ~ ., positions) #1 = N in indiv, 0 = no N, but includes loci that were and were not actually genotyped in that indiv
  unionNs <- outer( 1:ncol(positions.wide), 1:ncol(positions.wide), Vectorize(function(i, j) sum(pmax(positions.wide[,i],positions.wide[,j])))) 
  colnames(unionNs) <- colnames(positions.wide)
  rownames(unionNs) <- colnames(positions.wide)
  # build matrix of number of bp genotyped in both indivs (including Ns)
  #count number of bps (aka positions in each sequence - including Ns)
  df <- df %>% dplyr::mutate(n.bp = stringr::str_length(sequence))
  #get matrix of number of cogenotyped basepairs for each pair of individuals
  #make long data into wide and complete matrix (fill in missing data with zeros)
  df <- df %>% dplyr::select(clocus,sample,n.bp)
  df.wide <- stats::xtabs(n.bp ~ ., df)
  #take the crossproduct aka multiply number of basepairs by 1 if locus is cogenotyped and 0 if it's not, sum across all loci
  coGenowithNs <- crossprod(df.wide, df.wide>0) #total number of bps at loci genoed in both indivs (counting Ns)
  #subtract of Ns from intial coGeno
  coGeno = coGenowithNs - unionNs
  #add back number N's that were in indiv i but not genoed in j and vice versa
  for(i in 0:(Nindivs-2)){
  	for(j in (i+1):(Nindivs-1)){
  	  #print(paste("i is", i, "and j is", j))
  		inSampI <- which(positions$sample == paste0("sample",i))
  		missingBPinI <- unlist(lapply(strsplit(positions$N_ID[inSampI],":"),"[[",1))
  		missingLociInJ <- names(which(df.wide[,paste0("sample",j)]==0))
  		toAdd <- length(which(missingBPinI %in% missingLociInJ))
  		inSampJ <- which(positions$sample == paste0("sample",j))
  		missingBPinJ <- unlist(lapply(strsplit(positions$N_ID[inSampJ],":"),"[[",1))
  		missingLociInI <- names(which(df.wide[,paste0("sample",i)]==0))
  		toAdd <- toAdd + length(which(missingBPinJ %in% missingLociInI))
  		coGeno[paste0("sample",i),paste0("sample",j)] <- coGeno[paste0("sample",i),paste0("sample",j)] + toAdd
  		coGeno[paste0("sample",j),paste0("sample",i)] <- coGeno[paste0("sample",i),paste0("sample",j)]
  	}
    print(paste("done with sample", (i+1), "of", Nindivs))
  }
  #add (back) number of individuals genotyped per every locus
  df <- df %>% dplyr::add_count(clocus, name = "n_samps_genoed")
  #get number of loci in dataset
  #this number matches that reported in populations.log at least for test dataset (with no filtering), yay!
  Nloci = df %>% dplyr::distinct(clocus) %>% nrow()
  #get list of loci that we kept (so we can use it later to filter gt)
  lociIDskept <- df %>% dplyr::distinct(clocus) %>% dplyr::mutate(clocus = paste("^",clocus,":",sep=""))
  sampIDskept <- df %>% dplyr::distinct(sample)
  #get total number of unique base pairs in dataset
  #this number matches that reported in populations.log at least for test dataset (with no filtering), yay!
  Nbp = df %>% dplyr::select(clocus,n.bp) %>% 
    dplyr::distinct() %>% 
    dplyr::summarise(n=sum(n.bp))
  Nbp = Nbp$n
  
  #get number of individuals genotyped per every locus (and number of base pairs in locus)
  n_indiv_per_locus <- df %>% 
    dplyr::select(-sample) %>% dplyr::distinct()
  #for each position with at least one N, get number of samples with N 
  positions.summed <- positions %>% dplyr::group_by(N_ID) %>% dplyr::summarise(n_samps_with_N=n()) %>% 
    tidyfast::dt_separate(., N_ID, into = c("clocus","pos"), sep = ":", remove = T)
  #per locus, get total number of positions with at least one N
  total.n.per.locus <- positions.summed %>% dplyr::group_by(clocus) %>% summarise(n_pos_with_N=n())
  #get df of locus lengths and number of indivs they are genoed in, subtracting off any positions with N
  n_indiv_per_locus.noNs <- merge(n_indiv_per_locus, total.n.per.locus, by = "clocus", all.x = T) %>% 
    dplyr::mutate(n_pos_with_N = tidyr::replace_na(n_pos_with_N, 0)) %>% 
    dplyr::mutate(n.bp.noN = n.bp - n_pos_with_N) %>% 
    dplyr::select(clocus,n_samps_genoed,n.bp.noN) %>% 
    dplyr::rename("n.bp" = "n.bp.noN")
  #get df of number of indivs that ARE NOT N for each position with an N
  n_indiv_per_pos.withNs <- merge(n_indiv_per_locus %>% dplyr::select(clocus,n_samps_genoed), 
                                  positions.summed %>% dplyr::mutate(n.bp = 1), by = "clocus", all.y = T) %>% 
    dplyr::mutate(n_samps_really_genoed = n_samps_genoed - n_samps_with_N) %>% 
    dplyr::select(clocus,pos,n_samps_really_genoed,n.bp) %>% 
    dplyr::rename("n_samps_genoed" = "n_samps_really_genoed") %>% 
    dplyr::mutate(clocus = paste(clocus,pos,sep=":")) %>% 
    dplyr::select(clocus,n_samps_genoed,n.bp)
  #combine
  n_indiv_per_locus <- rbind(n_indiv_per_locus.noNs,n_indiv_per_pos.withNs)
  #get number of bp genotyped for each number of individuals from 1:N
  lociDistn <- sapply(1:Nindivs,
                      function(i){
                        sum(n_indiv_per_locus$n.bp[which(n_indiv_per_locus$n_samps_genoed==i)])
                      })
  
	# do checks on data - if we want
  if(checkforlowcovsamps) {
	alerts <- rep(0,3)
	# check for absolutely low genotyped bp
	if(any(diag(coGeno) < 250)){
		alerts[1] <- 1
	}
	# check for relatively low genotyped bp
	if(any(diag(coGeno) < median(diag(coGeno))*0.3)){
		alerts[2] <- 1
	}
	# check for relatively low mean coGeno
	if(any(rowMeans(coGeno) < median(rowMeans(coGeno))*0.3)){
		alerts[3] <- 1
	}
	# print alert file
	if(any(alerts==1)){
		alertFile <- paste0(outPath,"_ALERT")
		file.create(alertFile)
		if(alerts[1]){
			problemChildren <- names(diag(coGeno))[which(diag(coGeno)<250)]
			alertMsg <- sprintf("these sample(s) have fewer than 250 genotyped base-pairs:\n %s \n",
											paste0(problemChildren,collapse="\n"))
			cat(alertMsg,file=alertFile,append=TRUE)
		}
		if(alerts[2]){
			problemChildren <- names(diag(coGeno))[which(diag(coGeno) < median(diag(coGeno))*0.3)]
			alertMsg <- sprintf("these sample(s) have less than 30 percent of the median number of genotyped base-pairs:\n %s \n",
											paste0(problemChildren,collapse="\n"))
			cat(alertMsg,file=alertFile,append=TRUE)
		}
		if(alerts[3]){
			problemChildren <- names(diag(coGeno))[which(rowMeans(coGeno) < median(rowMeans(coGeno))*0.3)]
			alertMsg <- sprintf("these sample(s) have less than 30 percent of the median number of co-genotyped base-pairs:\n %s \n",
											paste0(problemChildren,collapse="\n"))
			cat(alertMsg,file=alertFile,append=TRUE)
		}
	}
  }
  
  #save relevant outputs
  BPstats <- list("lociDistn" = lociDistn,
                  "coGeno" = coGeno,
                  "nLoci" = Nloci,
                  "Nbp" = Nbp,
                  "Nindivs" = Nindivs,
                  "lociIDskept" = lociIDskept,
                  "sampIDskept" = sampIDskept)
  save(BPstats,file=paste0(outPath,"_BPstats.Robj"))
  return(invisible("BP stats generated!"))
}


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
vcf2R <- function(vcfFile,readDepth=FALSE,outPath=NULL,minPropIndivsScoredin){
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
	#filter to loci aka SNP positions in X % of indivs
	nindivsthreshold <- round(nrow(gt)*minPropIndivsScoredin,0)
	gt.long <- gt %>% as.data.frame() 
	gt.long$sampid <- rownames(gt.long)
	gt.long <- gt.long %>% dplyr::select(sampid,tidyr::everything())
	gt.long <- gt.long %>% tidyr::pivot_longer(.,names_to="SNPid",values_to="genotype",cols=2:ncol(gt.long))
	totoss <- gt.long %>% dplyr::filter(is.na(genotype)==T) %>% dplyr::group_by(SNPid) %>% dplyr::summarise(n_NAs = dplyr::n()) %>% dplyr::filter(n_NAs > (nrow(gt) - nindivsthreshold))
	gt.f <- gt.long %>% dplyr::filter(!SNPid %in% totoss$SNPid) %>% tidyfast::dt_pivot_wider(., names_from=SNPid, values_from=genotype) %>% as.data.frame()
	row.names(gt.f) <- gt.f$sampid
	gt.f <- gt.f %>% dplyr::select(-sampid)
	gt <- gt.f %>% as.matrix()
	#calc mean read depths - if we want
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
	dp.long <- dp.long %>% dplyr::group_by(SNPid) %>% dplyr::mutate(n_NAs = sum(is.na(read_depth)==T))
	#filter
	nindivsthreshold <- round(nrow(dp)*minPropIndivsScoredin,0)
	dp.long <- dp.long %>% dplyr::filter(n_NAs <= nrow(dp) - nindivsthreshold)
	dp.mean <- dp.long %>% stats::na.omit() %>% dplyr::group_by(sampid) %>% dplyr::summarise(mean.dp = mean(read_depth))
	return(dp.mean)
}

getallReadDepths <- function(vcfFile){
  sampid <- read_depth <- NULL
  `%>%` <- magrittr::`%>%`
  vcf <- vcfR::read.vcfR(vcfFile, verbose = FALSE)
  dp <- vcfR::extract.gt(vcf, "DP", as.numeric = T)
  dp <- as.data.frame(t(dp), stringsAsFactors = FALSE)
  dp$sampid <- row.names(dp)
  dp <- dp %>% dplyr::select(sampid, tidyr::everything())
  #change dataframe from wide matrix to long list
  dp.long <- tidyr::pivot_longer(data = dp, names_to = "SNPid", values_to = "read_depth", cols = 2:ncol(dp))
  dp.long <- dp.long %>% dplyr::group_by(SNPid) %>% dplyr::mutate(n_NAs = sum(is.na(read_depth)==T))
  return(dp.long)
}

