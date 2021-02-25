#' Calculate individual heterozygosity
#'
#' @param gt A genotype matrix with \emph{N} rows and 
#'			 \emph{L} columns, where \emph{N} is the number 
#'			 of individuals and \emph{L} is the number of 
#'			 loci, for which the \emph{i},\emph{j}th element 
#'			 gives the number of the counted allele in each 
#'			 individual at each locus.
#' @param nLoci An optional argument giving the number of genotyped 
#' 			base pairs in each individual. If not specified, 
#'			this is assumed to be the number of non-missing 
#'			genotypes in each row of \code{gt}. Default is 
#'			\code{NULL}.
#' @return A vector giving the proportion of loci that are 
#'			heterozygous in each individual.
#' @details The data matrix \code{gt} should consist of 
#'			0s, 1s, and 2s, with missing data indicated with 
#'			\code{NA}. This can be generated from a VCF file 
#'			using \code{\link{vcf2R}}. If \code{nLoci} is not 
#'			specified, and \code{gt} consists only of 
#'			polymorphic loci, this function will return 
#'			heterozygosity at polymorphic sites, which may not 
#'			be comparable across datasets.
#' @export
calcHet <- function(gt,nLoci=NULL){
	if(is.null(nLoci)){
		nLoci <- apply(gt,1,function(x){length(which(!is.na(x)))})
	}
	if(length(nLoci)!=nrow(gt)){
		stop("you must specify one value of nLoci for every row in your genotype matrix")
	}
	het <- unlist(lapply(1:nrow(gt),function(i){length(which(gt[i,]==1))/nLoci[i]}))
	return(het)
}


#' Run a principal components analysis
#'
#' @param gt A genotype matrix with \emph{N} rows and 
#'			 \emph{L} columns, where \emph{N} is the number 
#'			 of individuals and \emph{L} is the number of 
#'			 loci, for which the \emph{i},\emph{j}th element 
#'			 gives the number of the counted allele in each 
#'			 individual at each locus.
#' @param nPCs The number of principal component axes to retain.
#'			Default is 4.
#' @return A matrix with \emph{N} rows and \code{nPCs} columns 
#'			giving the position of each of the \emph{N} samples on 
#'			the first \code{nPCs} principal component axes.
#' @export
doPCA <- function(gt,nPCs = 4){
	sampleCov <- stats::cov(t(gt)/2,use="pairwise.complete.obs")
	pcAxes <- eigen(sampleCov)$vectors[,2:(nPCs+1)]
	row.names(pcAxes) <- row.names(gt)
	return(pcAxes)
}

detectPolymorphism <- function(genos){
	poly <- FALSE
	genos <- unique(genos[!is.na(genos)])
	if(length(genos) > 1 | any(genos==1)){
		poly <- TRUE    
	}
	return(poly)
}

calcThetaWsubSet <- function(gtSubset,n,L){
	k <- 2*n
	a_k <- sum(1/(1:(k-1)))
	S_k <- length(which(apply(gtSubset,2,detectPolymorphism)))
	thetaW <- S_k / (L*a_k)
	return(thetaW)
}

#' Calculate Wu and Watterson's Theta
#'
#' @param gt A genotype matrix with \emph{N} rows and 
#'			 \emph{L} columns, where \emph{N} is the number 
#'			 of individuals and \emph{L} is the number of 
#'			 loci, for which the \emph{i},\emph{j}th element 
#'			 gives the number of the counted allele in each 
#'			 individual at each locus.
#' @param lociDistn A vector of length \emph{N}, 
#'			  for which the \emph{i}th element 
#'			  gives the number of base pairs genotyped in exactly
#'			  \emph{i} individuals. E.g., the 3rd element of 
#'			  \code{lociDistn} gives the number of base pairs
#'			  genotyped in exactly 3 individuals. Default is 
#'			  \code{NULL}.
#' @return An estimate of Wu and Watterson's Theta.
#' @details The data matrix \code{gt} should consist of 
#'			0s, 1s, and 2s, with missing data indicated with 
#'			\code{NA}. This can be generated from a VCF file 
#'			using \code{\link{vcf2R}}, and the \code{lociDistn}
#'			can be generated using \code{\link{getBPstats}}.
#'			If \code{lociDistn} is not specified, 
#'			and \code{gt} consists only of polymorphic
#'			loci, this function will return Wu and 
#'			Watterson's Theta at polymorphic sites, 
#'			which may not be comparable across datasets.
#' @export
calcThetaW <- function(gt,lociDistn=NULL){
	N <- nrow(gt)
	if(is.null(lociDistn)){
		lociDistn <- rep(N,nrow(gt))
	}
	mdL <- apply(gt,2,function(x){N - length(which(is.na(x)))})
	thetaWs <- lapply(N:1,
					function(n){
						if(any(mdL == n)){
							gtSubset <- gt[,which(mdL == n),drop=FALSE]
							L <- lociDistn[n]
							tw <- calcThetaWsubSet(gtSubset,n,L)
						}
					})
	return(mean(unlist(thetaWs)))
}

#' Calculate pairwise pi between two diploid individuals
#'
#' @param ind1 A genotype vector for a single diploid individual
#'			 for which the \emph{i}th element 
#'			 gives the frequency of the counted 
#'			 allele in that individual at the \emph{i}th locus.
#' @param ind2 A genotype vector for a single diploid individual
#'			 for which the \emph{i}th element 
#'			 gives the frequency of the counted 
#'			 allele in that individual at the \emph{i}th locus.
#' @param L The number of loci genotyped in \emph{both} individuals 1 
#'			and 2.
#' @return An estimate of pairwise pi between individuals 1 and 2.
#' @details The genotype data should consist of {0,0.5,1}
#'			with missing data indicated with \code{NA}.
#'			If the genotype vectors \code{ind1} and \code{ind2} 
#'			consist only of polymorphic loci, and \code{L} is the 
#'			number of co-genotyped polymorphic loci, 
#'			this function will return pairwise pi at 
#'			polymorphic sites, which may not be 
#'			comparable across datasets.
#' @export
calcPWP <- function(ind1,ind2,L){
	diff.homs = sum(abs(ind1-ind2)==1,na.rm=TRUE)
	hets = sum(ind1==0.5 | ind2==0.5,na.rm=TRUE)
	return((diff.homs + hets/2)/L)
}

#' Calculate pairwise pi between all pairs of individuals in a dataset
#'
#' @param freqs  A genotype matrix with \emph{N} rows and 
#'			 \emph{L} columns, where \emph{N} is the number 
#'			 of individuals and \emph{L} is the number of 
#'			 loci, for which the \emph{i},\emph{j}th element 
#'			 gives the frequency of the counted allele in each 
#'			 individual at each locus.
#' @param coGeno A symmetric matrix (dimensions \emph{N} x \emph{N}) 
#'			for which the \emph{i},\emph{j}th element gives 
#'			the number of base pairs genotyped in \emph{both} 
#'			samples \emph{i} and \emph{j}. Default is \code{NULL}.
#'			If left unspecified, each cell will be assigned a value of 
#'			\code{L} (i.e., indicating no missing data).
#' @param quiet A \code{logical} value that controls whether 
#'			a progress par is displayed. Default is \code{FALSE}.
#' @return An estimate of pairwise pi between all individuals 
#'			in the data matrix.
#' @details The genotype data should consist of {0,0.5,1}
#'			with missing data indicated with \code{NA}.
#'			If the genotype matrix consists only of 
#'			polymorphic loci, and \code{coGeno} either 
#'			contains the number of co-genotyped polymorphic 
#'			loci for each sample pair or is not specified 
#'			this function will return pairwise pi at 
#'			polymorphic sites, which may not be 
#'			comparable across datasets.
#' @export
freqs2pairwisePi <- function(freqs,coGeno=NULL,quiet=FALSE){
	if(is.null(coGeno)){
		coGeno <- matrix(ncol(freqs),nrow(freqs),nrow(freqs))
	}
	n <- nrow(freqs)
	pwp <- matrix(NA,n,n)
	if(!quiet){
		prog <- utils::txtProgressBar(min=0,max=n+(n*(n-1))/2,char="*",style=3)
	}
	for(i in 1:n){
		for(j in i:n){
			if(!quiet){
				utils::setTxtProgressBar(prog,(i-1)*n+j-(i*(i-1)/2))
			}
			pwp[i,j] <- calcPWP(freqs[i,],freqs[j,],coGeno[i,j])
			if(i == j){
				pwp[i,i] <- 2 * pwp[i,i]
			} else if(i != j){
				pwp[j,i] <- pwp[i,j]
			}
		}
	}
	row.names(pwp) <- row.names(coGeno)
	colnames(pwp) <- row.names(coGeno)
	return(pwp)
}

#' Calculate a bunch of popgen stats from a VCF file
#'
#' @param vcfFile The filename (with necessary path) of the 
#'				  VCF file you want to read into R.
#' @param lociDistn A vector of length \emph{N}, where \emph{N} is the 
#'			  number of genotyped individuals, for which the
#'			  \emph{i}th element gives the number of base pairs
#'			  genotyped in exactly \emph{i} individuals.
#'			  E.g., the 3rd element of \code{lociDistn} gives 
#'			  the number of base pairs genotyped in exactly 3 
#'			  individuals. Default is \code{NULL}.
#' @param coGeno A symmetric matrix (dimensions \emph{N} x \emph{N}) 
#'			for which the \emph{i},\emph{j}th element gives 
#'			the number of base pairs genotyped in \emph{both} 
#'			samples \emph{i} and \emph{j}. The order of samples in 
#'			this matrix must be the same as the order of samples 
#'			in the VCF file. If left unspecified, each cell will 
#'			be assigned a value of \code{L} (i.e., indicating 
#'			no missing data). Default is \code{NULL}.
#' @param outPath The filename (with necessary path) of the 
#'					output object you want to generate.
#' @param nPCs The number of principal component axes to retain.
#'			Default is 4.
#' @return This function does not return a value. Instead, a 
#'			named list is saved as a .Robj file at the location 
#'			designated by the \code{outPath} argument. The elements 
#'			of the list are:
#'			\itemize{
#'				\item \code{thetaW} Wu and Watterson's theta, 
#'					as calculated by \code{\link{calcThetaW}}
#'				\item \code{pwp} a matrix of pairwise pi between 
#'					all samples, as calculated by \code{\link{freqs2pairwisePi}}
#'				\item \code{globalPi} global pi in the dataset 
#'					(the mean of the upper-triangle of \code{pwp})
#'				\item \code{pcs} the first \code{nPCs} principal 
#'					components, as calculated by \code{\link{doPCA}}
#'				\item \code{het} proportion of heterozygous loci in 
#'					each individual, as calculated by \code{\link{calcHet}}.
#'			}
#' @export
doAllTheThings <- function(vcfFile,lociDistn=NULL,coGeno=NULL,nPCs=4,outPath){
	gt <- vcf2R(vcfFile)
	thetaW <- calcThetaW(gt,lociDistn)
	pwp <- freqs2pairwisePi(freqs=gt/2,coGeno,quiet=FALSE)
	globalPi <- mean(pwp[upper.tri(pwp,diag=TRUE)])
	pcs <- doPCA(gt,nPCs)
	if(is.null(coGeno)){
		nLoci <- NULL
	} else {
		nLoci <- diag(coGeno)
	}
	het <- calcHet(gt,nLoci=nLoci)
	popgenstats <- list("thetaW" = thetaW,
						"pwp" = pwp,
						"globalPi" = globalPi,
						"pcs" = pcs,
						"het" = het)
	save(popgenstats,file=paste0(outPath,"_popgenstats.Robj"))
	return(invisible("popgen stuff done!"))
}