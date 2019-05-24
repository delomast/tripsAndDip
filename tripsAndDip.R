### tripsAndDip uses read counts for biallelic SNPS to determine if a sample is diploid or triploid
### Input
### counts is either a matrix or a dataframe with each row corresponding to a different sample
### the columns correspond to the read counts for each locus, in a two column per locus format
### So, column 1 is the read counts for locus1Allele1, column two is the read counts for locus1Allele2, locus2Allele1, locus2Allele2, ...
### the rownames of the matrix or dataframe should be the sample names
### h is a list of h values for each locus in the same order that the loci are ordered in counts
### eps is a list of epsilon (error rate per read) values for each locus in the same order that the loci are ordered in counts
### min_reads is the minimum number of reads needed to use a locus in the algorithm
### min_loci is the minimum number of loci needed to attempt to calculate an LLR
### Output
### a dataframe with column 1 containing sample names, column 2 containing calculated LLRs (larger means more likely triploid)
###    and column 3 containing the number of loci used to calculate the LLR
######
### note that constant values for h and epsilon can be easily implemented usign the rep() function in R, for example
###### constant_h <- rep(1, (ncol(allele_counts))/2)
###### constant_eps <- rep(.01, (ncol(allele_counts))/2)
###### results <- tripsAndDip(allele_counts, constant_h, constant_eps)

tripsAndDip <- function(counts, h, eps, min_reads = 30, min_loci = 15, binom_p_value = .05){
	### input error checking
	if(!is.matrix(counts) && !is.data.frame(counts)){
		cat("\nError. Counts must be either a matrix or a dataframe.\n")
		return()
	}
	if((ncol(counts)/2) != length(h)){
		cat("\nError. The number of columns of counts is not equal to twice the length of h.\n")
		return()
	}
	if(length(eps) != length(h)){
		cat("\nError. The length of h is not equal to the length of eps.\n")
		return()
	}
	if(min_reads < 1){
		cat("\nError. min_reads must be 1 or greater.\n")
		return()
	}
	if(min_loci < 1){
		cat("\nError. min_loci must be 1 or greater.\n")
		return()
	}
	if(binom_p_value < 0 || binom_p_value > 1){
		cat("\nError. binom_p_value must be between 0 and 1.\n")
		return()
	}
	
	### calculate llr for each sample
	num_counts <- ncol(counts)
	llr <- apply(counts, 1, function(x){
		#separate allele1 and allele2
		count1 <- x[seq(1, (num_counts - 1), 2)]
		count2 <- x[seq(2, num_counts, 2)]
		#caluculate total read count and which is larger
		n <- count1 + count2
		#determine whether to include loci based on read count
		count1 <- count1[n >= min_reads]
		count2 <- count2[n >= min_reads]
		h_corr <- h[n >= min_reads]
		eps_corr <- eps[n >= min_reads]
		n <- n[n >= min_reads]
		# determine if enough loci
		if (length(n) < min_loci){
			return(NA)
		}
		
		k <- mapply(max, count1, count2)
		#flip h as necessary
		h_corr[count2 > count1] <- 1/h_corr[count2 > count1]
		
		#calculate probabilities
		#based on model with error and allelic bias from Gerard et al. 2018
		p_temp <- (0.6666667*(1 - eps_corr) + (0.3333333)*eps_corr)
		prob_trip <- p_temp / (h_corr*(1 - p_temp) + p_temp)
		p_temp <- 0.5 # (0.5*(1 - eps_corr) + (0.5)*eps_corr) simplifies to 0.5
		prob_dip <- p_temp / (h_corr*(1 - p_temp) + p_temp)
		
		#determine which loci to include using binomial test
		binom_results <- mapply(function(x,y,z){
				return(binom.test(x, y, z, "greater")$p.value)
			}, k, n, prob_trip)
		binom_results <- binom_results > binom_p_value
		
		h_corr <- h_corr[binom_results]
		eps_corr <- eps_corr[binom_results]
		n <- n[binom_results]
		k <- k[binom_results]
		prob_trip <- prob_trip[binom_results]
		prob_dip <- prob_dip[binom_results]
		# determine if enough loci
		if (length(n) < min_loci){
			return(NA)
		}
		
		#calculate log likelihoods
		trip <- (k*log(prob_trip)) + ((n-k)*log(1-prob_trip))
		dip <- (k*log(prob_dip)) + ((n-k)*log(1-prob_dip))
		#sum accros loci
		trip <- sum(trip)
		dip <- sum(dip)
	
		return(c(trip - dip, length(n)))
			
	})
	# associate llr's with sample names and return
	return(data.frame(sample_name = rownames(counts),
			 LLR = llr[1,], loci_used=llr[2,], stringsAsFactors = F, row.names = NULL))
}