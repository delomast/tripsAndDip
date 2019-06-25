# tripsAndDip <br/>
differentiate triploid and diploid samples through amplicon sequencing <br/>
<br/><br/>
Instructions for using the python function   <br/>
tripsAndDip_mass.py <br/>
call from directory with .genos files containing GT-seq output for the samples you want to analyze <br/>
call as <br/>
python tripsAndDip_mass.py [options described below] <br/>
<br/><br/>
-h for h value to use across all markers <br/>
-eps for epsilon value to use across all markers <br/>
-marker_info for path to file containing h and epsilon values for all markers - tab-separated file with no header with: markername	h	epsilon <br/>
-pre for the prefix to use for the _ploidy_calls.csv output file, default is 'prefix' <br/>
-min_read minimum number of reads a marker must have to use it to determine ploidy (default 30) <br/>
-min_mark minimum number of markers required to attempt to determine ploidy	(default 15) <br/>
-llr_trip minimum llr to call a sample triploid (default 100) <br/>
-llr_dip maximum llr to call a sample diploid	(default -5) <br/>

 <br/> <br/> <br/> <br/>
 Instructions for using the R function   <br/>
 Note: an R package containing this function is in the repository "tripsAndDipR" <br/>
 It can be installed using: <br/>
      devtools::install_github("delomast/tripsAndDipR") <br/>
 Input <br/>
 counts is either a matrix or a dataframe with each row corresponding to a different sample <br/>
 the columns correspond to the read counts for each locus, in a two column per locus format <br/>
 So, column 1 is the read counts for locus1Allele1, column two is the read counts for locus1Allele2, locus2Allele1, locus2Allele2, ... <br/>
 the rownames of the matrix or dataframe should be the sample names <br/>
 h is a list of h values for each locus in the same order that the loci are ordered in counts <br/>
 eps is a list of epsilon (error rate per read) values for each locus in the same order that the loci are ordered in counts <br/>
 min_reads is the minimum number of reads needed to use a locus in the algorithm <br/>
 min_loci is the minimum number of loci needed to attempt to calculate an LLR <br/> <br/>
 Output <br/>
 a dataframe with column 1 containing sample names, column 2 containing calculated LLRs (larger means more likely triploid) <br/>
    and column 3 containing the number of loci used to calculate the LLR <br/>

 note that constant values for h and epsilon can be easily implemented usign the rep() function in R, for example <br/>
 constant_h <- rep(1, (ncol(allele_counts))/2) <br/>
 constant_eps <- rep(.01, (ncol(allele_counts))/2) <br/>
 results <- tripsAndDip(allele_counts, constant_h, constant_eps) <br/>
