# tripsAndDip <br/>
differentiate triploid and diploid samples through amplicon sequencing <br/>
<br/><br/>
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
