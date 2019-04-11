# tripsAndDip
differentiate triploid and diploid samples through amplicon sequencing

tripsAndDip_mass.py
call from directory with .genos files containing GT-seq output for the samples you want to analyze
call as:
python tripsAndDip_mass.py [options described below]

-h for h value to use across all markers
-eps for epsilon value to use across all markers
-marker_info for path to file containing h and epsilon values for all markers - tab-separated file with no header with: markername	h	epsilon
-pre for the prefix to use for the _ploidy_calls.csv output file, default is 'prefix'
-min_read minimum number of reads a marker must have to use it to determine ploidy (default 30)
-min_mark minimum number of markers required to attempt to determine ploidy	(default 15)
-llr_trip minimum llr to call a sample triploid (default 100)
-llr_dip maximum llr to call a sample diploid	(default -5)
