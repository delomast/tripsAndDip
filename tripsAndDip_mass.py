#!/usr/bin/env python
##### Tom Delomas 2019
## Call triploid vs diploid based on GT-seq results
## uses a likelihood-based method
##### This version is made to be flexible and work across multiple labs/panels/situations
##### There is another version made for use in EFGL with the panels used there and incorporates reprocessing of triploid samples
##### with adjusted cutoffs for calling heterozygotes
# -h for h value to use across all markers
# -eps for epsilon value to use across all markers
# -marker_info for path to file containing h and epsilon values for all markers - tab-separated file with no header with: markername	h	epsilon
# -pre for the prefix to use for the _ploidy_calls.csv output file, default is 'prefix'
# -min_read minimum number of reads a marker must have to use it to determine ploidy (default 30)
# -min_mark minimum number of markers required to attempt to determine ploidy	(default 15)
# -llr_trip minimum llr to call a sample triploid (default 100)
# -llr_dip maximum llr to call a sample diploid	(default -5)

import glob, re, sys
from math import log
from scipy.stats import binom_test

#function for calling triploid vs diploid using marker specific h and epsilon values
def marker_specific(m_info, files, markers_to_skip, min_read, min_mark, cut_llr_trip, cut_llr_dip):
	#using marker info from file
	#save marker info as two dictionaries with key as marker name and h or epsilon as value
	h_dict = {}
	eps_dict = {}
	with open(m_info, 'r') as info_file:
		line = info_file.readline()
		while line:
			sep = line.rstrip().split('\t')	#remove end of line and split
			h_dict[sep[0]] = float(sep[1])		#save h value
			eps_dict[sep[0]] = float(sep[2])		#save epsilon value
			line = info_file.readline()
	
	calls = []
	#calculate llr for all samples
	for sample in files:
		ll_t = 0	#log-likelihood triploid
		ll_d = 0	#log-likelihood diploid
		loci = 0	#number of loci used to calculate log-likelihoods
		with open(sample, 'r') as genos_file:
			next(genos_file) #skip header of .genos file
			line = genos_file.readline()
			while line:
				#get marker name and counts
				sep = line.split(',')	#marker name is sep[0]
				count1 = int(sep[1].split('=')[1])
				count2 = int(sep[2].split('=')[1])
				higher_count = max(count1, count2)
				total_reads = count1 + count2
				#ensure read counts meet the criteria for inclusion
				if total_reads >= min_read and sep[0] not in markers_to_skip:	#at least min_read reads, marker is not in the skip list
					#get h
					h_value = h_dict[sep[0]]
					#flip h if neccessary
					if count2 > count1:
						h_value = 1 / h_value
					#get epsilon
					eps_value = eps_dict[sep[0]]
					#calculate probability of read for triploid with h and epsilon
					p_binom = (0.6666667*(1 - eps_value)) + ((0.3333333)*eps_value)
					prob_trip = p_binom / ((h_value*(1 - p_binom)) + p_binom)
					#binomial test to determine if skew is greater than expected under triploid - indicates homozygosity or error and therefore not to include
					##no need to calculate binomial test if ratio is less than woudl be expected under triploidy
					if (higher_count / total_reads) <= prob_trip or binom_test(higher_count, total_reads, prob_trip, alternative = 'greater') > 0.05:
						loci += 1
						prob_dip = 0.5 / ((h_value*0.5) + 0.5)
						ll_t += (higher_count*log(prob_trip)) + ((total_reads-higher_count)*log(1-prob_trip))
						ll_d += (higher_count*log(prob_dip)) + ((total_reads-higher_count)*log(1-prob_dip))
				line = genos_file.readline()
		#########
		#### make the call
		llr = ll_t - ll_d
		if loci >= min_mark:		#if a large enough number of loci were useable in the log-likelihood calculation
			if llr >= cut_llr_trip: # triploid
				calls.append([sample, '3n', loci, llr])
			elif llr <= cut_llr_dip: # diploid
				calls.append([sample, '2n', loci, llr])
			else: #unknown
				calls.append([sample, 'U', loci, llr])
		else:	#unknown
			calls.append([sample, 'U', loci, llr])
	return calls
			

#function for calling triploid vs diploid using constant h and epsilon values
def constant(h, eps_value, files, markers_to_skip, min_read, min_mark, cut_llr_trip, cut_llr_dip):
	calls = []
	#calculate llr for all samples
	for sample in files:
		ll_t = 0	#log-likelihood triploid
		ll_d = 0	#log-likelihood diploid
		loci = 0	#number of loci used to calculate log-likelihoods
		with open(sample, 'r') as genos_file:
			next(genos_file) #skip header of .genos file
			line = genos_file.readline()
			while line:
				#get marker name and counts
				sep = line.split(',')	#marker name is sep[0]
				count1 = int(sep[1].split('=')[1])
				count2 = int(sep[2].split('=')[1])
				higher_count = max(count1, count2)
				total_reads = count1 + count2
				#ensure read counts meet the criteria for inclusion
				if total_reads >= min_read and sep[0] not in markers_to_skip:	#at least min_read reads, marker is not in the skip list
					#flip h if neccessary
					if count2 > count1:
						h_value = 1 / h
					else:
						h_value = h
					#calculate probability of read for triploid with h and epsilon
					p_binom = (0.6666667*(1 - eps_value)) + ((0.3333333)*eps_value)
					prob_trip = p_binom / ((h_value*(1 - p_binom)) + p_binom)
					#binomial test to determine if skew is greater than expected under triploid - indicates homozygosity or error and therefore not to include
					##no need to calculate binomial test if ratio is less than woudl be expected under triploidy
					if (higher_count / total_reads) <= prob_trip or binom_test(higher_count, total_reads, prob_trip, alternative = 'greater') > 0.05:
						loci += 1
						prob_dip = 0.5 / ((h_value*0.5) + 0.5)
						ll_t += (higher_count*log(prob_trip)) + ((total_reads-higher_count)*log(1-prob_trip))
						ll_d += (higher_count*log(prob_dip)) + ((total_reads-higher_count)*log(1-prob_dip))
				line = genos_file.readline()
		#########
		#### make the call
		llr = ll_t - ll_d
		if loci >= min_mark:		#if a large enough number of loci were useable in the log-likelihood calculation
			if llr >= cut_llr_trip: # triploid
				calls.append([sample, '3n', loci, llr])
			elif llr <= cut_llr_dip: # diploid
				calls.append([sample, '2n', loci, llr])
			else: #unknown
				calls.append([sample, 'U', loci, llr])
		else:	#unknown
			calls.append([sample, 'U', loci, llr])
	return calls

def Main():

	#list of markers to skip for determining if the sample is triploid or not
	##typically these are presence absence markers, such as a sex marker
	markers_to_skip = ['Ots_SEXY3-1', 'One_1b.75977-57']
	
	#default values for parameters
	panel = 'none'
	h = -9
	eps = -9
	m_info = 'none'
	prefix = 'prefix'
	min_read = 30
	min_mark = 15
	cut_llr_trip = 100
	cut_llr_dip = -5
		
	#read flags
	for flag in range(1, len(sys.argv), 1):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == '-marker_info':
			m_info = sys.argv[flag + 1]
		if sys.argv[flag] == '-h':
			h = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-eps':
			eps = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-t':
			threads = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-pre':
			prefix = sys.argv[flag + 1]
		if sys.argv[flag] == '-min_read':
			min_read = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-min_mark':
			min_mark = float(sys.argv[flag + 1])
		if sys.argv[flag] == '-llr_trip':
			cut_llr_trip = float(sys.argv[flag + 1])		
		if sys.argv[flag] == '-llr_dip':
			cut_llr_dip = float(sys.argv[flag + 1])
	
	if min_read < 1:
		print('\nmin_read cannot be less than 1. using a min_read of 1\n')
		min_read = 1
	if cut_llr_dip > cut_llr_trip:
		print('\nllr_dip cannot be less than llr_trip. Exiting.\n')
		return
	
	#get h and epsilon info, then make calls
	if m_info == 'none':
		if h == -9 or eps == -9:
			print('Did not find -marker_info, or (-h and -eps) specified. At least one of the two must be given. Exiting.')
			return
		else:
			# get .genos files
			files = glob.glob('*.genos')
			print('Calling triploidy/diploidy for ', str(len(files)), 'samples\n')
			#using h and eps as constant
			print('Using constant values of ', str(h), ' and ', str(eps), ' for h and epsilon, respectively\n')
			calls = constant(h, eps, files, markers_to_skip, min_read, min_mark, cut_llr_trip, cut_llr_dip)
	else:
		# get .genos files
		files = glob.glob('*.genos')
		print('Calling triploidy/diploidy for ', str(len(files)), 'samples\n')
		print('Using values of h and epsilon specified for each marker in ', m_info, '\n')
		calls = marker_specific(m_info, files, markers_to_skip, min_read, min_mark, cut_llr_trip, cut_llr_dip)
				
	### write output file for calling triploids
	count_trip = 0
	count_dip = 0
	count_unk = 0
	with open(prefix + '_ploidy_calls.csv', 'w') as outfile:
		outfile.write('Sample,Ploidy,Loci,Llr\n')
		for c in calls:
			if c[1] == '3n': count_trip += 1
			if c[1] == '2n': count_dip += 1
			if c[1] == 'U': count_unk += 1
			outfile.write(re.sub('^initial|^qc|^rr[0-9]|^f[0-9]|\.genos$', '', c[0]) + ',' + c[1] + ',' + str(c[2]) + ',' + str(c[3]) + '\n')
	print('Identified: \n\t', count_trip, ' triploids\n\t', count_dip, ' diploids \n\t', count_unk, ' samples that ploidy could not be determined\n\n')
	
Main()