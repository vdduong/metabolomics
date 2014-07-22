### ITERATIVE MATCHING
### includes :
### - a function listing all the BMRB database
### - a function computing the assignments between the theoretical pattern of the 
### metabolite and the experimental peak list
### - a function updating the chemical shifts of the metabolite from the previous assignment
### - a function repeating the previous steps for 7 times, returning the final assignment,
### the corresponding Hausdorff distance, the total and local flow.
### - a function writing out the output into a xml file to be read by TopSpin

from spectrum import *
from class_metabolite import *
from assignment_clustering import *
from final_assignment import *
from network_anchoring import *
import math
import re
import glob
import os
import time

tol_H = 0.1
tol_alignment = 0.04

def database_bmrb() :
	"""
	return the list of all files contained in the folder database
	"""
	path_local = os.getcwd()
	path_final = str(path_local) + '/database/*.str'
	list_file = glob.glob(path_final)
	return list_file

def assignment(metabolite, peak_list):
	"""
	return the dictionary of available assignments
	a dictionary in Python is the equivalent of a hash table
	"""
	global tol_H
	dict_assignment = dict()
	for key_tocsy in metabolite.tocsy_pattern.keys():
		item_tocsy = metabolite.tocsy_pattern[key_tocsy]
		cs_theor_1 = float(item_tocsy[0])
		cs_theor_2 = float(item_tocsy[1])
		dict_assignment[key_tocsy] = list()
		for key_peak_list in peak_list.keys():
			item_peak_list = peak_list[key_peak_list]
			cs_exp_1 = float(item_peak_list[0])
			cs_exp_2 = float(item_peak_list[1])
			if abs(cs_exp_1 - cs_theor_1) <= tol_H and \
				abs(cs_exp_2 - cs_theor_2) <= tol_H :
				proba = math.exp(-0.5*(cs_exp_1 - cs_theor_1)**2/ tol_H**2)*\
						math.exp(-0.5*(cs_exp_2 - cs_theor_2)**2/ tol_H**2)
				proba = float("%.4f"%proba)
				dict_assignment[key_tocsy].append((key_peak_list, proba))

	# sorting the assignment to show the best matching first
	for key_assignment in dict_assignment.keys() :
		item = dict_assignment[key_assignment]
		item.sort(key=lambda x:x[1], reverse= True)
	return dict_assignment

def update_shift(metabolite, peak_list, dict_assignment, best_pattern):
	"""
	update the shift list of the metabolite following the matching patterns
	the best pattern is created from the previous step called best_pattern_generation that
	take inputs from the ranking of different assignments
	we update the shift list based on this best pattern which is a self-consistent pattern
	"""
	global tol_H

	### shift_list is created to be a list of chemical shifts that has no duplicate 
	### within it, since the dictionary allows duplicates.
	shift_list = list()
	for key in metabolite.shift.keys() :
		if ("H" in key) :
			item = metabolite.shift[key]
			shift_list.append(item)
	shift_list = remove_duplicate(shift_list)

	shift_dict = {}
	for item_shift_list in shift_list :
		shift_dict[item_shift_list] = 0.0
	for key_shift_dict in shift_dict.keys() :
		sum_ponderation = 0.0
		# initialization of the sum of the contributions
		time_ponderation = 0.0
		# initialization of the sum of the corresponding probabilities
		for key_ranking in best_pattern.keys() :
			assignment = key_ranking
			peak_tocsy = metabolite.tocsy_pattern[assignment[0]]
			peak_exp = peak_list[assignment[1]]
			if (key_shift_dict in peak_tocsy[0]) or \
					(key_shift_dict in peak_tocsy[1]) :
					rank = best_pattern[key_ranking]
					if key_shift_dict in peak_tocsy[0] :
						sum_ponderation+= float(peak_exp[0])*rank
						# incrementation of the sum of the contributions by the rank of the assignment
						time_ponderation+=rank
					if key_shift_dict in peak_tocsy[1]:
						sum_ponderation+= float(peak_exp[1])*rank
						time_ponderation+=rank

		if time_ponderation > 0.0 :
			sum_ponderation = sum_ponderation/time_ponderation 
			# the sum of contributions is divided by the sum of probabilities
			shift_dict[key_shift_dict]= "%.4f"%(sum_ponderation)
		else :
			shift_dict[key_shift_dict] = "%.4f"%(sum_ponderation)
			# if the sum of probabilities is zero then return only the sum of contributions

	for key_shift in metabolite.shift.keys():
		item_shift = metabolite.shift[key_shift]
		for key_shift_dict in shift_dict :
			if key_shift_dict == item_shift and \
				float(shift_dict[key_shift_dict]) != 0.000 :
				metabolite.shift[key_shift]= shift_dict[key_shift_dict]
	return metabolite

def peaks_found(dict_assignment) :
	"""
	return the number of peaks found 
	"""
	nb_peaks_found = 0
	dict_set = dict()
	for key in dict_assignment.keys() :
		item = dict_assignment[key]
		dict_set[key] = set(item)
		#if len(item) > 0 : nb_peaks_found+=1

	big_set = []
	for key in dict_set.keys():
		if len(dict_set[key]) > 0:
			big_set.append(dict_set[key])
	try:
		intersection_set = set.intersection(*big_set)
	except TypeError as err:
		intersection_set = set()
	for key in dict_set.keys():
		if len(dict_set[key]) > 0:
			dict_set[key].difference_update(intersection_set)
		else: pass

	for key in dict_set.keys():
		if len(dict_set[key]) > 0:
			nb_peaks_found +=1
	nb_peaks_found += len(intersection_set)
	return nb_peaks_found

def iterative_matching(metabolite, peak_list):
	"""iterative matching between an experimental peak list and a renewing database"""
	global tol_alignment, tol_H
	nb_cycle = 0
	final_hausdorff = 100.0 # initialization
	fraction_hausdorff = len((metabolite.tocsy_pattern)) # initialization
	len_initial = len(metabolite.tocsy_pattern)
	while nb_cycle <= 7 :
		dict_assignment = assignment(metabolite, peak_list)
		print dict_assignment
		nb_peaks_found = peaks_found(dict_assignment)
		metabolite.tocsy_pattern = metabolite.tocsy()
		if len(metabolite.tocsy_pattern)==len_initial and len(metabolite.tocsy_pattern) > 1 and \
				float(nb_peaks_found)/float(len(metabolite.tocsy_pattern)) >= 0.7 :
				# more than 70% of peaks are assigned then we continue 
			time_initial=time.time()
			hausdorff, nb_peaks_hausdorff = hausdorffDistance(metabolite, peak_list)
			print "hausdorffDistance",time.time()-time_initial
			old_hausdorff = hausdorff
			time_initial=time.time()
			list_show_ranking = ranking(metabolite, peak_list,dict_assignment)
			print list_show_ranking
			# ranking the assignments
			print "ranking",time.time()-time_initial
			time_initial=time.time()
			sum_ranking, best_pattern = best_pattern_generation(list_show_ranking)
			# issuing the best pattern corresponding to the ranking table
			print "best_pattern_generation",time.time()-time_initial
			time_initial=time.time()
			metabolite = update_shift(metabolite, peak_list, dict_assignment, best_pattern)
			# updating the shift for metabolite object
			print "update_shift",time.time()-time_initial
			hausdorff, nb_peaks_hausdorff = hausdorffDistance(metabolite, peak_list)
			final_hausdorff = hausdorff
			final_nb_peaks_hausdorff = nb_peaks_hausdorff
			if hausdorff == old_hausdorff : # check if the Hausdorff distance is better 
				fraction_hausdorff = float(final_nb_peaks_hausdorff)/float(len(metabolite.tocsy_pattern))
				final_hausdorff = hausdorff
				 # stop the iteration since no more improvement achieved
		else : 
			break # stop the iteration since the number of peaks found less than a threshold of 
					# the total number
		nb_cycle+=1
	return metabolite, dict_assignment, final_hausdorff, fraction_hausdorff, sum_ranking, best_pattern

def write_out(metabolite) :
	"""
	write out the theoretical pattern of the metabolite into a xml file
	that can be read by TopSpin
	but we have to modify the name of the file into peaklist.xml
	"""
	name = metabolite.name
	f = open('/Users/James/Iterameta2.0/Theo_pattern/%s.xml'%name,'w')
	f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
	f.write('<PeakList>\n')
	f.write('  <PeakList2D>\n')
	f.write('    <PeakList2DHeader creator="ejobard@Ozzy" date="2010-07-07T14:49:00" expNo="43" name="IEM050110" \
						owner="cpontoiz" procNo="1" source="C:/Bruker/TOPSPIN">\n')
	f.write('      <PeakPickDetails># Manually picked peaks\n')
	f.write('</PeakPickDetails>\n')
	f.write('    </PeakList2DHeader>\n')
	for key_tocsy in metabolite.tocsy_pattern.keys() :
		cs_1 = float(metabolite.tocsy_pattern[key_tocsy][0])
		cs_2 = float(metabolite.tocsy_pattern[key_tocsy][1])
		f.write('    <Peak2D F1 ="%.3f" F2="%.3f" annotation="Hippuric" intensity="9105828.0" type="1"/>\n'%(cs_1,cs_2))
	f.write('  </PeakList2D>\n')
	f.write('</PeakList>\n')
	f.close()
	return None

##########################
## TEST for the module
if __name__ == '__main__':
	peak_list = {}
	path_local = str(os.getcwd())
	#name_spectrum = path_local + '/Exp_peak_list/Bru_tocsy.peaks'
	#name_spectrum = path_local + '/Exp_peak_list/first_spectrum.txt'
	#f = file(name_spectrum)
	#nb_peaks = 0
	#for line in f :
	#	line = line.rstrip('\n\r')
	#	data = re.split(' +', line)
	#	nb_peaks+=1
	#	peak_list[nb_peaks]=(data[0],data[1])
	#peak_list = first_spectrum()

	peak_list, volume_list = spectrum(path_local + "/Exp_peak_list/43-sn10.peaks")
	#name_file = path_local + "/database/bmse000051.str"
	nb_compounds = 0
	list_file = database_bmrb()
	for name_file in list_file :
	#if 1 :
		dict_assignment = {}
		c = Metabolite(name_file)
		shift = c.initial_shift(name_file)
		c = Metabolite(name_file, shift)
		c.tocsy_pattern = c.tocsy()
		tocsy_initial = c.tocsy_pattern
		#print c.tocsy_pattern
		len_initial = len(c.tocsy_pattern)
		dict_assignment = assignment(c, peak_list)
		nb_peaks_found = peaks_found(dict_assignment)
		if float(len(c.tocsy_pattern)) >1.0 and float(nb_peaks_found)/float(len(c.tocsy_pattern)) >=0.7 :
			c, dict_assignment, final_hausdorff, fraction_hausdorff,sum_ranking,best_pattern = iterative_matching(c, peak_list)
			print c.name,'%i over %i'%(peaks_found(dict_assignment), len_initial)
			#print final_hausdorff, fraction_hausdorff, sum_ranking
			#for key in dict_assignment.keys():
			#	print key, dict_assignment[key]
			print '***'
			for key in c.tocsy_pattern.keys():
				print key, c.tocsy_pattern[key], tocsy_initial[key], '%.4f'%(float(c.tocsy_pattern[key][0]) - float(tocsy_initial[key][0])),\
						'%.4f'%(float(c.tocsy_pattern[key][1])-float(tocsy_initial[key][1]))
			print '***'
			list_show_ranking = ranking(c, peak_list, dict_assignment)
			sum_ranking, best_pattern = best_pattern_generation(list_show_ranking)
			#print sum_ranking
			print [key for key in best_pattern.keys()]
			local_flow(c,peak_list,best_pattern)
			print '_'*20
			showAssignment(c, dict_assignment, peak_list, final_hausdorff, fraction_hausdorff)
			nb_compounds +=1
		else : pass
	print nb_compounds
