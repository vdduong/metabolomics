### FINAL ASSIGNMENT
### includes : 
### - a function computing Hausdorff distance between the theoretical pattern of 
### metabolite and the experimental peak list
### - a function printing the output assignment

from spectrum import *
from class_metabolite import *
from iterative_matching import *
from assignment_clustering import *
import os
import math
import itertools

# in the final assignment, we have a TOCSY pattern that is close to the experimental 
# pattern. The dictionary of final assignments is then easily obtained.
# The attribution of assignments will be completed by a matching of type 'stable marriage'
# The patterns are then computed to show Hausdorff distances to the experimental peak list
dict_assignment = {}
tol_H = 0.1
tol_alignment = 0.03
dict_ranking = {}


def hausdorffDistance(metabolite,peak_list):
	"""
	compute the Hausdorff distance between the TOCSY pattern of the metabolite
	against the experimental peak list
	return the fractional Hausdorff distance value and the fraction that we have to take
	out from the theoretical pattern in order to achieve this value
	"""
	global tol_H
	list_distance_h = []
	tocsy_pattern = metabolite.tocsy_pattern
	for key_tocsy in tocsy_pattern.keys():
		peak_tocsy = tocsy_pattern[key_tocsy]
		cs_theo_1 = float(peak_tocsy[0])
		cs_theo_2 = float(peak_tocsy[1])
		min_local = 10000.0 # default value for Hausdorff distance
		for key_exp in peak_list.keys():
			peak_exp = peak_list[key_exp]
			cs_exp_1 = float(peak_exp[0])
			cs_exp_2 = float(peak_exp[1])
			if abs(cs_theo_1 - cs_exp_1) <= tol_H and \
				abs(cs_theo_2 - cs_exp_2) <= tol_H :
					distance_h = (cs_theo_1-cs_exp_1)**2+(cs_theo_2-cs_exp_2)**2
					if distance_h <= min_local :
						min_local = distance_h
		list_distance_h.append(min_local)
	list_distance_h.sort(reverse=True) # sorting the distance set
	i = 0
	nb_peaks_hausdorff = 0
	while i <= len(list_distance_h)-1 :
		if list_distance_h[i] >= 100.0 :
			nb_peaks_hausdorff+=1
		i+=1

	return math.sqrt(list_distance_h[nb_peaks_hausdorff]), nb_peaks_hausdorff


def showAssignment(metabolite, dict_assignment, peak_list, hausdorff_distance, fraction, tocsy_initial):
	""" 
	show assignments whose matching probabilities are more than a threshold
	"""
	global tol_H
	def probability_threshold(h_distance):
		"""
		transform the upper-limit Hausdorff distance into a threshold probability
		"""
		return math.exp(-0.5*h_distance**2/tol_H**2)

	nb_peaks_fin = 0

	list_showing = []
	for key_assignment in dict_assignment.keys() :
		item_assignment = dict_assignment[key_assignment]
		peak_presence = False
		if len(item_assignment) > 0 :
			for assignment in item_assignment :
				peak_nb = assignment[0]
				peak = peak_list[peak_nb]
				if assignment[1] >= probability_threshold(hausdorff_distance):
					peak_presence = True
					list_showing.append((peak_nb, peak, key_assignment, \
						metabolite.tocsy_pattern[key_assignment], '%.4f'%(assignment[1])))
		if peak_presence == True : nb_peaks_fin+=1

	if nb_peaks_fin >= 0.7*len(tocsy_initial) :
		print metabolite.name, nb_peaks_fin,'over ', len(tocsy_initial)
		#print 'probability accepted : %.4f'%probability_threshold(hausdorff_distance)
		for item_showing in list_showing :
			print item_showing[0], item_showing[1],': ',item_showing[2], item_showing[3], \
				item_showing[4]
	else :
		pass
	return None

##############################
### TEST of the module

if __name__ == '__main__':
	path_local = str(os.getcwd())
	peak_list, volume_list = spectrum(path_local + "/Exp_peak_list/43-sn10.peaks")
	name_file = path_local + "/database/bmse000408.str"
	c = Metabolite(name_file)
	shift = c.initial_shift(name_file)
	c = Metabolite(name_file, shift)
	c.tocsy_pattern = c.tocsy()
	dict_assignment = assignment(c, peak_list)
	for key in dict_assignment :
		print key, dict_assignment[key]
	c, dict_assignment, final_hausdorff, fraction_hausdorff = iterative_matching(c, peak_list)
	print final_hausdorff, fraction_hausdorff
	#print peaks_found(dict_assignment)
	showAssignment(c, dict_assignment, peak_list, final_hausdorff, fraction_hausdorff)
	print peak_list[73]
