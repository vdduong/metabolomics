
from functions import *
from functions_txt import *
from functions_spectra import *
from hausdorff_distance import *
from class_metabolite import *

import math
import copy
import itertools
import thread
import datetime
import matplotlib.pyplot as plt

a = plt.gca()
a.set_xlim([-0.5,10.5])
a.set_ylim([-0.5,10.5])
a.set_ylim(a.get_ylim()[::-1])
a.set_xlim(a.get_xlim()[::-1])
plt.axhline(y=0.0)
plt.axvline(x=0.0)
plt.axhline(y=10.0)
plt.axvline(x=10.0)

t0 = datetime.datetime.now()

tol_H = 0.08 
cut_off_hausdorff = 4.0
tol_merge = 0.05
tol_hausdorff = 0.05

list_point = []
list_point_1 = []

name_file_spectrum = "/Users/James/Metabo_1/Exp_peak_list/43-sn10.peaks"
peak_list = spectrum(name_file_spectrum)
#volume_list = volume_list(name_file_spectrum)
dict_peak_list = dict()
for key in peak_list.keys():
  dict_peak_list[key] = list()


def proba_threshold(tol_H) :
	proba_base = math.exp(-0.5)
	coeff = (0.1/tol_H)**2
	proba_threshold = proba_base**coeff
	return proba_threshold
	
proba_threshold = proba_threshold(tol_H)
print "%.3f"%(proba_threshold)

def remove_duplicate_combinations(item) :
    	io = True
    	list_test = list()
    	for item_item in item :
    		list_test.append(item_item[1])
    	if len(list_test) == len(set(list_test)) :
    		return io
    	else : 
    		io = False
    		return io

def list_all_combinations(dict_assignment) :
	global proba_threshold
	list_all_combinations_new = list()
	list_all_combinations = list()
	list_all_assignments = list()
	for key_assignment in dict_assignment.keys() :
		list_all_assignments_local = list()
		list_assignment = dict_assignment[key_assignment]
		if len(list_assignment) > 0 :
			for item in list_assignment :
				if item[1] >= proba_threshold :
					list_all_assignments_local.append((key_assignment, item[0]))
			list_all_assignments.append(list_all_assignments_local)
	list_all_combinations = list(itertools.product(*list_all_assignments))
	for item_combinations in list_all_combinations :
		io = remove_duplicate_combinations(item_combinations) 
		if io == True :
			list_all_combinations_new.append(item_combinations)
	return list_all_combinations_new

def test_bilateral_alignment(tocsy_pattern, peak_list, tol_merge, item_list_all_combinations) :
	i = 0
	io_sum = 0
	while i < len(item_list_all_combinations) :
		item_i = item_list_all_combinations[i] 
		peak_theo_tocsy_i = item_i[0]
		cs_theo_1_i = float(tocsy_pattern[peak_theo_tocsy_i][0])
		cs_theo_2_i = float(tocsy_pattern[peak_theo_tocsy_i][1])
		peak_exp_i = item_i[1]
		cs_exp_1_i = float(peak_list[peak_exp_i][0])
		cs_exp_2_i = float(peak_list[peak_exp_i][1])
		
		j = i+1
		while j < len(item_list_all_combinations) :
			item_j = item_list_all_combinations[j]
			peak_theo_tocsy_j = item_j[0]
			cs_theo_1_j = float(tocsy_pattern[peak_theo_tocsy_j][0])
			cs_theo_2_j = float(tocsy_pattern[peak_theo_tocsy_j][1])
			peak_exp_j = item_j[1]
			cs_exp_1_j = float(peak_list[peak_exp_j][0])
			cs_exp_2_j = float(peak_list[peak_exp_j][1])
			
			io_1 = (abs(cs_theo_1_i - cs_theo_1_j) <= tol_merge) 
			io_2 = (abs(cs_theo_2_i - cs_theo_2_j) <= tol_merge)
			io_3 = (abs(cs_exp_1_i - cs_exp_1_j) <= tol_merge) 
			io_4 = (abs(cs_exp_2_i - cs_exp_2_j) <= tol_merge)
			
			
			#if (io_1 == True and io_3 == False) or \
			#	(io_2 == True and io_4 == False) : 
			#		io = False
			#else : io = True
			if (io_1 == True and io_3 == True) or \
				(io_2 == True and io_4 == True) :
					io_sum+=1
								
			j+=1
		i+=1
	return io_sum
	
    
def best_count(list_all_combinations, tocsy_pattern) :
	"""return the combinations giving the best horizontal and vertical alignment"""
	global peak_list, tol_merge
	list_best_combinations = list()
	io_max = 0
	for item in list_all_combinations :
		io_sum = test_bilateral_alignment(tocsy_pattern, peak_list, tol_merge, item)
		if io_sum > io_max :
			io_max = io_sum
		else :
			pass
	for item in list_all_combinations :
		io_sum = test_bilateral_alignment(tocsy_pattern, peak_list, tol_merge, item)
		if io_sum == io_max :
			list_best_combinations.append(item)
	#print io_max
	return list_best_combinations

	
#---------------Main function 


name_file_1 = '/Users/James/Metabo_1/database/bmse000408.str'
name_file_2 = '/Users/James/Metabo_1/database/bmse000028.str'
name_file_3 = '/Users/James/Metabo_1/database/bmse000120.str'
name_file_4 = '/Users/James/Metabo_1/database/bmse000052.str'
metabolite = name_metabolite(name_file_1)
m = matrix_adjacent(name_file_1)
shift_table = shift(name_file_1)
pattern_tocsy = pattern(shift_table, m) 





dict_assignment = assignment(pattern_tocsy, peak_list, tol_H)


def update_shift_table(shift_table,pattern_tocsy,peak_list, dict_best_combinations) :
	shift_list = list()
	for key in shift_table.keys() :
		if ("H" in key) == True :
			item = shift_table[key]
			shift_list.append(item)
	shift_list = remove_duplicate(shift_list)
	shift_dict = dict()
	for item_shift_list in shift_list :
		shift_dict[item_shift_list] = 0.0
	for key_shift_dict in shift_dict.keys() :
		sum_ponderation = 0.0
		time_ponderation = 0.0
		for item in dict_best_combinations.keys() : # attention : item means key for dict here
			for item_item in item :
				peak_1 = pattern_tocsy[item_item[0]]
				peak_2 = peak_list[item_item[1]]
				if (key_shift_dict in peak_1[0]) or \
					(key_shift_dict in peak_1[1]) :
					proba = proba_matching(peak_1, peak_2, tol_H)
					if proba >= 0.7 :
						if key_shift_dict in peak_1[0] :
							sum_ponderation += float(peak_2[0])*proba 
							time_ponderation +=proba 
						if key_shift_dict in peak_1[1] :
							sum_ponderation += float(peak_2[1])*proba 
							time_ponderation +=proba 
	
		if time_ponderation > 0.0 :
			sum_ponderation = sum_ponderation/time_ponderation
			shift_dict[key_shift_dict] = "%.4f"%(sum_ponderation) 
		else : 
			shift_dict[key_shift_dict] = "%.4f"%(sum_ponderation)
		
	for key_shift_table in shift_table.keys() :
		item_shift_table = shift_table[key_shift_table]
		for key_shift_dict in shift_dict :
			if key_shift_dict == item_shift_table :
				shift_table[key_shift_table] = shift_dict[key_shift_dict]
	return shift_table


def print_table(key_shift, dict_ref, peak_list, dict_assignment) :
	list_peak_local = dict_ref[key_shift]
	for item_peak_local in list_peak_local :
		key_pattern_tocsy = item_peak_local[0]
		position_peak = item_peak_local[1]
		for key_assignment in dict_assignment.keys() :
			if key_pattern_tocsy == key_assignment :
				list_assignment_local = dict_assignment[key_assignment]
				for item_assignment in list_assignment_local :
					proba = item_assignment[1]
					peak_assigned = item_assignment[0]
					cs_printed = peak_list[peak_assigned][position_peak]
					print key_shift, ": ", key_pattern_tocsy, ": ", \
							peak_assigned,":",cs_printed, ": ",pattern_tocsy[key_pattern_tocsy],": ",peak_list[peak_assigned],": ", proba
					


def iterative_matching(name_file, peak_list, tol_H, tol_merge, dtn) :
	list_point = []
	list_point_2 = []
	metabolite = name_metabolite(name_file)
	m = matrix_adjacent(name_file)
	shift_table = shift(name_file)
	pattern_tocsy = pattern(shift_table, m)
	
	len_initial = len(pattern_tocsy)
	cycle_nb = 0
	dict_assignment = assignment(pattern_tocsy, peak_list, tol_H)
	while cycle_nb <=20 :
		#print "CYCLE %i"%(cycle_nb)
		pattern_tocsy = pattern(shift_table,m)
		dict_assignment = assignment(pattern_tocsy, peak_list, tol_H)
		#for key in dict_assignment.keys() :
		#	print key, dict_assignment[key]
		list_all_combinations_new = list_all_combinations(dict_assignment)
		dict_best_combinations = dict()
		list_best_combinations = best_count(list_all_combinations_new, pattern_tocsy)
		for item in list_best_combinations :
			#print item
			distance = hausdorff_distance(pattern_tocsy, item, peak_list, tol_hausdorff)
			dict_best_combinations[item] = distance
		list_20_best,distance_min = twenty_best_matches(list_best_combinations, tol_hausdorff, pattern_tocsy, peak_list) 
		shift_table = update_shift_table(shift_table,pattern_tocsy,peak_list,dict_best_combinations)
		
		cycle_nb+=1
	print shift_table
		#print "************"
	if len_initial > 0.0 and float(len(list_20_best[0]))/float(len_initial) == 1.0 :
		print "_________________"
		print metabolite, "found : %i on %i"%(len(list_20_best[0]),len_initial)
		for key in pattern_tocsy.keys() :
			io = False
			for item in list_20_best[0]:
				if key == item[0] :
					io = True
					#print item[0],":", pattern_tocsy[item[0]], item[1], ":", peak_list[item[1]]
			if io == False :
				pass
				#print pattern_tocsy[key]
		#print list_20_best[0]
		
		print pattern_tocsy
		for item in list_20_best[0] :
			cs_exp_1 = peak_list[item[1]][0]
			cs_exp_2 = peak_list[item[1]][1]
			list_point.append((float(cs_exp_1), float(cs_exp_2)))
		for key in pattern_tocsy.keys() :
			item_tocsy = pattern_tocsy[key]
			cs_theo_1 = item_tocsy[0]
			cs_theo_2 = item_tocsy[1]
			list_point_2.append((float(cs_theo_1), float(cs_theo_2)))			
		plt.plot(*zip(*list_point_2), marker='o', color='y', ls='')
		#plt.plot(*zip(*list_point), marker= 'o', color='r', ls = '')
		#for point in list_point :
		#	plt.annotate(metabolite, xy = point, xytext = (-20, 20), textcoords = 'offset points', \
		#	ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.05', fc = 'yellow', alpha = 0.5),\
		#	arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

		
			
		#####dtn.write('<Peak2D F1 ="%s" F2 ="%s" annotation="%s" intensity = "%.3f" type="1"/>\n'%(cs_theo_1,cs_theo_2,metabolite, volume_list[item[1]]))
	else : pass		
		
	
print "#################"	
#---Execution

#list_file = data_file()
dtn = open("/Users/James/Metabo_1/tocsy_peaks/Final_assignment.xml", "w")
#dtn.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#dtn.write("<PeakList>\n")
#dtn.write("<PeakList2D>\n")
#dtn.write('    <PeakList2DHeader creator="ejobard@Ozzy" date="2010-07-07T14:49:00" expNo="43" name="IEM050110" owner="cpontoiz" procNo="1" source="C:/Bruker/TOPSPIN">\n')
#dtn.write("<PeakPickDetails># Manually picked peaks\n")
#dtn.write("</PeakPickDetails>\n")
#dtn.write("</PeakList2DHeader>\n")

#list_peak = []
#for key in peak_list.keys() :
#	list_peak.append((float(peak_list[key][0]), float(peak_list[key][1])))
#	plt.plot(*zip(*list_peak), marker = 'h', color='m', ls= '')
def main_metabolite() :
	list_file = [name_file_1, name_file_2,name_file_3]
	for name_file_new in list_file :
		list_point_1 = []
		iterative_matching(name_file_new, peak_list, tol_H, tol_merge,dtn)
		metabolite = name_metabolite(name_file_new)
		m = matrix_adjacent(name_file_new)
		shift_table = shift(name_file_new)
		pattern_tocsy = pattern(shift_table, m) 
		for key in pattern_tocsy.keys() :
			cs_theo_1 = pattern_tocsy[key][0]
			cs_theo_2 = pattern_tocsy[key][1]
			list_point_1.append((float(cs_theo_1), float(cs_theo_2)))
		plt.plot(*zip(*list_point_1), marker='o', color='b', ls='')
#iterative_matching(name_file_1, peak_list, tol_H, tol_merge, dtn)
#plt.plot(*zip(*list_point), marker='o', color='r', ls='')


main_metabolite()
plt.savefig("/Users/James/Metabo_1/spectrum_hippurate.png")

#print nb_components, " are found in the mixture"

dtn.write("</PeakList2D>\n")
dtn.write("</PeakList>\n")
dtn.close()

delta_t = datetime.datetime.now() - t0
print delta_t




  
    


