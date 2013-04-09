# functions : regrouping the whole set of functions used in the program
# graphical interface is not included

import re
import math
import glob
import itertools

def data_file() :
  """ load the database files into a list of file name"""
	list_file = glob.glob("/Users/James/Metabo_1/database/*.str")
	return list_file

def spectrum(name_file) :
    """ load the experimental peak list into a dictionary """
    spectrum = file(name_file)
    list = dict()
    nb_line = 13
    nb_count = 0
    while nb_count < nb_line :
          spectrum.readline()
          nb_count+=1
    for line_spectrum in spectrum :
        line_spectrum = line_spectrum.rstrip('\n\r')
        data_spectrum = re.split(" +", line_spectrum)
        peak_nb = int(data_spectrum[1])
        cs_1 = data_spectrum[2]
        cs_2 = data_spectrum[3]
        list[peak_nb] = (cs_1, cs_2)
    spectrum.close()      
    return list
    
def volume_list(name_file) :
    """ load the experimental peak list volume into a dictionary """ 
    spectrum = file(name_file)
    list = dict()
    nb_line = 13
    nb_count = 0
    while nb_count < nb_line :
          spectrum.readline()
          nb_count+=1
    for line_spectrum in spectrum :
        line_spectrum = line_spectrum.rstrip('\n\r')
        data_spectrum = re.split(" +", line_spectrum)
        peak_nb = int(data_spectrum[1])
        signal = data_spectrum[6]
        noise = data_spectrum[7]
        list[peak_nb] = float(signal)/ float(noise)
    spectrum.close()      
    return list
    
### some basic functions for type 'list':

def remove_duplicate(mylist) :
	""" remove duplicates from a list"""
	if mylist:
		mylist.sort()
		last = mylist[-1]
		for i in range(len(mylist)-2, -1, -1):
			if last == mylist[i]:
				del mylist[i]
			else:
				last = mylist[i]
	return mylist
    
def check_duplicate(mylist):
    """return True if duplicates found, False otherwise """
    return len(mylist)!=len(set(mylist))
    
### functions working on a particular database file

def name_metabolite(name_file) :
	""" return the name of the metabolite contained in the database file """
	data_file = file(name_file)
	found = False
	for line in data_file :
		if "_Entry.Title" in line :
			found = True
			line = line.rstrip('\n\r')
			data_line = line.split(" ")
			name_metabo = list()
			for item_data_line in data_line :
				if len(str(item_data_line)) > 0 and item_data_line != "_Entry.Title" :
					name_metabo.append(str(item_data_line))
	name_metabolite = ''
	for item in name_metabo :
		name_metabolite = name_metabolite + item
	return name_metabolite    

def matrix_adjacent(name_file) :
    """ return the adjacent matrix loaded from the database file
    	in the form of a 'list' """
    m_adjac = list()
    data_file = file(name_file)
    line_nb = 0
    found = False
    for line in data_file :
        line_nb +=1
        if "_Chem_comp_bond.Comp_ID" in line :
           #print line_nb
           break
    data_file.close()
    data_file = file(name_file)
    line_count = 0
    while line_count <= line_nb :
          line_count+=1
          data_file.readline()
    
    for line in data_file :
        line_nb+=1
        io = "stop_" in line
        if io == False : 
           line = line.rstrip("\n\r")
           data_line = re.split(" +", line)
           m_adjac.append((str(data_line[4]),str(data_line[5])))
        
        else : break
    data_file.close()
    return m_adjac
    
def shift(name_file) :
    """ load the chemical shifts of the metabolite into a dictionary of 'list'
    	since one atom can have different chemical shifts (bmrb file) """
    shift_list = dict()
    data_file = file(name_file)
    line_nb = 0
    found = False
    for line in data_file :
        line_nb+=1
        if "_Atom_chem_shift.Assigned_chem_shift_list_ID" in line :
           #print line_nb
           break
    data_file = file(name_file) 
    line_count = 0
    list_shift = list()
    while line_count <= line_nb :
          line_count+=1
          data_file.readline()
    for line in data_file :
        line_nb+=1
        io = "stop_" in line 
        if io == False :
          try :
            line = line.rstrip("\n\r")
            data_line = re.split(" +", line)
            key = str(data_line[5])
            if key in list_shift :
                item = str(data_line[8])
            	shift_list[key].append(item)
            else :
            	shift_list[key] = list()
            	item = str(data_line[8])
            	shift_list[key].append(item)
          except IndexError as err :
            break
        else : break
           
    return shift_list
    
def shortest_path(top_a, top_b, m) : 
    """ return the bonding distance that separates two atoms """
    found = False
    distance = 0
    
    list_initial = list()
    list_initial.append(top_a)
    
    list_1 = list()
    for item in m :
        if item[0] == top_a :
           list_1.append(item[1])
        elif item[1] == top_a :
           list_1.append(item[0])
        else : pass
    if top_b in list_1 :
       distance+=1
       found = True
    else :
       distance+=1
    
    while (found == False) and (distance <= 4):
          list_2 = list()
          for item_1 in list_1 :
              for item in m :
                  if item[0] == item_1 :
                     io = (item[0] in list_initial)
                     if io == False  : list_2.append(item[1])
                     else : pass
                  elif item[1] == item_1 :
                     io = (item[1] in list_initial) 
                     if io == False  : list_2.append(item[0])
                     else : pass
                  else : pass
          
          if top_b in list_2 :
             distance+=1
             found = True
          else :   
             list_initial = list_1
             list_1 = list_2
             distance+=1
    return distance

def pattern(shift_list, m_adjac) :
    """ load the theoretical Tocsy pattern of the metabolite into a dictionary"""

    pattern_tocsy = dict()
    peak_nb = 0
    pattern_tocsy_list = list()
    shift_list_proton = list()
    for key in shift_list :
        if ('H' in key) == True : # we take only the protons
           shift_list_proton.append(key)
       
    pattern_peak = list()
    index_proton_1 = 0
    while index_proton_1 < len(shift_list_proton) :
          top_a = shift_list_proton[index_proton_1]
          pattern_peak.append((top_a, top_a))
          index_proton_2 = index_proton_1 + 1
          while index_proton_2 < len(shift_list_proton) :
                top_b = shift_list_proton[index_proton_2]
                distance = shortest_path(top_a, top_b, m_adjac)
                if distance <=4 : 
                   pattern_peak.append((top_a, top_b))
                   pattern_peak.append((top_b, top_a)) # add symmetrical crosspeak
                index_proton_2+=1
          index_proton_1+=1
       
    for item_peak in pattern_peak :
        cs_a = 0.0
        cs_b = 0.0
        top_a = item_peak[0]
        top_b = item_peak[1]
        for key in shift_list :
            if key == top_a : list_cs_a = shift_list[key]
        for key in shift_list :
            if key == top_b : list_cs_b = shift_list[key]
        for item_a in list_cs_a :
        	cs_a = item_a
        	for item_b in list_cs_b :
        		cs_b = item_b
        		cross_peak = (cs_a, cs_b)
        pattern_tocsy_list.append(cross_peak)
    
    pattern_tocsy_list = remove_duplicate(pattern_tocsy_list)
    for item in pattern_tocsy_list :
        peak_nb+=1
        pattern_tocsy[peak_nb] = item
    return pattern_tocsy
    
### matching functions

def proba_2_peaks(peak_1, peak_2) :
	""" return the probability of matching between two peaks """
	cs_1_1 = float(peak_1[0])
	cs_1_2 = float(peak_1[1])
	cs_2_1 = float(peak_2[0])
	cs_2_2 = float(peak_2[1])
	proba_2_peaks = math.exp(-0.5*(cs_1_1 - cs_2_1)**2/0.05**2)*\
					math.exp(-0.5*(cs_1_2 - cs_2_2)**2/0.05**2)
	return proba_2_peaks
	
def proba_matching(peak_1, peak_2, tol_H):
	"""define the probability of matching
		peak_1, peak_2 in form of (cs_a, cs_b)
		given tolerance of matching"""
	proba_matching = 0.0
	cs_1_a = float(peak_1[0])
	cs_1_b = float(peak_1[1])
	cs_2_a = float(peak_2[0])
	cs_2_b = float(peak_2[1])
	proba_matching = math.exp(- 0.5*(cs_1_a -cs_2_a)**2/tol_H**2)*math.exp(- 0.5*(cs_1_b - cs_2_b)**2/tol_H**2)
	return proba_matching

def distance(peak_1, peak_2, tol_hausdorff):
	distance = (float(peak_1[0]) - float(peak_2[0]))**2 + (float(peak_1[1]) - float(peak_2[1]))**2
	distance = distance/tol_hausdorff**2
	return math.sqrt(distance)

def hausdorff_distance(pattern_tocsy, pattern_assigned, peak_list, tol_hausdorff):
	"""compute the hausdorff distance between the two patterns
		experimental and theoretical"""
	hausdorff_distance = 0.0
	for item_assigned in pattern_assigned :
		nb_peak_1 = item_assigned[0]
		nb_peak_2 = item_assigned[1]
		peak_1 = pattern_tocsy[nb_peak_1]
		peak_2 = peak_list[nb_peak_2]
		distance_local = distance(peak_1, peak_2, tol_hausdorff)
		if distance_local > hausdorff_distance :
			hausdorff_distance = distance_local
		else : pass
	return hausdorff_distance
	
def twenty_best_matches(list_best_combinations, tol_hausdorff, pattern_tocsy, peak_list) :
	"""compute 20 best_matches in the list of best combinations"""
	list_20_best = list()
	list_20_best_local = list()
	i = 0
	while i < len(list_best_combinations) :
		pattern_assigned = list_best_combinations[i]
		distance = hausdorff_distance(pattern_tocsy, pattern_assigned, peak_list, tol_hausdorff)
		element = (i, distance)
		list_20_best_local.append(element)
		i+=1
	list_20_best_local.sort(key=lambda element: element[1])
	if len(list_20_best_local) > 0 :
		distance_min = list_20_best_local[0][1]
	else :
		distance_min = 1000.0
	j = 0
	if len(list_best_combinations) >= 20 :
		while j < 20 :
			pattern_assigned_nb = list_20_best_local[j][0]
			item = list_best_combinations[pattern_assigned_nb]
			list_20_best.append(item)
			j+=1
	else :
		while j < len(list_best_combinations) :
			pattern_assigned_nb = list_20_best_local[j][0]
			item = list_best_combinations[pattern_assigned_nb]
			list_20_best.append(item)
			j+=1
	return list_20_best, distance_min
    
def assignment(pattern_tocsy, peak_list, tol_H) :
    """ return a dictionary of assignment available for every metabolite """
    #tol_H = 0.05
    dict_assignment = dict()
    
    for key_pattern_tocsy in pattern_tocsy.keys() :
        item_pattern_tocsy = pattern_tocsy[key_pattern_tocsy]
        cs_theor_1 = float(item_pattern_tocsy[0])
        cs_theor_2 = float(item_pattern_tocsy[1])
        
        cs_str_1 = '%.3f' %(cs_theor_1)
        cs_str_2 = '%.3f' %(cs_theor_2)
        
        dict_assignment[key_pattern_tocsy] = list()
        for key_peak_list in peak_list.keys() :
            item_peak_list = peak_list[key_peak_list]
            cs_exp_1 = float(item_peak_list[0])
            cs_exp_2 = float(item_peak_list[1])
            dict_assignment_proba = dict()
            if abs(cs_exp_1 - cs_theor_1) <= tol_H and \
               abs(cs_exp_2 - cs_theor_2) <= tol_H :
               
               proba = math.exp(-0.5*(cs_exp_1 - cs_theor_1)**2/ tol_H**2) * \
                       math.exp(-0.5*(cs_exp_2 - cs_theor_2)**2/ tol_H**2)
                       
               proba = float("%.3f"%(proba))
               dict_assignment[key_pattern_tocsy].append((key_peak_list, proba))
        
    return dict_assignment

	
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
	""" load all the non-duplicated combinations of assigned peaks into a list """
	list_all_combinations_new = list()
	list_all_combinations = list()
	list_all_assignments = list()
	for key_assignment in dict_assignment.keys() :
		list_all_assignments_local = list()
		list_assignment = dict_assignment[key_assignment]
		if len(list_assignment) > 0 :
			for item in list_assignment :
				if item[1] > 0.7 :	
					list_all_assignments_local.append((key_assignment, item[0]))
			list_all_assignments.append(list_all_assignments_local)
	list_all_combinations = list(itertools.product(*list_all_assignments))
	
	for item_combinations in list_all_combinations :
		io = remove_duplicate_combinations(item_combinations) 
		if io == True :
			list_all_combinations_new.append(item_combinations)
	return list_all_combinations_new


def test_bilateral_alignment(tocsy_pattern, peak_list, tol_merge, item_list_all_combinations) :
	""" return the best score of alignment within a combination of assigned peaks 
		we suppose that the perfectly assigned pattern must have the best score
		of alignment """
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

			if (io_1 == True and io_3 == True) or \
				(io_2 == True and io_4 == True) :
					io_sum+=1
								
			j+=1
		i+=1
	return io_sum
	
def best_count(list_all_combinations, tocsy_pattern, peak_list, tol_merge) :
	""" return the combinations giving the best horizontal and vertical alignment 
		in the form of a 'list' """
	#global peak_list, tol_merge
	#global tol_merge
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
	return list_best_combinations
	
	
def update_shift_table(shift_table,pattern_tocsy,peak_list, dict_best_combinations, tol_H) :
	""" compute the shift of the theoretical pattern in comparison with the experimental 
		one
		and return the new shift table for the metabolite """

	shift_list = list()
	for key in shift_table.keys() :
		if ("H" in key) == True :
			item = shift_table[key]
			for item_item in item :
				shift_list.append(item_item)
	shift_list = remove_duplicate(shift_list)
	shift_dict = dict()
	for item_shift_list in shift_list :
		shift_dict[item_shift_list] = []
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
					if key_shift_dict in peak_1[0] :
						sum_ponderation += float(peak_2[0])*proba 
						time_ponderation +=proba 
					if key_shift_dict in peak_1[1] :
						sum_ponderation += float(peak_2[1])*proba 
						time_ponderation +=proba 
		if time_ponderation > 0.0 :
			sum_ponderation = sum_ponderation/time_ponderation
			shift_dict[key_shift_dict].append("%.4f"%(sum_ponderation)) 
		else : 
			shift_dict[key_shift_dict].append("%.4f"%(sum_ponderation))
		
	for key_shift_table in shift_table.keys() :
		item_shift_table = shift_table[key_shift_table]
		for key_shift_dict in shift_dict :
			for item_item_shift_table in item_shift_table :
				if key_shift_dict == item_item_shift_table :
					item_item_shift_table = shift_dict[key_shift_dict]
	return shift_table


def iterative_matching(name_file, peak_list, tol_H, tol_merge, nb_components, list_point, list_point_1) :    
																	#, dtn, dtn_theo) :# depends on whenever people 
																	#	want to load into bruker topspin
																	
	""" matching the experimental peak list against the continuingly shifted theoretical 
	pattern in order to find the best matching """
	
	tol_hausdorff = tol_merge
	metabolite = name_metabolite(name_file)
	m = matrix_adjacent(name_file)
	shift_table = shift(name_file)
	pattern_tocsy = pattern(shift_table, m)
	len_initial = len(pattern_tocsy)
	cycle_nb = 0
	while cycle_nb <=20 :
		pattern_tocsy = pattern(shift_table,m)
		dict_assignment = assignment(pattern_tocsy, peak_list, tol_H)
		
		nb_found = 0
		
		for key in dict_assignment.keys() :
			item = dict_assignment[key]
			if len(item) > 0 : nb_found+=1
		if len_initial > 0 and float(nb_found)/float(len_initial) < 0.8 :
			break
		else : pass 
		
		list_all_combinations_new = list_all_combinations(dict_assignment)
		dict_best_combinations = dict()
		list_best_combinations = best_count(list_all_combinations_new, pattern_tocsy, peak_list, tol_merge)
		for item in list_best_combinations :
			distance = hausdorff_distance(pattern_tocsy, item, peak_list, tol_hausdorff)
			dict_best_combinations[item] = distance
			#list_20_best,distance_min = twenty_best_matches(list_best_combinations, tol_hausdorff, pattern_tocsy, peak_list) 
			shift_table = update_shift_table(shift_table,pattern_tocsy,peak_list,dict_best_combinations, tol_H)
		cycle_nb+=1
	print shift_table
	nb_found_final = 0
	for key in dict_assignment.keys() :
		item = dict_assignment[key]
		if len(item) > 0 : nb_found_final+=1
		
	if len_initial >= 1.0 and float(nb_found_final)/float(len_initial) == 1.0 :
		score_max, score_min = 0.0, 1.0
		# condition of approving a metabolite
	
		list_20_best,distance_min = twenty_best_matches(list_best_combinations, tol_hausdorff, pattern_tocsy, peak_list) 
		if len(list_20_best) > 0 :
			print metabolite, "found : %i on %i"%(len(list_20_best[0]), len(pattern_tocsy))
			if len(list_20_best[0]) > 0 :
				for item in list_20_best[0]:
					peak_1 = pattern_tocsy[item[0]]
					peak_2 = peak_list[item[1]]
					proba = proba_2_peaks(peak_1, peak_2)
					print pattern_tocsy[item[0]], peak_list[item[1]], item[1], ':', proba
					if score_max < proba : score_max = proba
					if score_min > proba : score_min = proba
				
				nb_components+=1
				print nb_components
					#peaks_nb+=1
					#dtn.write("	%i	%s	%s	0 U          0.000e+00  0.00e+00 -   0    0    0 0 # %s\n"%(peaks_nb, cs_exp_1, cs_exp_2, metabolite))			
				#dtn.write("</PeakList2D>\n")
				#dtn.write("</PeakList>\n")				
				
				m = matrix_adjacent(name_file)
				shift_table = shift(name_file)
				pattern_tocsy = pattern(shift_table, m)
				for key in pattern_tocsy.keys() :
					cs_theo_1 = pattern_tocsy[key][0]
					cs_theo_2 = pattern_tocsy[key][1]
					list_point_1.append((float(cs_theo_1), float(cs_theo_2)))
					
					#dtn_theo.write('<Peak2D F1 ="%s" F2 ="%s" annotation="%s" intensity = "9105828.0" type="1"/>\n'%(cs_theo_1,cs_theo_2,metabolite))
				
				for item in list_20_best[0] :
					peak_2 = peak_list[item[1]]
					cs_exp_1 = peak_2[0]
					cs_exp_2 = peak_2[1]
					list_point.append((float(cs_exp_1), float(cs_exp_2)))
					
					#dtn.write('<Peak2D F1 ="%s" F2 ="%s" annotation="%s" intensity = "9105828.0" type="1"/>\n'%(cs_exp_1,cs_exp_2,metabolite))
				score_final = math.sqrt(score_max*score_min)
				print metabolite, ": ", score_final
			
			print "_________________"
	
	return nb_components, list_point, list_point_1


	
def main_run(peak_list, list_file, tol_H, tol_merge):
	nb_components = 0
	list_point = []
	list_point_1 = []
	for name_file_new in list_file :
		nb_components, list_point, list_point_1 = iterative_matching(name_file_new, peak_list, tol_H, tol_merge,nb_components, list_point, list_point_1) #,dtn, dtn_theo)
		
	return nb_components, list_point, list_point_1
	
	
	
