### NETWORK ANCHORING
### includes :
### a complex function ranking the assignments issued from the data comparison
### - functions defining alignment and ranking probabilities
### - ranking each assignment due to its relative strength compared to others
### - good assignments support each other
### a function generating the best matching pattern from the ranking table
### a function computing the local maximum flow 

from class_metabolite import *
from iterative_matching import *
from spectrum import *
import math

tol_H = 0.08
tol_alignment = 0.04

def ranking(metabolite, peak_list, dict_assignment) :
	network = {}
	for key_assignment in dict_assignment.keys() :
		item_assignment = dict_assignment[key_assignment]
		sum_proba_local = 0.0
		if len(item_assignment) > 0 : 
			for assignment in item_assignment :
				sum_proba_local+=assignment[1]
			for assignment in item_assignment :
				network[(key_assignment, assignment[0])] = assignment[1]#float('%.6f'%(assignment[1]/sum_proba_local))
	print network
	def proba_ranking(assignment_a, assignment_b) :
		"""
		ranking probability
		"""
		def proba_alignment(assignment_a, assignment_b) :
			"""
			alignment probability
			"""
			global tol_alignment
			cs_theo_a_1 = float(metabolite.tocsy_pattern[assignment_a[0]][0])
			cs_theo_a_2 = float(metabolite.tocsy_pattern[assignment_a[0]][1])
			cs_theo_b_1 = float(metabolite.tocsy_pattern[assignment_b[0]][0])
			cs_theo_b_2 = float(metabolite.tocsy_pattern[assignment_b[0]][1])

			cs_exp_a_1 = float(peak_list[assignment_a[1]][0])
			cs_exp_a_2 = float(peak_list[assignment_a[1]][1])
			cs_exp_b_1 = float(peak_list[assignment_b[1]][0])
			cs_exp_b_2 = float(peak_list[assignment_b[1]][1])

			distance_min = min(abs(cs_theo_a_1-cs_theo_b_1-cs_exp_a_1+cs_exp_b_1), \
								abs(cs_theo_a_2-cs_theo_b_2-cs_exp_a_2+cs_exp_b_2))
			# return the minimum deviation value within the two coordinates
			proba_alignment = math.exp(-0.5*distance_min**2/tol_alignment**2)
			return proba_alignment

		if 1 : 
			proba_ranking = 0.0
			time_ranking = 0.0
			for key_ranking in network.keys() :
				assignment_c = key_ranking
				time_ranking+=network[assignment_c]
				proba_ranking+=proba_alignment(assignment_a, assignment_c)*\
							proba_alignment(assignment_c, assignment_b)*\
							network[assignment_c]
			if time_ranking > 0.0 :
				proba_ranking = proba_ranking/time_ranking
				return proba_ranking
			else :
				return proba_ranking

	list_show_ranking = [] # initialization
	for key_a in network.keys() :
		assignment_a = key_a
		sum_ranking = 0.0
		time_ranking = 0.0
		for key_b in network.keys() :
			assignment_b = key_b
			sum_ranking+=proba_ranking(assignment_a, assignment_b)*\
						network[assignment_b]
			time_ranking+=network[assignment_b]
		if time_ranking > 0.0 :
			sum_ranking = sum_ranking/time_ranking
			list_show_ranking.append((assignment_a, sum_ranking))
		else :
			list_show_ranking.append((assignment_a, sum_ranking))
	list_show_ranking.sort(key = lambda x:x[1], reverse = True)
	# list sorting to show biggest elements first 
	return list_show_ranking

def best_pattern_generation(list_show_ranking) :
	"""
	return the best feasible pattern obtained from experimental data
	and the total flow based on the total assignment
	"""
	dict_ranking = {}
	for item_ranking in list_show_ranking :
		dict_ranking[item_ranking[0]] = item_ranking[1]
	sum_ranking = 0.0
	pattern = {}
	for item_ranking in list_show_ranking :
		if (item_ranking[0][0] not in [x[0] for x in pattern]) and \
			(item_ranking[0][1] not in [x[1] for x in pattern]) :
			pattern[item_ranking[0]] = item_ranking[1]
			sum_ranking+=item_ranking[1]
		else : pass
	return sum_ranking, pattern

def local_flow(metabolite, peak_list, best_pattern):
	"""
	print the local (best) flow corresponding to the best pattern
	"""
	global tol_H
	dict_local_assignment = {}
	for key_tocsy in metabolite.tocsy_pattern.keys():
		dict_local_assignment[key_tocsy] = list()
		for key_pattern in best_pattern.keys():
			if key_pattern[0] == key_tocsy :
				peak_tocsy = metabolite.tocsy_pattern[key_tocsy]
				peak_exp = peak_list[key_pattern[1]]
				proba_matching = math.exp(-0.5*(float(peak_tocsy[0])-float(peak_exp[0]))**2/tol_H**2)*\
							math.exp(-0.5*(float(peak_tocsy[1])-float(peak_exp[1]))**2/tol_H**2)
				dict_local_assignment[key_tocsy].append((key_pattern[1], proba_matching))
	list_show_ranking_local = ranking(metabolite, peak_list,dict_local_assignment)
	sum_ranking_local, pattern = best_pattern_generation(list_show_ranking_local)
	#print sum_ranking_local
	return None

#########################
#### TEST for the module
if __name__ == '__main__':
	path_local = str(os.getcwd())
	peak_list, volume_list = spectrum(path_local + "/Exp_peak_list/43-sn10.peaks")
	#peak_list[5000]=('7.808','7.808')
	peak_list[5001]=('7.53','7.808')
	#peak_list[73] = ('7.655','7.653')
	name_file = path_local + "/database/bmse000408.str"
	c = Metabolite(name_file)
	shift = c.initial_shift(name_file)
	c = Metabolite(name_file, shift)
	c.tocsy_pattern = c.tocsy()
	dict_assignment = assignment(c, peak_list)
	for key in dict_assignment.keys() :
		print key, dict_assignment[key]
	print '_'*10
	list_show_ranking = ranking(c,peak_list,dict_assignment)
	sum_ranking, best_pattern = best_pattern_generation(list_show_ranking)
	print sum_ranking
	print c.tocsy_pattern
	for item in list_show_ranking :
		print item
	c, dict_assignment, final_hausdorff,fraction_hausdorff,sum_ranking,best_pattern = iterative_matching(c, peak_list)
	print sum_ranking
	print c.tocsy_pattern
	#print '_'*10
	#for key in dict_assignment.keys():
	#	print key, dict_assignment[key]
	print '_'*10
	#list_show_ranking = ranking(c, peak_list, dict_assignment)
	#for item in list_show_ranking :
	#	print item, '*'
	local_flow(c,peak_list,best_pattern)
	print '_'*10
	for key_tocsy in c.tocsy_pattern.keys():
		print key_tocsy, c.tocsy_pattern[key_tocsy]
