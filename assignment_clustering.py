### ASSIGNMENT CLUSTERING
### assignment clustering using two main methods : k-means clustering and 
### quality threshold clustering



from spectrum import *
from class_metabolite import *
from iterative_matching import *
from final_assignment import *
import math
import re
import os

tol_H = 0.09
tol_alignment = 0.02
peak_list = {}

def list_point_append(dict_assignment, peak_list, metabolite) :
	"""
	create the list of secondary chemical shifts
	scattering (X,Y) in order to take the k-mean algo in count"""
	list_point = []
	for key_assignment in dict_assignment.keys():
		peak_tocsy = metabolite.tocsy_pattern[key_assignment] # name the peak in tocsy
		cs_tocsy_1 = float(peak_tocsy[0]) # assigned chemical shifts
		cs_tocsy_2 = float(peak_tocsy[1])
		item_assignment = dict_assignment[key_assignment]
		if len(item_assignment) > 0 :
			for item in item_assignment :
				cs_exp_1 = float(peak_list[item[0]][0])
				cs_exp_2 = float(peak_list[item[0]][1])
				coord_1 = cs_tocsy_1 - cs_exp_1
				coord_2 = cs_tocsy_2 - cs_exp_2
				coord_1 = float(('%.4f')%(coord_1))
				coord_2 = float(('%.4f')%(coord_2))
				list_point.append((coord_1,coord_2,key_assignment,item[0]))
		else : 
			pass
	return list_point

#####################
## K-MEANS CLUSTERING

def distanceBetweenPoints(peak_1, peak_2):
	global tol_H
	distance = math.sqrt((peak_1[0] - peak_2[0])**2 + \
					(peak_1[1] - peak_2[1])**2)/ tol_H
	return distance

def assignment_centre(list_centre, list_point) :
	"""
	classify the elements of list_point into a dictionary
	whose keys are centres"""
	global tol_H
	dict_centre = dict()
	for item in list_centre :
		dict_centre[item] = list() # initialization of dictionary of centres

	for peak_1 in list_point :
		min_local = 100.0 
		listDistanceLocal = []
		for key_point in dict_centre.keys() :
			listDistanceLocal.append((distanceBetweenPoints(peak_1, key_point,tol_H), \
														key_point))
			if distanceBetweenPoints(peak_1, key_point, tol_H) < min_local :
				min_local = distanceBetweenPoints(peak_1, key_point, tol_H)
		for item in listDistanceLocal :
			if item[0] == min_local :
				dict_centre[item[1]].append(peak_1)
	return dict_centre

def update_centre(dict_centre) :
	"""update the centres until a stationary state"""
	list_centre = list()
	for key_centre in dict_centre.keys() :
		new_coord_1 = 0.0
		new_coord_2 = 0.0
		item_centre = dict_centre[key_centre] #it's a cluster of points
		for point in item_centre :
			new_coord_1+=point[0]
			new_coord_2+=point[1]
		if len(item_centre) > 0 :
			new_coord_1 = '%.4f'%(new_coord_1/float(len(item_centre)))
			new_coord_2 = '%.4f'%(new_coord_2/float(len(item_centre)))
		new_coord_1 = float(new_coord_1)
		new_coord_2 = float(new_coord_2)
		new_key_centre = (new_coord_1, new_coord_2)
		list_centre.append(new_key_centre)
	return list_centre

def k_means_clustering(list_point) :
	"""
	randomly chose the centres from list_point
	to fit the number k in k_means
	then 
	make the k_means algorithm
	assignment_centre
	update_centre
	"""
	global tol_H
	nb_stab_limit = 100 # set a limit for ones that never converge
	filter_clustering = list() # to regroup the best clustering group of points
	filter_clustering_complet = list() # to regroup over the different clusters
                          # to form the best combinations
	k = round(math.sqrt(float(len(list_point))/2.0)) # thumb number formula
	list_centre = list()
	i = 1
	while i <=k :
          chococo = choice(list_point)
          if chococo not in list_centre :
            list_centre.append(chococo)
            i+=1
	nb_time_clustering = 1
	while 1 :
		dict_centre = assignment_centre(list_centre, list_point, tol_H)
		list_centre_1 = update_centre(dict_centre)
		if list_centre_1 == list_centre : 
			break
			#print 'situation stablized after %i times around %i centre'%(nb_time_clustering, \
											#len(list_centre))
		else : 
			list_centre = list_centre_1

		nb_time_clustering+=1
		if nb_time_clustering <= nb_stab_limit : pass
		else : 
			break
			#print 'situation has no convergent pattern'
	longest_element = 0
	for key in dict_centre.keys():
		if len(dict_centre[key]) > longest_element :
			longest_element = len(dict_centre[key])
	for key in dict_centre.keys():
		if len(dict_centre[key]) == longest_element :
			item_centre = dict_centre[key]
			for item_local in item_centre :
				filter_clustering.append((item_local[2], item_local[3]))
        for item in filter_clustering :
          filter_clustering_complet.append(item) # initialization by including all best elements
                                  # in the total group
        for key in dict_centre.keys():
          item_centre = dict_centre[key]
          for item_local in item_centre :
            if item_local[2] not in [item[0] for  item in filter_clustering] :
              filter_clustering_complet.append((item_local[2], item_local[3]))

        return filter_clustering, filter_clustering_complet

###############################
## QUALITY THRESHOLD CLUSTERING


def distanceBetweenPoints_qt(peak_1, peak_2):
	"""
	define the distance between two points in quality threshold clustering
	"""
	#distance_min = min(abs(peak_1[0]-peak_2[0]), abs(peak_1[1]-peak_2[1]))
	#distance = distance_min
	distance = math.sqrt( (peak_1[0] - peak_2[0])**2 + \
							(peak_1[1]-peak_2[1])**2)
	return distance

def distanceBetweenPoint_List_qt(peak_1, list_1) :
	list_distance = []
	for item in list_1 :
		distance = distanceBetweenPoints_qt(peak_1, item)
		list_distance.append(distance)
	list_distance.sort(reverse=True)
	return list_distance[0]


def qt_clustering(list_point):
	"""
	the items in list_point are initially ranked False
	every time an item is found to be ranked in a cluster due to its distance below
	the quality threshold tol alignment, it will be ranked True
	Whenever an item is ranked True, we no longer try to cluster this item
	"""
	global tol_alignment
	i = 0
	dict_clustering = {}
	dict_ranked = {}
	for item_point in list_point :
		dict_ranked[item_point] = False # all items are not ranked. initialization
	nb_clusters = 0
	while i <= len(list_point)-1 :
		point_i = list_point[i]
		if dict_ranked[point_i] == False :
			nb_clusters+=1
			dict_clustering[nb_clusters] = []
			for item_point in list_point :
				if len(dict_clustering[nb_clusters]) > 0 :
					list_local = dict_clustering[nb_clusters]
					if distanceBetweenPoint_List_qt(item_point, list_local) <= tol_alignment:
						dict_clustering[nb_clusters].append(item_point)
						dict_ranked[item_point]= True
				else :
					if distanceBetweenPoints_qt(point_i, item_point) <= tol_alignment :
						dict_clustering[nb_clusters].append(item_point)
						dict_ranked[item_point]= True
		else : pass
		i+=1
	return dict_clustering

def pattern_production(dict_assignment, peak_list, c, tocsy_initial):
	"""
	from the results of qt_clustering, we retain only the patterns that contain over 
	50 percent of the total peak numbers
	"""
	global tol_alignment
	dict_patterns = {}
	list_point = list_point_append(dict_assignment, peak_list, c)
	dict_clustering = qt_clustering(list_point)
	len_tocsy = len(tocsy_initial)
	nb_patterns = 0
	for key in dict_clustering.keys():
		if len(dict_clustering[key]) >= 0.7*float(len_tocsy) :
		 # the consistent pattern should contain more than 70 percent of the total pattern
			nb_patterns+=1
			dict_patterns[nb_patterns]= dict_clustering[key]
		else : pass
	return dict_patterns


#########################
## TEST of the module
if __name__ == '__main__':
	path_local = str(os.getcwd())
	peak_list, volume_list = spectrum(path_local + "/Exp_peak_list/43-sn10.peaks")
	#peak_list[5001]=('7.53','7.808')
	#name_file = path_local + "/database/bmse000408.str"
	list_file = database_bmrb()
	nb_compounds = 0
	for name_file in list_file :
		c = Metabolite(name_file)
		shift = c.initial_shift(name_file)
		c = Metabolite(name_file, shift)
		c.tocsy_pattern = c.tocsy()
		dict_assignment = assignment(c, peak_list)
		#for key in dict_assignment.keys():
		#	print key, dict_assignment[key]
		dict_patterns = pattern_production(dict_assignment,peak_list,c)
		if len(dict_patterns) > 0 and len(c.tocsy_pattern) >=4.0 :
			print c.name
			nb_compounds +=1
			for key in dict_patterns.keys():
				list_local = [(item[2], item[3]) for item in dict_patterns[key]]
				print 'quality-threshold clustering %i : '%key, list_local
				print '%.4f'%(float(len(dict_patterns[key]))/float(len(c.tocsy_pattern)))
				print '*'*5
			print '_'*30
	print nb_compounds

