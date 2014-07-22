### CLASS METABOLITE
### Metabolite is defined as an object with attributs and functions.
### Attributs : name, matrix, shift, tocsy_pattern
### - name : name of the metabolite
### - shift : the chemical shifts of the metabolite
### - matrix : the adjacent matrix of the metabolite
### - tocsy_pattern : the TOCSY pattern of the metabolite
### Functions : name_metabolite, matrix_adjacent, initial_shift, tocsy
###
### at the very end of this module, we define two functions : shortest_path using the 
### algorithm of Dijkstra in order to define the pattern of tocsy ; and remove_duplicate
### in order to remove duplicates in a list

import re # module to split space
import os

class Metabolite():
	def __init__(self, name_file):
		self.name = self.name_metabolite(name_file)
		self.matrix = self.matrix_adjacent(name_file)
		self.shift = self.initial_shift(name_file)
		self.tocsy_pattern = self.tocsy()

	def __str__(self):
		print self.name, self.matrix, self.shift, self.tocsy_pattern

	def name_metabolite(self,name_file):
		"""return the name of the metabolite"""
		data_file = file(name_file)
		for line in data_file:
			if "_Entry.Title" in line :
				line = line.rstrip('\n\r')
				data_line = re.split(" ", line)
				name_metabo = list()
				for item_data_line in data_line :
					if len(str(item_data_line)) > 0 and item_data_line != "_Entry.Title":
						name_metabo.append(str(item_data_line))
		name_metabolite = ''
		for item in name_metabo :
			name_metabolite = name_metabolite + item
		return name_metabolite

	def matrix_adjacent(self, name_file):
		"""
		return the connection/adjacent matrix 
		""" 
		m_adjac = list()
		try:
			data_file=file(name_file)
			line=data_file.readline()
			while not "_Chem_comp_bond.Comp_ID" in line :
				line=data_file.readline()
			data_file.readline()
			line=data_file.readline()
			while not "stop_" in line :
				line = line.rstrip('\n\r')
				data_line = re.split(' +', line)
				m_adjac.append((str(data_line[4]), str(data_line[5])))
				line=data_file.readline()
		finally:
			data_file.close()
		return m_adjac

	def initial_shift(self, name_file) :
		"""
		load chemical shift into a dictionary 
		"""
		shift_list = dict()
		try :
			data_file = file(name_file)
			line=data_file.readline()
			while line != '':
				while not '_Atom_chem_shift.ID\n' in line and line != '':
					line=data_file.readline()
				nb_column=1
				nb_column_key=0
				nb_column_key=0
				while not '_Atom_chem_shift.Assigned_chem_shift_list_ID\n' in line and line != '':
					line=data_file.readline()
					nb_column+=1
					if '_Atom_chem_shift.Atom_ID\n' in line : nb_column_key=nb_column
					if '_Atom_chem_shift.Val\n' in line : nb_column_item=nb_column
				data_file.readline()
				line=data_file.readline()
				while not "stop_" in line and line != '':
					try :
						line = line.rstrip('\n\r')
						data_line = re.split(" +", line)
						key = str(data_line[nb_column_key])
						item = str(data_line[nb_column_item])
						shift_list[key] = item
						line=data_file.readline()
					except IndexError as err :
							break
		finally:
			data_file.close()
		return shift_list

	def tocsy(self):
		"""
		return the TOCSY pattern in the form of dictionary
		"""
#		pattern_tocsy = dict()
		pattern_tocsy=list()
		peak_nb = 0
		pattern_tocsy_list = list()
		shift_list_proton = list()
		for key in self.shift.keys() :
			if 'H' in key :
				shift_list_proton.append(key) # we take only the protons
		pattern_peak = list()
		index_proton_1 = 0
		while index_proton_1 < len(shift_list_proton):
			top_a = shift_list_proton[index_proton_1]
			pattern_peak.append((top_a,top_a))
			index_proton_2 = index_proton_1 + 1
			while index_proton_2 < len(shift_list_proton):
				top_b = shift_list_proton[index_proton_2]
				distance = shortest_path(top_a, top_b, self.matrix)
				if distance <= 4 :
					pattern_peak.append((top_a, top_b))
					pattern_peak.append((top_b, top_a)) # add symmetrical peak
				index_proton_2+=1
			index_proton_1+=1
		for item_peak in pattern_peak:
			cs_a = 0.0
			cs_b = 0.0
			top_a = item_peak[0]
			top_b = item_peak[1]
			for key in self.shift.keys():
				if key == top_a : cs_a = self.shift[key]
				if key == top_b : cs_b = self.shift[key]
			cross_peak = (cs_a, cs_b)
			pattern_tocsy_list.append(cross_peak)
		pattern_tocsy = remove_duplicate(pattern_tocsy_list)
#		for item in pattern_tocsy_list :
#			peak_nb+=1
#			pattern_tocsy[peak_nb] = item
		return pattern_tocsy

#### Outside functions :

def remove_duplicate(mylist) :
	if mylist :
		mylist.sort()
		last = mylist[-1]
		for i in range(len(mylist)-2, -1, -1):
			if last == mylist[i]:
				del mylist[i]
			else :
				last = mylist[i]
	return mylist

def shortest_path(top_a, top_b, m) : 
    """ 
    return the bonding distance that separates two atoms 
    """
    found = (top_a == top_b)
    distance = 0
    
    set_checked = set() # set of atoms already visited

    set_test=set_checked # set of atoms which are going to be visited during this iteration
    set_test.add(top_a) # the first atom to test is the one given
    
    while (found == False) and (distance <= 4):
          set_next = set() # set of atoms which will be visited for the next iteration
          for item in m : # lisiting all bonds
              if (item[0] in set_test) : # if the atom visited must be tested
                 set_next.add(item[1]) # its neighbour will be tested next time
              elif (item[1] in set_test) : # idem
                 set_next.add(item[0])
              else : pass    
          distance+=1      
          if top_b in set_next :
             found = True
          else :
             set_checked=set_checked.union(set_test)
             set_next=set_next.difference(set_checked)
             set_heteroatomes=set()
             for item in set_next :
                 if not 'C' in item:
                    set_heteroatomes.add(item)
             set_next=set_next.difference(set_heteroatomes)
             set_test=set_next
    return distance

########################
### TEST of the module

if __name__ == '__main__':
	path_local = str(os.getcwd())
	name_file = path_local + "/database/bmse000408.str"
	c = Metabolite(name_file)
	shift = c.initial_shift(name_file)
	c = Metabolite(name_file, shift)
	c.tocsy_pattern = c.tocsy()
	for i in range(1,len(c.tocsy_pattern)+1) :
		print i, c.tocsy_pattern[i][0], c.tocsy_pattern[i][1]



