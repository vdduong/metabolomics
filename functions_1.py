def shift(name_file) :
    """ return the list of chemical shifts 
    the input is the database file name
    the output is the shift list in the form of a dictionary :
        the keys are the names of atom
        the items are the shifts """
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
              shift_list[key] = item
            else :
            	item = str(data_line[8])
            	shift_list[key] = item
          except IndexError as err :
            break
        else : break
           
    return shift_list
          

				

def shortest_path(top_a, top_b, m) : 
    """ return the distance between two protons
    the input : the numbers corresponding to the two protons and the adjacent matrix
    the ouput is in the form of an integer """
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
    """ return the tocsy pattern if the covalent distance between two protons
        is less than 4 bonds.
        The output is in the form of a LIST """

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
            if key == top_a : cs_a = shift_list[key]
        for key in shift_list :
            if key == top_b : cs_b = shift_list[key]
        cross_peak = (cs_a, cs_b)
        pattern_tocsy_list.append(cross_peak)
    
    pattern_tocsy_list = remove_duplicate(pattern_tocsy_list)
    for item in pattern_tocsy_list :
        peak_nb+=1
        pattern_tocsy[peak_nb] = item
    return pattern_tocsy
