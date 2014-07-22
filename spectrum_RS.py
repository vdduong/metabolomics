### SPECTRUM
### spectrum routine to extract peak list from ATNOS automated peak-picking output

import re
import os

def spectrum(name_file) :
    """ load the experimental peak list into two dictionaries
    one for peak-list
    other for volume-list 
    """
    spectrum = file(name_file)
    peak_list = list() # list the peak list
    volume_list = list() # list the volume list
    nb_line = 13
    nb_count = 0
    while nb_count < nb_line :
          spectrum.readline()
          nb_count+=1
    for line_spectrum in spectrum :
        line_spectrum = line_spectrum.rstrip('\n\r')
        data_spectrum = re.split(" +", line_spectrum)
        cs_1 = data_spectrum[2]
        cs_2 = data_spectrum[3]
        signal = data_spectrum[6]
        peak_list.append((cs_1, cs_2))
        volume_list.append(float(signal))
    spectrum.close()      
    return peak_list, volume_list

#############################
### TEST of the module
if __name__ == '__main__':
	path_local = str(os.getcwd())
	peak_list, volume_list = spectrum(path_local + "/Exp_peak_list/43-sn10.peaks")
	for key in peak_list.keys():
		print key, peak_list[key], volume_list[key]
