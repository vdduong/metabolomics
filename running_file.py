import matplotlib.pyplot as plt
from functions import *
import datetime

#src = open("/Users/James/Metabo_1/Exp_peak_list/tocsy1GHz_QC.peaks","r")
#src.readline()

tol_H = 0.09
tol_merge = 0.05
tol_hausdorff = tol_H
peak_list = dict()

t0 = datetime.datetime.now()


#for line in src :
#  line = re.split(" +", line)
#	nb_peak = int(line[1])
#	cs_1 = line[2]
#	cs_2 = line[3]
#	peak_list[nb_peak] = (cs_1, cs_2)

name_file_spectrum = "/Users/James/Metabo_1/Exp_peak_list/43-sn10.peaks"
peak_list = spectrum(name_file_spectrum)
list_peak = []
for key in peak_list.keys():
	#print key, peak_list[key]
	cs_1 = float(peak_list[key][0])
	cs_2 = float(peak_list[key][1])
	list_peak.append((cs_1,cs_2))


#list_file = data_file()
name_file_1 = '/Users/James/Metabo_1/database/bmse000408.str'

list_file = []
list_file.append(name_file_1)
#nb_components, list_point, list_point_1 = main_run(peak_list, list_file, tol_H, tol_merge)
nb_components = 0
list_point, list_point_1 = [],[]
nb_components, list_point, list_point_1 = iterative_matching(name_file_1, peak_list, tol_H, tol_merge, nb_components, list_point, list_point_1)
plt.figure()
a = plt.gca()
a.set_xlim([-0.5,10.5])
a.set_ylim([-0.5,10.5])
a.set_ylim(a.get_ylim()[::-1])
a.set_xlim(a.get_xlim()[::-1])
plt.axhline(y=0.0)
plt.axvline(x=0.0)
plt.axhline(y=10.0)
plt.axvline(x=10.0)
#plt.plot(*zip(*list_point_1), marker='o', color='b', ls='')
plt.plot(*zip(*list_peak), marker='o', color='g', ls='')
plt.plot(*zip(*list_point), marker='o', color='r', ls='')
plt.savefig("/Users/James/Iterameta/spectrum.png")	


print nb_components, " are found in the mixture"

	
delta_t = datetime.datetime.now() - t0
print delta_t
