# comparison between BMRB and TOCSY pattern from Madison
# tyrosine, tryptophan, and hippuric acid
from class_metabolite_RS_modif import *
from iterative_matching_RS_modif import *
from final_assignment_RS_modif import *
from assignment_clustering_RS_modif import *

import re
import os
import glob
import threading
import time

def list_all_file(path_local) :
  """ load the metabolites into a list """
  list_file = glob.glob(path_local + "/database/*.str")
  return list_file

tol_H = 0.1 # matching tolerance, defined by user


path_local = str(os.getcwd()) # to know where we are

name_file_spectrum = path_local + "/Exp_peak_list/PeakList.peaks"
peak_list, volume_list = spectrum(name_file_spectrum) # experimental peak-list, volume list

list_file = list_all_file(path_local)

class thread_matching(threading.Thread):
    def __init__(self,list_file,peak_list,nb_threads,index_thread):
        threading.Thread.__init__(self)
        self.list_file=list_file
        self.peak_list=peak_list
        self.nb_threads=nb_threads
        self.index_thread=index_thread
        self.nb_compounds=0
        self.nb_checked=0
        self.halt=False
    def run(self):
        for file_bmrb in self.list_file[self.index_thread:50:self.nb_threads]:
            #print "Starting file : "+file_bmrb
            dict_assignment = {}
            c = Metabolite(file_bmrb)
            shift = c.initial_shift(file_bmrb)
            c = Metabolite(file_bmrb, shift)
            c.tocsy_pattern = c.tocsy()
            
            tocsy_initial = c.tocsy_pattern
            len_initial = len(c.tocsy_pattern)
            
            dict_assignment = assignment(c, self.peak_list)
            nb_peaks_found = peaks_found(dict_assignment)
            
            if float(len(c.tocsy_pattern)) >1.0 and float(nb_peaks_found)/float(len(tocsy_initial)) >=0.7 :
               c, dict_assignment, final_hausdorff, fraction_hausdorff,sum_ranking,best_pattern = iterative_matching(c, self.peak_list)
               showAssignment(c, dict_assignment, self.peak_list, final_hausdorff, fraction_hausdorff, tocsy_initial)
               self.nb_compounds+=1
               print '_'*20
            else : pass
            #print "Exiting : " + file_bmrb
            self.nb_checked+=1
            if self.halt :
               break
    def halt(self):
        self.halt=True

threads=[]
nb_threads_initial=threading.activeCount() # Main_thread is theoretically always active
nb_threads=10 # number of threads, defined by user

for index_thread in range(nb_threads):
    threads.append(thread_matching(list_file,peak_list,nb_threads,index_thread))
    threads[-1].start()

initial_time=time.time()
time.sleep(1)

nb_compounds=sum([thread.nb_compounds for thread in threads])
while threading.activeCount() > nb_threads_initial:
      print "_"*30
      print "Time enlapsed : ",int(time.time()-initial_time)," seconds. Threads active : ", threading.activeCount()-nb_threads_initial
      nb_compounds=sum([thread.nb_compounds for thread in threads])
      nb_checked=sum([thread.nb_checked for thread in threads])
      print nb_compounds," compounds found on ",nb_checked," checked"
      time.sleep(1)

nb_compounds=sum([thread.nb_compounds for thread in threads])

print '%i compounds found'%(nb_compounds)
