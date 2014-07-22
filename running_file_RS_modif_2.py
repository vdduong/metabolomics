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

lockfile=threading.Lock()

class thread_matching(threading.Thread):
    def __init__(self,list_file,peak_list):
        threading.Thread.__init__(self)
        self.list_file=list_file
        self.peak_list=peak_list
        self.nb_compounds=0
        self.nb_checked=0
        self.halt=False
    def run(self):
        global index_file
        lockfile.acquire()
        try:
           index=index_file
           index_file+=1
        finally:
           lockfile.release()
        while index < len(self.list_file):
            file_bmrb=list_file[index]
            print "Starting file : "+file_bmrb
            list_assignment = {}
            time_initial=time.time()
            c = Metabolite(file_bmrb)
            if float(len(c.tocsy_pattern)) >1.0:
            print "Metabolite",time.time()-time_initial
            #len_initial = len(c.tocsy_pattern)
            time_initial=time.time()
            list_assignment = assignment(c, self.peak_list)
            print "dict_assignment",time.time()-time_initial
            time_initial=time.time()
            nb_peaks_found = peaks_found(list_assignment)
            print "peaks_found",time.time()-time_initial
            
            if float(len(c.tocsy_pattern)) >1.0 and float(nb_peaks_found)/float(len(c.tocsy_pattern)) >=0.7 :
               time_initial=time.time()
               c, dict_assignment, final_hausdorff, fraction_hausdorff,sum_ranking,best_pattern = iterative_matching(c, self.peak_list)
               print "iterative_matching",time.time()-time_initial
               time_initial=time.time()
               showAssignment(c, dict_assignment, self.peak_list, final_hausdorff, fraction_hausdorff, c.tocsy_pattern)
               print "showAssignment",time.time()-time_initial
               self.nb_compounds+=1
               print '_'*20
            else : pass
            print "Exiting : " + file_bmrb
            self.nb_checked+=1
            if self.halt :
               break
            lockfile.acquire()
            try:
               index=index_file
               index_file+=1
            finally:
               lockfile.release()
    def halt(self):
        self.halt=True

threads=[]
nb_threads_initial=threading.activeCount() # Main_thread is theoretically always active
nb_threads=1 # number of threads, defined by user
index_file=0

for index_thread in range(nb_threads):
    threads.append(thread_matching(list_file[0:10:],peak_list))
    threads[-1].start()

initial_time=time.time()
time.sleep(1)

while threading.activeCount() > nb_threads_initial:
      print "#"*30
      print "Time enlapsed : ",int(time.time()-initial_time)," seconds. Threads active : ", threading.activeCount()-nb_threads_initial
      nb_compounds=sum([thread.nb_compounds for thread in threads])
      nb_checked=sum([thread.nb_checked for thread in threads])
      print nb_compounds," compounds found on ",nb_checked," checked"
      time.sleep(1)

nb_compounds=sum([thread.nb_compounds for thread in threads])

print '%i compounds found'%(nb_compounds)
