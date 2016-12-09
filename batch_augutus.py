#!/usr/bin/env python

import os
import datetime
from threading import Thread
from Queue import Queue


class RUN:
    def __init__(self, file):
        self.file = file
    def do(self):
        #print "current running %s..." %(self.file)
        st = datetime.datetime.now()
        # main command here
        os.system("augustus --species=PhelliusNoxius --extrinsicCfgFile=extrinsic.E.cfg --singlestrand=true --alternatives-from-evidence=true --hintsfile=~/augustus/rna_hints/ALL/hints.gff --protein=on --introns=on --start=off --stop=off --cds=on --codingseq=on --uniqueGeneId=true %s > %s.gtf" %(self.file, self.file.replace(".fasta","")))
        td = datetime.datetime.now() - st
        print "finish %s run in %s" %(self.file, format(td))

# start time
total_start = datetime.datetime.now()

# create a queue
que = Queue()

# put jobs into queue
for file in sorted(os.listdir("./")):
    if file.endswith(".fasta"): # check filename suffix
        que.put(RUN(file))

# run job function
def doJob(*args):
    queue = args[0] 
    while queue.qsize() > 0: # same to: while not queue.empty()
          job = queue.get()
          job.do()

thread_list = []

# open 8 theads using for loop
for i in range(8):
    thd = Thread(target=doJob, args=(que,)) # args requires a tuple, so (que,)
    thread_list.append(thd)

# start 8 threads
for thread in thread_list:
    thread.start()

# wait for all threads are finish
for thread in thread_list:
    thread.join()

# old version
''' 
# open 3 threads
thd1 = Thread(target=doJob, name='Thd1', args=(que,))
thd2 = Thread(target=doJob, name='Thd2', args=(que,))
thd3 = Thread(target=doJob, name='Thd3', args=(que,))

thd1.start()
thd2.start()
thd3.start()

thd1.join()
thd2.join()
thd3.join()
'''

total_time = datetime.datetime.now() - total_start
print "\nAll jobs finish in %s!" %(format(total_time))
