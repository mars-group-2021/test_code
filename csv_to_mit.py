# Convert csv file to wfdb
# MARS Reverse-Engineering
# University of Maryland Global Campus
# Authors: Maura Franz

##Before running script, make sure you have downloaded wfdb
##pip install wfdb
##if it still does not run, make sure it is the current version (3.3.0 I think)

import fileinput
import sys
import pandas as pd
import wfdb
import re

from wfdb.io import _header
from wfdb.io import _signal
from wfdb.io import download

from wfdb.io.record import (Record, MultiRecord, rdheader, rdrecord, rdsamp, wrsamp,
                            dl_database, edf2mit, mit2edf, wav2mit, mit2wav, wfdb2mat,
                            csv2mit, sampfreq, signame, SIGNAL_CLASSES)
from wfdb.io._signal import est_res, wr_dat_file
from wfdb.io.annotation import (Annotation, rdann, wrann, show_ann_labels,
                                show_ann_classes, ann2rr)
from wfdb.io.download import get_dbs, get_record_list, dl_files, set_db_index_url
from wfdb.io.tff import rdtff

csv_name= input("Enter csv name: ")
hdr_name= input("Enter hdr name: ")

#function
#name new csv file
position = csv_name.find('.')
strip_pos=len(csv_name)-position
if (position != -1):
    new_csv_name = csv_name[:-strip_pos] +"_new" + ".csv"
else:
    new_csv_name = csv_name + ".csv"


#Function
csv=open(csv_name)
#find the date and time
#we might not need it, but it's good to have it just in case

for inline in csv:
    if (re.match ("^([0-9]{4})", inline)):
        timeline = inline
    else:
        break
    
#split it into a list using comma as delineator

firstline = timeline.split(',')
starttimestring = firstline[0]
starttimeparts = re.search("(.+) [\S]*$",starttimestring)
starttime = starttimeparts.group(1)
#split time and date into list using space as delineator
time_parts=starttime.split(' ')
time=time_parts[1]
date=time_parts[0]
print("time: ", time, "date: ", date)
csv.close()

#extract information from hdr file
ids=[]
labels=[]
units=[]
frequency=[]
hdr=open(hdr_name)
for line in hdr:
    parts = line.split(",")
    ids.append(parts[0].strip("{id: "))
    labels.append(parts[1].strip("label: "))
    units.append(parts[2].strip("unit: "))
    frequency.append(parts[3].strip("period: ").strip(" ms}\n"))

#remove first column from csv, add headers
df=pd.read_csv(csv_name, header=None, names=labels, low_memory=False)
df_mod=df.apply(pd.to_numeric, errors='coerce')
df_mod.to_csv(new_csv_name, index=False)

#to find data types, run the next two lines.
##datatypes=df.dtypes
##print(datatypes)
#datatypes with our csv and with the model csv were the same. no idea what's wrong

#determine sample frequency:
frequency = [int(i) for i in frequency]
sample_frequency = 0
for fs in frequency:
    if fs > sample_frequency:
        sample_frequency = fs
hz= 1000//sample_frequency

#convert csv to hea & dat
#this works on the practice csv, but not ours.   
wfdb.io.csv2mit(new_csv_name,fs=hz, units=units)
