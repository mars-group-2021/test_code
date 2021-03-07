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

#name new csv file
position = csv_name.find('.')
strip_pos=len(csv_name)-position
if (position != -1):
    new_csv_name = csv_name[:-strip_pos] +"_new" + ".csv"
else:
    new_csv_name = csv_name + ".csv"

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
print(starttime)
csv.close()

#NOTE:I've tried doing this portion with csv reader, but it doesn't like to work with wfdb.io for some reason

#remove first column
df=pd.read_csv(csv_name, index_col=None, usecols=[1,2], low_memory=False)
df_mod=df.apply(pd.to_numeric, errors='coerce')
df_mod.to_csv(new_csv_name)
#some issues with the output:
##it adds an index column. I can't figure out how to get rid of it
###i've tried removing the index_col parameter and that doesn't help
##it keeps changing the first value in the second column to -40.46.1 

#to find data types, run the next two lines.
##datatypes=df.dtypes
##print(datatypes)
#datatypes with our csv and with the model csv were the same. no idea what's wrong


#convert csv to hea & dat
#can parse out units from hdr file
#this works on the practice csv, but not ours.
#fs is sample frequency. Not provided. 360 is the most common i've seen.
#wfdb.csv2mit(new_csv_name,fs=360, units='mV')

