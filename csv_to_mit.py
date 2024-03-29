# Convert csv file to wfdb
# MARS Reverse-Engineering
# University of Maryland Global Campus
# Authors: Maura Franz

##Before running script, make sure you have downloaded wfdb
##pip install wfdb
##if it still does not run, make sure it is the current version (3.3.0)

import pandas as pd
import re
import datetime
import wfdb
from wfdb.io.record import (csv2mit)

#parse out info from hdr file
def hdrparse(file_name):    
    ids=[]
    labels=[]
    units=[]
    frequency=[]
    hdr=open(file_name)
    for line in hdr:
        parts = line.split(",")
        #added "ID" before the id number to appease csv2mit
        ids.append("ID"+parts[0].strip("{id: "))
        labels.append(parts[1].strip("label: "))
        units.append(parts[2].strip("unit: "))
        frequency.append(int(parts[3].strip("period: ").strip(" ms}\n")))
    hdr.close()
    return(ids,units,frequency)

#name new csv
def namecsv(file_name):
    position = file_name.find('.')
    strip_pos=len(file_name)-position
    if (position != -1):
        new_csv_name = file_name[:-strip_pos] +"_new" + ".csv"
    else:
        new_csv_name = file_name + ".csv"
    return new_csv_name

#find date and time
def findtimedate(file_name):
    csv=open(file_name)
    #find lines that start with four numbers
    for inline in csv:
        if (re.match ("^([0-9]{4})", inline)):
            timeline = inline
        else:
            break
    #split it into a list using comma as delineator (because data is in the first line too)
    firstline = timeline.split(',')
    #first portion of the first line is the date and time
    starttimestring = firstline[0]
    #split it into list
    starttimeparts = re.search("(.+) [\S]*$",starttimestring)
    #single out the first part of the list (date and time)
    starttime = starttimeparts.group(1)
    #split time and date into list using space as delineator & assign to variables
    time_parts=starttime.split(' ')
    time=time_parts[1]
    date=datetime.datetime.strptime(time_parts[0],"%Y-%m-%d").date().strftime("%d/%m/%Y")
    csv.close()
    return(time, date)

#find next valid value in the column
def findnextvalue(column, position):
    value = column[position]
    while value  == ' ':
        position += 1
        value = column[position]
    return value
        
#turn dataframe columns to list, replace empty cells with average of last value and next value
def repaircolumn(column):
    newcolumn=[]
    position = 0
    for value in column:
        if value == ' ':
            position+=1
            nextvalue = float(findnextvalue(column, position))
            #find average of last value and next valid value
            newvalue = (currentvalue + nextvalue)/2
            #append the new value to the list
            newcolumn.append(newvalue)
        elif type(value)==float:
            newcolumn.append(value)
            currentvalue=value
            position+=1
        else:
            newcolumn.append(float(value.strip()))
            currentvalue=float(value)
            position+=1
    return newcolumn

#determine sample frequency:
def samplefrequency(frequency_list):
    frequency_list = [int(i) for i in frequency_list]
    sample_frequency = 2000000
    for fs in frequency:
        if fs < sample_frequency:
            sample_frequency = fs
    hz= 1000//sample_frequency
    return hz



##main
csv_name = input("Enter csv name: ")
hdr_name = input("Enter hdr name: ")

#extract info from hdr
hdrcolumns=hdrparse(hdr_name)

#unpack tuple
(ids, units, frequency)=hdrcolumns
print("IDS: ", ids)
print("Units: ", units)
print("frequency: ", frequency)

#find new csv name
new_csv_name = namecsv(csv_name)

#extract date and time from csv
timedate = findtimedate(csv_name)

#unpack tuple
(time, date) = timedate

#read original csv into dataframe
numbers=[]
for n in range(1,(len(ids)+1)):
    numbers.append(n)
print("Columns to be used: ", numbers)
try:
    df=pd.read_csv(csv_name, header=None, names=ids, usecols=numbers, low_memory=False)
except:
    print("An error occured reading the csv file.")
    
#create new dataframe
df_mod=pd.DataFrame()

#cycle columns from data frame through repair function
for columnname in ids:
    print("Working on column ", columnname)
    #extract column from dataframe
    columnlist=df[columnname].tolist()
    #new column with no missing values
    repairedcolumn = repaircolumn(columnlist)
    #add repaired column to new dataframe
    df_mod[columnname] = repairedcolumn

#find sample frequency
sfreq = samplefrequency(frequency)

#turn dataframe into a csv
df_mod.to_csv(new_csv_name, index=False)

#convert csv file to hea and dat
print("Converting to MIT format")
try:
    wfdb.io.csv2mit(new_csv_name,fs=sfreq, units=units, samps_per_frame=1, base_time=time, base_date=date)
except:
    print("An error occurred when converting to MIT format")
