
import pandas as pd
import wfdb
from wfdb.io.record import (csv2mit)


def hdrparse(hdr_name):    
    ids=[]
    labels=[]
    units=[]
    frequency=[]
    hdr=open(hdr_name)
    for line in hdr:
        parts = line.split(",")
        ids.append("ID"+parts[0].strip("{id: "))
        labels.append(parts[1].strip("label: "))
        units.append(parts[2].strip("unit: "))
        frequency.append(int(parts[3].strip("period: ").strip(" ms}\n")))
    hdr.close()
    return(ids,units,frequency)

#name new csv
def namecsv(csv_name):
    position = csv_name.find('.')
    strip_pos=len(csv_name)-position
    if (position != -1):
        new_csv_name = csv_name[:-strip_pos] +"_new" + ".csv"
    else:
        new_csv_name = csv_name + ".csv"
    return new_csv_name

#turn dataframe columns to list, replace empty cells with 0.1
def repaircolumn(column):
    newcolumn=[]
    for value in column:
        if value == ' ':
            newcolumn.append(0.1)
        elif value==float:
            newcolumn.append(value)
        else:
            newcolumn.append(value.strip())
    return newcolumn

#determine sample frequency:
def samplefrequency(frequency_list):
    frequency_list = [int(i) for i in frequency_list]
    sample_frequency = 0
    for fs in frequency:
        if fs > sample_frequency:
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
#find new csv name
new_csv_name = namecsv(csv_name)
#read original csv into dataframe
df=pd.read_csv(csv_name, header=None, names=ids, low_memory=False)
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
wfdb.io.csv2mit(new_csv_name,fs=sfreq, units=units)

