# MARS Reverse-Engineering
# University of Maryland Global Campus
# Authors: Alex Mancera, Stephen Panossian, Analia Treviño-Flitton

import datetime as dt
import os
from decimal import Decimal
import numpy as np
import heartpy as hp
from scipy.signal import butter, filtfilt, iirnotch, savgol_filter
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objs as go






def cat_file_parser(cat_file):
    del_num=0
    hdr =[]
    intro = open(cat_file, 'r').readlines()
    ## Split the hdr contents into own list
    for i in range(len(intro)):
        if intro[i].startswith('{id:'):
            hdr.append(intro[i])
            if i > del_num:
                del_num = i

## delete hdr contents from csv list
    for j in range(0,(del_num+1)):
        (intro.pop(0))

## Pack lists into tuples & return
    contents_tup = (hdr, intro)
    return contents_tup




## Opens files & returns them to global lists
def file_opener(patient):

    exists = 0
    cur_dir= os.listdir(os.getcwd())

    while exists == 0:
        if (patient+'.csv' and patient + '.hdr') in cur_dir:
            print('Both csv and hdr files found')
            hdr = open(patient + '.hdr', 'r')
            csv = open(patient + '.csv', 'r')
            contents = (hdr, csv)
            exists=1

        elif (patient+'.cat') in cur_dir:
            print("The csv and hdr file have been found in a concatenated file")
            contents = cat_file_parser(patient + '.cat')
            exists=1


        elif (patient+'.csv') in cur_dir:
           contents = cat_file_parser(patient + '.csv')
           exists=1

        else:
            print('not found')
            #reprompt?

    return contents




def hdr_data(hdr):
    # parses hdr data and saves in a nested dictionary
    hdr_dat = {}
    node = 1
    nodes_4ms = []
    for line in hdr:
        l_split = line.strip('{').strip('\n').strip('}').split(', ')
        hdr_dat[node] = {} #creates new nested dictionary with each line read in .hdr file
        for item in l_split:
            hdr_dat[node][item.split()[0].strip(':')] = item.split()[1]
            if item.split()[0].strip(':') == 'period' and item.split()[1] == '4ms':
                nodes_4ms.append(node)
        node += 1

    hdr_info=(hdr_dat, node, nodes_4ms)

    return  hdr_info




## Program begins
if __name__ == '__main__':
    patient = input('Enter patient name (for now a,b, or c): ')

    # Call file opener, returns tuple of the hdr & csv lists
    contents = file_opener(patient)
    # Unpack tuple
    (hdr, csv) = contents

    ## Call hdr_data- dictionary builder & node info
    hdr_info = hdr_data(hdr)
    (hdr_dat, node, nodes_4ms) = hdr_info



# determines if csv is pre-aligned from the beginning
aligned=''
num_nodes = len(hdr_dat) # number of line in .hdr file
num_cols_start = len(csv.readline().split(',')) # number of columns at the beginning of .csv file
num_cols_mis = num_nodes + 1 - num_cols_start # missing number of columns at the beginning of .csv file
if  num_cols_mis == 0:  # if no columns missing 
    print('csv file pre-aligned')
    aligned='yes'
else:
    print('csv file not pre-aligned')
    aligned='no'


print('Collecting data')

file_dat = {} 
file_dat['Time'] = []
for key in hdr_dat:
    if hdr_dat[key]['label'] in file_dat:
        file_dat[hdr_dat[key]['label']+'(2)']=[]
    else:
        file_dat[hdr_dat[key]['label']]=[]
       
last_line_time=''
two_ms = dt.timedelta(milliseconds=2)
csv.seek(0)
for line in csv:
    if not line.startswith(','): # if timestamp already given
        curr_line_time = dt.datetime.strptime(line.split(', ')[0], '%Y-%m-%d %H:%M:%S.%f %z') # current time is kept as provided
        if last_line_time != '' and curr_line_time != last_line_time + two_ms: # checks for time gap
            print('Time gap at',line.split(', ')[0])
        last_line_time = curr_line_time # saves current time for reference on next line
    else: # if timestamp not already given
        last_line_time += two_ms # adds 2ms to saved reference time
    file_dat['Time'].append(last_line_time)

    # generates printed lines to outfile after header line
    for n in range(1,num_cols_start):
        file_dat[list(file_dat)[n]].append(float(line.split(', ')[n].strip() or 0)) # fills in blanks with recorded values or 0 where data was not recorded
    if aligned == 'no':
        for n in range(num_cols_mis):
            if len(line.split(', ')) < num_nodes +1: # if empty columns exist
                file_dat[list(file_dat)[num_cols_start+n]].append(0) # fills in blanks with 0 where columns are empty
            else:
                file_dat[list(file_dat)[num_cols_start+n]].append(float((line.split(', ')[num_cols_start+n].strip() or 0))) # fills in blanks with recorded values or 0 where data was not recorded

    if len(nodes_4ms) > 0: #if there are any 4ms nodes
        for node in nodes_4ms:
            if len(file_dat[list(file_dat)[node]])>3 and file_dat[list(file_dat)[node]][-2] == 0 and file_dat[list(file_dat)[node]][-1] != 0:
                file_dat[list(file_dat)[node]][-2] = (file_dat[list(file_dat)[node]][-3]+file_dat[list(file_dat)[node]][-1])/2
                    
print('Done pulling data')




# direct copy from HeartPy source code (start)------------------

__all__ = ['filter_signal',
           'hampel_filter',
           'hampel_correcter',
           'smooth_signal']

def butter_lowpass(cutoff, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a
def butter_highpass(cutoff, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_bandpass(lowcut, highcut, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def filter_signal(data, cutoff, sample_rate, order=2, filtertype='highpass', return_top = False):       # changed 'lowpass' to 'highpass' for more accuracy
    if filtertype.lower() == 'lowpass':
        b, a = butter_lowpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'highpass':
        b, a = butter_highpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'bandpass':
        assert type(cutoff) == tuple or list or np.array, 'if bandpass filter is specified, \
cutoff needs to be array or tuple specifying lower and upper bound: [lower, upper].'
        b, a = butter_bandpass(cutoff[0], cutoff[1], sample_rate, order=order)
    elif filtertype.lower() == 'notch':
        b, a = iirnotch(cutoff, Q = 0.005, fs = sample_rate)
    else:
        raise ValueError('filtertype: %s is unknown, available are: \
lowpass, highpass, bandpass, and notch' %filtertype)

    filtered_data = filtfilt(b, a, data)
    
    if return_top:
        return np.clip(filtered_data, a_min = 0, a_max = None)
    else:
        return filtered_data
    
def remove_baseline_wander(data, sample_rate, cutoff=0.05):
    return filter_signal(data = data, cutoff = cutoff, sample_rate = sample_rate, filtertype='notch')


# direct copy from HeartPy source code (end)------------------






print('Plotting points')

plots = len(file_dat)-1
fig = make_subplots(rows=len(file_dat)-1, cols=1, shared_xaxes=True, vertical_spacing=0.02)

for key, items in file_dat.items():
    if key != 'Time':
        fig.add_trace(go.Scatter(x=file_dat['Time'][0::4],y=file_dat[key][0::4],name='Unedited '+key, line=dict(color='rgb(251,180,174)')),row=plots, col=1)
        fig.add_trace(go.Scatter(x=file_dat['Time'][0::4],y=remove_baseline_wander(file_dat[key], 2000.0)[0::4],name='Corrected '+key, line=dict(color='royalblue')),row=plots, col=1)
        fig.update_yaxes(title_text=key, row=plots, col=1)
        plots += (-1)

fig.update_layout(title_text='Data from '+patient+'.csv and '+patient+'.hdr')

print('Printing graph')
plotly.offline.plot(fig, filename=patient+'_blc.html')

print('Done')
