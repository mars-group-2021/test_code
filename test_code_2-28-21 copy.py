################################################################################
#
# Script:           MARS Reverse-Engineering
#
# Programmers:      Alex Mancera
#                   Analia Treviño-Flitton
#                   Stephen Panossian
#
# Class:            BIOT 670 9040 Spring 2021
#                   University of Maryland Global Campus
#
# Assignment:       Capstone Project - Stand alone Application
#
# Purpose:          Perform the following tasks to view data imported into
#                   the GE MARS ECG System:
#
#                   1) 
#
# Date:             March 2, 2021
#
# Version:          0.0
#
# Notes:            1. 
#
# References:       
#
################################################################################

import numpy as np
import heartpy as hp
from heartpy import hampel_filter, hampel_correcter, smooth_signal
from heartpy.peakdetection import fit_peaks
from heartpy.datautils import rolling_mean, _sliding_window
from scipy.signal import butter, filtfilt, iirnotch, savgol_filter

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


def filter_signal(data, cutoff, sample_rate, order=2, filtertype='highpass',
                  return_top=False):  # changed 'lowpass' to 'highpass' for more accuracy
    if filtertype.lower() == 'lowpass':
        b, a = butter_lowpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'highpass':
        b, a = butter_highpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'bandpass':
        assert type(cutoff) == tuple or list or np.array, 'if bandpass filter is specified, \
cutoff needs to be array or tuple specifying lower and upper bound: [lower, upper].'
        b, a = butter_bandpass(cutoff[0], cutoff[1], sample_rate, order=order)
    elif filtertype.lower() == 'notch':
        b, a = iirnotch(cutoff, Q=0.005, fs=sample_rate)
    else:
        raise ValueError('filtertype: %s is unknown, available are: \
lowpass, highpass, bandpass, and notch' % filtertype)

    filtered_data = filtfilt(b, a, data)

    if return_top:
        return np.clip(filtered_data, a_min=0, a_max=None)
    else:
        return filtered_data


def remove_baseline_wander(data, sample_rate, cutoff=0.05):
    return filter_signal(data=data, cutoff=cutoff, sample_rate=sample_rate, filtertype='notch')

def hampel_filter(data, filtsize=6):
    #generate second list to prevent overwriting first
    #cast as array to be sure, in case list is passed
    output = np.copy(np.asarray(data)) 
    onesided_filt = filtsize // 2
    for i in range(onesided_filt, len(data) - onesided_filt - 1):
        dataslice = output[i - onesided_filt : i + onesided_filt]
        mad = MAD(dataslice)
        median = np.median(dataslice)
        if output[i] > median + (3 * mad):
            output[i] = median
    return output



def hampel_correcter(data, sample_rate):
    return data - hampel_filter(data, filtsize=int(sample_rate))



def quotient_filter(RR_list, RR_list_mask = [], iterations=2):
    if len(RR_list_mask) == 0:
        RR_list_mask = np.zeros((len(RR_list)))
    else:
        assert len(RR_list) == len(RR_list_mask), \
        'error: RR_list and RR_list_mask should be same length if RR_list_mask is specified'

    for iteration in range(iterations):
        for i in range(len(RR_list) - 1):
            if RR_list_mask[i] + RR_list_mask[i + 1] != 0:
                pass #skip if one of both intervals is already rejected
            elif 0.8 <= RR_list[i] / RR_list[i + 1] <= 1.2:
                pass #if R-R pair seems ok, do noting
            else: #update mask
                RR_list_mask[i] = 1
                #RR_list_mask[i + 1] = 1

    return np.asarray(RR_list_mask)


def smooth_signal(data, sample_rate, window_length=None, polyorder=3):
    if window_length == None:
        window_length = sample_rate // 10
        
    if window_length % 2 == 0 or window_length == 0: window_length += 1

    smoothed = savgol_filter(data, window_length = window_length,
                             polyorder = polyorder)

    return smoothed
# direct copy from HeartPy source code (end)------------------










## Program functions begin

############################################################################
#
# Function:     cat_file_parser()
#
# Purpose:      Splits concatenated files & returns them as lists
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        cat_file = file whose data and header are concatenated
#
# Output:       contents_tup = data lists returned as tuples
#
# Notes:        1. 
#
############################################################################

## Need to add an error catch if it is just a plain csv file but hdr is missing
def cat_file_parser(cat_file):
    del_num = 0
    hdr = []
    with open(cat_file, 'r') as f:
        intro = f.readlines()

    ## Split the hdr contents into own list
    for i in range(len(intro)):
        if intro[i].startswith('{id:'):
            hdr.append(intro[i])
            if i > del_num:
                del_num = i

    ## delete hdr contents from csv list
    for j in range(0, (del_num + 1)):
        (intro.pop(0))

    ## Pack lists into tuples & return
    contents_tup = (hdr, intro)
    # (cat_file + '.cat').close()
    return contents_tup


############################################################################
#
# Function:     new_folder()
#
# Purpose:      Creates patient folder.
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        None.
#
# Output:       None.
#
# Notes:        None. 
#
############################################################################

def new_folder():
    import os
    path = os.path.join(os.getcwd(),patient)
    try: 
        os.makedirs(path, exist_ok = True) 
        print('Patient folder created successfully') 
    except OSError as error: 
        print('Patient folder cannot be created')
    os.chdir(path)


############################################################################
#
# Function:     file_opener()
#
# Purpose:      Opens files & returns them to lists
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        patient  = input string indicating patient file name
#
# Output:       contents = patient data as tuples.
#
# Notes:        None. 
#
############################################################################

def file_opener(patient):
    import os
    
    exists = 0
    cur_dir = os.listdir(os.getcwd())

    # Reads in hdr & csv to lists
    while exists == 0:
        if (patient + '.csv' and patient + '.hdr') in cur_dir:
            print('Both csv and hdr files found')
            with open(patient + '.hdr', 'r') as h:
                hdr = h.readlines()
            with open(patient + '.csv', 'r') as v:
                csv = v.readlines()
            contents = (hdr, csv)
            exists = 1
            return contents

        elif (patient + '.cat') in cur_dir:
            print("The csv and hdr file have been found in a concatenated file")
            contents = cat_file_parser(patient + '.cat')
            exists = 1
            return contents


        elif (patient + '.csv') in cur_dir:
            contents = cat_file_parser(patient + '.csv')
            exists = 1
            return contents

        else:
            print('not found')
            contents = "end"
            return contents

  
############################################################################
#
# Function:     hdr_data()
#
# Purpose:      Parses hdr data line by line and saves in a nested
#               dictionary
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        hdr      = input header file
#
# Output:       hdr_info = nested dictionary
#
# Notes:        None. 
#
############################################################################

def hdr_data(hdr):
    hdr_dat = {}
    node = 1
    nodes_not_2ms = {}
    for line in hdr:
        l_split = line.strip('{').strip('\n').strip('}').split(', ')
        hdr_dat[node] = {}  # creates new nested dictionary with each line read in .hdr file
        for item in l_split:
            hdr_dat[node][item.split()[0].strip(':')] = item.split()[1]
            if item.split()[0].strip(':') == 'period' and item.split()[1] != '2ms':
                nodes_not_2ms[int(item.split()[1].strip('ms'))] = node
        node += 1

    hdr_info = (hdr_dat, nodes_not_2ms)

    return hdr_info


############################################################################
#
# Function:     align_plus_dat()
#
# Purpose:      Determines if csv is pre-aligned from the beginning
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        hdr_dat = data from data file header
#               csv     = comma separated data file
#
# Output:       dat_plus_align = aligned data file as a dictionary
#
# Notes:        None. 
#
############################################################################

def align_plus_dat(hdr_dat, csv):
    aligned = ''
    num_nodes = len(hdr_dat)  # number of lines in .hdr file
    num_cols_start = len((csv[0]).split(','))  # number of columns at the beginning of .csv file
    num_cols_mis = num_nodes + 1 - num_cols_start  # missing number of columns at the beginning of .csv file
    if num_cols_mis == 0:  # if no columns missing
        print('csv file pre-aligned')
        aligned = 'yes'
    else:
        print('csv file not pre-aligned')
        aligned = 'no'


    file_dat = {}
    file_dat['Time'] = []
    for key in hdr_dat:
        if hdr_dat[key]['label'] in file_dat:
            file_dat[hdr_dat[key]['label'] + '(2)'] = []
        else:
            file_dat[hdr_dat[key]['label']] = []
    hdr_dat.clear()


    align_contents = (num_nodes, num_cols_start, num_cols_mis, aligned)
    dat_plus_align = (align_contents, file_dat)

    return dat_plus_align


############################################################################
#
# Function:     fpeaks()
#
# Purpose:      Fit with varying peak detection thresholds given a heart
#               rate signal.
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        patient = patient name string
#               hrdata  = patient data from .csv file
#
# Output:       fp = dictionary object of patient data with fitted peaks
#
# Notes:        1. Adapted from HeartPy Quickstart Guide and
#                  Docs > API Reference > Peakdetection > fit_peaks
#
############################################################################

def fpeaks(patient, hrdata):

    # load patient data without header info
    #hrdata = hp.get_data(patient + '.csv', column_name = 'datetime')

    # estimate sample rate
    sample_rate = hp.get_samplerate_datetime(hrdata, timeformat='%Y-%m-%d %H:%M:%S.%f %z')

    # need to compute rolling mean before peak fitting
    windowsize = 0.75 # rolling_mean default
    rol_mean   = rolling_mean(hrdata, windowsize, sample_rate)

    working_data, measures = hp.process(hrdata, sample_rate)

    # find min, max beats per minute
    # max heart rate = 220 - adult (18+ y. o.) patient's age [bpm] per CDC
    bpmmin =  40.0 # default
    bpmmax = 180.0 # default

    # fit_peaks: HeartPy function fitting with varying peak detection
    #            thresholds given a heart rate signal
    fp = fit_peaks(hrdata, rol_mean, sample_rate, bpmmin, bpmmax)
    
    return fp


############################################################################
#
# Function:     negdata_flip()
#
# Purpose:      Call HeartPy flip_signal to invert any signals with negative
#               mV peaks, bringing them back to normal ECG signal form.
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        data_section = section to be evaluated (1d list/numpy array)
#
# Output:       out = 1d array of normal ECG signal
#
# Notes:        1. Adapted from HeartPy Quickstart Guide and
#                  Docs > API Reference > Peakdetection > fit_peaks
#
############################################################################

def negdata_flip(data_section):
    # invert any raw negative mV peaks to normal ECG
    # input:  data = section to be evaluated (1d list or numpy array)
    # output: out  = 1d array

    # defaults: enhance_peaks = F, keep_range = T
    keep_range = T
    out_array  = flip_signal (data_section, keep_range)

    return out


############################################################################
#
# Function:     time_out()
#
# Purpose:      
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        csv           =
#               file_dat      =
#               nodes_not_2ms =
#               align_cont    =  tuple of aligned data in dictionary
#
# Output:       file_dat = 1d array of normal ECG signal
#
# Notes:        1. Adapted from HeartPy Quickstart Guide and
#                  Docs > API Reference > Peakdetection > fit_peaks
#
############################################################################

def time_out(csv, file_dat, nodes_not_2ms, align_cont):
    import datetime as dt

    # Unpack tuple from align_plus_dat
    (num_nodes, num_cols_start, num_cols_mis, aligned) = align_cont


    last_line_time = ''
    two_ms = dt.timedelta(milliseconds=2)

    for line in csv:
        if not line.startswith(','): # if timestamp already given
            curr_line_time = dt.datetime.strptime(line.split(', ')[0], '%Y-%m-%d %H:%M:%S.%f %z') # current time is kept as provided
            if last_line_time != '' and curr_line_time != last_line_time + two_ms: # checks for time gap
                time_dif = (curr_line_time - last_line_time).total_seconds()
                if time_dif < 180:
                    print(' >',time_dif,'second gap starting at',last_line_time)
                else:
                    print(' >',round(time_dif/60,3),'minute gap starting at',last_line_time)
            last_line_time = curr_line_time # saves current time for reference on next line
        else: # if timestamp not already given
            last_line_time += two_ms # adds 2ms to saved reference time
        file_dat['Time'].append(last_line_time)

        # generates printed lines to outfile after header line
        for n in range(1, num_cols_start):
            file_dat[list(file_dat)[n]].append(float(line.split(', ')[n].strip() or 0.5))  # fills in blanks with recorded values or 0.5 where data was not recorded
        if aligned == 'no':
            for n in range(num_cols_mis):
                if len(line.split(', ')) < num_nodes + 1:  # if empty columns exist
                    file_dat[list(file_dat)[num_cols_start + n]].append(0.5)  # fills in blanks with 0.5 where columns are empty
                else:
                    file_dat[list(file_dat)[num_cols_start + n]].append(float((line.split(', ')[num_cols_start + n].strip() or 0.5)))  # fills in blanks with recorded values or 0.5 where data was not recorded


        if len(nodes_not_2ms) > 0:  # if there are any non-2ms nodes
            for interval , node in nodes_not_2ms.items():
                node_list = file_dat[list(file_dat)[node]]
                min_start = (interval//2)+1
                empty_vals = node_list[-(interval//2):-1]
                if len(node_list) > min_start and set(empty_vals) == {0.5} and node_list[-1] != 0.5:
                    mis_val_int = (node_list[-min_start] - node_list[-1])/(len(empty_vals)+1)
                    mult = 1
                    for n in range(len(empty_vals)):
                        node_list[-2-n] = node_list[-1] + (mult * mis_val_int)
                        mult += 1
        
    print('Done pulling data')
    return file_dat


############################################################################
#
# Function:     print_dict()
#
# Purpose:      
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        dictionart =
#               outfile    =
#
# Output:       file_dat = 1d array of normal ECG signal
#
# Notes:        Called from & returns to print_out.
#
############################################################################

def print_dict(dictionary, outfile):
    header = []
    for key in dictionary:
        header.append(key)
    outfile.write(", ".join(header) + '\n')

    print_line = []
    for key, items in dictionary.items():
        for i in range(len(dictionary[key])):
            for k in range(len(dictionary)):
                print_line.append(str(dictionary[list(dictionary)[k]][i]))
            outfile.write(", ".join(print_line) + '\n')
            print_line.clear()


############################################################################
#
# Function:     flat_peak_reduct()
#
# Purpose:      
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        lst =
#
# Output:       lst = 
#
# Notes:        
#
############################################################################

def flat_peak_reduct(lst):
    rep_val = None
    for i, v in enumerate(lst):
        if rep_val is None and v != 0.5 and v == lst[i - 1]:
            rep_val = v
            rv_start = i-1
        if rep_val is not None and v != rep_val:
            rv_end = i
            len_rep_val = rv_end - rv_start
            if len_rep_val >= 5:
                lst[rv_start:rv_end] = [0.5]*len_rep_val
            rep_val = None
    return lst


############################################################################
#
# Function:     corrected_dict()
#
# Purpose:      
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        dictionary =
#
# Output:        = 
#
# Notes:        1. dictionary format:
#               dict_dat = {'Time' : [datetime, datetime, datetime...],
#                           'node' : [value, value, value...],
#                           'node2' : [value, value, value...]  }
#
############################################################################

## Called from & returns to print_out
def corrected_dict(dictionary):
    corrected = {}
    for key, items in dictionary.items():
        if key == 'Time':
            corrected['Time'] = tuple(dictionary['Time'])
        else:
            sample_rate = 2000.0
            corrected[key] = tuple(remove_baseline_wander(
                                        smooth_signal(
                                            flat_peak_reduct(
                                                dictionary[key]
                                            )
                                        ,sample_rate, window_length=14, polyorder=3)
                                    ,sample_rate)
                                   )
            #corrected[key] = tuple(fpeaks(corrected[key], sample_rate))

    return corrected


############################################################################
#
# Function:     plot_out()
#
# Purpose:      Provide printing options to user.
#
# Programmer:   
#
# Class:        BIOT 670 9040 Spring 2021
#               University of Maryland Global Campus
#
# Assignment:   Capstone Project - Stand alone Application
#
# Date:         March 2, 2021
#
# Version:      1.0
#
# Input:        final_dat =
#               patient   = string of patient file name.
#
# Output:       None.
#
# Notes:        1. dictionary format:
#               dict_dat = {'Time' : [datetime, datetime, datetime...],
#                           'node' : [value, value, value...],
#                           'node2' : [value, value, value...]  }
#
############################################################################

## Printing options
def plot_out(final_dat, patient):
    import plotly
    import plotly.express as px

    new_folder()
    
    print(
        '\nWhich set of data would you like to use?\n 1) Original, unaltered data \n 2) Data with baseline correction, noise reduction, and flat peak reduction')
    while True:
        op1 = int(input('Choice: '))
        if op1 in (1, 2):
            break
        print('Invalid selection')

    if op1 == 1:
        data_pref = final_dat
        title = 'Unaltered'
    if op1 == 2:
        negdata_flip(data_section)
        corrected = corrected_dict(final_dat)
        #sample_rate=2000.0
        #corrected = fpeaks(corrected, sample_rate)
        data_pref = corrected
        final_dat.clear()
        title = 'Corrected'
   
    print('\nWould you like to print '+title.lower()+' data?')
    while True:
        op2 = input('Y/N: ').upper()
        if op2 in ('Y','N'):
            break
        print('Invalid selection')
        
    if op2 == 'Y':
        with open(patient + '_concat_'+title.lower()+'.csv', 'w') as out:
            print_dict(data_pref, out)
        print('\nData printed to outfile')


    print('\nPlotting points')    
    for key, items in data_pref.items():
        if key != 'Time':
            fig = px.line(x=data_pref['Time'][0::4],
                          y=data_pref[key][0::4],
                          labels={'x':'', 'y':key})
            items=tuple()
            fig.update_layout(title_text=title + ' data from ' + patient+', node '+key,showlegend=False)
            plotly.offline.plot(fig, filename= patient+'_'+key+'_'+title.lower()+'.html')
            print('node',key,'chart printed')

    data_pref.clear()



## Program begins
if __name__ == '__main__':
    patient = input('Enter patient name (for now a,b, or c): ')
    # Call file opener
    contents = file_opener(patient)


    if contents != 'end':
        # If file present, then unpack tuple
        (hdr, csv) = contents

        ## Call hdr_data- dictionary builder & node info
        hdr_info = hdr_data(hdr)
        (hdr_dat, nodes_not_2ms) = hdr_info

        # Call dat_file generation
        dat_align = align_plus_dat(hdr_dat, csv)
        (align_cont, file_dat) = dat_align

        # call fit peaks
        file_dat = fpeaks(patient, file_dat)

        # Call alignment checker & outfile writer
        final_dat = time_out(csv, file_dat, nodes_not_2ms, align_cont)

        # Call printing options
        plot_out(final_dat, patient)

    else:
        print("Please check the file name/directory & try again")

    print('Done')
