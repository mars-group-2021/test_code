# MARS Reverse-Engineering
# University of Maryland Global Campus
# Authors: Alex Mancera, Stephen Panossian, Analia Treviño-Flitton

import numpy as np
from heartpy import smooth_signal
from scipy.signal import butter, filtfilt, iirnotch, savgol_filter

# direct copy from HeartPy source code (start)------------------


__all__ = ['filter_signal',
           'smooth_signal',
           'get_samplerate_datetime']

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


def smooth_signal(data, sample_rate, window_length=None, polyorder=3):
    if window_length == None:
        window_length = sample_rate // 10
        
    if window_length % 2 == 0 or window_length == 0: window_length += 1

    smoothed = savgol_filter(data, window_length = window_length,
                             polyorder = polyorder)
    return smoothed

def get_samplerate_datetime(datetimedata, timeformat='%H:%M:%S.%f'):
    from datetime import datetime
    datetimedata = np.asarray(datetimedata, dtype='str') #cast as str in case of np.bytes type
    elapsed = ((datetime.strptime(datetimedata[-1], timeformat) -
                datetime.strptime(datetimedata[0], timeformat)).total_seconds())
    sample_rate = (len(datetimedata) / elapsed)
    return sample_rate
    
# direct copy from HeartPy source code (end)------------------










## Program functions begin
## Splits concatenated files & returns them as lists
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

## Creates patient folder
def new_folder():
    import os
    path = os.path.join(os.getcwd(),patient)
    try: 
        os.makedirs(path, exist_ok = True) 
        print('Patient folder created successfully') 
    except OSError as error: 
        print('Patient folder cannot be created')
    os.chdir(path)


## Opens files & returns them to lists
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
        
# parses hdr data and saves in a nested dictionary
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


# determines if csv is pre-aligned from the beginning
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
                while curr_line_time != file_dat['Time'][-1] + two_ms:
                    file_dat['Time'].append(file_dat['Time'][-1] + two_ms)
                    for key, items in file_dat.items():
                        if key != 'Time':
                            items.append(0.5)
                    if curr_line_time - two_ms <= file_dat['Time'][-1]:
                        break
                print(' Gap filled')
            last_line_time = curr_line_time # saves current time for reference on next line
        else: # if timestamp not already given
            last_line_time += two_ms # adds 2ms to saved reference time
        file_dat['Time'].append(last_line_time)

        # records values from csv and fills in blank columns/values
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

## Called from & returns to print_out
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

def flat_peak_reduct(lst):
    rep_val = None

    steep_slope = 0.5
    not_flat = 0.05
    i_start = None
    v_start = None
    i_end = None
    v_end = None
    ss_check1 = 'N'
    ss_check2 = 'N'
    scale_max = 3
    ss_1 = 0
    
    for i, v in enumerate(lst):
        
        if i > 0 and rep_val is None and v != 0.5 and v == lst[i - 1]:
            rep_val = v
            rv_start = i-1
        if rep_val is not None and v != rep_val:
            rv_end = i
            len_rep_val = rv_end - rv_start
            if len_rep_val >= 10 and ss_check1 == 'N':
                lst[rv_start:rv_end] = [0.5]*len_rep_val
            rep_val = None
        
        if i > 0:
            slope = v - lst[i-1]
            if i_start == None and abs(slope) >= not_flat:
                i_start = i-1
                v_start = lst[i-1]
            if abs(slope) >= steep_slope:
                ss_check1 = 'Y'
                ss_1 = slope
                if i_start == None:
                    i_start = i-1
                    v_start = lst[i-1]
            if slope != 0 and ss_1 != 0 and ss_check1 == 'Y' and abs(slope) >= not_flat and int(slope/abs(slope)) == int(ss_1/abs(ss_1)*(-1)):
                ss_check2 = 'Y'
            if i_start != None and ss_check1 == 'N' and i - i_start >= 50:
                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'
            if i_start != None and ss_check2 == 'Y' and slope != 0 and abs(slope) <= not_flat and rep_val == None: 
                if i_start == 0:
                    v_start = v
                i_end = i
                v_end = v
                len_ss = i_end - i_start
                fill = []
                v_dif = v_end - v_start
                step = v_dif / len_ss
                val = v_start

                vals = lst[i_start:i_end]
                most_dif_val = max(map(lambda x: abs(x-v_start),vals))
                if  most_dif_val >= scale_max:
                    for gap in range(len_ss):
                        fill.append(val)
                        val += step
                    lst[i_start:i_end] = fill

                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'
                ss_1 =  0
                
    return lst

def tall_peak_reduct(lst):
    steep_slope = 0.3
    not_flat = 0.01
    i_start = None
    v_start = None
    i_end = None
    v_end = None
    ss_check1 = 'N'
    ss_check2 = 'N'
    scale_max = 2
    ss_1 = 0
    
    for i, v in enumerate(lst):
        if i > 0:
            slope = v - lst[i-1]
            if i_start == None and abs(slope) >= not_flat:
                i_start = i-1
                v_start = lst[i-1]
            if abs(slope) >= steep_slope:
                ss_check1 = 'Y'
                ss_1 = slope
                if i_start == None:
                    i_start = i-1
                    v_start = lst[i-1]
            if slope != 0 and ss_1 != 0 and ss_check1 == 'Y' and abs(slope) >= not_flat and int(slope/abs(slope)) == int(ss_1/abs(ss_1)*(-1)):
                ss_check2 = 'Y'
            if i_start != None and ss_check1 == 'N' and i - i_start >= 50:
                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'
            if i_start != None and ss_check2 == 'Y' and abs(slope) <= not_flat:
                i_end = i
                v_end = v
                len_ss = i_end - i_start
                vals = lst[i_start:i_end]

                
                fill = []
                v_dif = v_end - v_start
                step = v_dif / len_ss
                val = v_start
                most_dif_val = max(map(lambda x: abs(x-v_start),vals))
                if  most_dif_val >= scale_max:
                    for gap in range(len_ss):
                        fill.append(val)
                        val += step
                    lst[i_start:i_end] = fill
                
                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'
                ss_1 =  0

    return lst

## Called from & returns to print_out
def corrected_dict(dictionary):
    corrected = {}
    for key, items in dictionary.items():
        if key == 'Time':
            corrected['Time'] = tuple(dictionary['Time'])
        else:
            sample_rate = round(get_samplerate_datetime(dictionary['Time'], timeformat = '%Y-%m-%d %H:%M:%S.%f%z'), 3)
            corrected[key] = tuple(smooth_signal(
                                    tall_peak_reduct(
                                        remove_baseline_wander(
                                                flat_peak_reduct(
                                                                dictionary[key]
                                                                 ) # for flat_peak_reduct
                                                  ,sample_rate) # for remove_baseline_wander
                                                    ) # for tall_peak_reduct
                                                ,sample_rate,window_length=16, polyorder=3) # for smooth_signal
                                   )
            print('node',key,'corrected')                       
    return corrected

## Printing options
def plot_out(final_dat, patient):
    import plotly
    import plotly.express as px

    new_folder()
    
    print('\nWhich set of data would you like to use?')
    print(' 1) Original, unaltered data')
    print(' 2) Data with baseline correction, noise reduction,')
    print('    and flat/steep peak reduction')
    while True:
        op1 = int(input('Choice: '))
        if op1 in (1, 2):
            break
        print('Invalid selection')
    print()
    if op1 == 1:
        data_pref = final_dat
        title = 'Unaltered'
    if op1 == 2:
        corrected = corrected_dict(final_dat)
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

        # Call alignment checker & outfile writer
        final_dat = time_out(csv, file_dat, nodes_not_2ms, align_cont)

        # Call printing options
        plot_out(final_dat, patient)

    else:
        print("Please check the file name/directory & try again")

    print('Done')
