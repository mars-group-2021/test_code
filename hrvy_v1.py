## MARS Group Reverse-Engineering
## MS: Biotechnology- Bioinformatics Capstone Spring 2021
## University of Maryland Global Campus
## Authors: Alex Mancera, Stephen Panossian, Analia Treviño-Flitton
## HRVY- Heart Rate Viewer in PYthon Version 1.0

import numpy as np
from heartpy import smooth_signal
from scipy.signal import butter, filtfilt, iirnotch, savgol_filter

# --------------------- HeartPy Source Code Start --------------------------------------------------
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

    smoothed = savgol_filter(data, window_length=window_length,
                             polyorder=polyorder)
    return smoothed


def get_samplerate_datetime(datetimedata, timeformat='%H:%M:%S.%f'):
    from datetime import datetime
    datetimedata = np.asarray(datetimedata, dtype='str')  # cast as str in case of np.bytes type
    elapsed = ((datetime.strptime(datetimedata[-1], timeformat) -
                datetime.strptime(datetimedata[0], timeformat)).total_seconds())
    sample_rate = (len(datetimedata) / elapsed)
    return sample_rate


def scale_data(data, lower=0, upper=1024):
    rng = np.max(data) - np.min(data)
    minimum = np.min(data)
    data = (upper - lower) * ((data - minimum) / rng) + lower
    return data


def enhance_peaks(hrdata, iterations=2):
    scale_data(hrdata)
    for i in range(iterations):
        hrdata = np.power(hrdata, 2)
        hrdata = scale_data(hrdata)
    return hrdata


def flip_signal(data, enhancepeaks=False, keep_range=True):
    data_mean = np.mean(data)
    data_min = np.min(data)
    data_max = np.max(data)

    # invert signal
    data = (data_mean - data) + data_mean

    if keep_range:
        # scale data so original range is maintained
        data = scale_data(data, lower=data_min, upper=data_max)
    if enhancepeaks:
        data = enhance_peaks(data)
    return data


# --------------------- HeartPy Source Code End ---------------------------------------------------


''' HRVY- Heart Rate Viewer in PYthon Version 1.0 '''


## Opens files & returns them as lists
def file_opener(dataset):
    import os
    contents = 'end'
    cur_dir = os.listdir(os.getcwd())
    hdr_found = False

    ## Check files in cwd for approved file types
    for file in cur_dir:
        ## Look for csv first
        if file.startswith(dataset) and file.endswith('.csv'):
            ## Checks for matching hdr file
            for second_file in cur_dir:
                if second_file.startswith(dataset) and second_file.endswith('.hdr'):
                    try:
                        with open(dataset + '.hdr', 'r') as h:
                            hdr = h.readlines()
                        with open(dataset + '.csv', 'r') as v:
                            csv = v.readlines()
                        print('Both csv and hdr files found')
                        contents = (hdr, csv)
                        hdr_found = True
                        break
                    except FileNotFoundError:
                        print('Not a matching csv and hdr file')

            ## If hdr wasn't found, check csv for hdr data
            if hdr_found != True:
                try:
                    contents = cat_file_parser(dataset + '.csv')
                    break
                except FileNotFoundError:
                    print('Not a csv file nested with hdr data')
                break

        ## If not one of the above, check for cat file
        elif file.startswith(dataset) and file.endswith('.cat'):
            try:
                print("The cat file has been found in the directory")
                contents = cat_file_parser(dataset + '.cat')
            except FileNotFoundError:
                print('Not a cat file')
            break

    ## If no file found
    if contents == 'end':
        print(
            'Cannot locate a viable file, please make sure data set name, the file type, and the directory are correct')

    return contents


  
## Splits concatenated files & returns them as lists to main
def cat_file_parser(cat_file):
    del_num = 0
    hdr = []

    with open(cat_file, 'r') as f:
        intro = f.readlines()

    ## Split the hdr contents into own list
    for i in range(len(intro)):
        if intro[i].startswith('{id:'):
            hdr.append(intro[i])
            del_num += 1


        ## If no hdr file contents found
        elif not intro[i].startswith('{id:') and del_num == 0:
            print('This is not a concatenated file, no hdr data found')
            contents = "end"
            return contents

    ## If it is a cat file, hdr contents will be deleted from csv list
    if del_num != 0:
        print('Concatenated data found')
        for j in range(0, del_num):
            del intro[0]

        ## Pack lists into tuples & return
        contents_tup = (hdr, intro)

        return contents_tup


      
## Parses hdr data and saves in a nested dictionary
def hdr_data(hdr):
    hdr_dat = {}
    node = 1
    nodes_not_2ms = {}

    ## Node information is formatted
    for line in hdr:
        l_split = line.strip('{').strip('\n').strip('}').split(', ')

        ## Creates new nested dictionary with each line read in .hdr file
        hdr_dat[node] = {}
        for item in l_split:
            hdr_dat[node][item.split()[0].strip(':')] = item.split()[1]
            if item.split()[0].strip(':') == 'period' and item.split()[1] != '2ms':
                nodes_not_2ms[int(item.split()[1].strip('ms'))] = node
        node += 1

    ## Packed into tuple and returned
    hdr_info = (hdr_dat, nodes_not_2ms)

    return hdr_info


  
## Determines if csv is pre-aligned
def align_plus_dat(hdr_dat, csv):
    aligned = ''

    ## Number of lines in .hdr file
    num_nodes = len(hdr_dat)

    ## Number of columns at the beginning of .csv file
    num_cols_start = len((csv[0]).split(','))

    ## Missing number of columns at the beginning of .csv file
    num_cols_mis = num_nodes + 1 - num_cols_start

    ## If no columns are missing
    if num_cols_mis == 0:
        print('No missing columns found, csv file pre-aligned')
        aligned = 'yes'
    else:
        print('Missing columns found, csv file not pre-aligned')
        aligned = 'no'

    ## Create new dictionary with the time for keys
    file_dat = {}
    file_dat['Time'] = []
    for key in hdr_dat:
        if hdr_dat[key]['label'] in file_dat:
            file_dat[hdr_dat[key]['label'] + '(2)'] = []
        else:
            file_dat[hdr_dat[key]['label']] = []
    hdr_dat.clear()

    ## Pack into tuples & return
    align_contents = (num_nodes, num_cols_start, num_cols_mis, aligned)
    dat_plus_align = (align_contents, file_dat)

    return dat_plus_align

  

## Checks for time gaps & fills is missing
def time_out(csv, file_dat, nodes_not_2ms, align_cont):
    import datetime as dt

    ## Unpack tuple from align_plus_dat
    (num_nodes, num_cols_start, num_cols_mis, aligned) = align_cont

    ## Set time object
    last_line_time = ''
    two_ms = dt.timedelta(milliseconds=2)

    for line in csv:
        ## If there is a timestamp, current time is kept as provided
        if not line.startswith(','):
            curr_line_time = dt.datetime.strptime(line.split(', ')[0], '%Y-%m-%d %H:%M:%S.%f %z')

            ## Checks for time gap
            if last_line_time != '' and curr_line_time != last_line_time + two_ms:
                time_dif = (curr_line_time - last_line_time).total_seconds()
                if time_dif < 180:
                    print(' >', time_dif, 'second gap starting at', last_line_time)

                ## If time gap is over a minute long
                else:
                    print(' >', round(time_dif / 60, 3), 'minute gap starting at', last_line_time)
                while curr_line_time != file_dat['Time'][-1] + two_ms:
                    file_dat['Time'].append(file_dat['Time'][-1] + two_ms)
                    for key, items in file_dat.items():
                        if key != 'Time':
                            items.append(0.5)
                    if curr_line_time - two_ms <= file_dat['Time'][-1]:
                        break
                print(' Gap filled')
            ## Saves current time for reference on next line
            last_line_time = curr_line_time

        ## If timestamp not already given adds 2ms to saved reference time
        else:
            last_line_time += two_ms
        file_dat['Time'].append(last_line_time)

        ## Records values from csv and fills in blank columns/values with recorded values
        ## Fills in blanks with recorded values or 0.5 (baseline) where data was not recorded
        for n in range(1, num_cols_start):
            file_dat[list(file_dat)[n]].append(float(line.split(', ')[n].strip() or 0.5))

        ## If empty columns exist blanks are filled with 0.5 where columns are empty
        if aligned == 'no':
            for n in range(num_cols_mis):
                if len(line.split(', ')) < num_nodes + 1:
                    file_dat[list(file_dat)[num_cols_start + n]].append(
                        0.5)
                ## Fills in blanks with recorded values or 0.5 where data was not recorded
                else:
                    file_dat[list(file_dat)[num_cols_start + n]].append(float((line.split(', ')[
                                                                                   num_cols_start + n].strip() or 0.5)))
        ## If there are any non-2ms nodes, derive values for 2ms interval
        if len(nodes_not_2ms) > 0:
            for interval, node in nodes_not_2ms.items():
                node_list = file_dat[list(file_dat)[node]]
                min_start = (interval // 2) + 1
                empty_vals = node_list[-(interval // 2):-1]
                if len(node_list) > min_start and set(empty_vals) == {0.5} and node_list[-1] != 0.5:
                    mis_val_int = (node_list[-min_start] - node_list[-1]) / (len(empty_vals) + 1)
                    mult = 1
                    for n in range(len(empty_vals)):
                        node_list[-2 - n] = node_list[-1] + (mult * mis_val_int)
                        mult += 1

    print('Done pulling data')
    return file_dat


  
## Data processing, printing, and plotting options
def plot_out(final_dat, dataset):
    import plotly
    import plotly.express as px

    ## Create data set folder in cwd
    new_folder()

    ## Present user with options
    print('\nWhich set of data would you like to use?')
    print(' 1) Original, unaltered data')
    print(' 2) Data with baseline correction, noise reduction,')
    print('    and flat/steep peak reduction')

    ## If baseline corrections etc selected, corrected_dict is called
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

    print('\nWould you like to print ' + title.lower() + ' data?')
    while True:
        op2 = input('Y/N: ').upper()
        if op2 in ('Y', 'N'):
            break
        print('Invalid selection')

    ## If saving locally for printing, print_dict is called
    if op2 == 'Y':
        with open(dataset + '_concat_' + title.lower() + '.csv', 'w') as out:
            print_dict(data_pref, out)
        print('\nData has been printed to the data set folder')

    ## Plotting data begins
    print('\nPlotting points')
    for key, items in data_pref.items():
        if key != 'Time':
            fig = px.line(x=data_pref['Time'][0::4],
                          y=data_pref[key][0::4],
                          labels={'x': '', 'y': key})
            items = tuple()
            fig.update_layout(title_text=title + ' data from ' + dataset + ', node ' + key, showlegend=False)

            ## Offline plot file named here
            plotly.offline.plot(fig, filename=dataset + '_' + key + '_' + title.lower() + '.html')
            print('node', key, 'chart printed')

    data_pref.clear()

    

## Creates data set folder in the current directory
def new_folder():
    import os

    ## Get path for cwd
    path = os.path.join(os.getcwd(), dataset)

    ## Make new folder named after data set
    try:
        os.makedirs(path, exist_ok=True)
        print('Data set folder created successfully')
        os.chdir(path)
    except OSError:
        print('Data set folder cannot be created')

        

## Data corrections occur here
def corrected_dict(dictionary):
    corrected = {}

    ## Calls HeartPy functions to remove baseline & HRVY peak reductions & inversion
    for key, items in dictionary.items():
        if key == 'Time':
            corrected['Time'] = tuple(dictionary['Time'])
        else:
            sample_rate = round(get_samplerate_datetime(dictionary['Time'], timeformat='%Y-%m-%d %H:%M:%S.%f%z'), 3)
            corrected[key] = tuple(smooth_signal(
                tall_peak_reduct(
                    remove_baseline_wander(
                        flat_peak_reduct(
                            dictionary[key]
                        )  # for flat_peak_reduct
                        , sample_rate)  # for remove_baseline_wander
                )  # for tall_peak_reduct
                , sample_rate, window_length=16, polyorder=3)  # for smooth_signal
            )
            print('node', key, 'corrected')
    return corrected


  
## Reduces long flat peaks
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

        ## If segement of repeated values over 20ms is found, reduce to baseline 0.5
        if i > 0 and rep_val is None and v != 0.5 and v == lst[i - 1]:
            rep_val = v
            rv_start = i - 1
        if rep_val is not None and v != rep_val:
            rv_end = i
            len_rep_val = rv_end - rv_start
            if len_rep_val >= 10 and ss_check1 == 'N':
                lst[rv_start:rv_end] = [0.5] * len_rep_val
            rep_val = None

        ## If unrealistic tall peaks or steep slopes found, reduce
        if i > 0:
            slope = v - lst[i - 1]
            if i_start == None and abs(slope) >= not_flat:
                i_start = i - 1
                v_start = lst[i - 1]
            if abs(slope) >= steep_slope:
                ss_check1 = 'Y'
                ss_1 = slope
                if i_start == None:
                    i_start = i - 1
                    v_start = lst[i - 1]
            if slope != 0 and ss_1 != 0 and ss_check1 == 'Y' and abs(slope) >= not_flat and int(
                    slope / abs(slope)) == int(ss_1 / abs(ss_1) * (-1)):
                ss_check2 = 'Y'

            ## If not a steep slope, reset values for next iteration
            if i_start != None and ss_check1 == 'N' and i - i_start >= 50:
                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'

            ## If it is a steep slope, reduce segment
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
                most_dif_val = max(map(lambda x: abs(x - v_start), vals))
                if most_dif_val >= scale_max:
                    for gap in range(len_ss):
                        fill.append(val)
                        val += step
                    lst[i_start:i_end] = fill

                ## Clears values for next iteration
                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'
                ss_1 = 0

    return lst


  
## Removes unrealistic point and checks for data inversion
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

    ## Placeholder for the starting index of sample selection
    samp_start = None
    ## Placeholder for the ending index of sample selection
    samp_end = None
    ## List of guesses for each selection, whether selection is thought to be inverted, normal, or unknown
    inv_guesses = []

    for i, v in enumerate(lst):
        ## Start of the segment index & value taken,if unrealistic tall peaks created in flat_peak_reduct, reduce
        if i > 0:
            slope = v - lst[i - 1]
            if i_start == None and abs(slope) >= not_flat:
                i_start = i - 1
                v_start = lst[i - 1]
            if abs(slope) >= steep_slope:
                ss_check1 = 'Y'
                ss_1 = slope
                if i_start == None:
                    i_start = i - 1
                    v_start = lst[i - 1]
            if slope != 0 and ss_1 != 0 and ss_check1 == 'Y' and abs(slope) >= not_flat and int(
                    slope / abs(slope)) == int(ss_1 / abs(ss_1) * (-1)):
                ss_check2 = 'Y'

            ## If not a steep slope, reset values for next iteration
            if i_start != None and ss_check1 == 'N' and i - i_start >= 50:
                i_start = None
                v_start = None
                i_end = None
                v_end = None
                ss_check1 = 'N'
                ss_check2 = 'N'

            ## If it is a steep slope, reduce segment
            if i_start != None and ss_check2 == 'Y' and abs(slope) <= not_flat:
                i_end = i
                v_end = v
                len_ss = i_end - i_start
                vals = lst[i_start:i_end]
                fill = []
                v_dif = v_end - v_start
                step = v_dif / len_ss
                val = v_start
                most_dif_val = max(map(lambda x: abs(x - v_start), vals))
                if most_dif_val >= scale_max:
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
                ss_1 = 0

            ##Detection of inverted signal by taking samples (1000 datapoints) of normal sections of data
            ## After slope is not flat, index & value taken for start of the sample selection
            if samp_start == None and slope <= not_flat:
                samp_start = i - 1
                ss_v = v

            ## If sample selection has started and either a steep slope or tall peak is detected, restart selection
            if samp_start != None and (slope >= steep_slope or abs(
                    v - ss_v) > scale_max):
                samp_start = None

            ## Once clean selection of 1000 data points is found, take last index of sample selection
            ## Create a list of values from start to end of sample
            if samp_start != None and i - samp_start == 1000:
                samp_end = i
                samp_vals = lst[samp_start:samp_end]

                ## Max, min, & midpoint values of sample
                max_sv = max(samp_vals)
                min_sv = min(samp_vals)
                mid_sv = (max_sv + min_sv) / 2

                ## Find the percentage of values within the sample that fit the top & bottom half of the range
                top = sum(1 for x in samp_vals if mid_sv <= x <= max_sv) / len(
                    samp_vals) * 100
                bottom = sum(1 for x in samp_vals if min_sv <= x <= mid_sv) / len(
                    samp_vals) * 100

                ## If 75% or more of the data points are in the top half of the sample range, label guess- 'inverted'
                ## If data points are in bottom half of the range, label guess- 'normal'
                ## Everything else is labeled- 'unknown', then restart sample selection
                if top >= 75:
                    inv_guesses.append('inv')
                elif bottom >= 75:
                    inv_guesses.append('norm')
                else:
                    inv_guesses.append('unk')
                samp_start = None

    ## For the list of guesses, find highest occurrence & probability
    if len(inv_guesses) > 0:
        guess = max(inv_guesses, key=inv_guesses.count)
        guess_perc = round(inv_guesses.count(guess) / len(inv_guesses) * 100, 2)

        ## If best guess is 'inverted' and probability is higher than 60%, call negdata_flip to invert signal
        if guess == 'inv' and guess_perc > 60:
            print(' Inverted signal detected')
            print(' Un-inverting data')
            lst[:] = negdata_flip(lst)

    return lst

  

## Inverts any raw negative mV peaks data to positive, normal ECG
def negdata_flip(data_section):
    ## HeartPy Defaults: enhance_peaks = F, keep_range = T
    enhance_peaks = False
    keep_range = True
    out_array = flip_signal(data_section, enhance_peaks, keep_range)

    return out_array


  
## Provides user option to save or print processed data locally
def print_dict(dictionary, outfile):
    ## Write the header
    header = []
    for key in dictionary:
        header.append(key)
    outfile.write(", ".join(header) + '\n')

    ## Write data set values
    print_line = []
    for key, items in dictionary.items():
        for i in range(len(dictionary[key])):
            for k in range(len(dictionary)):
                print_line.append(str(dictionary[list(dictionary)[k]][i]))
            outfile.write(", ".join(print_line) + '\n')
            print_line.clear()


            
## Program begins
if __name__ == '__main__':

    ## Ask user for data set or file name
    dataset = input('HRVY Begin\nEnter the name of the data set or the file name: ')

    ## Call file opener
    contents = file_opener(dataset)

    if contents != 'end':
        ## If file present, then unpack tuple
        (hdr, csv) = contents

        ## Call hdr_data- dictionary builder & node info
        hdr_info = hdr_data(hdr)
        (hdr_dat, nodes_not_2ms) = hdr_info

        ## Call dat_file generation
        dat_align = align_plus_dat(hdr_dat, csv)
        (align_cont, file_dat) = dat_align

        ## Call alignment checker & outfile writer
        final_dat = time_out(csv, file_dat, nodes_not_2ms, align_cont)

        ## Call printing options
        plot_out(final_dat, dataset)

    else:
        print("Please try again")

    print('HRVY complete')
