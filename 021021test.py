# MARS Reverse-Engineering
# University of Maryland Global Campus
# Authors: Alex Mancera, Stephen Panossian, Analia Treviño-Flitton

import datetime as dt
import os
from decimal import Decimal







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
            hdr = open(patient + '.hdr', 'r').readlines()
            csv = open(patient + '.csv', 'r').readlines()
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
num_nodes = len(hdr) # number of line in .hdr file
num_cols_start = len((csv[0]).split(',')) # number of columns at the beginning of .csv file
num_cols_mis = num_nodes + 1 - num_cols_start # missing number of columns at the beginning of .csv file
if  num_cols_mis == 0:  # if no columns missing 
    print('csv file pre-aligned')
    aligned='yes'
else:
    print('csv file not pre-aligned')
    aligned='no'


print('Processing output file')
with open(patient+'_comb_dat.csv','w') as output:
    
    # writes header line (column names) to out file
    output.write('Time,')
    cols = ''
    for key in hdr_dat:
        cols+=hdr_dat[key]['label']+','
    output.write(cols[:-1]+'\n')
    
    #generates timestamp data
    last_line_time=''
    two_ms = dt.timedelta(milliseconds=2)
    temp_line_count = 0
    prev_line = []
    temp_line = []
    for line in csv:
        print_line = [] 
        if not line.startswith(','): # if timestamp already given
            curr_line_time = dt.datetime.strptime(line.split(', ')[0], '%Y-%m-%d %H:%M:%S.%f %z') # current time is kept as provided
            if last_line_time != '' and curr_line_time != last_line_time + two_ms: # checks for time gap
                print('Time gap at',line.split(', ')[0])
            last_line_time = curr_line_time # saves current time for reference on next line
        else: # if timestamp not already given
            last_line_time += two_ms # adds 2ms to saved reference time
        out_date = dt.datetime.strftime(last_line_time,'%Y-%m-%d %H:%M:%S.%f %z')
        print_line.append(out_date[:-2].replace('000 ',' ')+':'+out_date[-2:]) # saves timestamp data to print_line list
            
        # generates printed lines to outfile after header line
        for n in range(1,num_cols_start):
            print_line.append(str(line.split(', ')[n].strip() or 0)) # fills in blanks with recorded values or 0 where data was not recorded
        if aligned == 'no':
            for n in range(num_cols_mis):
                if len(line.split(', ')) < num_nodes +1: # if empty columns exist
                    print_line.append('0') # fills in blanks with 0 where columns are empty
                else:
                    print_line.append(str(line.split(', ')[num_cols_start+n].strip() or 0)) # fills in blanks with recorded values or 0 where data was not recorded
        
        if len(nodes_4ms) > 0: #if there are any 4ms nodes
            if temp_line_count == 1: # if blank line in 4ms node exists
                for node in nodes_4ms:
                    if len(prev_line) > 0 :
                        avg_val = round(Decimal((float(prev_line[node])+float(print_line[node]))/2),3)
                        temp_line[node] = str(avg_val)
                    else:
                        temp_line[node] = print_line[node]
                    
                if len(prev_line)>0:
                    output.write(', '.join(prev_line)+'\n')
                    prev_line = []
                output.write(', '.join(temp_line)+'\n')
                
                temp_line_count = 0
                
            for node in nodes_4ms:
                if print_line[node] == '0':  # if blank line in 4ms node exists
                    temp_line = print_line
                    temp_line_count +=1
                    break
                else:
                    prev_line = print_line
                    
        else:
            output.write(', '.join(print_line)+'\n') #prints out lines
    output.write(', '.join(print_line)+'\n')

print('done')
