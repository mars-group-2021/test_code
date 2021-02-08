# MARS Reverse-Engineering
# University of Maryland GLobal Campus
# Authors: Alex Mancera, Stephen Panossian (avg2)

import datetime as dt
import os
from decimal import Decimal


exists = 0
cur_dir= os.listdir(os.getcwd())
while exists == 0:
    patient=input('Enter patient name (for now a,b, or c): ')
    if (patient+'.csv' and patient+'.hdr') in cur_dir:
        print('Both csv and hdr files found')
        exists=1
    else:
        print('not found')

# determines if csv is pre-aligned from the beginning
aligned='' #straightforward
with open(patient+'.hdr','r') as hdr, open(patient+'.csv','r') as csv:
    num_nodes = len(hdr.readlines())
    num_cols_start = len(csv.readline().split(','))
    num_cols_mis = num_nodes + 1 - num_cols_start
    if  num_cols_mis == 0:
        print('csv file pre-aligned')
        aligned='yes'
    else:
        print('csv file not pre-aligned')
        aligned='no'

# parces hdr data and saves in a nested dictionary
hdr_dat={}
node=1
nodes_4ms=[]
with open(patient+'.hdr','r') as hdr:
    for line in hdr:
        l_split=line.strip('{').strip('\n').strip('}').split(', ')
        hdr_dat[node]={}
        for item in l_split:
            hdr_dat[node][item.split()[0].strip(':')]=item.split()[1]
            if item.split()[0].strip(':') == 'period' and item.split()[1] == '4ms':
                nodes_4ms.append(node)
        node+=1
        
print('Processing output file')

with open(patient+'.csv','r') as csv, open(patient+'_comb_dat.csv','w') as output:

    # writes first line (column names) to out file
    output.write('Time,')
    cols = ''
    for key in hdr_dat:
        cols+=hdr_dat[key]['label']+','
    output.write(cols[:-1]+'\n')

    
    last_line_time=''
    two_ms = dt.timedelta(milliseconds=2)
    temp_line_count = 0
    prev_line = []
    temp_line = []
    for line in csv:
        if not line.startswith(','): # if timestamp already given
            curr_line_time = dt.datetime.strptime(line.split(', ')[0], '%Y-%m-%d %H:%M:%S.%f %z')
            if last_line_time != '' and curr_line_time != last_line_time + two_ms: # checks for missing time
                print('Time gap at',line.split(', ')[0])
            last_line_time = curr_line_time
        else:
            last_line_time += two_ms
            
        print_line = []
        for n in range(1,num_cols_start):
            print_line.append(str(line.split(', ')[n].strip() or 0)) # fills in blanks with 0 where not data was recorded
        if aligned == 'no':
            for n in range(num_cols_mis):
                if len(line.split(', ')) < num_nodes +1: # fills in empty columns
                    print_line.append('0') # fills in blanks with 0 where columns are empty
                else:
                    print_line.append(str(line.split(', ')[num_cols_start+n].strip() or 0))
        out_date = dt.datetime.strftime(last_line_time,'%Y-%m-%d %H:%M:%S.%f %z')
        print_line.insert(0,out_date[:-2].replace('000 ',' ')+':'+out_date[-2:])
        
        if len(nodes_4ms) > 0:
            if not line.strip():
                print('here')
            if temp_line_count == 1:
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
                if print_line[node] == '0':
                    temp_line = print_line
                    temp_line_count +=1
                    break
                else:
                    prev_line = print_line
                    
        else:
            output.write(', '.join(print_line)+'\n') #prints out lines
    output.write(', '.join(print_line)+'\n')

print('done')
