# MARS Reverse-Engineering
# University of Maryland GLobal Campus
# Authors: Alex Mancera, Stephen Panossian (avg2)

import datetime as dt
import os

def avg2(val2back,current_val): 
    return (val2back + current_val) / 2.0

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

# parses hdr data and saves in a nested dictionary
hdr_dat={}
node=1
with open(patient+'.hdr','r') as hdr:
    for line in hdr:
        l_split=line.strip('{').strip('\n').strip('}').split(', ')
        hdr_dat[node]={}
        for item in l_split:
            hdr_dat[node][item.split()[0].strip(':')]=item.split()[1]
        node+=1
        


print('Processing output file')           
with open(patient+'.csv','r') as csv, open(patient+'_comb_dat.csv','w') as output:

    # writes first line (column names) to out file
    output.write('Time,')
    cols = ''
    for key in hdr_dat:
        cols+=hdr_dat[key]['label']+','
    output.write(cols[:-1]+'\n')

    # values 2 lines and 1 line behind current line values, respectively
    vals2back = []
    vals1back = []
    cur_vals  = []
    
    count = 0 # loop counter
    lnum  = 0 # integer line number
    
    last_line_time=''
    two_ms = dt.timedelta(milliseconds=2)
    
    # process entire csv file line by line
    for line in csv:
        # find line number in output file
        lnum = output.tell()
        count += 1
        
        # checks for timestamp, send user alert if missing
        if not line.startswith(','): # if timestamp not already given
            curr_line_time = dt.datetime.strptime(line.split(', ')[0], '%Y-%m-%d %H:%M:%S.%f %z')
            print(curr_line_time)
            if last_line_time != '' and curr_line_time != last_line_time + two_ms: # checks for missing time
                print('Time discrepancy at',line.split(', ')[0])
            last_line_time = curr_line_time
            output.write(line)
        else:
            # time stamp exists, proceed with current line and process empty value column(s)
            last_line_time += two_ms # advance timestamp
            rest_line = ''           # rest of line past timestamp, value(s) for 1+ columns

            for n in range(1,num_cols_start):
                # fill in value or 0 for each column past timestamp
                rest_line += ', '+str(line.split(', ')[n].strip() or 0) # fills in blanks with 0 where not data was recorded

            if aligned == 'no' and len(line.split(', ')) < num_nodes+1: # fills in empty columns
                for n in range(num_cols_mis):
                    rest_line += ', 0' # fills in blanks with 0 where columns are empty

            # Check for value=0 in previous line
            elif (count >= 3):
                # Check and replace 0 values if any across line
                for n in range(1,num_cols_start-1):
                    cur_vals += str(line.split(', ')[n].strip())
                    if (abs(vals1back[n]) < 0.00000001):
                        # average/linear interpolate for 0 value
                        vals1back[n] = avg2(vals2back[n],cur_vals[n])
                        output.seek(lnum-1) # overwrite last line's 0 value
                    # end if average calls
                    output.write(str(last_line_time)+str(vals1back[n])+'\n')
                # end for search for 0 values through every column
            elif (count == 1):
                vals2back = rest_line
                vals1back = rest_line
            elif (count >= 2):
                vals2back = vals1back
                vals1back = cur_vals
            # endif - checks for timestamp and zero value(s)
        
            out_date = dt.datetime.strftime(last_line_time,'%Y-%m-%d %H:%M:%S.%f %z')
            output.write(out_date[:-2].replace('000 ',' ')+':'+out_date[-2:]+rest_line+'\n') #prints out lines
        # end if checks
    # end for "line in csv"
            
print('done')

