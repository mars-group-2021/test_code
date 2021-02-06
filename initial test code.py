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

    # values 2 lines and 1 line behind current value, respectively
    val2back = 0.0
    val1back = 0.0
    cur_val  = 0.0
    count    = 0
    lnum     = 0
    
    last_line_time=''
    two_ms = dt.timedelta(milliseconds=2)
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
            #time stamp exists, but check for value in previous line
            if ((count >= 3) and (abs(val1back) < 0.00000001)):
                # average/linear interpolate for 0 value
                val1back = avg2(val2back,cur_val)

            output.seek(lnum-1)
            output.write(last_line_time+str(val1back)+'\n')

            val2back = val1back
            val1back = cur_val

            # proceed with current line
            last_line_time += two_ms
            rest_line = ''
            for n in range(1,num_cols_start):
                rest_line += ', '+str(line.split(', ')[n].strip() or 0) # fills in blanks with 0 where not data was recorded
            if aligned == 'no' and len(line.split(', ')) < num_nodes +1: # fills in empty columns
                for n in range(num_cols_mis):
                    rest_line += ', 0' # fills in blanks with 0 where columns are empty
                    
            # thank you Stack Overflow:
            # https://stackoverflow.com/questions/30112357/
            #    typeerror-descriptor-strftime-requires-a-datetime-date-object-but-received
            print(last_line_time, val2back, val1back, cur_val)
            #out_date_obj = dt.datetime.strptime(last_line_time, '%Y-%m-%d %H:%M:%S.%f %z')
            #out_date_str = dt.datetime.strftime(out_date_obj,   '%Y-%m-%d %H:%M:%S.%f %z')
            #output.write(out_date_str[:-2].replace('000 ',' ')+':'+out_date_str[-2:]+rest_line+'\n') #prints out lines
            #out_date = dt.datetime.strptime(last_line_time,'%Y-%m-%d %H:%M:%S.%f')

            out_date = dt.datetime.strftime(last_line_time,'%Y-%m-%d %H:%M:%S.%f %z')
            output.write(out_date[:-2].replace('000 ',' ')+':'+out_date[-2:]+rest_line+'\n') #prints out lines

            
print('done')

