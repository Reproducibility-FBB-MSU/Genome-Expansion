# kinda parser 4 our needs, make dataframes of regions

from Bio import AlignIO
import pandas as pd
import re

a = AlignIO.parse(str(" "), "maf") # MAF-file adress
final_list = []
for i in a:
    result = pd.DataFrame()
    for j in range(len(i[0].seq)):
        result[j] = list(i[:,j])
        result[j][1:] = result[j][1:].map(lambda x:'-' if i[:,j][0] == '-' else(0 if x == '-' else 1)) 
    idlist = [[i[k].id for k in range(len(result[0]))], [i[0].annotations["start"] for k in range(len(result[0]))]] 

    result.index = idlist
    final_list.append(result)
	
# technical piece of code, for check, DO NOT look
# pd.set_option('display.max_columns', 500)
# for i,j in enumerate(test_list):
#     print('___________________________________________',i+1,'___________________________________________')
#     print(j)


test_list = final_list
damn = []
print(len(test_list))			# technical, 4 progress
for df in range(len(test_list)):             # go through the entire list of dataframes
	print(df)					# technical, 4 progress
    reg = {'length': 0, 'start': 0, 'name':'nan', 'pattern':[]} # in reg we write current claster/region
    for col in test_list[df]:                # take one dataframe and go through its columns
        switch, first = 0, 0 
        if len(test_list[df]) == 1:         # check 4 human-only (single line) alignment
            damn.append([df + 1, (test_list[df].index[0][0], test_list[df].index[0][1]), 0, test_list[df].columns.values[-1]])
            break                          
        for i,j in enumerate(test_list[df][col]):      # go through its one specific column
            if j == '-':                             # if it is '-'column close cluster/end region
                if reg['length'] > 1:            # if its length was 2 columns or more - write in the output file
                    damn.append([df + 1, reg['name'], reg['start'], col - 1])
                reg = {'length': 0, 'start': 0, 'name':'nan', 'pattern':[]}
                break
            if i > 1 and j != test_list[df][col][i-1]: # check quantity of vector-changes, minding the number of first change ('0' to '1' 4 example, if there was only 0s in the column and than we meet 1
                if switch == 0:
                    first = test_list[df].index[i-1]
                switch +=1
        if test_list[df][col].tolist()[1:] == reg['pattern']: # if the region continues(we dont close it) - add 1 to his length
            reg['length'] +=1
            # write down in the output vector, minding its type (nan or not) 
        elif test_list[df][col].tolist()[1:] != reg['pattern'] and switch > 1: # if switch > 1, we have no age 4 this one  
            if reg['length'] > 1:
                damn.append([df + 1, reg['name'], reg['start'], col - 1])
            reg = {'length': 1, 'start': col, 'name': ['NaN',test_list[df].index[1][1]] , 'pattern':test_list[df][col].tolist()[1:]}
        elif test_list[df][col].tolist()[1:] != reg['pattern'] and switch == 1: # if switch == 1, we have vector with age 
            if reg['length'] > 1:
                damn.append([df + 1, reg['name'], reg['start'], col - 1])
            reg = {'length': 1, 'start': col, 'name': first, 'pattern':test_list[df][col].tolist()[1:]}
        elif test_list[df][col].tolist()[1:] != reg['pattern'] and switch == 0: # if switch == 0, we have '1'-vector
            if reg['length'] > 1:
                damn.append([df + 1, reg['name'], reg['start'], col - 1])
            reg = {'length': 1, 'start': col, 'name': test_list[df].index[len(test_list[df][col].tolist())-1], 'pattern':test_list[df][col].tolist()[1:]}
    if reg['length'] > 1: # check 4 end-region
        damn.append([df + 1, reg['name'], reg['start'], col]) 
# print(damn)

# since there are alignments in the files without real alignment in that (human only), you should delete them (they interfere with the normal operation of the code)
for i,j in enumerate(damn):
    if j[1] == 'nan':
        damn.pop(i)

# write down in the file-output
nanlen = 0
with open('project_output', 'w') as file: # adress of output MIND IT we ned THE SAME in "transposon_overlap" 4 input
    file.write('№\tname\tstart\tstop\n')
    for i in damn:
        if i[1][0] == 'NaN':
            nanlen += i[3] - i[2]
        file.write(str(str(i[0]) + '\t' + str(i[1][0]) + '\t' + str(i[1][1] + i[2]) + '\t' + str(i[1][1] + i[3]) + '\n' ))
    file.write(str('#total length of the "NaN" regions is ' + str(nanlen) + '\n'))
file.close()




