# check for overlap with transposons
# we just pass through the entire list and check each element for correspondence with any part of the transposons from the file of our chromosome

proj_out = pd.read_table('project_output', comment = '#')  # here we put the adress (name) of the "read_and_write"-output 
regions = pd.read_table('or')         # here we put the adress (name) of the transposons_of our chr
regions = regions.sort_values(by=['genoStart'])
regions.reset_index(drop = True, inplace = True)
flag = 0
regions

# gen = ([regions['genoStart'][i], regions['genoEnd'][i], i] for i in range(flag, regions.index[-1])) 
# for i in gen:
#     print(i)

def overlap(reg, flag, file = regions):
    output = []
    gen = ([file['genoStart'][i], file['genoEnd'][i], i, file['repClass'][i]] for i in range(flag,  file.index[-1]+1)) 
    for j in gen:
        if reg['start'] > j[1]:
            flag +=1
        if j[0] > reg['stop']:
            if output != []:
                return [output, flag]
            else:
                return
        inter = sorted(set(list(range(int(reg['start']), int(reg['stop'])+1))).intersection(list(range(int(j[0]), int(j[1])+1))))
        if inter:
            output.append([[list(inter)[0], list(inter)[-1]], j[3]])
    if output != []:
        return [output, 0]
        
with open('test_output','w') as out:  # output adress
    out.write('№\tname\tstart\tstop\tintersect\n')
    for i in proj_out.index:
        temp = overlap(proj_out.loc[i], flag)
        if temp != None:    
            flag = temp[1]
            a = (str(proj_out.loc[i]['№']) + '\t' + str(proj_out.loc[i]['name']) + '\t' + 
                    str(proj_out.loc[i]['start']) + '\t' + str(proj_out.loc[i]['stop']) + '\t' + str(temp[0]))
            out.write(a + '\n')
