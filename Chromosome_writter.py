# This part of the code is needed if there are no separate files of transposons on the chromosomes, and there is only one large

table = pd.read_table('repeats') # here  put the address of the large file with transposons (RepClass-file)
for i in list(range(1,23)) + ['M', 'Un', 'X', 'Y']:
    with open(str('chrs/chr' + str(i)), 'w') as file: # path to chromosome-directory
        file.write('#bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd	repLeft	id\n')
for i in table.index:
    a = ''
    for j in table.loc[i]:
        j = re.sub(r'_.*', '', str(j))
        a += str(j) + '\t'
    with open(str('chrs/' + re.match(r'chr[^_]{1,2}',table['genoName'][i]).group()), 'a') as file: # here we put path to chromosome directory
        file.write(a[:-1] + '\n')
other_table = pd.read_table('other_repeats') # here put the address of file with Simple Repeats
for i in other_table.index:
    string = str('nan\tnan\tnan\tnan\tnan\t' + str(other_table.loc[i]['genoName']) + '\t' + str(other_table.loc[i]['genoStart']) 
            + '\t' + str(other_table.loc[i]['genoEnd']) + '\tnan\tnan\tnan\t' 
            + str(other_table.loc[i]['repClass']) + '\tnan\tnan\tnan\tnan\tnan' + '\n')
    string = re.sub(r'_\w*\s', '\t', string)
    with open ('chrs/' + re.search(r'chr[^_\s]{1,2}',string).group(), 'a') as f: # here we put path to chromosome directory
        f.write(string)
