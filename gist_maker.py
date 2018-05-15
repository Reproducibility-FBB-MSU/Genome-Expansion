summ_Y = transpos.append(no_transpos)
summ_Y = summ_Y[['N', 'name', 'start', 'stop', 'intersect']].reset_index()
summ_Y

def lenall(i):
    len_df = pd.read_csv('read_and_write_chr' + str(i)+'/read_and_write_chr'+str(i), sep="\t")
    len_df.columns = ['nom', 'ids', 'start', 'stop']
    len_s = len_df['stop'] - len_df['start']
    return len_s.sum(axis = 0)

chr22 = lenall(22)
chr21 = lenall(21)
chrY = lenall('Y')
chr19 = lenall(19)
lenchr = chr19+chr21+chr22+chrY
lenchr

def totalnans(i):
    nanslen = 0
    with open('total_nans_chr' + str(i)) as f:
        data = f.readlines()
    for i in data:
        nanlen = int(re.search(r'is\s\d+', i).group()[3:])
        nanslen += nanlen
    return nanslen

chr22 = totalnans(22)
chr21 = totalnans(21)
chrY = totalnans('Y')
chr19 = totalnans(19)
lennan = chr19+chr21+chr22+chrY
lennan

import pandas as pd
import re
import collections
import matplotlib.pyplot as plt

dic = collections.defaultdict(list)
with open('/home/pig/Рабочий стол/I_HATE_THIS_WORLD', 'r') as file:
    for i in file:
        if re.search(':', i):
            temp = re.sub(r'\W', '', i)
        elif re.search(r'\S',i):
            dic[temp] += [re.sub(r'  +', '', i[:-1])]
			
temp_dic = collections.defaultdict(list)
with open('/home/pig/Рабочий стол/chromodict', 'r') as file:
    for i in file:
        if not re.match('[=|\n]', i):
            a = re.sub('  .*', '', re.sub(r'  [^/]*/', '!@!', i))[:-1].split('!@!')
            temp_dic[a[0]] = a[1]
for key, value in dic.items():
    a = [temp_dic[i] for i in value]
    dic[key] = a
	
def animal_group(ids, dic = dic):
    for key, value in dic.items():
        for i in value:
            if i == ids:
                return key 

# file_input = pd.read_table('/home/pig/Рабочий стол/test_inut') # here's outpute file address 
file_input = summ_Y
ag = file_input.groupby('name')
file_input.head()
totallen = 0
for i in file_input.index:
    totallen += file_input.loc[i]['stop'] - file_input.loc[i]['start']
# totallen


# sum te length + sum species length
summ = 0
transposons = collections.defaultdict(int)
try:
    for i in file_input['intersect'].map(lambda x:re.sub(r'[\[\]\',]', '', re.sub(r'\], \[', '|', x)).split('|')):
        for j in i:
            transposons[j.split(' ')[2]] = int(transposons[j.split(' ')[2]]) + int(j.split(' ')[1]) - int(j.split(' ')[0])
    for i in transposons.values():summ += i
except TypeError:
    pass

# sum amount of alignmets for each age
total_len = {}
temp = []
for i in ag:
    for j in i:
        if type(j) == str:
            temp = [j]
        else:
            temp.append(j.shape[0])
            total_len[animal_group(re.search(r'[^\.]+', temp[0]).group())] = temp[1]
            string = []
total_len

summary = dic.copy() #creating 3 histogram subgroups for each column
for i in summary.keys():
    summary[i] = {'repeats': 0, 'other_repeats': 0, 'nothing': 0}
summary

nandict = {'repeats': 0, 'other_repeats': 0, 'nothing': 0}
for j in file_input.index:
    if type(file_input.loc[j][4]) != float:
        sim_tr, tr = 0, 0
        lenght = file_input.loc[j]['stop'] - file_input.loc[j]['start']
        for i in re.sub(r'[\[\]\',]', '', re.sub(r'\], \[', '|', file_input.loc[j]['intersect'])).split('|'):
            print(i)
            if i.split(' ')[2]:
                tr += int(i.split(' ')[1]) - int(i.split(' ')[0])
            else:
                sim_tr += int(i.split(' ')[1]) - int(i.split(' ')[0])
        if str(file_input['name'][j]) != 'nan' and str(file_input['name'][j]) != 'NaN':
            ids = re.match(r'[^\.]+', file_input['name'][j]).group()
            for key, value in dic.items():
                if ids in value:
                    ids = key
            summary[ids]['repeats'] += tr
            summary[ids]['other_repeats'] += sim_tr
            summary[ids]['nothing'] += lenght - tr - sim_tr
    else:
        summary[ids]['nothing'] += file_input.loc[j]['stop'] - file_input.loc[j]['start']

for i in summary:
    if i == 'Human':
        dfs = pd.DataFrame({i : summary[i]}).transpose()
    else:
        dfs = dfs.append(pd.DataFrame({i : summary[i]}).transpose())
dfs = dfs.apply(lambda x: round(x/lenall,4))
dfs = dfs.drop(['nothing'])
dfs
 

fig, ax = plt.subplots()
fig.set_size_inches(6,6)

dfs.plot.bar(stacked=True, ax=ax);

ax.set_title("Nucleotides distribution")
ax.legend(loc='upper left')
plt.show()