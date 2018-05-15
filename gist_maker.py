import pandas as pd
import re
import collections
import matplotlib.pyplot as plt

def inter(list1, list2):
    if len(list1) == 2 and len(list2) == 2:
        if set(list(range(list1[0], list1[1]))).intersection(list(range(list2[0], list2[1]))):
            templist = list1 + list2
            return max(templist) - min(templist) + 1
        else:
            if list1[0] != list1[1] and list2[0] != list2[1]:
                return list1[1] - list1[0] + list2[1] - list2[0] + 2
            else:
                return list1[1] - list1[0] + list2[1] - list2[0] + 1
    if list1 == []:
        return list2[1] - list2[0] + 1
    if list2 == []:
        return list1[1] - list1[0] + 1


 no_transpos.columns = ['N', 'name', 'start', 'stop']
transpos = pd.read_table('need_for_transfer/transposon_chrY')
summ_Y = transpos.append(no_transpos)
summ_Y = summ_Y[['N', 'name', 'start', 'stop', 'intersect']].reset_index()
def lenall(i):
#     len_df = pd.read_csv('read_and_write_chr' + str(i)+'/read_and_write_chr'+str(i), sep="\t")
    len_df = pd.read_csv(i, sep="\t")
    len_df.columns = ['nom', 'ids', 'start', 'stop']
    len_s = len_df['stop'] - len_df['start']
    return len_s.sum(axis = 0)

chrY = lenall(r'read_and_write_chrY')
chr19 = lenall(r'read_and_write_chr19')
lenchr = chr19 + chrY

def totalnans(i):
    nanslen = 0
#     with open('total_nans_chr' + str(i)) as f:
    with open(r'total_nans_chr' + str(i)) as f:
        data = f.readlines()
    for i in data:
        nanlen = int(re.search(r'is\s\d+', i).group()[3:])
        nanslen += nanlen
    return nanslen

chrY = totalnans('Y')
chr19 = totalnans(19)
lennan = chr19 + chrY
total_nan_percent = lennan/lenchr

dic = collections.defaultdict(list)
with open(r'group_related_taxonomy', 'r') as file:
    for i in file:
        if re.search(':', i):
            temp = re.sub(r'\W', '', i)
        elif re.search(r'\S',i):
            dic[temp] += [re.sub(r'  +', '', i[:-1])]
            
def animal_group(ids, dic = dic):
    for key, value in dic.items():
        for i in value:
            if i == ids:
                return key 
            
            
temp_dic = collections.defaultdict(list)      
with open(r'chromodict', 'r') as file:
    for i in file:
        if not re.match('[=|\n]', i):
            a = re.sub('  .*', '', re.sub(r'  [^/]*/', '!@!', i))[:-1].split('!@!')
            temp_dic[a[0]] = a[1]
for key, value in dic.items():
    a = [temp_dic[i] for i in value]
    dic[key] = a

file_input = summ_Y
ag = file_input.groupby('name')
file_input.head()
totallen = 0
for i in file_input.index:
    totallen += file_input.loc[i]['stop'] - file_input.loc[i]['start']

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

def inter(list1, list2):
    if len(list1) == 2 and len(list2) == 2:
        if set(list(range(list1[0], list1[1]))).intersection(list(range(list2[0], list2[1]))):
            templist = list1 + list2
            return max(templist) - min(templist) + 1
        else:
            if list1[0] != list1[1] and list2[0] != list2[1]:
                return list1[1] - list1[0] + list2[1] - list2[0] + 2
            else:
                return list1[1] - list1[0] + list2[1] - list2[0] + 1
    if list1 == []:
        return list2[1] - list2[0] + 1
    if list2 == []:
        return list1[1] - list1[0] + 1
    
file_input = summ_Y
for j in file_input.index:
    if type(file_input.loc[j][4]) != float:
        sim_tr, tr = 0, 0
        trs, sim_trs = [], []
        lenght = file_input.loc[j]['stop'] - file_input.loc[j]['start'] +1
        for i in re.sub(r'[\[\]\',]', '', re.sub(r'\], \[', '|', file_input.loc[j]['intersect'])).split('|'):
            if i.split(' ')[2] != 'trf':
                tr += int(i.split(' ')[1]) - int(i.split(' ')[0])
                trs = [int(i.split(' ')[0]), int(i.split(' ')[1])]
            else:
                sim_tr += int(i.split(' ')[1]) - int(i.split(' ')[0])
                sim_trs = [int(i.split(' ')[0]), int(i.split(' ')[1])]
        if str(file_input['name'][j]) != 'nan' and str(file_input['name'][j]) != 'NaN':
            ids = re.match(r'[^\.]+', file_input['name'][j]).group()
            for key, value in dic.items():
                if ids in value:
                    ids = key
            summary[ids]['repeats'] += tr
            summary[ids]['other_repeats'] += sim_tr
            summary[ids]['nothing'] += lenght - inter(trs, sim_trs)
            if lenght - inter(trs, sim_trs) < 0:
                print("AAAA", list(file_input.loc[j]), trs, sim_trs, lenght, lenght - inter(trs, sim_trs))
    else:
        summary[ids]['nothing'] += file_input.loc[j]['stop'] - file_input.loc[j]['start']

for i in summary:
    if i == 'Human':
        dfs = pd.DataFrame({i : summary[i]}).transpose()
    else:
        dfs = dfs.append(pd.DataFrame({i : summary[i]}).transpose())
dfs = dfs.apply(lambda x: round(x/lenchr,4))
fig, ax = plt.subplots()
fig.set_size_inches(6,6)

dfs.plot.bar(stacked=True, ax=ax);

ax.set_title("Nucleotides distribution")
ax.legend(loc='upper left')
plt.show()
print(total_nan_percent, 'out of all nucleotides belongs to NaNs')
