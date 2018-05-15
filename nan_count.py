files = []
inp = input('Input file names: ')
while inp != '':
    file_name, head = inp.split()
    files.append([file_name, head])
    inp = input()

print('File name\tPercentage of NaNs')
for unit in files:
    file_name, head = unit[0], unit[1]
    with open(file_name, 'r') as file:
        count_line = 0
        count_nans = 0
        
        if head == 'y':
            header = file.readline()
        
        for line in file:
            if 'nan' in line or 'NaN' in line:
                count_nans += 1
            count_line += 1
    
        ratio = count_nans / count_line
        print(file_name + '\t{:.2%}'.format(ratio))
