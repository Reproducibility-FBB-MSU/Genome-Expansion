# here we got THE CODE for finding NOT-transposons and NOT-simple-repeats overlapped regions 

ourlist, output = [], []
start = 0
mainfile = open('/home/pig/Загрузки/temp/!!!!test_output', 'r')    # here put read_and_write-output, probably with name read_and_write_chr_№
sidefile = open('/home/pig/Загрузки/temp/!!test_output', 'r')  # here put transposal-output, probably with name transposal_chr_№
endfile = open('no_tranposons', 'w')     # outputfile
gen = (i for i in mainfile)
for i in sidefile:
    ourlist.append(i.split('\t')[2:4])#[:-1])
for i in gen:
    flag = 0
    for j in ourlist[start:]:
        print([i.split('\t')[2], i.split('\t')[3][:-1]], j)
        if [i.split('\t')[2], i.split('\t')[3][:-1]] == j:
            start += 1
            flag = 1
            break
    if flag == 0:
#         print('A')
        endfile.write(i)
mainfile.close()
sidefile.close()
endfile.close()
