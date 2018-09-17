

def Divide_TkList(fmax, rez, slices_q, client, rank, begin=0,): # BLOGAI VEIKIA
    intervalsList=[] # sorted by all program proceses (with SuperPC up to 48...)
    #if (fmax != 1): rez += 1 # !@@@ Jei fmax == 10**0, ty bus fj==1, iskreips grafika
    print ('fmax', fmax,'rez', rez,'slices_q', slices_q,'client', client,'begin', begin)

    sveikojiDalis =  int(rez / slices_q)
    liekana = rez % slices_q
    if (client == 'workerCPU' and rank != 1): end = begin + sveikojiDalis-1
    else: end = sveikojiDalis-1
    for processNumber in range (1, slices_q+1): # Rezoliucijos intervalas padalinamas i intervaliukus procesams
        if (liekana > 0): 
            end += 1
            liekana -= 1
        #if (client == 'workerCPU'): intervalsList.insert(processNumber*self.rank, [begin, end])
        intervalsList.insert(processNumber-1, [begin, end])  #tk_IntervalsList.append( tkList[begin:end+1] ) 
        interval = [begin, end] # tik printinimui
        if (rank !=1): TheNumber = processNumber+(slices_q*(rank-1))
        else: TheNumber = processNumber
        begin = 1 + end # koreguojam nes tai index, kuris rodys i list[i] nuo 0 ir imtinai
        end = end + sveikojiDalis 
        if (client == 'mainCPU'):     print ('CPU:',rank,' Process nr:', processNumber, 'sk. intervalas:', interval)
        elif (client == 'workerCPU'): print ('CPU:',rank,' Process nr:',TheNumber, 'sk. intervalas:', interval)
        elif (client == 'Thread'):    print(' CPU:',rank,'Thread nr:', TheNumber, 'skaiciavimu intervalas:', interval) 
    return intervalsList # grazinamas list[[]], vidiniu listu - kiek procesu uzsakyta metodui

#fake data:
fmax = 10**-1
rez = 300
slices_q = 3
client = 'mainCPU'
begin=0
rank=0
countCPU = 4
    
interval = Divide_TkList(fmax, rez, slices_q, client, rank, begin)
for CPU in range (1,countCPU):
    client = 'workerCPU'
    slices_q = 12
    inner_rez = int(rez / (countCPU-1))
    rank = CPU
    begin = interval[CPU-1][0]
    Divide_TkList(fmax, inner_rez, slices_q, client, rank, begin)
