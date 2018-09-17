import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pylab as plt
import math
import random
import time
import threading
import multiprocessing as mltp
from multiprocessing import Process, JoinableQueue
import itertools
import sys
from mpi4py import MPI 

# By Joris Mykolas Medeisis

class Sprendimas:
    tkList = []
    xList = []
    yList = []
    N = 10**5
    lmd = 2.9
    t0 = 0
    tN = None
    rez = 1
    threads_quantity = None
    countCPU = None
    threadsList = []
    threads_intervalList = [[]]
    queuePoints = None
    safeQueue = JoinableQueue()
    fmin = 10**-6 # UZDUOTIES GRAFIKAS
    fmax = 10**0 

    #Konstruktorius
    def __init__(self, rank, fmin, fmax, lamda, N, rezoliucija, threads=None, countCPU=1):
        self.fmin = fmin
        self.fmax = fmax
        self.lmd = lamda
        self.N = N
        self.rez = rezoliucija
        self.countCPU = countCPU
        if (threads != None or threads != 1 or threads != 0):
            self.threads_quantity = threads
            self.queuePoints = JoinableQueue()
        self.FindTn()
        print('Paleidau Sprendimo objk')
        self.FindPointsXY()
             
    def PoissonRandom(self, lmd): # Poisson taskinio proceso pasiskirstymas
        a = np.e**(-1.0 * lmd)
        r = 1
        n = -1
        while r > a:
            u = np.random.random(1)
            r *= u 
            n += 1
        return n
    # Generuojamas Puasso list ir tN
    def FindTn (self): # Pagal: Tk = Tk-1 + P(lmd)k;  P(lmd)-Poiss random nr; T0 = 0
        lmd = self.lmd
        N = self.N
        tk = self.t0
        for k in range(1, N+1):
            #tk += np.random.poisson(lmd) # padaryti rankom
            tk += self.PoissonRandom(lmd)
            self.tkList.append(tk)
            continue
        self.tN = tk # N ir lmd - ivedami; Gaunam: tList ir Tn(tList pask. nr.)
        return

    # Apskaiciuojam funkcija S(f) arba kitaip y pagal duota formule
    def FindSf (self, f, t0, tN, tkList):
        c = 0.0
        s = 0.0
        for k in range(0, len(tkList)): # Sigma susumavimas sin ir cos reiskiniu; N - koks bus tkList size
            c += np.cos(2.0 * np.pi * f * tkList[k])
            s += np.sin(2.0 * np.pi * f * tkList[k])
        c*=c # Keliam kvadratu susumuota reiskini  
        s*=s
        Sf = (2.0/(tN - t0))*(c+s)
        return Sf
        #Apskaiciuojamas daznio f kordinates (x, y)
    def FindPointsXY(self):
        tN = self.tN
        rez = self.rez
        fmin = self.fmin #Tikrieji sk. 10**-6
        fmax = self.fmax
        t0 = self.t0
        tkList = self.tkList
        step = (math.log10(fmax) - math.log10(fmin)) / rez # surandame zingsnio ilgi h, bet skaiciuojam pagal log10 skale
        start_time = time.time()
        # Pasirenkamas budas skaiciavimui atlikti:
        if (self.countCPU <= 1 and self.threads_quantity == None): # Single CPU, one_thread
            print('Patekau i SingleCPU, SingleThread')
            self.xList, self.yList = self.CalculatePointsXY_simple(t0, tN, tkList, rez, fmin, step)

        elif (self.countCPU <= 1 and self.threads_quantity != None): # Single CPU, multi_thread
            print('Patekau i elif SingleCPU, Thread count:', self.threads_quantity)
            worker = Worker(rank)
            if (fmax != 1): rez += 1 # !@@@ Jei fmax == 10**0, ty bus fj==1, iskreips grafika
            tk_Intervals = worker.Divide_TkList( fmax, rez, self.threads_quantity, 'Thread' ) #grazinamas list[[beg, end]...] nuorodos
            self.xList, self.yList = worker.CalcPointsXY_toThread( t0, tN, step, fmin, rez, self.threads_quantity, tkList, tk_Intervals ) #Paziureti kiek arg priima def...

        elif (self.countCPU > 1): # Multi CPU
            print('Patekau i MultiCPU')
            worker = Worker(rank)
            tk_Intervals = worker.Divide_TkList( fmax, rez, self.countCPU-1, 'mainCPU') #  #grazinamas list[[beg, end]...] nuorodos
            print ('Padalinau tkList to CPUs')
            xyList = [[]]
            #Kad butu fmax imtinai, kompnesuojame rez+1 (padaryti paskutini zingsniuka iki fmax)
            if (fmax != 1): rez += 1 # !@@@ Jei fmax == 10**0, ty bus fj==1, iskreips grafika
            liekana = rez % (self.countCPU-1)
            for CPU in range(1, self.countCPU): # Paskirstom rez(zingsnelius po workerCPUs)
                if (liekana > 0): 
                    inner_rez = 1+int(rez / (self.countCPU-1))
                    liekana -= 1
                else: inner_rez = int(rez / (self.countCPU-1)) # kadandi dar dalinsime tkList uzduotis threadams CPU viduje
                print('CPU:', CPU)
                data = t0, tN, inner_rez, fmin, fmax, step, tk_Intervals[CPU-1], tkList #  N galima pakeisti len(tkList)
                comm.send( data, dest=CPU )
                print ('Issiunciau data i CPU:', CPU)

            for CPU in range(1, self.countCPU):
                data = comm.recv( source=CPU ) 
                self.xList.append(list(data[0]))
                self.yList.append(list(data[1]))
                print ('Gavau data is CPU:', CPU, 'xList Size:', len(list(data[0])), 'yList Size:', len(list(data[1])) )

            self.xList = list(itertools.chain.from_iterable(self.xList)) # is list[[]] i list[]
            self.yList = list(itertools.chain.from_iterable(self.yList))
            print('xList size:',len(self.xList), 'yList size:', len(self.yList))     
        print("Skaiciavimu trukme: %s s" % (time.time() - start_time))
        return

    def CalculatePointsXY_simple(self, t0, tN, tkList, rez, fmin, step):
        if (self.fmax == 1): rez -= 1 # !@@@ fmax == 10**0, ty fj==fmax, iskreips grafika
        for j in range(0, rez+1): # k - zinsniu kiekis, kiek h reikia padaryti fmin -> fmax; 
            fj = fmin*10**(j*step) # formule: gaunam f pagal duota rezoliucija, kitaip x koord
            Sfj = self.FindSf(fj, t0, tN, tkList) # S(f) - y koord
            self.xList.append(fj)
            self.yList.append(Sfj)
        return

class Atlikimas:
    minX = 0 
    maxX = 0
    minY = 0 
    maxY = 0
        #Konstruktorius
    def __init__(self, minX, maxX, xList, yList, savefile):
        self.minX = minX
        self.maxX = maxX
        self.minY, self.maxY = self.Find_minY_maxY(yList)
        self.DrawGraphic(self.minX, self.maxX, self.minY, self.maxY, xList, yList, savefile)

    def Find_minY_maxY(self, yList):
        minY = yList[0]
        maxY = yList[0]
        for y in range(len(yList)):
            if (yList[y] < minY): minY = yList[y]
            if (yList[y] > maxY): maxY = yList[y]
        return minY, maxY

    def FindAlfa(self, minY, maxY, xList, yList):
        #atkarpa tarp 25% ir 75%(nuo didziausios y reiksmes) Y asies(S(f)) reiksmiu
        sliceListX=[]
        sliceListY=[]
        delta = np.log10(maxY) - np.log10(minY)
        minY_lim = delta*0.25 + np.log10(minY)
        maxY_lim = delta*0.75 + np.log10(minY) # Cia abu turi buti +minY arba abu -maxY
        for y in range(len(yList)):
            if(maxY_lim > np.log10(yList[y]) and np.log10(yList[y]) > minY_lim): #Velgi tikrinam log-log mastelyje
                sliceListX.append(xList[y])
                sliceListY.append(yList[y])
        logX = np.log10(sliceListX) #Pasiverciam i log skale, kad butu tiesinis alfa
        logY = np.log10(sliceListY)
        alfa, const = np.polyfit(logX, logY, 1) # fit log(S(f)) = alfa*log(f) + const
        FitSf = 10**(alfa * np.log10(sliceListX) + const) # apskaiciuojam Fit y ir ATSIVERCIAM is log skales atgal i norml
        return sliceListX, FitSf, sliceListY

    def Format_GraphLim(self, minX, maxX, minY, maxY):
        minX = 10**(np.log10(minX) - 0.5)
        maxX = 10**(np.log10(maxX) + 0.5)
        minY = 10**(np.log10(minY) - 0.5)
        maxY = 10**(np.log10(maxY) + 0.5)
        return minX, maxX, minY, maxY

    def DrawGraphic(self, minX, maxX, minY, maxY, xList, yList, savefile):
        #MatPlotLib grafiko nustatymai:
        plt.loglog (basex=10, basey=10) #(baseX, baseY) - automatiskai visas paduodamas reiksmes konvertuoja i log10xlog
        g_minX, g_maxX, g_minY, g_maxY = self.Format_GraphLim(minX, maxX, minY, maxY)
        plt.axis([g_minX, g_maxX, g_minY, g_maxY]) #[Xmin, Xmax], [Ymin, Ymax]
        plt.grid (True)
        # Taskinio proceso grafikas
        plt.plot (xList, yList, label='NMM1', color='red') # grafikas lauzte
        #plt.scatter (xList, yList, 5, label='NMM1', color='red') # grafikas taskais
        
        slice_f_List, fit_Sf_list, slice_Sf_List = self.FindAlfa(minY, maxY, xList, yList) # Fitinimas pagal maxY
        plt.plot(slice_f_List, fit_Sf_list, color="green") 
        plt.scatter(slice_f_List, slice_Sf_List, 5, label='NMM1', color='blue')
        # Draw or Imaging
        plt.savefig(savefile) #issaugo plot obj faile
        plt.show ()

class Worker():

    rank = None
    queuePoints = None

    def __init__(self, rank=0):
        self.rank = rank
        self.queuePoints = JoinableQueue()

    def Divide_TkList(self, fmax, rez, slices_q, client, begin=0): # nereikia f max
        intervalsList=[] # sorted by all program proceses (with SuperPC up to 48...)
        sveikojiDalis =  int(rez / slices_q)
        liekana = rez % slices_q
        if (client == 'workerCPU' and rank != 1): end = begin + sveikojiDalis-1
        else: end = sveikojiDalis-1
        for processNumber in range (1, slices_q+1): # Rezoliucijos intervalas padalinamas i intervaliukus procesams
            if (liekana > 0):
                end += 1
                liekana -= 1
            intervalsList.insert(processNumber-1, [begin, end])  #tk_IntervalsList.append( tkList[begin:end+1] ) 
            interval = [begin, end] # tik printinimui
            begin = 1 + end # koreguojam nes tai index, kuris rodys i list[i] nuo 0 ir imtinai
            end = end + sveikojiDalis 

            if (rank !=1): TheNumber = processNumber+(slices_q*(self.rank-1))
            else: TheNumber = processNumber
            if (client == 'mainCPU'):     print ('CPU:',self.rank,' Process nr:', processNumber, 'sk. intervalas:', interval)
            elif (client == 'workerCPU'): print ('CPU:',self.rank,' Process nr:', TheNumber, 'sk. intervalas:', interval)
            elif (client == 'Thread'):    print (' CPU:',self.rank,'Thread nr:', TheNumber, 'skaiciavimu intervalas:', interval) 
        return intervalsList # grazinamas list[[]], vidiniu listu - kiek procesu uzsakyta metodui
                
    def CalcPointsXY_toThread(self, t0, tN, step, fmin, rez, threads_q, tk_Intervals, tkList): # Naudoja kiekvienas CPU (isskyrus MAIN, nebent jis vienintelis)
    # Funkcija: paskirstyti darba (tk Intervalus) po Hipergijas (viename CPU)
        threadsList = []
        for threadNr in range (0, threads_q): # Rezoliucijos intervalas padalinamas i intervaliukus procesam
            thread = Process(target = self.ThreadJob, args=(t0, tN, step, fmin, threadNr, tk_Intervals[threadNr], tkList))
            threadsList.append(thread)
            threadsList[ threadNr ].start() # Startuojamas def ThreadJob()
        time.sleep(1) # main thread miega 5sec, procesai dirba
        self.WaitUntil_Join(self.queuePoints, threads_q)
        xList, yList = self.FormatData_toList(self.queuePoints)
        return xList, yList

    def ThreadJob(self, t0, tN, step, fmin, threadNumber, tkInterval, tkList):
        xList_Inner = []
        yList_Inner = []
        xList_Inner.append(threadNumber) # !!!!@@@@ xList_Inner[0] uzkoduojamas threadNr
        for j in range(tkInterval[0], tkInterval[1]+1): # j - zinsniu kiekis, kiek h reikia padaryti fmin -> fmax;
            fj = fmin*10**(j*step) #  formule: gaunam f pagal duota rezoliucija, kitaip x koord
            ###### S(f) apskaiciavimas #######
            c = 0.0
            s = 0.0
            for k in range(len(tkList)): # Sigma susumavimas sin ir cos reiskiniu
                c += np.cos(2.0 * np.pi * fj * tkList[k])
                s += np.sin(2.0 * np.pi * fj * tkList[k])
            c*=c # Keliam kvadratu susumuota reiskini  
            s*=s
            Sfj = (2.0/(tN - t0))*(c+s) # S(f) - y koord
            xList_Inner.append(fj)
            yList_Inner.append(Sfj)
            time.sleep(0.001) # Kad procesoriaus neperkaistu
            continue
        pList=[xList_Inner, yList_Inner]
        self.queuePoints.put(pList) # Kadangi xList[0] ihardkodintas eiles numeris, del eiles tvarkos nesibaiminti
        self.queuePoints.task_done()
        sys.exit() # Paprograme (Process) pats save susinaikina sistemoje "End Task"

    def WaitUntil_Join(self, queuePoints, threads_q):
        working = True
        while (working):
            if (queuePoints.qsize() == threads_q):
                working = False
                break;
            else:
                time.sleep(1)
        return

    def FormatData_toList(self, queuePoints):
        # Is Sync Queue issitraukiam data ir paverciam i obj array
        xArr = np.empty(queuePoints.qsize(), dtype=object)
        yArr = np.empty(queuePoints.qsize(), dtype=object)
        for i in range(queuePoints.qsize()):
            pList = queuePoints.get()
            listX = pList[0]
            listY = pList[1]
            threadNr = listX.pop(0)
            xArr[threadNr] = listX
            yArr[threadNr] = listY
        # Sujungiam lista (is list[[]] i list[])
        xList = list(itertools.chain.from_iterable( list(xArr) ))
        yList = list(itertools.chain.from_iterable( list(yArr) ))
        return xList, yList

#Valdymo skydelis:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
countCPU = comm.Get_size()
host = MPI.Get_processor_name()
if (countCPU > 1): threads_quantity = 12 #Kiek SUPERKOMPAS TURI HIPERGIJU

#Viskas jau ir taip CPU for loop (kazkodel) lyg for( rank in MPI.COMM_WORLD.Get_size() ):
if (rank != 0): # for rank in (1, CPUcount+1): Arba kitaip visiems darbiniams CPU
    data = comm.recv(source=0) # wait to get data from MAIN
    t0, tN, inner_rez, fmin, fmax, step, tk_Interval, tkList = data # prisideda fmax, nes reikia dar karta divideList (palyginus su 1xCPU atveju) 
    worker = Worker(rank)
    tk_Intervals2D = worker.Divide_TkList( fmax, inner_rez, threads_quantity, 'workerCPU', tk_Interval[0] ) # NZN kaip del sito ka siusti kaip skirstyti
    xList, yList = worker.CalcPointsXY_toThread(t0, tN, step, fmin, inner_rez, threads_quantity, tk_Intervals2D, tkList)
    data = [xList, yList]
    comm.send(data, dest=0) # send data to MAIN CPU

if (rank == 0):
    print ('@@@@@@@@@@@@@@@@@ START @@@@@@@@@@@@@@@@@')
    if __name__ == '__main__':
        print ('Main Program:', ' main CPU:', rank, ' CPU kiekis:', countCPU, ' Socket:', host)
        fmin = 10**-7
        fmax = 10**-1
        lmd = 1.1
        N = 10**5
        rez = 301
        threads = None # (Galioja tik su 1xCPU) if None or 0 or 1 - procesas be paprogramiu
        savefile = "NMMgraph.pdf"

        spr = Sprendimas(rank, fmin, fmax, lmd, N, rez, threads, countCPU) # (lmd, N, rez, threads_quntity (none if not multithreaded))
        Atlikimas(fmin, fmax, spr.xList, spr.yList, savefile)
        print('@@@@@@@@@@@@@@@@@  Darbas atliktas sekmingai! Geros dienos! @@@@@@@@@@@@@@@@@')

# By Joris Mykolas Medeisis
'''PASTABOS:
+ Kazkodel FindAlfa aptinkama, jog PointList size virsija: +1
+ Galima perkelti send ir revc i Worker klase
+ Isplesti skydeli: nustatyti grafiko failo pavadinima, kur ji saugoti, nustatyti vienoje vietoje spalvas, grafiko formatavima jau automatizavau
+ ar tN galima pakeisti i tkList[-1] ?
'''