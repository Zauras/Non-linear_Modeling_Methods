import numpy as np
import matplotlib.pylab as plt
#matplotlib.use('pdf') bet reikia import matplotlib
import math
import random
import time
import threading
from multiprocessing import Process, JoinableQueue
import itertools
import sys

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
    threadsList = []
    threads_intervalList = [[]]
    queuePoints = None
    safeQueue = JoinableQueue()
    fmin = 10**-6 # UZDUOTIES GRAFIKAS
    fmax = 10**0 

    #Konstruktorius
    def __init__(self, fmin, fmax, lamda, N, rezoliucija, threads=None):
        self.fmin = fmin
        self.fmax = fmax
        self.lmd = lamda
        self.N = N
        self.rez = rezoliucija
        self.FindTn()
        if (threads != None):
            self.threads_quantity = threads
            self.queuePoints = JoinableQueue()
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
            #tk += np.random.poisson(lmd) #rankom
            tk += self.PoissonRandom(lmd)
            self.tkList.append(tk)
            continue
        self.tN = tk # N ir lmd - ivedami; Gaunam: tList ir Tn(tList pask. nr.)
        return

    # Apskaiciuojam funkcija S(f) arba kitaip y pagal duota formule
    def FindSf (self, f, t0, tN, N, tkList):
        c = 0.0
        s = 0.0
        for k in range(0, N): # Sigma susumavimas sin ir cos reiskiniu
            c += np.cos(2.0 * np.pi * f * tkList[k])
            s += np.sin(2.0 * np.pi * f * tkList[k])
        c*=c # Keliam kvadratu susumuota reiskini  
        s*=s
        Sf = (2.0/(tN - t0))*(c+s)
        return Sf
        #Apskaiciuojamas daznio f kordinates (x, y)
    def FindPointsXY(self):
        tN = self.tN
        N = self.N
        rez = self.rez
        fmin = self.fmin #Tikrieji sk. 10**-6
        fmax = self.fmax
        t0 = self.t0
        tkList = self.tkList
        step = (math.log10(fmax) - math.log10(fmin)) / rez # surandame zingsnio ilgi h, bet skaiciuojam pagal log10 skale
        start_time = time.time()
        # Pasirenkamas budas skaiciavimui atlikti:
        if (self.threads_quantity == None):
            self.CalculatePointsXY_simple(t0, tN, tkList, N, rez, fmin, step)
        else: 
            self.CalculatePointsXY_multithreading(t0, tN, tkList, N, rez, fmin, step)
        print("Skaiciavimu trukme: %s s" % (time.time() - start_time))
        return

    def CalculatePointsXY_simple(self, t0, tN, tkList, N, rez, fmin, step):
        if (self.fmax == 1): rez -= 1 # !@@@ fmax == 10**0, ty fj==fmax, iskreips grafika
        for j in range(0, rez+1): # k - zinsniu kiekis, kiek h reikia padaryti fmin -> fmax; 
            fj = fmin*10**(j*step) # formule: gaunam f pagal duota rezoliucija, kitaip x koord
            Sfj = self.FindSf(fj, t0, tN, N, tkList) # S(f) - y koord
            self.xList.append(fj)
            self.yList.append(Sfj)
        return

    def CalculatePointsXY_multithreading(self, t0, tN, tkList, N, rez, fmin, step):
        intervalsList = self.threads_intervalList
        threads_q = self.threads_quantity
        queuePoints = self.queuePoints
        if (self.fmax != 1): rez += 1 # !@@@ Jei fmax == 10**0, ty bus fj==1, iskreips grafika
        sveikojiDalis =  int(rez / threads_q)
        liekana = rez % threads_q
        begin=0
        end = sveikojiDalis
        # Veiksmai su Procesais
        for threadNumber in range (0, threads_q): # Rezoliucijos intervalas padalinamas i intervaliukus procesams
            if (liekana > 0): 
                end += 1
                liekana -= 1
            intervalsList.insert(threadNumber, [begin, end])
            begin = end
            end = end + sveikojiDalis 
            # Proceso kurimas
            thread = Process(target = self.ThreadJob, args=(threadNumber, intervalsList, step, fmin, t0, tN, tkList, N))
            self.threadsList.append(thread)
            self.threadsList[threadNumber].start() # Startuojamas def ThreadJob()
        #Threadu intervaliniu listu sujungimas i baigtinius List X, Y 
        time.sleep(5) # main thread miega 5sec, procesai dirba
        self.WaitUntil_ProcessOver(queuePoints, threads_q) # Kadangi join() ir terminate() neveikia, MECHANISKAS join()
        self.FormatData_toList(queuePoints) # Is Sync Queue issitraukiam data ir paverciam i listus
        return

    def ThreadJob(self, threadNumber, intervalsList, step, fmin, t0, tN, tkList, N):
        xList_Inner = []
        yList_Inner = []
        interval = intervalsList[threadNumber]  
        print ('Process nr:'+str(threadNumber)+' skaiciavimu intervalas: '+ str(interval[1]-1))
        xList_Inner.append(threadNumber) # !@@@@ xList_Inner[0] uzkoduojamas threadNr
        for j in range(interval[0], interval[1]): # j - zinsniu kiekis, kiek h reikia padaryti fmin -> fmax;
            fj = fmin*10**(j*step) #  formule: gaunam f pagal duota rezoliucija, kitaip x koord
            ###### S(f) apskaiciavimas #######
            c = 0.0
            s = 0.0
            for k in range(0, N): # Sigma susumavimas sin ir cos reiskiniu
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
        self.queuePoints.put(pList)
        self.queuePoints.task_done()
        sys.exit() # Paprograme (Process) pats save susinaikina sistemoje "End Task"
    
    def WaitUntil_ProcessOver(self, queuePoints, threads_q):
        working = True
        while (working):
            if (queuePoints.qsize() == threads_q):
                working = False
                break;
            else:
                time.sleep(1)
        return
    def FormatData_toList(self, queuePoints):
        # Is Sync Queue issitraukiam data ir paverciam i listus
        for i in range(queuePoints.qsize()): # Gali reikti que.size naudoti
            pList = queuePoints.get()
            listX = pList[0]
            listY = pList[1]
            threadNr = listX[0]
            listX.remove(listX[0])
            self.xList.insert(threadNr, listX)
            self.yList.insert(threadNr, listY)
        # Sujungiam lista (is list[[]] i list[])
        self.xList = list(itertools.chain.from_iterable(self.xList)) 
        self.yList = list(itertools.chain.from_iterable(self.yList))
        return

class Atlikimas:
    minX = 0 
    maxX = 0
    minY = 0 
    maxY = 0
        #Konstruktorius
    def __init__(self, minX, maxX, xList, yList):
        self.minX = minX
        self.maxX = maxX
        self.minY, self.maxY = self.Find_minY_maxY(yList)
        self.DrawGraphic(self.minX, self.maxX, self.minY, self.maxY, xList, yList)

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
        maxY_lim = delta*0.75 + np.log10(minY)

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

    def DrawGraphic(self, minX, maxX, minY, maxY, xList, yList):
        #MatPlotLib grafiko nustatymai:
        image = plt.figure()
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
        image.savefig("NMMgraph.pdf") #issaugo "figura"(plot) faile
        plt.show ()

if __name__ == '__main__':
    fmin = 10**-7
    fmax = 10**-1
    lmd = 1.1
    N = 10**5
    rez = 300
    threads = 4 # if None - procesas be paprogramiu

    spr = Sprendimas(fmin, fmax, lmd, N, rez, threads) # (lmd, N, rez, threads_quntity (none if not multithreaded))
    Atlikimas(fmin, fmax, spr.xList, spr.yList)