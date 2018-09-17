import numpy as np
import matplotlib.pylab as plt
#matplotlib.use('pdf') bet reikia import matplotlib
import math
import random
import time
import threading
#from processing import process as worker
from multiprocessing import Process
#import multiprocessing
import itertools

class Sprendimas:

    tkList = []
    xList = []
    yList = []
    N = pow(10,5)
    lmd = 1.1
    t0 = 0
    tN = None
    rez = 1
    threads_quantity = None
    threadsList = []
    threads_intervalList = [[]]

    #fmax = 10**-0 
    #fmin = 10**-8 # PILNO VAIZDO GRAFIKAS
    fmax = 10**-1 
    fmin = 10**-7 # UŽDUOTIES GRAFIKAS
    #fmax = 10**-3.5 
    #fmin = 10**-5.5 # GRAFIKAS FITINIMUI (geriausias rezultatas) - alfos(slope) paieškoms

    #Konstruktorius
    def __init__(self, lamda, N, rezoliucija, threads=None):
        self.lmd = lamda
        self.N = N
        self.rez = rezoliucija
        self.FindTn()
        if (threads != None):
            self.threads_quantity = threads
        self.FindPointsXY()
        


    def PoissonRandom(self, lmd):
        a = np.e**(-1.0 * lmd)
        r = 1
        n = -1
        while r > a:
            u = np.random.random(1)
            r *= u 
            n += 1
        return n

    #Generuojamas Puasso list ir tN
    #Pagal: Tk = Tk-1 + P(lmd)k;  P(lmd)-Poiss random nr; T0 = 0
    # N ir lmd - ivedami; Gaunam: tList ir Tn(tList pask. nr.)
    def FindTn (self):
        lmd = self.lmd
        N = self.N
        tk = self.t0

        for k in range(1, N+1):
            #tk += np.random.poisson(lmd) #rankom
            tk += self.PoissonRandom(lmd)
            self.tkList.append(tk)
            continue
        self.tN = tk
        return

    #Apskaiciuojam funkcija S(f) arba kitaip y pagal duota formule
    def FindSf (self, f, tN, N):
        t0 = self.t0
        c = 0.0
        s = 0.0
        for k in range(0, N): # Sigma susumavimas sin ir cos reiskiniu
            c += np.cos(2.0 * np.pi * f * self.tkList[k])
            s += np.sin(2.0 * np.pi * f * self.tkList[k])
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
        step = (math.log10(fmax) - math.log10(fmin)) / rez # surandame zingsnio ilgi h, bet skaiciuojam pagal log10 skale
        start_time = time.time()
        # Pasirenkamas budas skaiciavimui atlikti:
        if (self.threads_quantity == None):
            self.CalculatePointsXY_simple(tN, N, rez, fmin, step)
        else: 
        	self.CalculatePointsXY_multithreading(tN, N, rez, fmin, step)
        print("Skaiciavimu trukme: %s s" % (time.time() - start_time))
        return

    def CalculatePointsXY_simple(self, tN, N, rez, fmin, step):   
        for j in range(0, rez): # k - zinsniu kiekis, kiek h reikia padaryti fmin -> fmax; 
                                #fmax != 1, ty fj!=fmax, nes iskreips grafika
            fj = fmin*10**(j*step) #  formule: gaunam f pagal duota rezoliucija, kitaip x koord
            Sfj = self.FindSf(fj, tN, N) # S(f) - y koord
            self.xList.append(fj)
            self.yList.append(Sfj)
        return

    def ThreadJob(self, threadNumber, intervalsList, step, fmin, tN, N):
        while True:
            xList_Inner = []
            yList_Inner = []
            interval = intervalsList[threadNumber]
            print ('Thread nr:'+str(threadNumber)+' skaiciavimu intervalas: '+ str(interval))
            for j in range(interval[0], interval[1]): # k - zinsniu kiekis, kiek h reikia padaryti fmin -> fmax; 
                                                            #fmax != 1, ty fj!=fmax, nes iskreips grafika
                fj = fmin*10**(j*step) #  formule: gaunam f pagal duota rezoliucija, kitaip x koord
                Sfj = self.FindSf(fj, tN, N) # S(f) - y koord
                xList_Inner.append(fj)
                yList_Inner.append(Sfj)

            self.xList.insert(threadNumber, xList_Inner)
            self.yList.insert(threadNumber, yList_Inner)
        return

    def CalculatePointsXY_multithreading(self, tN, N, rez, fmin, step):
        intervalsList = self.threads_intervalList
        threads_q = self.threads_quantity

        sveikojiDalis =  int(rez / threads_q)
        liekana = rez % threads_q
        begin=0
        end = sveikojiDalis
        # Veiksmai su Threadais
        for threadNumber in range (0, threads_q):
            # Rezoliucijos intervalas padalinamas i intervaliukus gijoms
            if (liekana > 0): 
                end += 1
                liekana -= 1
            intervalsList.insert(threadNumber, [begin, end])
            begin = end
            end = end + sveikojiDalis
            # Threadu sukurimas su eiles numeriu
            
            thread = threading.Thread(target = self.ThreadJob, args=(threadNumber, intervalsList, step, fmin, tN, N,))
            self.threadsList.append(thread) #idedam thread, bet pakeiciau pavadinima
            self.threadsList[threadNumber].start() # Startuojamas def ThreadJob()

        #Threadu intervaliniu listu sujungimas i baigtinius List X, Y
        [t.join() for t in self.threadsList] # Palaukiam, kol gijos baigs darba
        self.xList = list(itertools.chain.from_iterable(self.xList)) # Sujungiam lista (is list[[]] i list[])
        self.yList = list(itertools.chain.from_iterable(self.yList))
        return

    def FindAlfa_byY(self):
        xList = self.xList
        yList = self.yList
        #atkarpa tarp 25% ir 75%(nuo didziausios y reiksmes) Y asies(S(f)) reiksmiu
        maxY = yList[0]
        for y in range(len(yList)):
            if (yList[y] > maxY): maxY = yList[y]
        sliceListX=[]
        sliceListY=[]
        for y in range(len(yList)):
            if(yList[y] < maxY*0.75 or yList[y] > maxY*0.25):
                sliceListX.append(xList[y])
                sliceListY.append(yList[y])

        logX = np.log10(sliceListX) #Pasiverciam i log skale, kad butu tiesinis alfa
        logY = np.log10(sliceListY) #
        alfa, const = np.polyfit(logX, logY, 1) # fit log(S(f)) = alfa*log(f) + const

        #FitSf = 10**(alfa * logX + const) # apskaiciuojam Fit y ir ATSIVERCIAM is log skales atgal i norml
        #return sliceListX, FitSf
        FitSf_full = 10**(alfa * np.log10(xList) + const) # Fit breziama per visa grafika
        return xList, FitSf_full

    # Prie jau gauto grafiko, prifitinama y=1/f**alfa (log skaleje: y=-alfa*x+const) pasirinktam diapozone
    def FindAlfa_byPoints(self):
        xList = self.xList
        yList = self.yList
        #atkarpa tarp 25% ir 75% grafiko tasku
        yBeg_index = int(len(yList)*0.25)
        yEnd_index = int(len(yList)*0.75)
        sliceListX = xList [yBeg_index : yEnd_index]
        sliceListY = yList[yBeg_index : yEnd_index]

        logX = np.log10(sliceListX) #Pasiverciam i log skale, kad butu tiesinis alfa
        logY = np.log10(sliceListY) #
        alfa, const = np.polyfit(logX, logY, 1) # fit log(S(f)) = alfa*log(f) + const

        #FitSf = 10**(alfa * logX + const) # apskaiciuojam Fit y ir ATSIVERCIAM is log skales atgal i norml
        #return sliceListX, FitSf
        FitSf_full = 10**(alfa * np.log10(xList) + const) # Fit breziama per visa grafika
        return xList, FitSf_full

class Atlikimas:
    spr = Sprendimas(2.9, 100000, 50, 3) # (lmd, N, rez, threads_quntity (none if not multithreaded))
    xList = spr.xList
    yList = spr.yList

    #Grafikas:
    plt.loglog (basex=10, basey=10) #(baseX, baseY) - automatiskai visas paduodamas reiksmes konvertuoja i log10xlog
    plt.axis([10**-8, 10**1, 10**-3, 10**5,]) #[Xmin, Xmax], [Ymin, Ymax]
    plt.grid (True)

    plt.plot (xList, yList, label='NMM1', color='red' )
    #plt.scatter (xList, yList, 5, label='NMM1', color='red' )
    
    #Fitinimas pagal Point List
    fit_f_list, fit_Sf_list = spr.FindAlfa_byPoints()
    plt.plot(fit_f_list, fit_Sf_list, color="blue") #nepamirsti, jog grafikas auto-loglog, todel reikia norml reiksmiu

    #Fitinimas pagal maxY
    fit_f_list, fit_Sf_list = spr.FindAlfa_byY()
    plt.plot(fit_f_list, fit_Sf_list, color="violet") 
  
    plt.show ()

'''
    import multiprocessing

class Worker(multiprocessing.Process):

    def run(self):
        print 'In %s' % self.name
        return

if __name__ == '__main__':
    jobs = []
    for i in range(5):
        p = Worker()
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()
'''