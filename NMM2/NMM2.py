# Uzduotis: 22.  Lygtis F, algoritmas F1, kraštinės sąlygos II.
#Netiesinė Kuramoto-Cuzuki lygtis
#Lygits: ∂u/∂t = (a^2 + i)* ∂^2*u / ∂*x^2 + i*c*u + id|u|^2*u + f(x,t); 0 < x < 1; t > 0
#Algoritmas: (unj-uj)/r = (a^2 + i)*(1/2)*((un[j+1] - 2unj + un[j-1]) / h2 + (u[j+1] - 2uj + u[j-1]) / h2) +
#                 ic*((unj+uj)/2) + id*|(unj + uj)/2|^2 * ((unj+uj)/2) + ((fnj + fj)/2)
#                j = [1,2,..,N-1]

# Krastines salygos: (∂u/∂x|x=0) = (∂u/∂x|x=1) = 0; t>0 (Grafiko galuose y = 0)

#u(x,0) = u0(x) # 0<=x<=1
#For Wolfram: (uNj - uj )/ Tau = (a^2+i)*0.5*((uNjp - 2*uNj + uNjm)/h^2 + (ujp - 2*uj + ujm)/h^2) + i*c((uNj +uj)/2) + i*d*|(uOj + uj) /2|^2 * ((uOj+uj)/2) + ((fNj +fj) / 2)
#F = 2*h**2/(a**2+i)*Tau*uj + (ujP - 2*uj + ujM) + i*c*h**2 / (a**2+i) * (Uj_old + uj) + i*d*h**2 / (a**2+i) * (abs((Uj_old+uj)/2)**2)*Uj_old+uj + h**2*((Fj+fj)/(a**2+i))
#C = 2+ 2*h**2/(a**2+i)*Tau
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pylab as plt
import math
import numpy as np
import random
import cmath

pi = cmath.pi
i = 1j
# Algoritmui:
# koeficiantai gali buti max int(3), kitaip overflow!
a=2.2;
c=0.6;
d=0.1;
# Aproksimacijai:
N = 100
h = 1/N
T = 0.5
t = 0.0
Tau = 0.01;
# Testams:
delta = 10**-6;
test_x = 0.74;
test_t = 3.45;

# Darbiniai metodai:
def sin(x):
    return cmath.sin(x)

def cos(x):
    return cmath.cos(x)

def Get_MaxValue_List(uList):
    maxU = uList[0]
    for u in uList:
        if (u > maxU): maxU = u
    return maxU

# Pirmasis sprendimas: (f(x,t) ir u0(x,0) radimas)
def u_tikslus_Point(x, t):
    #uPoint = (1-i*t)*(sin(pi*x)**2) #is paskaitos
    #uPoint = 3.0*(2.0-i*t)*((x**3)/3.0 - (x**2)/2.0)
    #uPoint = (1 - i*t)*(-x**2/2 + x**3/3)
    uPoint =  (2 - i*t**2)*(-x + x**3/3)
    #print ('uPoint = ', uPoint)
    return uPoint

def f_Point(x,t):
    u = u_tikslus_Point(x, t)
    #udt = -i*pi*x*(sin(x)**2)
    #uddxx = 2*pi*(1 - i*t)*(x*cos(2*x) + sin(2*x))
    #udt = 0.5*i*(3.0-2.0*x)*x**2
    #uddxx = 3.0*(2.0-i*t)*(2.0*x-1.0)
    #udt = -1/6*i*x**2*(-3 + 2*x)
    #uddxx = (1 - i*t)*(-1 + 2*x)
    udt = -2/3*i*t*x*(-3 + x**2)
    uddxx = 2*(2 - i*t**2)*x
    f = udt - (a**2+i)*uddxx - i*c*u - i*d*(abs(u)**2)*u
    return f

def f_List(t):
    pi = math.pi
    i = 1j
    fList = []
    for j in range (0,N+1): #[1,N-1]
        f = f_Point(j*h, t)
        fList.append(f)
    return fList

def u_tikslus_List(t):
    uList = []
    for j in range (1, N): #[1,N-1]
        u = u_tikslus_Point(j*h, t)
        uList.append(u)
    uList.insert(0, uList[0])
    uList.append(uList[-1])
    return uList

def Init_uf_Points(x, t, h, Tau):
    uj = u_tikslus_Point(x, t)
    ujP = u_tikslus_Point(x+h, t)
    ujM = u_tikslus_Point(x-h, t)

    uj_St = u_tikslus_Point(x, t+Tau)
    ujP_St = u_tikslus_Point(x+h, t+Tau)
    ujM_St = u_tikslus_Point(x-h, t+Tau)

    fj = f_Point(x, t)
    fj_St = f_Point(x, t+Tau)

    return uj, ujP, ujM, uj_St, ujP_St, ujM_St, fj, fj_St

def Alg_leftSide(uj, uj_St, Tau):
    return (uj_St-uj)/Tau

def Alg_right_A(uj, ujP, ujM, uj_St, ujP_St, ujM_St, h, Tau):
    return (a**2 + i)*0.5*(((ujP_St - 2.0*uj_St + ujM_St)  +  (ujP - 2.0*uj + ujM)) / (h**2))

def Alg_right_C(uj, uj_St, kof=1):
    return kof*i*c*((uj_St+uj)/2.0)

def Alg_right_D(uj, uj_St, kof=1):
    return kof*i*d*(abs((uj_St+uj)/2.0)**2) * ((uj_St+uj)/2.0)

def Alg_right_f(fj, fj_St, kof=1):
    return kof*((fj_St + fj)/2.0)

def F_Point (u, uP, uM, U_old, f, fSt, h, Tau):
    kof = 2.0*h**2/(a**2+i)
    F = kof/Tau*u + (uP - 2.0*u + uM) + Alg_right_C(u, U_old, kof) + Alg_right_D(u, U_old, kof) + Alg_right_f(f, fSt, kof)
    return F

def Fj_radimas(t, u, U_old):
    Fj = []
    f = f_List(t)
    fSt = f_List(t+Tau)
    for j in range(1, N): # [1; N]
        F = F_Point(u[j+1], u[j], u[j-1], U_old[j], f[j], fSt[j], h, Tau)
        #F = (2*h**2/(a**2+i)*Tau)*u[j] + (u[j+1] - 2*u[j] + u[j-1]) + (i*c*h**2 / (a**2+i)) * (U_old[j] + u[j]) + (i*d*h**2 / (a**2+i)) * (abs((U_old[j]+u[j])/2)**2)*U_old[j]+u[j] + h**2*((fSt[j]+f[j])/(a**2+i))
        Fj.append(F)
    #print(Fj)
    return Fj

def C_radimas(h, Tau):
    C = 2.0 + 2.0*h**2/((a**2+i)*Tau)
    return C


####### TESTAS 1 - f(x,t) ##################################################
def Testas_1(Tau, h):
    print ('=============== Testas 1-Nr. ================')
    errorList = []
    #uj, ujP, ujM, uj_St, ujP_St, ujM_St, fj, fj_St = Init_uf_Points(test_x, test_t, h, Tau)
    for tikslumas in range (1,5):
        uj, ujP, ujM, uj_St, ujP_St, ujM_St, fj, fj_St = Init_uf_Points(test_x, test_t, h, Tau)
        leftSide = Alg_leftSide(uj, uj_St, Tau)
        right_A = Alg_right_A(uj, ujP, ujM, uj_St, ujP_St, ujM_St, h, Tau)
        right_C = Alg_right_C(uj, uj_St)
        right_D = Alg_right_D(uj, uj_St)
        right_f = Alg_right_f(fj, fj_St)
        neatitiktis = abs(leftSide - right_A - right_C - right_D - right_f)
        errorList.append(neatitiktis)
        h = h/10.0
        Tau = Tau/10.0
    for k in range(0, len(errorList)-1):
        error = errorList[k] / errorList[k+1]
        if (70.0 < error):
            print ("Paklaida:", "%.2f"%error, 'kartu', 'Testas 1 - SUCCESS')
        else:
            print ("Paklaida:", "%.2f"%error,'kartu',  'Testas 1 - FAIL')
    return

####### TESTAS 2 - F(u, uSt_old, fj, fj_St) ################################
def Testas_2 (Tau, h):
    print ('=============== Testas 2-Nr. ================')

    errorList = []
    for tikslumas in range (1,5):
        uj, ujP, ujM, uj_St, ujP_St, ujM_St, fj, fj_St = Init_uf_Points(test_x, test_t, h, Tau)
        C = C_radimas(h, Tau)
        F = F_Point (uj, ujP, ujM, uj_St, fj, fj_St, h, Tau)
        neatitiktis = abs(ujP_St - C*uj_St + ujM_St + F)
        errorList.append(neatitiktis)
        h = h/10.0
        Tau = Tau/10.0

    for k in range(0, len(errorList)-1):
        error = errorList[k] / errorList[k+1]
        if (9500 < error and error < 12000):
            print ("Paklaida:", "%.2f"%error, 'kartu', 'Testas 1 - SUCCESS')
        else:
            print ("Paklaida:", "%.2f"%error,'kartu',  'Testas 1 - FAIL')
    return

####### TESTAS 3 - Thomas Alg. #############################################
def Testas_3():
    print ('=============== Testas 3-Nr. ================')
    # Pasigaminam test y List len: N+1
    y_test = []
    for j in range(1, N):
        y_test.append(complex(random.uniform(1000.0, -1000.0), random.uniform(1000.0, -1000.0))) # Test y
    y_test.insert(0, y_test[0])
    y_test.insert(-1, y_test[-1])
    # Pasigaminam test Fj List len: N-1
    C = C_radimas(h, Tau)
    Fj = []
    for j in range(1, N):
        Fj.append( C*y_test[j] - y_test[j+1] -y_test[j-1])

    # Testuojam Thomas_Alg su testiniais C ir Fj, gaunam yList
    y_temp = Thomas_Algoritmas(C, Fj)
    errorList = []
    for j in range (0, N+1):
        errorList.append(abs( y_test[j] - y_temp[j] ))

    maxError = Get_MaxValue_List(errorList)
    if(maxError < delta): print ('Paklaida:', maxError,'Testas 3 - SUCCESS')
    else: print('Paklaida:', maxError,'Testas 3 - FAIL')
    return

################### Algoritmo sudedamosios dalys ###########################

#Thomas Alg
def Thomas_Algoritmas(C, Fj):
    # Pirmas zingsnis pagal 3 punkta - rasti alfa1 ir beta1
    alfa_j = kapa1 = kapa2 = 1.0
    alfaList = [alfa_j] # nes alfa0=alfa1
    beta_j = gama1 = gama2 = 0.0
    betaList = [beta_j] # nes beta0=beta1

    # Antras zingsnis pagal 2 punkta - rasti afla_List ir beta_List

    for j in range(1, len(Fj)+1): # j = [1; N-1], is viso alfj = [0;N], nes negalime suskaiciuoti  Fj[N]
        alfa_j = 1.0 / (C - alfa_j)
        beta_j = (Fj[j-1] + beta_j)*alfa_j
        alfaList.append(alfa_j)
        betaList.append(beta_j)
    #print ('F len', len(Fj), 'alfa len', len(alfaList))
    # Trecias zingsnis pagal 4 punkta - rasti yN
    yN = (kapa2 * betaList[-1] + gama2) / (1 - kapa2 * alfaList[-1])  # yN = (kapa2 * beta[])
    yList = [yN]

    # Ketvirtas zingsnis pagal 1 punkta - atgaline tvarka rasti yList nuo y[N-1; 0]
    for j in range(N-1, 0, -1): # nes betaList [1;N-1]
        y = alfaList[j] * yList[-1] + betaList[j]
        yList.append(y)
    # Kadangi II krastiniu salygos ty y0=y1 ir yN-1 = yN, todel idedam i List galus, kad uzbaigtume lygciu sistema len(N+1) 
    yList.append(alfaList[0]*yList[-1] + betaList[0]) # prideda y0, kurio anksciau nepridejome
    yList.reverse()
    return yList

def uNew_Tikslinimas (ujSt_Old, ujSt_New):
    errorList = []
    utk = u_tikslus_List(t+Tau)
    for j in range (0, N+1):
        errorList.append(abs( ujSt_Old[j] - ujSt_New[j] ))
    maxError = Get_MaxValue_List(errorList)
    #print (maxError)
    if(maxError < delta):
        print (maxError, "Konvergacija - SUCCESS")
        return True
    else: 
        return False

###################### ALGORITMAS ############################
def Algoritmas(t, T):
    Testas_1(Tau, h)
    Testas_2(Tau, h)
    Testas_3()
    print ('================= Atlikimas =================')
    uj = ujSt_Old = u_tikslus_List(0) # randamas u0
    print('u0_List ilgis:', len(uj)) # j = [0;N] ty len(N+1)
    #DrawGraphic(uj)

    while (t <= T):
        #rasti ujNew
        Fj = Fj_radimas(t, uj, ujSt_Old)
        C = C_radimas(h, Tau)
        ujSt_New = Thomas_Algoritmas(C, Fj)
        iterCount = 0 # apsauga nuo while(true)

        #Testuojam uj^ : Paklaida = max|ujNew - ujOld)|
        while (not uNew_Tikslinimas(ujSt_Old, ujSt_New)):
            iterCount +=1
            ujSt_Old = ujSt_New
            C = C_radimas(h, Tau)
            Fj = Fj_radimas(t, uj, ujSt_New)
            ujSt_New = Thomas_Algoritmas(C, Fj)
            if iterCount > 20:
                print('Konvergacija - FAIL')
                break
        #DrawGraphic(ujSt_New)
        t = t+Tau
        uj = u_tikslus_List(t)
        ujSt_Old = ujSt_New


def DrawGraphic(yList):
    #MatPlotLib grafiko nustatymai:
    #g_minX, g_maxX, g_minY, g_maxY = self.Format_GraphLim(minX, maxX, minY, maxY)
    #plt.axis([g_minX, g_maxX, g_minY, g_maxY]) #[Xmin, Xmax], [Ymin, Ymax]
    plt.grid (True)
    xList = [0]
    for j in range(1, len(yList)):
        xList.append(j*h)
    # Taskinio proceso grafikas
    xEndings = [xList[0], xList[1], xList[-2], xList[-1]]
    yEndings = [yList[0], yList[1], yList[-2], yList[-1]]
    plt.scatter (xEndings, yEndings, 10, label='NMM1', color='blue') # grafikas taskais
    plt.plot (xList, yList, label='NMM2', color='red') # grafikas lauzte

    # Draw or Imaging
    #plt.savefig("NMMgraph.pdf") #issaugo "figura"(plot) faile
    plt.show ()

print ('############# SPRENDIMAS ####################')
print ('N =', N, '| h =', h, '| Tau =', Tau, '| a =', a, '| c =', c,'| d =', d)
print ('Test x =', test_x, ' Test t = ', test_t)
Algoritmas(t, T)