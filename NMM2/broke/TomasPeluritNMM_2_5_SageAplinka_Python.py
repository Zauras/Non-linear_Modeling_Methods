SageThe Sage Notebook Version 5.11
Tomas_Peluritis Home Published Log Settings Help Report a Problem Sign out
NMM_2_5
last edited Dec 16, 2013 7:47:33 AM by Tomas_Peluritis
sage: N=10000
sage: a=0.01
sage: tau=0.00001
sage: delta=1.0/1000000
sage: Time=2
sage: h=(1.0/N)
sage: xt=0.7
sage: tt=0.2
sage: c=(2*h*h/(a*a*tau))+2
sage: t=0
sage: delta
sage: def ThomasSolver(d0, d1, d2, fj):
...       n = len(fj) 
...       for i in xrange(n-1):
...           fj[i+1] -= fj[i] * d0[i] / d1[i]
...           d1[i+1] -= d2[i] * d0[i] / d1[i]
...       for i in reversed(xrange(n-1)):
...           fj[i] -= fj[i+1] * d2[i] / d1[i+1]
...       return [fj[i] / d1[i] for i in xrange(n)]
sage: def uxt(x,t):
...     #ats=complex(x*x-x,(x*x - x)*t)
...     tmp=complex(1,t)
...     ats=tmp*(x*x*x-x*x)
...     return(ats)
sage: uxt(3,0)
sage: def fxt(x,t):
...       #re=0
...       #im=0
...       #re=-2*a*a-(2*x-1)*(x*x-x)*(1-t*t)
...       #im=-2*a*a*t + x*x-x - 2*t*(2*x-1)*(x*x-x)
...       tmp=complex(1,t)
...       ats=-tmp*a*a*(6*x-2)+complex(0,x*x*x-x*x)-tmp*tmp*(x*x*x-x*x)*(3*x*x-2*x)
...       return(ats)
sage: def Fj(x,t,uom1,uo,uop1):
...       u_xt=uxt(x,t)
...       u_xp1t=uxt(x+h,t)
...       u_xm1t=uxt(x-h,t)
...       F=complex(0,0)
...       F=(fxt(x,t+tau)-fxt(x,t))*(h*h)/(a*a)+uxt(x+h,t)-2*uxt(x,t)+uxt(x-h,t)
...       F=F+(h)/(a*a)*((uo-uxt(x,t))/8*(uop1-uom1+uxt(x+h,t)-uxt(x-h,t))+(1/4)*(uo*(uop1-uom1)+uxt(x,t)*(uxt(x+h,t)-uxt(x-h,t))))+2*h*h*uxt(x,t)/(a*a)
...       return(F)
sage: diag= [0] * (N)
sage: diag1=[0] * (N+1)
sage: diag2=[0] * (N)
sage: uold=[0] * (N+1)
sage: fj_di=[0] * (N+1)
sage: diag[N-1]=0
sage: diag2[0]=0
sage: diag1[0]=1
sage: diag1[N]=1
sage: for i in range(1,N):
...       diag[i-1]=1
...       diag1[i]=-c
...       diag2[i]=1
sage: fj_di[0]=0
sage: for i in range(1,N):
...       fj_di[i]=-Fj(i*h,t,uold[i-1],uold[i],uold[i+1])
sage: ABS1=complex(0,0)
sage: #ABS1=[complex(0,0)]*(N+1)
sage: for i in range(0,(N+1)):
...       uold[i]=uxt(i*h,t)
...
sage: paklaida=0
sage: stoti=0
sage: t=0
sage: k=0
sage: while (t<Time):
...       while (stoti<>1):
...           fj_di[0]=0
...           for i in range(1,N):
...               fj_di[i]=-Fj(i*h,t,uold[i-1],uold[i],uold[i+1])
...           un=ThomasSolver(diag,diag1,diag2,fj_di)
...           progr=0
...           k=k+1
...           for i in range(0,N):
...               ABS1=uold[i]-un[i]
...               progr=max(abs(ABS1),progr)
...           if (progr < delta):
...               stoti=1
...           else:
...               uold=un
...       for i in range(0,N):
...           paklaida=max(abs(uold[i]-uxt(i*h,t)),paklaida)
...       k=0 
...       t=t+tau
...       stoti=0
sage: paklaida
sage: N=1000
sage: h=1.0/N
sage: tau=0.001
sage: xt=0.7
sage: tt=0.2
sage: a=10
sage: ujt=uxt(xt,tt)
sage: ujt_stog=uxt(xt,tt+tau)
sage: ujt_p1=uxt(xt+h,tt)
sage: ujt_m1=uxt(xt-h,tt)
sage: ujt_stog_m1=uxt(xt-h,tt+tau)
sage: ujt_stog_p1=uxt(xt+h,tt+tau)
sage: abs((ujt_stog-ujt)/tau - (a*a)/(2*h*h) * ((ujt_stog_p1-2*ujt_stog+ujt_stog_m1)+(ujt_p1-2*ujt+ujt_m1)) -1/(8*h) * ((ujt_stog+ujt)/2 *(ujt_stog_p1-ujt_stog_m1+ujt_p1-ujt_m1) + (ujt_stog *(ujt_stog_p1-ujt_stog_m1)+ujt*(ujt_p1-ujt_m1)))-(fxt(xt,tt+tau)+fxt(xt,tt))/2)
sage: N=1000
sage: h=1.0/N
sage: tau=0.001
sage: xt=0.7
sage: tt=0.2
sage: a=10
sage: c=(2*h*h/(a*a*tau))+2
sage: ujt=uxt(xt,tt)
sage: ujt_stog=uxt(xt,tt+tau)
sage: ujt_p1=uxt(xt+h,tt)
sage: ujt_m1=uxt(xt-h,tt)
sage: ujt_stog_m1=uxt(xt-h,tt+tau)
sage: ujt_stog_p1=uxt(xt+h,tt+tau)
sage: abs(ujt_stog_p1-c*ujt_stog+ujt_stog_m1+Fj(xt,tt,ujt_stog_m1,ujt_stog,ujt_stog_p1))
sage: di0=[1,1,1,0]
sage: di1=[1,-4,-4,-4,1]
sage: di2=[0,1,1,1]
sage: fjd=[0,-11,1,-7,0]
sage: ThomasSolver(di0,di1,di2,fjd)
