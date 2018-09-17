rm(list=ls())
N<-10
a=0.01
tau=0.0001
delta=10^(-8)
Time=2
h=1/N
xt=0.7
tt=0.2
c=(2*h*h/(a*a*tau))-2

N<-4
h<-1/N
c<-3
fj=c(0,1,2,3,0)
thomas(fj)
for (i in 2:N+1){
  alfa[i]=1/(c-alfa[i-1])}
  beta[i]=(fj_21[i-1]+beta[i-1])*alfa[i]
}
ats[N+1]<-0
for (i in N:2){
  ats[i]<-alfa[i+1]*ats[i+1]+beta[i+1]
}

uxt<-function(x,t)
{
  ats<-complex(real=(x*x-x),imaginary=((x*x - x)*t))
  return(ats)
}
fxt<-function(x,t){
  re<-0
  im<-0
  re=-2*a*a-(2*x-1)*(x*x-x)*(1-t*t)
  im=-2*a*a*t + x*x-x - 2*t*(2*x-1)*(x*x-x)
return(complex(real=re,imaginary=im))
}
Fj<-function(x,t,uold){
  u_xt=uxt(x,t)
  u_xp1t=uxt(x+h,t)
  u_xm1t=uxt(x-h,t)
  xas=as.integer(x*N)+1
  xm1=as.integer((x-h)*N)+1
  xp1=as.integer((x+h)*N)+1
  Fj=complex(real=0,imaginary=0)
  Fj=(fxt(x,t+tau)-fxt(x,t))*(h*h)/(a*a)+uxt(x+h,t)-2*uxt(x,t)+uxt(x-h,t)
  Fj=Fj+(h)/(a*a)*((uold[xas]-uxt(x,t))/8*(uold[xp1]-uold[xm1]+uxt(x+h,t)-uxt(x-h,t))+(1/4)*(uold[xas]*(uold[xp1]-uold[xm1])+uxt(x,t)*(uxt(x+h,t)-uxt(x-h,t))))
  return(Fj)
}

t=0
un<-complex(real=0,imaginary=0)
uold<-complex(real=0,imaginary=0)



s=0
fj_di=0
for(i in 1:N+1){
  s=(i-1)*h
fj_di[i]=-Fj(s,t,uold)
}



x<-rep(0,(N+1)^2)

A<-matrix(x,N+1,byrow=T)
x<-0
i=3
b=0
for (i in 2:(N)){
  b=0
  b=i-2
    x<-c(rep(0,b),1,-c,1,rep(0,N-i))
#    print(x)
    A[i,]=x
}

A[1,]=c(1,0,rep(0,N-1))
A[N+1,]=c(rep(0,N-1),0,1)
d0<-c(rep(1,N-1),0)
d2<-c(0,rep(1,N-1))
d1<-c(1,rep(-c,N-1),1)
progr=0

stoti=0
k=0
while (stoti!=1){
k=k+1
un<-thomas(d0,d1,d2,fj_di)
un[is.na(un)] <- 0
progr<-max(abs(un))
if (progr < delta) {stoti=1}
else{uold=un 
     t=t+tau}
s=0
for(i in 1:N+1){
  uold[i]<-uxt((i-1)*h,t)
  
}
for(i in 1:N+1){
  s=(i-1)*h
  fj_di[i]=-Fj(s,t,uold)
}
#print(un)
}

# 1 - funkcijos f(x,t) parinkimo testas
# N<-10000
# h<-1/N
# tau=1/N*10
# ujt=uxt(xt,tt)
# ujt_stog=uxt(xt,tt+tau)
# ujt_p1=uxt(xt+h,tt)
# ujt_m1=uxt(xt-h,tt)
# ujt_stog_m1=uxt(xt-h,tt+tau)
# ujt_stog_p1=uxt(xt+h,tt+tau)
# max(abs((ujt_stog-ujt)/tau - (a^2)/(2*h^2) * ((ujt_stog_p1-2*ujt_stog+ujt_stog_m1)+(ujt_p1-2*ujt+ujt_m1)) -1/(8*h) * ((ujt_stog+ujt)/2 *(ujt_stog_p1-ujt_stog_m1+ujt_p1-ujt_m1) + (ujt_stog *(ujt_stog_p1-ujt_stog_m1)+ujt*(ujt_p1-ujt_m1)))-(fxt(xt,tt+tau)+fxt(xt,tt))/2) )
# 2 - testas
# tau<-tau/10
# N<-N/10
# h<-h/10
# c=(2*h*h/(a*a*tau))-2
# pro<-0
# abs(uxt(xt+h,tt)-c*uxt(xt,tt)+uxt(xt-h,tt)+Fj(xt,tt,uold))
# 3 - testas
