/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nmm_2_5;

/**
 *
 * @author Tomas
 */
public class NMM_2_5 {

    /**
     * @param args the command line arguments
     */
    /* Konstantu apsibrezimas*/
    final static double a=100;
    final static double tau=0.01;
    final static double delta=0.0001;
    final static double T=2;
    final static int N=100;
    final static double h=1/N;
    final static double xt=0.56;
    final static double tt=0.3;
    static Complex y[] = new Complex[N+1];
    /*Kompleksiniu skaiciu klase naudota ir pirmoje dalyje (rasta internete)*/
    public static class Complex {
    private final double re;   // the real part
    private final double im;   // the imaginary part

    // create a new object with the given real and imaginary parts
    public Complex(double real, double imag) {
        re = real;
        im = imag;
    }

    // return a string representation of the invoking Complex object
    public String toString() {
        if (im == 0) return re + "";
        if (re == 0) return im + "i";
        if (im <  0) return re + " - " + (-im) + "i";
        return re + " + " + im + "i";
    }

    // return abs/modulus/magnitude and angle/phase/argument
    public double abs()   { return Math.hypot(re, im); }  // Math.sqrt(re*re + im*im)
    public double phase() { return Math.atan2(im, re); }  // between -pi and pi
    // return a new Complex object whose value is (this + b)
    public Complex plus(Complex b) {
        Complex a = this;             // invoking object
        double real = a.re + b.re;
        double imag = a.im + b.im;
        return new Complex(real, imag);
    }
    public Complex minus(Complex b) {
        Complex a = this;             // invoking object
        double real = a.re - b.re;
        double imag = a.im - b.im;
        return new Complex(real, imag);
    }

    // return a new Complex object whose value is (this - b)

    // return a new Complex object whose value is (this * b)
    public Complex times(Complex b) {
        Complex a = this;
        double real = a.re * b.re - a.im * b.im;
        double imag = a.re * b.im + a.im * b.re;
        return new Complex(real, imag);
    }

    // scalar multiplication
    // return a new object whose value is (this * alpha)
    public Complex times(double alpha) {
        return new Complex(alpha * re, alpha * im);
    }

     // return the real or imaginary part
    public double re() { return re; }
    public double im() { return im; }

    // return a / b
    public Complex div(Complex b) {
        Complex a = this;
        double den = (b.abs())*(b.abs());
        Complex c = new Complex(0,0);
        c=(a.times(b));
        return new Complex(c.re()/den,c.im()/den);
    }
    }
    public static double solve(Complex[] cd, Complex[] rd, Complex[] ld, Complex[] fj, Complex[] y)
		{
		
		//nauja centrinė
                Complex[] alpha = new Complex[N + 1];
                
		//dešinė dalis
                Complex[] z = new Complex [N + 1];
                //D centrinė įstrižainė cd
                //f kairė įstrižainė ld
                //e dešinė rd
                //b - Fj
                alpha[0] = cd[0];
		z[0] = fj[0];
                for (int i = 1; i < N+1; i++)
			{
			Complex t = new Complex(0,0);
                        t=ld[i-1].div(alpha[i-1]);
			alpha[i] = cd[i].minus(rd[i-1].times(t));
			z[i] = fj[i].minus(z[i-1].times(t));
			if (alpha[i].abs() == 0.0)
				{
				return(-1);
				}
			}

		// Back substitution.
		int Nminus1 = N;
		y[Nminus1] = z[Nminus1].div(alpha[Nminus1]);
		for (int i = N-1; i >= 0; i--)
			{
			y[i] = z[i].minus(rd[i].times(y[i+1])).div(alpha[i]);
			}
                return(1);
		}
    
    public static double fsmalltest(){
        return(0);
    }
    public static double fbigtest(){
        return(0);
    }
    public static double ThomasTest(){
        return(0);
    }
    /*Funkciju u(x,t) (ir ju isvestiniu) ir f(x,t) skaiciavimas, formules gautos is Wolfram Mathematica*/

    public static Complex uxt(double x, double t){
        double reu=x*x-x;
        double imu=(x*x - x)*t;
        Complex b=new Complex(reu,imu);
        return(b);
    }
    public static Complex uxt_x(double x, double t){
        double reu=2*x-1;
        double imu=(2*x - 1)*t;
        Complex b=new Complex(reu,imu);
        return(b);
    }
    public static Complex uxt_x2(double x, double t){
        double reu=2;
        double imu=2*t;
        Complex b=new Complex(reu,imu);
        return(b);
    }
    public static Complex uxt_t(double x, double t){
        double reu=0;
        double imu=x*x-x;
        Complex b=new Complex(reu,imu);
        return(b);
    }
    public static Complex fxt(double x, double t){
        double reu=-2*a*a+2*(x*x*x)-2*(x*x*x)*t*t-2*x*x+2*(x*t)*(x*t)-x*x+x-x*t*t;
        double imu=x*x-x-2*a*a*t+4*(x*x*x)*t-4*x*x*t-2*t*x*x+2*x*t;
        Complex b=new Complex(reu,imu);
        return(b);
    }
    public static Complex Fj_fun(double x, double t, Complex u_old[]){
        Complex nFj= new Complex(0,0);
        Complex oFj=new Complex(0,0);
        Complex ooFj=new Complex(0,0);
        Complex u_xt=uxt(x,t);
        // u (j+1,t)
        Complex u_xp1t=uxt(x+h,t);
        //u(j -1, t) 
        Complex u_xm1t=uxt(x-h,t);
        nFj=u_xt.times((2*h*h)/(tau*a*a));
        oFj=(u_xp1t.minus(u_xt.times(2))).plus(u_xm1t);
        nFj=nFj.plus(oFj);
        oFj=u_old[(int)(x*N)].plus(u_xt);
        ooFj=u_old[(int)((x+h)*N)].minus(u_old[(int)((x-h)*N)]);
        ooFj=ooFj.plus(u_xp1t).minus(u_xm1t);
        oFj=oFj.times(ooFj).times(h/(8*a*a));
        nFj=nFj.plus(oFj);
        oFj=u_old[(int)(x*N)].times(u_old[(int)((x+h)*N)].minus(u_old[(int)((x-h)*N)]));
        ooFj=u_xt.times(u_xp1t.minus(u_xm1t));
        nFj=nFj.plus(oFj.plus(ooFj).times(h/(4*a*a)));
        oFj=(fxt(x+h,t+tau).plus(fxt(x,t))).times((h*h)/(a*a));
        nFj=nFj.plus(oFj);
        return(nFj);
    }
    
    /*=================================================================================================*/
    public static void main(String[] args) {
    Complex u[][] = new Complex [N][(int)(T/tau)];
    Complex u_o[] = new Complex[N];
    Complex cdiag[] = new Complex[N+1];
    Complex ldiag[] = new Complex[N+1];
    Complex rdiag[] = new Complex[N+1];
    Complex fj_di[] = new Complex[N+1];
    Complex c= new Complex(2+(2*h*h)/(a*a*tau),0);
    //double tau1=0;
    for (int i=0;i<N;i++){
        u[i][0]=uxt(i,0);
        u_o[i]=u[i][0];
    }
    for (int i=0;i<=N;i++){
        cdiag[i]=c.times(-1);
        ldiag[i]=new Complex(1,0);
        rdiag[i]=new Complex(1,0);
        fj_di[i]=Fj_fun(i*h,0,u_o);
    }
    ldiag[0]=new Complex(0,0);
    rdiag[N]=new Complex(0,0);
    cdiag[0]=new Complex(1,0);
    cdiag[N]= new Complex(1,0);
    int stop=0;
    double progresas=0;
    while (stop!=1)
        {
           if( solve(cdiag,ldiag,rdiag,fj_di,y)==1) {
            for(int i=0;i<N;i++){
                if((y[i].minus(u_o[i])).abs()>progresas)
                    progresas=(y[i].minus(u_o[i])).abs();
            }
            if (progresas<delta)
                stop=1;
            }
           u_o=y;
        }
    
    }
}
