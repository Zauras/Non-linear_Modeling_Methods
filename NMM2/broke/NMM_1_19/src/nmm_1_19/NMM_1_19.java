/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nmm_1_19;
import java.awt.*;
import javax.swing.*;        
/**
 *
 * @author Tomas
 */
public class NMM_1_19 {
    
	/*
	Resolution parameters
	*/
    final static int rezx=600;
    final static int rezy=600;
    /*
	Amount of steps
	*/
	final static int nmax=1000;
    /*
	Artificial infinity
	*/
	final static int infinity=10000;
    
	/* 
	Minimal x value
	*/
    final static double xmin=-2;
	/* 
	Maximal x value
	*/
    final static double xmax=2;
	/* 
	Minimal x value
	*/
    final static double ymin=-2;
	/* 
	Maximal y value
	*/
    final static double ymax=2;
    /*
	Complex numbers class from internet.
	*/
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
 
    
    }
/*
Result matrix for visual representation
*/
public static int rezult[][] = new int[rezx][rezy];
/*
Class for visual representing
*/
public static class Frame extends JFrame  
{  
public Frame()  
{  
 super("Atvaizdavimas");
 setContentPane(new DrawPanel());
 setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);  
 setSize(rezx,rezy);  
 setVisible(true);
}
class DrawPanel extends JPanel{
        @Override
        public void paintComponent(Graphics g){
           for (int i =0;i<rezx;i++){
               for(int j = 0;j<rezy;j++){
                   /*
				   Change color depending on count |Z_n|<|Z_(n+1)|
				   */
				   g.setColor(Color.getHSBColor(rezult[i][j] / 256.0f, 1.0f, 1.0f));
                   g.fillRect(i,j, i+1, j+1);
               }
           }
           
         }
     }
}
   
    public static void main(String[] args) {
      
       /*
	   Step size on x axis
	   */
	   double hx=0d;
       hx=(xmax-xmin);
       hx=hx/rezx;
       /*
	   Step size on y axis
	   */
	   double hy=0d;
       hy=(ymax-ymin);
       hy=hy/rezy;
       double x=-hx+xmin;
       double y=hy+ymax;
       Complex c = new Complex(0,0);
       Complex z[]=new Complex[nmax];
       for (int i=0; i<rezx;i++){
            //x=x+hx;
            for (int j=0;j<rezy;j++){
              //  y=y+hy;
                x=xmin+i*hx;
                y=ymin + j*hy;
                z[0]=new Complex(0.15,0);
                c=new Complex(x,y);
                rezult[i][j]=0;
                for (int n=1;n<nmax;n++){
                    z[n]=new Complex(0,0);
                    z[n]=z[n-1].times(z[n-1]).plus(c);
                    if (z[n].abs()<infinity)
                    {
                        if (z[n].abs()>z[n-1].abs())
                            rezult[i][j]=rezult[i][j]+1;
                    }
                    if (z[n].abs()>infinity)
                    {
                            rezult[i][j]+=nmax-n;
                            break;
                    }
                }
           
           }
       }
    Frame Remas = new Frame();
    }

}