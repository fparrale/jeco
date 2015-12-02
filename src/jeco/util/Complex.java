/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jeco.util;

/**
 *
 * @author José Luis Risco Martín <jlrisco at ucm.es>
 * @author Josué Pagán Ortiz <jpagan at ucm.es>
 */
public class Complex {

    protected double real;
    protected double imag;

    public Complex(double real, double imag) {
        this.real = real;
        this.imag = imag;
    }
    
    @Override
    public Complex clone() {
        Complex clone = new Complex(this.real, this.imag);
        return clone;
    }

    // return a string representation of the invoking Complex object
    @Override
    public String toString() {
        if (imag == 0) {
            return real + "";
        }
        if (real == 0) {
            return imag + "i";
        }
        if (imag < 0) {
            return real + " - " + (-imag) + "i";
        }
        return real + " + " + imag + "i";
    }

    // return abs/modulus/magnitude and angle/phase/argument
    public double abs() {
        return Math.hypot(real, imag);
    }  // Math.sqrt(re*re + im*im)

    public double[] abs(Complex[] c) {
        double[] a = new double[c.length];
        
        for (int i=0; i<c.length; i++){
            a[i] = c[i].abs();
        }
        return a;
    }  
     
    public double phase() {
        return Math.atan2(imag, real);
    }  // between -pi and pi

    // return a new Complex object whose value is (this + b)
    public Complex plus(Complex b) {
        Complex a = this;             // invoking object
        double re = a.real + b.real;
        double im = a.imag + b.imag;
        return new Complex(re, im);
    }

    // return a new Complex object whose value is (this - b)
    public Complex minus(Complex b) {
        Complex a = this;
        double re = a.real - b.real;
        double im = a.imag - b.imag;
        return new Complex(re, im);
    }

    // return a new Complex object whose value is (this * b)
    public Complex times(Complex b) {
        Complex a = this;
        double re = a.real * b.real - a.imag * b.imag;
        double im = a.real * b.imag + a.imag * b.real;
        return new Complex(re, im);
    }

    // scalar multiplication
    // return a new object whose value is (this * alpha)
    public Complex times(double alpha) {
        return new Complex(alpha * real, alpha * imag);
    }

    // return a new Complex object whose value is the conjugate of this
    public Complex conjugate() {
        return new Complex(real, -imag);
    }

    // return a new Complex object whose value is the reciprocal of this
    public Complex reciprocal() {
        double scale = real * real + imag * imag;
        return new Complex(real / scale, -imag / scale);
    }

    // return the real or imaginary part
    public double getReal() {
        return real;
    }

    public double getImag() {
        return imag;
    }

    // return a / b
    public Complex divides(Complex b) {
        Complex a = this;
        return a.times(b.reciprocal());
    }

    // return a new Complex object whose value is the complex exponential of this
    public Complex exp() {
        return new Complex(Math.exp(real) * Math.cos(imag), Math.exp(real) * Math.sin(imag));
    }

    // return a new Complex object whose value is the complex sine of this
    public Complex sin() {
        return new Complex(Math.sin(real) * Math.cosh(imag), Math.cos(real) * Math.sinh(imag));
    }

    // return a new Complex object whose value is the complex cosine of this
    public Complex cos() {
        return new Complex(Math.cos(real) * Math.cosh(imag), -Math.sin(real) * Math.sinh(imag));
    }

    // return a new Complex object whose value is the complex tangent of this
    public Complex tan() {
        return sin().divides(cos());
    }

    // a static version of plus
    public static Complex plus(Complex a, Complex b) {
        double re = a.real + b.real;
        double im = a.imag + b.imag;
        Complex sum = new Complex(re, im);
        return sum;
    }
}
