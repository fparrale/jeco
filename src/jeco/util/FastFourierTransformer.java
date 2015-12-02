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
public class FastFourierTransformer {

    public FastFourierTransformer() {
    }

    public Complex[] completeWithZero(Complex[] x) {
        int powerOfTwo = 1;
        long maxPowerOfTo = 2147483648L;
        while (powerOfTwo < x.length && powerOfTwo < maxPowerOfTo) {
            powerOfTwo *= 2;
        }
        Complex[] xx = new Complex[powerOfTwo];
        for (int i = 0; i < x.length; ++i) {
            xx[i] = x[i];
        }
        for (int i = x.length; i < powerOfTwo; ++i) {
            xx[i] = new Complex(0, 0);
        }
        return xx;
    }

    public Complex[] fft(Complex[] xx) {
        Complex[] x = completeWithZero(xx);
        int N = x.length;

        // base case
        if (N == 1) {
            return new Complex[]{x[0]};
        }

        // radix 2 Cooley-Tukey FFT
        if (N % 2 != 0) {
            throw new RuntimeException("N is not a power of 2");
        }

        // fft of even terms
        Complex[] even = new Complex[N / 2];
        for (int k = 0; k < N / 2; k++) {
            even[k] = x[2 * k];
        }
        Complex[] q = fft(even);

        // fft of odd terms
        Complex[] odd = even;  // reuse the array
        for (int k = 0; k < N / 2; k++) {
            odd[k] = x[2 * k + 1];
        }
        Complex[] r = fft(odd);

        // combine
        Complex[] y = new Complex[N];
        for (int k = 0; k < N / 2; k++) {
            double kth = -2 * k * Math.PI / N;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k] = q[k].plus(wk.times(r[k]));
            y[k + N / 2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }

    public Complex[] fft(double[] x) {
        Complex[] s = new Complex[x.length];
        
        for (int i=0; i<x.length; i++) {
            s[i] = new Complex(x[i], 0);
        }
        Complex[] fft = fft(s);
        return fft;
    }
    
    
    // compute the inverse FFT of x[], assuming its length is a power of 2
    public Complex[] ifft(Complex[] x) {
        int N = x.length;
        Complex[] y = new Complex[N];

        // take conjugate
        for (int i = 0; i < N; i++) {
            y[i] = x[i].conjugate();
        }

        // compute forward FFT
        y = fft(y);

        // take conjugate again
        for (int i = 0; i < N; i++) {
            y[i] = y[i].conjugate();
        }

        // divide by N
        for (int i = 0; i < N; i++) {
            y[i] = y[i].times(1.0 / N);
        }

        return y;

    }

    public Complex[] ifft(double[] x) {
        Complex[] s = new Complex[x.length];
        
        for (int i=0; i<x.length; i++) {
            s[i] = new Complex(x[i], 0);
        }
        return ifft(s);       
    }
    
    // compute the circular convolution of x and y
    public Complex[] cconvolve(Complex[] x, Complex[] y) {

        // should probably pad x and y with 0s so that they have same length
        // and are powers of 2
        if (x.length != y.length) {
            throw new RuntimeException("Dimensions don't agree");
        }

        int N = x.length;

        // compute FFT of each sequence
        Complex[] a = fft(x);
        Complex[] b = fft(y);

        // point-wise multiply
        Complex[] c = new Complex[N];
        for (int i = 0; i < N; i++) {
            c[i] = a[i].times(b[i]);
        }

        // compute inverse FFT
        return ifft(c);
    }

    // compute the linear convolution of x and y
    public Complex[] convolve(Complex[] x, Complex[] y) {
        Complex ZERO = new Complex(0, 0);

        Complex[] a = new Complex[2 * x.length];
        for (int i = 0; i < x.length; i++) {
            a[i] = x[i];
        }
        for (int i = x.length; i < 2 * x.length; i++) {
            a[i] = ZERO;
        }

        Complex[] b = new Complex[2 * y.length];
        for (int i = 0; i < y.length; i++) {
            b[i] = y[i];
        }
        for (int i = y.length; i < 2 * y.length; i++) {
            b[i] = ZERO;
        }

        return cconvolve(a, b);
    }

    // display an array of Complex numbers to standard output
    public void show(Complex[] x, String title) {
        System.out.println(title);
        System.out.println("-------------------");
        for (Complex x1 : x) {
            System.out.println(x1);
        }
        System.out.println();
    }
}
