/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jeco.util;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author jlrisco
 */
public class Maths {

  public static double sum(List<Double> numbers) {
      double res = 0;
      for (Double number : numbers) {
          res += number;
      }
      return res;
  }
  
  public static double sum(double[] numbers) {
      double res = 0;
      if (!Double.isNaN(numbers[0])){
          for (Double number : numbers) {
              res += number;
          }
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }

  public static double mean(List<Double> numbers) {
    if (numbers.isEmpty()) {
      return 0;
    }
    double res = sum(numbers) / numbers.size();
    return res;
  }
  
  public static double mean(double[] numbers) { 
      if (!Double.isNaN(numbers[0])){
          double res = sum(numbers)/numbers.length;
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }

  public static double median(List<Double> numbers) {
    Collections.sort(numbers);
    int middle = numbers.size() / 2;
    if (numbers.size() % 2 == 1) {
      return numbers.get(middle);
    } else {
      return (numbers.get(middle - 1) + numbers.get(middle)) / 2.0;
    }
  }

  public static double std(List<Double> numbers) {
    double res = 0;
    double avg = mean(numbers);
    for(Double number : numbers) {
      res += Math.pow(number-avg, 2);
    }
    res = Math.sqrt(res/(numbers.size()-1));
    return res;
  }
  
  public static double std(double[] numbers) {
      if (!Double.isNaN(numbers[0])){
          double res = 0;
          double avg = mean(numbers);
          for (int i = 0; i <= numbers.length-1; i++) {
              res += Math.pow(numbers[i] - avg, 2);
          }
          res = Math.sqrt(res/numbers.length);
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }
  
  
  public static double min(double[] numbers) {
      if (!Double.isNaN(numbers[0])){
          double res = Double.POSITIVE_INFINITY;
          for(int i=0; i<=numbers.length-1; i++){
              if (numbers[i] < res) {
                  res = numbers[i];
              }
          }
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }
  
  public static double max(double[] numbers) {
      if (!Double.isNaN(numbers[0])){
          double res = Double.NEGATIVE_INFINITY;
          for(int i=0; i<=numbers.length-1; i++){
              if (numbers[i] > res) {
                  res = numbers[i];
              }
          }
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }
  
  public static double totalVar(double[] numbers) {      
      if (!Double.isNaN(numbers[0])){
          double res = 0.0;
          double[] derv = new double[numbers.length-1];
          for(int i=0; i<=numbers.length-2; i++){
              derv[i] = numbers[i+1]-numbers[i];
          }
          for (int i = 0; i <= derv.length-1; i++) {
              res += Math.abs(derv[i]);
          }
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }
  
  public static double pod(double[] numbers) {
      if (!Double.isNaN(numbers[0])){
          double res = 1.0;
          for (int i = 0; i <= numbers.length-1; i++) {
              res *= numbers[i];
          }
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }
  
  public static double geoMean(double[] numbers) {
      if (!Double.isNaN(numbers[0])){
          double res = Math.pow(pod(numbers), 1/(numbers.length));
          return res;
      }
      else {
          return Double.POSITIVE_INFINITY;
      }
  }
  
  // Temporal functions
    public static double[] pow(double[] numbers, double exp) {
        if (!Double.isNaN(numbers[0])){
            for (int i = 0; i <= numbers.length-1; i++) {
                numbers[i] = Math.pow(numbers[i], exp);
            }
            return numbers;
        }
        else {
            return new double[] {Double.NaN};
        }
    }
    
    public static double[] conv(double[] x, double[] h) {
        if ((!Double.isNaN(x[0])) && (!Double.isNaN(h[0]))){
            double[] res = new double[x.length + h.length - 1];
            for (int i = 0; i <= res.length -1; i++ )	{
                res[i] = 0;                       // set to zero before sum
                for (int j = 0; j <= h.length-1; j++ ) {
                    if ((j <= i) && ((i-j) < x.length)) {
                        res[i] += x[i - j] * h[j];    // convolve: multiply and accumulate
                    }
                }
            }
            return res;
        }
        else {
            return new double[] {Double.NaN};
        }
    }
    
    public static double[] diff(double[] numbers) {
        if (!Double.isNaN(numbers[0])){
            double[] res = new double[numbers.length-1];
            for (int i = 1; i <= numbers.length-1; i++) {
                res[i-1] = numbers[i] - numbers[i-1];
            }
            return res;
        }
        else {
            return new double[] {Double.NaN};
        }
    }
    
    public static double[] abs(double[] numbers) {
        if (!Double.isNaN(numbers[0])){
            for (int i = 0; i <= numbers.length-1; i++) {
                numbers[i] = Math.abs(numbers[i]);
            }
            return numbers;
        }
        else {
            return new double[] {Double.NaN};
        }
    }
}
