/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jeco.util.classifier;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import test.parkinson.ParkinsonAdaBoostClassifier;

/**
 *
 * @author josueportiz
 *  TODO: hardcoded to Parkinson GE classifier. Fix to general classifier
 */


public class AdaBoost {
    protected ParkinsonAdaBoostClassifier geClassifier;
    protected ClassifierEvaluator confMatrix;
    protected ArrayList<double[]> originalData; // Original dataset
    protected ArrayList<double[]> weightedData; // Current weighted dataset
    
    protected int depth; // Depth of algortihm
    protected int numSamples; // Num of elements in the dataset
    protected int numFeatures; // Number of features
    
    protected double[][] d;

    
    // Constructor:
    //public AdaBoost(int k, int n, int f) {
    public AdaBoost(ParkinsonAdaBoostClassifier classifier, ClassifierEvaluator cm, ArrayList<double[]> data, int k) {
        this.geClassifier = classifier;
        this.confMatrix = cm;
        this.originalData = data;
        this.weightedData = data;
        this.depth = k;
        this.numSamples = data.size();
        this.numFeatures = data.get(0).length;

        // Init parameters
        d = new double[numSamples][numFeatures]; // Initialize uniform importance to each example
        for (int n=0; n<numSamples; n++){
            for (int f=0; f<numFeatures; f++){
                d[n][f] = 1.0/(double)(numSamples);
            }            
        }
    }
    
    public ArrayList<double[]> weightDataSet(double[][] w, ArrayList<double[]> data){
        ArrayList<double[]> newData = new ArrayList<>();
        Iterator<double[]> itr = data.iterator();

        int n = 0;
        while (itr.hasNext()){
            double[] row = itr.next();
            for (int f=0; f<data.get(0).length; f++){
                row[f] = row[f]*w[n][f];            
            } 
            newData.add(n, row);
            n++;
        }
        return newData;
    }
    
    public double computeError(double[][] d, double[] x, double[] y){
        int[] booleanResult = new int[y.length];
        double error = 0.0;
        
        // Compute boolean result of the classification
        // Compute the error
        
        for (int i=0; i<=y.length; i++){
            if(x[i]==y[i]){
                booleanResult[i] = 1;
            } else {
                booleanResult[i] = -1;
            }
            error += error
        }
        
        
        return error;        
    }
    
    public void run() throws IOException{
        ArrayList<double[]> wData = weightedData;
        
        for (int k=1; k<=depth; k++){
            wData = weightDataSet(d, wData); // Weight the data
            geClassifier.runLearner(wData); // Train kth classifier on weighted data and make predictions on training data
            computeError(d, confMatrix.getOriginal(), confMatrix.getPredicted()); // Compute weighted training error
            weightedData = wData; // Update values
        }
    }
    
    
}
