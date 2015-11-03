package jeco.util;

/**
 * Class to manage the confusion matrix and metric parameters of a classifier.
 *
 * @author José Luis Risco Martín
 * @author Josué Pagán Ortiz
 */
public class classifierEvaluator {
    protected int[][] confusionMatrix;
    
    public classifierEvaluator(int C) {
        this.confusionMatrix = new int[C][C];
    }
    
    public void setValue(int originalClass, int classifiedClass, int V) {
        confusionMatrix[classifiedClass][originalClass] = V;
    }
    
    public void setValue(int[][] CM) {
        confusionMatrix = CM;
    }
     
    public int[][] getConfusionMatrix() {
        return confusionMatrix;
    }
    
    // Calculate the OSR (Overall Success Rate)
    public double getClassificationRate() {
        double osr = 0.0;
        for (int i=0; i < confusionMatrix.length; i++) {
            osr += confusionMatrix[i][i];
        }        
        return osr/confusionMatrix.length;
    }
    
    // Marginal rates functions
    public double getSensitivity(int classC) { // Or TPR
        return getTruePositives(classC)/(getTruePositives(classC)+getFalseNegatives(classC));
    }
    
     public double getSpecificity(int classC) {  // Or TNR
        return getTrueNegatives(classC)/(getTrueNegatives(classC)+getFalsePositives(classC));
    }
     
      public double getPrecission(int classC) { // Or PpV
        return getTruePositives(classC)/(getTruePositives(classC)+getFalsePositives(classC));
    }
    
    // Marginal rates
    // Marginal true positives
    public double getTruePositives(int classC) {
        return confusionMatrix[classC][classC];
    }
    
    // Marginal false positives
    public double getFalsePositives(int classC) {
        double niC = 0.0;
        for (int i=0; i < confusionMatrix.length; i++) {
            niC += confusionMatrix[classC][i];
        }  
        return niC-confusionMatrix[classC][classC];
    }
    
    // Marginal false negatives
    public double getFalseNegatives(int classC) {
        double nCi = 0.0;
        for (int i=0; i < confusionMatrix.length; i++) {
            nCi += confusionMatrix[i][classC];
        }  
        return nCi-confusionMatrix[classC][classC];
    }
    
    // Marginal true negatives
    public double getTrueNegatives(int classC) {
        double n = 0.0;
        for (int i=0; i < confusionMatrix.length; i++) {
            for (int j=0; j < confusionMatrix.length; j++) {
                n += confusionMatrix[i][j];
            }
        }  
        return n-getTruePositives(classC)-getFalsePositives(classC)-getFalseNegatives(classC);
    }
}