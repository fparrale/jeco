package jeco.util;

/**
 * Class to manage the confusion matrix and metric parameters of a classifier.
 *
 * @author José Luis Risco Martín
 * @author Josué Pagán Ortiz
 */
public class ClassifierEvaluator {
    protected int[][] confusionMatrix;
    
    public ClassifierEvaluator(int c) {
        this.confusionMatrix = new int[c][c];
    }
    
    public void resetConfusionMatrix() {
        for (int i=0; i<confusionMatrix.length; i++) {
            for (int j=0; j<confusionMatrix.length; j++) {
                confusionMatrix[i][j] = 0;
            }
        }
    }
    
    public void setConfusionMatrix(int[][] cm) {
        confusionMatrix = cm;
    }

    public void setValue(int originalClass, int classifiedClass, int v) {
        confusionMatrix[classifiedClass][originalClass] += v;
    }
    
    public int[][] getConfusionMatrix() {
        return confusionMatrix;
    }

    public int getN() {
        int n = 0;
        for (int i=0; i < confusionMatrix.length; i++) {
            for (int j=0; j < confusionMatrix.length; j++) {
                n += confusionMatrix[i][j];
            }
        }
        return n;
    }
    
    // Calculate the OSR (Overall Success Rate)
    public double getClassificationRate() {
        double osr = 0.0;
        for (int i=0; i < confusionMatrix.length; i++) {
            osr += confusionMatrix[i][i];
        }        
        return osr/getN();
    }
    
    // Marginal rates functions
    public double getSensitivity(int classC) { // Or TPR
        return getTruePositives(classC)/(getTruePositives(classC)+getFalseNegatives(classC));
    }
    
     public double getSpecificity(int classC) {  // Or TNR
        return getTrueNegatives(classC)/(getTrueNegatives(classC)+getFalsePositives(classC));
    }
     
      public double getPrecision(int classC) { // Or PPV
        return getTruePositives(classC)/(getTruePositives(classC)+getFalsePositives(classC));
    }
    
    // Marginal rates
    // Marginal true positives
    public int getTruePositives(int classC) {
        return confusionMatrix[classC][classC];
    }
    
    // Marginal false positives
    public int getFalsePositives(int classC) {
        int niC = 0;
        for (int i=0; i < confusionMatrix.length; i++) {
            niC += confusionMatrix[classC][i];
        }  
        return niC-confusionMatrix[classC][classC];
    }
    
    // Marginal false negatives
    public int getFalseNegatives(int classC) {
        int nCi = 0;
        for (int i=0; i < confusionMatrix.length; i++) {
            nCi += confusionMatrix[i][classC];
        }  
        return nCi-confusionMatrix[classC][classC];
    }
    
    // Marginal true negatives
    public int getTrueNegatives(int classC) {
        return getN()-getTruePositives(classC)-getFalsePositives(classC)-getFalseNegatives(classC);
    }
}