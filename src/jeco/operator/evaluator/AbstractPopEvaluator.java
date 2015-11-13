/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/
package jeco.operator.evaluator;

import java.util.ArrayList;

/**
 * @author José Luis Risco Martín
 * @author Josué Pagán Ortiz
 */
public abstract class AbstractPopEvaluator {
    
    protected ArrayList<double[]> dataTable;
    protected int[] dataLimits;

    public abstract void evaluateExpression(int idxExpr);
    public abstract double evaluate(int idxExpr, int k);
    
    public void setDataTable(ArrayList<double[]> dataTable) {
        this.dataTable = dataTable;
    }
    
    public void setDataLimits(int[] dataLimits) {
        this.dataLimits = dataLimits;
    }
    
    public ArrayList<double[]> getDataTable() {
        return dataTable;
    }
    
    public int[] getDataLimits(int ex, int f) {
        int[] limits = new int[2];
        limits[0] = dataLimits[4*ex+2*f];
        limits[1] = dataLimits[4*ex+2*f+1];
        return limits;
    }
    
    public double getDataTable(int idxVar, int k){
        if (k < 0) {            
           return dataTable.get(0)[idxVar];
        }
        else if (k >= dataTable.size()) {
            return dataTable.get(dataTable.size()-1)[idxVar];
        }
        else {
            return dataTable.get(k)[idxVar];
        }
    }
}
