/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/
package jeco.operator.evaluator;

import java.util.ArrayList;

/**
 *
 * @author José Luis Risco Martín
 */
public abstract class AbstractPopEvaluator {
    
    protected ArrayList<double[]> table;
    
    public abstract void evaluateExpression(int idxExpr);
    public abstract double evaluate(int idxExpr, int k);
    
    public void setDataTable(ArrayList<double[]> dataTable) {
        this.table = dataTable;
    }
    
    public ArrayList<double[]> getDataTable() {
        return table;
    }
    
    public double getDataTable(int idxVar, int k){
        if (k < 0) {
            return table.get(0)[idxVar];
        }
        else if (k > table.size() ) {
            return table.get(table.size())[idxVar];
        }
        else {
            return table.get(k)[idxVar];
        }
    }
}
