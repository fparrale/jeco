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

    protected ArrayList<double[]> dataTable;

    public abstract void evaluateExpression(int idxExpr);
    public abstract double evaluate(int idxExpr, int k);

    public void setDataTable(ArrayList<double[]> dataTable) {
        this.dataTable = dataTable;
    }

    public ArrayList<double[]> getDataTable() {
        return dataTable;
    }
}
