/*
 * Copyright (C) 2010-2015 José Luis Risco Martín <jlrisco@ucm.es> and 
 * José Manuel Colmenar Verdugo <josemanuel.colmenar@urjc.es>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *  - Josué Pagán Ortíz
 *  - José Luis Risco Martín
 */
package test.parkinson;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Logger;
import jeco.algorithm.moge.AbstractProblemGE;
import jeco.operator.evaluator.AbstractPopEvaluator;
import jeco.problem.Solution;
import jeco.problem.Variable;

/**
 * Class to manage a normalized data table. Originally, the data table is passed
 * to this class as a regular data table. After the constructor, the data table
 * is normalized in the interval [1,2].
 *
 * @author José Luis Risco Martín
 */
public class DataTable {

    private static final Logger logger = Logger.getLogger(DataTable.class.getName());

    protected ParkinsonClassifier problem;
    protected String trainingPath = null;
    protected ArrayList<double[]> trainingTable = new ArrayList<>();
    protected int idxBegin = -1;
    protected int idxEnd = -1;
    protected int numInputColumns = 0;
    protected int numTotalColumns = 0;
    protected double[] xLs = null;
    protected double[] xHs = null;

    protected double bestFitness = Double.POSITIVE_INFINITY;

    public DataTable(ParkinsonClassifier problem, String trainingPath, int idxBegin, int idxEnd) throws IOException {
        this.problem = problem;
        this.trainingPath = trainingPath;
        logger.info("Reading data file ...");
        fillDataTable(trainingPath, trainingTable);
        this.idxBegin = (idxBegin == -1) ? 0 : idxBegin;
        this.idxEnd = (idxEnd == -1) ? trainingTable.size() : idxEnd;
        logger.info("Evaluation interval: [" + this.idxBegin + "," + this.idxEnd + ")");
        logger.info("... done.");
    }

    public DataTable(ParkinsonClassifier problem, String trainingPath) throws IOException {
        this(problem, trainingPath, -1, -1);
    }

    public final void fillDataTable(String dataPath, ArrayList<double[]> dataTable) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(new File(dataPath)));
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] parts = line.split(";");
            if (parts.length > numInputColumns) {
                numInputColumns = parts.length;
                numTotalColumns = numInputColumns + 1;
            }
            double[] dataLine = new double[numTotalColumns];
            for (int j = 0; j < numInputColumns; ++j) {
                dataLine[j] = Double.valueOf(parts[j]);
            }
            dataTable.add(dataLine);
        }
        reader.close();
    }

    public double evaluate(AbstractPopEvaluator evaluator, Solution<Variable<Integer>> solution, int idx) {
        String functionAsString = problem.generatePhenotype(solution).toString();
        double fitness = computeFitness(evaluator, idx);
        if (fitness < bestFitness) {
            bestFitness = fitness;
            for (int i = 0; i < numTotalColumns; ++i) {
                if (i == 0) {
                    functionAsString = functionAsString.replaceAll("getVariable\\(" + i + ",", "yr\\(");
                } else if (i == numTotalColumns - 1) {
                    functionAsString = functionAsString.replaceAll("getVariable\\(" + i + ",", "yp\\(");
                } else {
                    functionAsString = functionAsString.replaceAll("getVariable\\(" + i + ",", "u" + i + "\\(");
                }
            }
            logger.info("Best FIT=" + (100 * (1 - bestFitness)) + "; Expresion=" + functionAsString);
        }
        return fitness;
    }

    public double computeFitness(AbstractPopEvaluator evaluator, int idx) {
        double resultGE =  evaluator.evaluate(idx, -1);
        
        double qResult = quantizer(resultGE); 
        
        double fitness = Math.abs(qResult-Double.parseDouble(problem.properties.getProperty("ParkinsonLevel")));
        return fitness;
    }

    public ArrayList<double[]> getDataTable() {
        return trainingTable;
    }
    
    public double quantizer(double currFitness) {
        // Hardcode H&Y Parkinson Scale 0 to 3 (0 means no PD)
        double qFitness = 0;
        
        if (currFitness >= 2.5) {
            qFitness = 3;
        } else if (currFitness >= 1.5){
            qFitness = 2;
        } else if (currFitness >= 0.5){
            qFitness = 2;
        } else {
            qFitness = 0;
        }
        return qFitness;
    }

}
