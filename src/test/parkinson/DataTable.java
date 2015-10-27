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
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import jeco.operator.evaluator.AbstractPopEvaluator;
import jeco.problem.Solution;
import jeco.problem.Variable;

/**
 * Class to manage a normalized data table. Originally, the data table is passed
 * to this class as a regular data table. After the constructor, the data table
 * is normalized in the interval [1,2].
 *
 * @author José Luis Risco Martín
 * @author Josué Pagán Ortiz
 */
public class DataTable {
    
    private static final Logger logger = Logger.getLogger(DataTable.class.getName());
    
    protected ParkinsonClassifier problem;
    protected String baseTrainingPath = null;
    protected ArrayList<double[]> trainingTable = new ArrayList<>();
    protected ArrayList<double[]> subTrainingTable = new ArrayList<>();
    protected int idxBegin = -1;
    protected int idxEnd = -1;
    protected int numInputColumns = 0;
    protected int numTotalColumns = 0;
    protected double[] xLs = null;
    protected double[] xHs = null;
    
    protected double bestFitness = Double.POSITIVE_INFINITY;
    
    protected String foot = null;
    protected Double patientPDLevel = null;
    protected ArrayList<double[]> clinicalTable = new ArrayList<>();
    protected int lengthIni = 0;
    protected int lengthEnd = 0;
    protected int[][] patientsIdXs = new int[clinicalTable.size()][2];
    protected double[] error = new double[92];
    protected double cumulatedFitness = 0.0;
    
    protected FileHandler fh;  

 
    
    
    public DataTable(ParkinsonClassifier problem, String baseTrainingPath, int idxBegin, int idxEnd) throws IOException {
        this.problem = problem;
        this.baseTrainingPath = baseTrainingPath;
        logger.info("Reading data file ...");
        
        readData(problem.properties.getProperty("ClinicalPath"), clinicalTable, false);
        
        fillTrainingDataTable(trainingTable);
        this.idxBegin = (idxBegin == -1) ? 0 : idxBegin;
        this.idxEnd = (idxEnd == -1) ? trainingTable.size() : idxEnd;
        logger.info("Evaluation interval: [" + this.idxBegin + "," + this.idxEnd + ")");
        logger.info("... done.");
        
        //try {
            // This block configure the logger with handler and formatter
        //    fh = new FileHandler("/tmp/MyLogFile.log");
        //    logger.addHandler(fh);
        //   SimpleFormatter formatter = new SimpleFormatter();
        //    fh.setFormatter(formatter);
        //} catch (SecurityException | IOException e) {
        //}
    }
    
    public DataTable(ParkinsonClassifier problem, String baseTrainingPath) throws IOException {
        this(problem, baseTrainingPath, -1, -1);
    }
    
    
    public final void fillTrainingDataTable(ArrayList<double[]> dataTable) throws IOException {
// exercises : exercises = {walk, cycling, hoolToe};
        String exercises = problem.properties.getProperty("Exercises");
        String[] exercisesTrunc = exercises.split(",");
        
        numInputColumns = 0;
        numTotalColumns = 0;
        patientsIdXs = new int[clinicalTable.size()][2];
        
        for (int p = 0; p < clinicalTable.size(); p++) {
            String patientID = String.valueOf((int)clinicalTable.get(p)[1]);               // Get the code GAxxxxxx
            patientPDLevel = clinicalTable.get(p)[8];              // Get the level scale H&Y
            logger.info("PatientID: GA" + patientID + ", PDlevel: " + patientPDLevel);
            
            String absoluteBasePath = problem.properties.getProperty("DataPathBase");
            
            lengthIni = trainingTable.size();
            for (int f = 0; f<=1; f++){	// For each foot
                foot = (f == 0) ? "RightFoot_" : "LeftFoot_";
                
                for (String ex : exercisesTrunc) {
                    String absoluteDataPath = absoluteBasePath + "/GA" + patientID + "/" + foot + ex + ".csv";
                    //logger.info("Data: " + absoluteDataPath);
                    readData(absoluteDataPath, trainingTable, true);
                }
            }
            if (lengthIni < trainingTable.size()-1){
                patientsIdXs[p][0] = lengthIni;
                patientsIdXs[p][1] = trainingTable.size()-1;
            }
            else {
                System.out.println("No data for Patient: GA" + patientID);
                patientsIdXs[p][0] = lengthIni;
                patientsIdXs[p][1] = lengthIni;
            }
        }
    }
    
    
    public final void readData(String dataPath, ArrayList<double[]> dataTable, Boolean addOutputLine) throws IOException {
        File file = new File(dataPath);
        if (file.exists()){
            
            try (BufferedReader reader = new BufferedReader(new FileReader(new File(dataPath)))) {
                String line;
                
                while ((line = reader.readLine()) != null) {
                    if (line.isEmpty() || line.startsWith("#")) {
                        continue;
                    }
                    String[] parts = line.split(";");
                    if (parts.length == 1){
                        parts = line.split(",");
                    }
                    if (parts.length > numInputColumns) {
                        numInputColumns = parts.length;
                        numTotalColumns = numInputColumns + 1;
                    }
                    
                    double[] dataLine = new double[numTotalColumns];
                    for (int j = 0; j < numInputColumns; j++) {
                        dataLine[j] = Double.valueOf(parts[j]);
                    }
                    if (addOutputLine) {
                        dataLine[numTotalColumns-1] = patientPDLevel;
                    }
                    dataTable.add(dataLine);
                }
                reader.close();
            }
        }
        else {
            logger.info("File: " + dataPath + " DOES NOT EXIST");
        }
    }
        
    
    public ArrayList<double[]> getDataTable(String  type) {
        switch (type) {
            case "training":
                return trainingTable;
            case "clinical":
                return clinicalTable;
            default:
                return trainingTable;
        }
    }
    
    public ArrayList<double[]> getDataTable(String  type, int idx1, int idx2) {
        switch (type) {
            case "training":
                subTrainingTable = new ArrayList(trainingTable.subList(idx1, idx2));
                return new ArrayList(trainingTable.subList(idx1, idx2));
            case "clinical":
                return new ArrayList(clinicalTable.subList(idx1, idx2));
            default:                
                return trainingTable;
        }
    }

    public int[][] getPatientsIdXs(){
        return patientsIdXs;
    }
    
    public double evaluate(AbstractPopEvaluator evaluator, Solution<Variable<Integer>> solution, int patientNo, int idx) {
        String functionAsString = problem.generatePhenotype(solution).toString();
        if (patientNo == 0.0) { 
            cumulatedFitness = 0.0;
        }
        error[patientNo] = computeError(evaluator, patientNo, idx);
        cumulatedFitness += error[patientNo]/(Integer.valueOf(problem.properties.getProperty("MaxPDLevel"))*clinicalTable.size());

        if (patientNo == clinicalTable.size()-1) {            
            if (cumulatedFitness < bestFitness) {
                bestFitness = cumulatedFitness;
                logger.info("Best FIT=" + (100 * (1 - bestFitness)) + "; Expresion=" + functionAsString);
            }   
        }
        return cumulatedFitness;
    }
    
    public double computeError(AbstractPopEvaluator evaluator, int patientNo, int idx) {
        double resultGE =  evaluator.evaluate(idx, -1);
        
        double qResult = quantizer(resultGE);
        //logger.info("Solution, " + idx + ", patient, " + patientNo + ", ResultGE, " + resultGE + ", qResult, " + qResult);

        // Get the PD H&Y level and compute fitness
        double dif = Math.abs(qResult-clinicalTable.get(patientNo)[Integer.valueOf(problem.properties.getProperty("PDLevelCol"))]);
        return dif;
    }
    public double quantizer(Double currFitness) {
        // Hardcode H&Y Parkinson Scale 0 to 3 (0 means no PD)
        double qFitness = 0.0;
        
        if (currFitness >= 4.5) {
            qFitness = 5.0;
        } else if (currFitness >= 3.5) {
            qFitness = 4.0;
        } else if (currFitness >= 2.5) {
            qFitness = 3.0;
        } else if (currFitness >= 1.5){
            qFitness = 2.0;
        } else if (currFitness >= 0.5){
            qFitness = 2.0;
        } else {
            qFitness = 0.0;
        }
        return qFitness;
    }
    
}