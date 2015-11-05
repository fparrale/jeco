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
import java.util.Random;
import java.util.logging.Logger;

/**
 * Class to manage a data table. The data table is passed
 * to this class as a regular data table.
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
    
    protected String foot = null;
    protected Double patientPDLevel = null;
    protected ArrayList<double[]> clinicalTable = new ArrayList<>();
    protected int lengthIni = 0;
    protected int lengthEnd = 0;
    protected int[][] limitMarkers;
    protected int[][] patientsIdxs;
    protected String exercises;
    protected String[] exercisesTrunc;
    
    public DataTable(ParkinsonClassifier problem, String baseTrainingPath, int idxBegin, int idxEnd) throws IOException {
        this.problem = problem;
        this.baseTrainingPath = baseTrainingPath;
        logger.info("Reading data file ...");
        
        readData(problem.properties.getProperty("ClinicalPath"), clinicalTable, false);

        this.exercises = problem.properties.getProperty("Exercises");
        this.exercisesTrunc = exercises.split(",");
        this.limitMarkers = new int[clinicalTable.size()][2*exercisesTrunc.length];
     
        fillTrainingDataTable(trainingTable);
        this.idxBegin = (idxBegin == -1) ? 0 : idxBegin;
        this.idxEnd = (idxEnd == -1) ? trainingTable.size() : idxEnd;
        logger.info("Evaluation interval: [" + this.idxBegin + "," + this.idxEnd + ")");
        logger.info("... done.");
    }
    
    public DataTable(ParkinsonClassifier problem, String baseTrainingPath) throws IOException {
        this(problem, baseTrainingPath, -1, -1);
    }
    
    public final void fillTrainingDataTable(ArrayList<double[]> dataTable) throws IOException {
        numInputColumns = 0;
        numTotalColumns = 0;
        
        String absoluteBasePath = problem.properties.getProperty("DataPathBase");
        
        for (int p = 0; p < clinicalTable.size(); p++) {
            String patientID = String.valueOf((int)clinicalTable.get(p)[1]);               // Get the code GAxxxxxx
            patientPDLevel = clinicalTable.get(p)[8];              // Get the level scale H&Y
            logger.info("PatientID: GA" + patientID + ", PDlevel: " + patientPDLevel);
            
            for (int ex = 0; ex < exercisesTrunc.length; ex++) { // For each exercise
                lengthIni = trainingTable.size();
                
                for (int f = 0; f<=1; f++){	// For each foot
                    foot = (f == 0) ? "RightFoot_" : "LeftFoot_";
                    
                    String absoluteDataPath = absoluteBasePath + "/GA" + patientID + "/" + foot + exercisesTrunc[ex] + ".csv";
                    //logger.info("Data: " + absoluteDataPath);
                    readData(absoluteDataPath, trainingTable, true);
                }
                
                // Store indexes: from-to for each exercise
                if (lengthIni < trainingTable.size()-1){
                    limitMarkers[p][2*ex] = lengthIni;
                    limitMarkers[p][2*ex+1] = trainingTable.size()-1;
                }
                else {
                    System.out.println("No data for Patient: GA" + patientID);
                    limitMarkers[p][2*ex] = -1;
                    limitMarkers[p][2*ex+1] = -1;
                }
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

    public int[][] getLimitMarkers(){
        return limitMarkers;
    }

    public int[][] getPatientsIdXs(boolean crossVal){
        // Check N fold cross-validation
        if (crossVal) {
            patientsIdxs = randomizeDataSelection(clinicalTable.size(), Integer.valueOf(problem.properties.getProperty("N")), true);
        } else {
            patientsIdxs = randomizeDataSelection(clinicalTable.size(), 1, false);
        }        
        return patientsIdxs;
    }
    
    public int[][] randomizeDataSelection(int elements, int groups, boolean randomize){
        int[][] randomTable = new int[groups][elements/groups];
        int[] elementsAvailable = new int[elements];
        
        for (int i=0; i<= elements-1; i++){
            elementsAvailable[i] = i;
        }
        if (randomize) {
            // Implementing Fisher–Yates shuffle
            Random rnd = new Random();
            for (int i = elementsAvailable.length - 1; i > 0; i--) {
                int index = rnd.nextInt(i + 1);
                // Simple swap
                int a = elementsAvailable[index];
                elementsAvailable[index] = elementsAvailable[i];
                elementsAvailable[i] = a;
            }
        }
        
        for (int i=0; i<= groups-1; i++){
            for (int j=0; j<= (elements/groups)-1; j++){
                randomTable[i][j] = elementsAvailable[j+(i*elements/groups)];
            }
        }
        return randomTable;
    }
}