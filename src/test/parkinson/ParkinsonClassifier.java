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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;
import jeco.algorithm.ga.SimpleGeneticAlgorithm;
import jeco.algorithm.moge.AbstractProblemGE;
import jeco.algorithm.moge.Phenotype;
import jeco.operator.comparator.SimpleDominance;
import jeco.operator.crossover.SinglePointCrossover;
import jeco.operator.evaluator.AbstractPopEvaluator;
import jeco.operator.evaluator.AbstractPopEvaluator;
import jeco.operator.mutation.IntegerFlipMutation;
import jeco.operator.selection.BinaryTournament;
import jeco.problem.Solution;
import jeco.problem.Solutions;
import jeco.problem.Variable;
import jeco.util.compiler.MyCompiler;
import jeco.util.compiler.MyLoader;
import jeco.util.logger.JecoLogger;

public class ParkinsonClassifier extends AbstractProblemGE {
    
    private static final Logger logger = Logger.getLogger(ParkinsonClassifier.class.getName());
    
    protected int threadId;
    protected MyCompiler compiler;
    protected DataTable dataTable;
    protected Properties properties;
    protected AbstractPopEvaluator evaluator;
    
    public ParkinsonClassifier(Properties properties, int threadId) throws IOException {
        super(properties.getProperty("BnfPathFile"), 1);
        this.properties = properties;
        this.threadId = threadId;
        compiler = new MyCompiler(properties);
        dataTable = new DataTable(this, properties.getProperty("DataPath"), Integer.valueOf(properties.getProperty("IdxBegin", "-1")), Integer.valueOf(properties.getProperty("IdxEnd", "-1")));
    }
    
    
    @Override
    public void evaluate(Solutions<Variable<Integer>> solutions) {
        StringBuilder currentJavaFile = new StringBuilder();
        int[][] patientsIdXs = dataTable.getPatientsIdXs();
        
        currentJavaFile.append("public class PopEvaluator").append(threadId).append(" extends jeco.operator.evaluator.AbstractPopEvaluator {\n\n");
        currentJavaFile.append("\tdouble[] var0 ={0.0};\n");
        currentJavaFile.append("\tdouble[] var1 ={1.0};\n");
        currentJavaFile.append("\tdouble[] var2 ={2.0};\n");
        currentJavaFile.append("\tdouble[] var3 ={3.0};\n");
        currentJavaFile.append("\tdouble[] var4 ={4.0};\n");
        currentJavaFile.append("\tdouble[] var5 ={5.0};\n");
        currentJavaFile.append("\n");
        
        currentJavaFile.append("public double MyDrv(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\t\treturn (data[allIndexes[1]] - data[allIndexes[0]])/(allIndexes[1] - allIndexes[0] + 1);\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double MySum(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\t\tdouble res = 0.0;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\t\tfor (int i = allIndexes[0]; i <= allIndexes[1]; ++i) {\n");
        currentJavaFile.append("\t\t\tres += data[i];\n");
        currentJavaFile.append("\t\t}\n");
        currentJavaFile.append("\t\treturn res;\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double MyPod(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble mypod = 1.0;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tfor (int i = allIndexes[0]; i <= allIndexes[1]; i++) {\n");
        currentJavaFile.append("\tmypod *= data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn mypod;\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double MyAvg(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble res = 0.0;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tres = MySum(allIndexes[0], allIndexes[1], array) / (allIndexes[1] - allIndexes[0] + 1);\n");
        currentJavaFile.append("\treturn res;\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double MyGeoAvg(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble mygeoavg = 0.0;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tmygeoavg = Math.pow(MyPod(allIndexes[0], allIndexes[1], array), 1/(allIndexes[1]-allIndexes[0]+1));\n");
        currentJavaFile.append("\treturn mygeoavg;\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double MyStd(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble mystd;\n");
        currentJavaFile.append("\tdouble res = 0.0;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tdouble[] varPow = new double[allIndexes[1]-allIndexes[0]+1];\n");
        currentJavaFile.append("\tfor(int i=allIndexes[0]; i<=allIndexes[1]; i++){\n");
        currentJavaFile.append("\tvarPow[i] = Math.pow(data[i], 2);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\tfor (int i = allIndexes[0]; i <= allIndexes[1]; i++) {\n");
        currentJavaFile.append("\tres += varPow[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t\n");
        currentJavaFile.append("\tmystd = (res/(allIndexes[1]-allIndexes[0]))-MyAvg(allIndexes[0], allIndexes[1], array);\n");
        currentJavaFile.append("\treturn mystd;\n");
        currentJavaFile.append("\t}\n");
        
        
        currentJavaFile.append("public double[] MyConv(int idx1_1, int idx2_1, int idx1_2, int idx2_2, double[] array1, double[] array2) {\n");
        currentJavaFile.append("\tint[] allIndexes1 = calculateIndexes(idx1_1, idx1_2, array1);\n");
        currentJavaFile.append("\tint[] allIndexes2 = calculateIndexes(idx2_1, idx2_2, array2);\n");
        currentJavaFile.append("\tdouble[] x = (array1.length > 1) ? array1 : getData((int)array1[0], allIndexes1[0], allIndexes1[1]);\n");
        currentJavaFile.append("\tdouble[] h = (array2.length > 1) ? array2 : getData((int)array2[0], allIndexes2[0], allIndexes2[1]);\n");
        currentJavaFile.append("\tdouble[] myconv = new double[(allIndexes1[1]-allIndexes1[0])+(allIndexes2[1]-allIndexes2[0])-1];\n");
        currentJavaFile.append("\tfor (int i = 0; i < (allIndexes1[1]-allIndexes1[0])+(allIndexes2[1]-allIndexes2[0])-1; i++ )	{\n");
        currentJavaFile.append("\tmyconv[i] = 0;                       // set to zero before sum\n");
        currentJavaFile.append("\tfor (int j = 0; j < (allIndexes2[1]-allIndexes2[0]); j++ ) {\n");
        currentJavaFile.append("\tif ((j < i) && ((i-j) < (allIndexes1[1]-allIndexes1[0]))) {\n");
        currentJavaFile.append("\tmyconv[i] += x[i - j] * h[j];    // convolve: multiply and accumulate\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn myconv;\n");
        currentJavaFile.append("\t}\n");
        
        
        currentJavaFile.append("public double MyMax(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble mymax = Double.NEGATIVE_INFINITY;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tfor(int i=allIndexes[0]; i<=allIndexes[1]; i++){\n");
        currentJavaFile.append("\tif (data[i] > mymax) {\n");
        currentJavaFile.append("\tmymax = data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("}\n");
        currentJavaFile.append("\treturn mymax;\n");
        currentJavaFile.append("}\n");
        
        currentJavaFile.append("\tpublic double MyTotalVar(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble mytotalvar = 0.0;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tdouble[] derv = new double[allIndexes[1]-allIndexes[0]-1];\n");
        currentJavaFile.append("\tfor(int i=allIndexes[0]; i<=allIndexes[1]-1; i++){\n");
        currentJavaFile.append("\tderv[i] = data[i+1]-data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\tfor (int i = allIndexes[0]; i <= allIndexes[1]; i++) {\n");
        currentJavaFile.append("\tmytotalvar += Math.abs(derv[i]);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn mytotalvar;\n");
        currentJavaFile.append("\t}\n");
        
        
        currentJavaFile.append("public double MyMin(int idx1, int idx2, double[] array) {\n");
        currentJavaFile.append("\tdouble mymin = Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\tint[] allIndexes = calculateIndexes(idx1, idx2, array);\n");
        currentJavaFile.append("\tdouble[] data = (array.length > 1) ? array : getData((int)array[0], allIndexes[0], allIndexes[1]);\n");
        currentJavaFile.append("\tfor(int i=allIndexes[0]; i<=allIndexes[1]; i++){\n");
        currentJavaFile.append("\tif (data[i] < mymin) {\n");
        currentJavaFile.append("\tmymin = data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn mymin;\n");
        currentJavaFile.append("}\n");
        
        
        
        currentJavaFile.append("\tpublic void evaluateExpression(int idxExpr) {\n");
        currentJavaFile.append("\t\treturn;\n");
        currentJavaFile.append("\t}\n\n");
        
        currentJavaFile.append("public double[] getData(int idx, int from, int to) {\n");
        currentJavaFile.append("\tdouble[] data = null;\n");
        currentJavaFile.append("\tfor(int i=from; i<=to; i++){\n");
        currentJavaFile.append("\tdata[i] = getDataTable(idx, i);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn data;\n");
        currentJavaFile.append("}\n");
        
        currentJavaFile.append("public int[] calculateIndexes(int idx1, int idx2, double[] array){\n");
        currentJavaFile.append("\tint[] idxs = null;\n");
        currentJavaFile.append("\tif (idx1 < idx2){\n");
        currentJavaFile.append("\tidxs[0] = idx1;\n");
        currentJavaFile.append("\tidxs[1] = idx2;\n");
        currentJavaFile.append("\t} else {\n");
        currentJavaFile.append("\tidxs[0] = idx2;\n");
        currentJavaFile.append("\tidxs[1] = idx1;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\tif (array.length > 1){\n");
        currentJavaFile.append("\tidxs[0] = (int)Math.round(idxs[0]*array.length/100);\n");
        currentJavaFile.append("\tidxs[1] = (int)Math.round(idxs[1]*array.length/100);\n");
        currentJavaFile.append("\tidxs[2] = (int)array[0];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\tidxs[0] = (int)Math.round(idxs[0]*table.size()/100);\n");
        currentJavaFile.append("\tidxs[1] = (int)Math.round(idxs[1]*table.size()/100);\n");
        currentJavaFile.append("\tidxs[2] = 0;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn idxs;\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("\tpublic double evaluate(int idxExpr, int k) {\n");
        currentJavaFile.append("\t\tdouble result = 0.0;\n");
        currentJavaFile.append("\t\ttry {\n");
        
        currentJavaFile.append("\t\t\tswitch(idxExpr) {\n");
        for (int i = 0; i < solutions.size(); ++i) {
            currentJavaFile.append("\t\t\t\tcase ").append(i).append(":\n");
            Solution<Variable<Integer>> solution = solutions.get(i);
            Phenotype phenotype = generatePhenotype(solution);
            if (correctSol) {
                currentJavaFile.append("\t\t\t\t\tresult = ").append(phenotype.toString()).append(";\n");
            } else {
                currentJavaFile.append("\t\t\t\t\tresult = Double.POSITIVE_INFINITY;\n");
            }
            currentJavaFile.append("\t\t\t\t\tbreak;\n");
        }
        currentJavaFile.append("\t\t\t\tdefault:\n");
        currentJavaFile.append("\t\t\t\t\tresult = Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t\t\t}\n"); // End switch
        
        currentJavaFile.append("\t\t}\n"); // End try
        currentJavaFile.append("\t\tcatch (Exception ee) {\n");
        currentJavaFile.append("\t\t\t// System.err.println(ee.getLocalizedMessage());\n");
        currentJavaFile.append("\t\t\tresult = Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t\t}\n"); // End catch
        currentJavaFile.append("\t\tif(Double.isNaN(result)) {\n");
        currentJavaFile.append("\t\t\tresult = Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t\t}\n");
        currentJavaFile.append("\t\treturn result;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("}\n");
        // Compilation process:
        try {
            File file = new File(compiler.getWorkDir() + File.separator + "PopEvaluator" + threadId + ".java");
            BufferedWriter writer = new BufferedWriter(new FileWriter(file));
            writer.write(currentJavaFile.toString());
            writer.flush();
            writer.close();
            LinkedList<String> filePaths = new LinkedList<>();
            filePaths.add(file.getAbsolutePath());
            boolean sucess = compiler.compile(filePaths);
            if (!sucess) {
                logger.severe("Unable to compile, with errors:");
                logger.severe(compiler.getOutput());
            }
        } catch (Exception ex) {
            logger.severe(ex.getLocalizedMessage());
        }
        // Get clinical information
        ArrayList<double[]> clinicalTable = dataTable.getDataTable("clinical");
        //ArrayList<double[]> rawDataTable = dataTable.getDataTable("training");
        
        // And now we evaluate all the solutions with the compiled file:
        evaluator = null;
        double fitness;
        
        try {
            evaluator = (AbstractPopEvaluator) (new MyLoader(compiler.getWorkDir())).loadClass("PopEvaluator" + threadId).newInstance();
        } catch (Exception ex) {
            logger.severe(ex.getLocalizedMessage());
        }
        
        for (int i = 0; i < solutions.size(); ++i) {
            // For every patient
            for (int p = 0; p < clinicalTable.size(); p++) {
                logger.info("Solution: " + i + ", Patient: " + p);
                logger.info("Id1: " + patientsIdXs[p][0] + ", Id2: " + patientsIdXs[p][1]);
                evaluator.setDataTable((ArrayList<double[]>) dataTable.getDataTable("training", patientsIdXs[p][0], patientsIdXs[p][1]));
                
                Solution<Variable<Integer>> solution = solutions.get(i);
                fitness = dataTable.evaluate(evaluator, solution, p, i);
                if (Double.isNaN(fitness)) {
                    logger.info("I have a NaN number here");
                }
                solution.getObjectives().set(0, fitness);
            }
        }
    }
    
    
    @Override
    public void evaluate(Solution<Variable<Integer>> solution) {
        logger.severe("The solutions should be already evaluated. You should not see this message.");
    }
    
    @Override
    public void evaluate(Solution<Variable<Integer>> solution, Phenotype phenotype) {
        logger.severe("The solutions should be already evaluated. You should not see this message.");
    }
    
    @Override
    public ParkinsonClassifier clone() {
        ParkinsonClassifier clone = null;
        try {
            clone = new ParkinsonClassifier(properties, threadId + 1);
        } catch (IOException ex) {
            logger.severe(ex.getLocalizedMessage());
        }
        return clone;
    }
    
    public static Properties loadProperties(String propertiesFilePath) {
        Properties properties = new Properties();
        try {
            properties.load(new BufferedReader(new FileReader(new File(propertiesFilePath))));
            File clsDir = new File(properties.getProperty("WorkDir"));
            URLClassLoader sysloader = (URLClassLoader) ClassLoader.getSystemClassLoader();
            Class<URLClassLoader> sysclass = URLClassLoader.class;
            Method method = sysclass.getDeclaredMethod("addURL", new Class[]{URL.class});
            method.setAccessible(true);
            method.invoke(sysloader, new Object[]{clsDir.toURI().toURL()});
        } catch (Exception ex) {
            logger.severe(ex.getLocalizedMessage());
        }
        return properties;
    }
    
    public static void runGE(Properties properties, int threadId) {
        JecoLogger.setup(properties.getProperty("LoggerBasePath") + "_" + threadId + ".log", Level.parse(properties.getProperty("LoggerLevel")));
        
        ParkinsonClassifier problem = null;
        try {
            problem = new ParkinsonClassifier(properties, threadId);
        } catch (IOException ex) {
            logger.severe(ex.getLocalizedMessage());
        }
        // Second create the algorithm
        IntegerFlipMutation<Variable<Integer>> mutationOperator = new IntegerFlipMutation<>(problem, 1.0 / problem.reader.getRules().size());
        SinglePointCrossover<Variable<Integer>> crossoverOperator = new SinglePointCrossover<>(problem, SinglePointCrossover.DEFAULT_FIXED_CROSSOVER_POINT, SinglePointCrossover.DEFAULT_PROBABILITY, SinglePointCrossover.AVOID_REPETITION_IN_FRONT);
        SimpleDominance<Variable<Integer>> comparator = new SimpleDominance<>();
        BinaryTournament<Variable<Integer>> selectionOp = new BinaryTournament<>(comparator);
        SimpleGeneticAlgorithm<Variable<Integer>> algorithm = new SimpleGeneticAlgorithm<>(problem, Integer.valueOf(properties.getProperty("NumIndividuals")), Integer.valueOf(properties.getProperty("NumGenerations")), true, mutationOperator, crossoverOperator, selectionOp);
        algorithm.initialize();
        algorithm.execute();
    }
    
    public static void main(String[] args) {
        String propertiesFilePath = "test" + File.separator + ParkinsonClassifier.class.getSimpleName() + ".properties";
        int threadId = 1;
        if (args.length == 1) {
            propertiesFilePath = args[0];
        } else if (args.length >= 2) {
            propertiesFilePath = args[0];
            threadId = Integer.valueOf(args[1]);
        }
        Properties properties = loadProperties(propertiesFilePath);
        runGE(properties, threadId);
    }
}
