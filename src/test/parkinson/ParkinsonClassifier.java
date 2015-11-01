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
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;
import jeco.algorithm.ga.SimpleGeneticAlgorithm;
import jeco.algorithm.moge.AbstractProblemGE;
import jeco.algorithm.moge.Phenotype;
import jeco.operator.comparator.SimpleDominance;
import jeco.operator.crossover.SinglePointCrossover;
import jeco.operator.evaluator.AbstractPopEvaluator;
import jeco.operator.mutation.IntegerFlipMutation;
import jeco.operator.selection.BinaryTournament;
import jeco.optimization.threads.MasterWorkerThreads;
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
        currentJavaFile.append("\t\n");
        currentJavaFile.append("\tint[] ex0 ={0};\n");
        currentJavaFile.append("\tint[] ex1 ={1};\n");
        currentJavaFile.append("\tint[] ex2 ={2};\n");
        currentJavaFile.append("\tint[] ex3 ={3};\n");
        currentJavaFile.append("\tint[] ex4 ={4};\n");
        currentJavaFile.append("\tint[] noEx ={-1};\n");


        /**
         * Implementation of functions of the Grammar
         * */
               
        currentJavaFile.append("public double MyAvg(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble res = 0.0;\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tint[] limits = {0, data.length-1};\n");
        currentJavaFile.append("\tres = MySum(data, limits)/data.length;\n");
        currentJavaFile.append("\treturn res;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");   
    
        currentJavaFile.append("public double MySum(double[] array, int[] ex) {\n");
        currentJavaFile.append("\t\tdouble res = 0.0;\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\t\tfor (int i = 0; i <= data.length-1; i++) {\n");
        currentJavaFile.append("\t\t\tres += data[i];\n");
        currentJavaFile.append("\t\t}\n");
        currentJavaFile.append("\t\treturn res;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");        
        currentJavaFile.append("\t}\n");
               
       
        currentJavaFile.append("public double MyMax(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tdouble mymax = Double.NEGATIVE_INFINITY;\n");
        currentJavaFile.append("\tfor(int i=0; i<=data.length-1; i++){\n");
        currentJavaFile.append("\tif (data[i] > mymax) {\n");
        currentJavaFile.append("\tmymax = data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("}\n");
        currentJavaFile.append("\treturn mymax;\n");
        currentJavaFile.append("}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("}\n");
        
        currentJavaFile.append("public double MyMin(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tdouble mymin = Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\tfor(int i=0; i<=data.length-1; i++){\n");
        currentJavaFile.append("\tif (data[i] < mymin) {\n");
        currentJavaFile.append("\tmymin = data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn mymin;\n");
        currentJavaFile.append("}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("\tpublic double MyStd(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tint[] limits = {0, data.length-1};\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tdouble mystd;\n");
        currentJavaFile.append("\tdouble[] res = new double[data.length];\n");
        currentJavaFile.append("\tdouble avg = MyAvg(data, limits);\n");
        currentJavaFile.append("\tfor (int i = 0; i <= data.length-1; i++) {\n");
        currentJavaFile.append("\tres[i] = Math.pow(data[i] - avg, 2);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\tmystd = Math.pow(MyAvg(res, limits), 0.5);\n");
        currentJavaFile.append("\treturn mystd;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        
        currentJavaFile.append("\tpublic double MyTotalVar(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tdouble mytotalvar = 0.0;\n");
        currentJavaFile.append("\tdouble[] derv = new double[data.length-1];\n");
        currentJavaFile.append("\tfor(int i=0; i<=data.length-2; i++){\n");
        currentJavaFile.append("\tderv[i] = data[i+1]-data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\tfor (int i = 0; i <= derv.length-1; i++) {\n");
        currentJavaFile.append("\tmytotalvar += Math.abs(derv[i]);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn mytotalvar;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        
        currentJavaFile.append("public double MyPod(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");        
        currentJavaFile.append("\tdouble mypod = 1.0;\n");
        currentJavaFile.append("\tfor (int i = 0; i <= data.length-1; i++) {\n");
        currentJavaFile.append("\tmypod *= data[i];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn mypod;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double MyGeoAvg(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");   
        currentJavaFile.append("\tdouble mygeoavg = 0.0;\n");
        currentJavaFile.append("\tint[] limits = {0, data.length-1};\n");
        currentJavaFile.append("\tmygeoavg = Math.pow(MyPod(data, limits), 1/(data.length));\n");
        currentJavaFile.append("\treturn mygeoavg;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        
        currentJavaFile.append("\tpublic double[] MyPow(double[] array, int[] ex, double pow) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");           
        currentJavaFile.append("\tfor (int i = 0; i <= data.length-1; i++) {\n");
        currentJavaFile.append("\tdata[i] = Math.pow(data[i], pow);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn data;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn new double[] {Double.NaN};\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("public double[] MyConv(double[] array1, double[] array2, int[] ex1, int[] ex2) {\n");
        currentJavaFile.append("\tdouble[] x = getData(array1, ex1);\n");
        currentJavaFile.append("\tdouble[] h = getData(array2, ex2);\n");
        currentJavaFile.append("\tif ((!Double.isNaN(x[0])) && (!Double.isNaN(h[0]))){\n");           
        currentJavaFile.append("\tdouble[] myconv = new double[x.length + h.length - 1];\n");
        currentJavaFile.append("\tfor (int i = 0; i <= myconv.length -1; i++ )	{\n");
        currentJavaFile.append("\tmyconv[i] = 0;                       // set to zero before sum\n");
        currentJavaFile.append("\tfor (int j = 0; j <= h.length-1; j++ ) {\n");
        currentJavaFile.append("\tif ((j <= i) && ((i-j) < x.length)) {\n");
        currentJavaFile.append("\tmyconv[i] += x[i - j] * h[j];    // convolve: multiply and accumulate\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn myconv;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn new double[] {Double.NaN};\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("\tpublic double[] MyDiff(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tdouble[] diff = new double[data.length-1];\n");
        currentJavaFile.append("\tfor (int i = 1; i <= data.length-1; i++) {\n");
        currentJavaFile.append("\tdiff[i-1] = data[i] - data[i-1];\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn diff;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn new double[] {Double.NaN};\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t\n");

        currentJavaFile.append("\tpublic double[] MyAbs(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tdouble[] data = getData(array, ex);\n");
        currentJavaFile.append("\tif (!Double.isNaN(data[0])){\n");
        currentJavaFile.append("\tfor (int i = 0; i <= data.length-1; i++) {\n");
        currentJavaFile.append("\tdata[i] = Math.abs(data[i]);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn data;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn new double[] {Double.NaN};\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t\n");
        /**
         * Utils
         * */
        currentJavaFile.append("\tpublic double[] getData(double[] array, int[] ex) {\n");
        currentJavaFile.append("\tif (array.length > 1) {\n");
        currentJavaFile.append("\tdouble[] data = array;\n");
        currentJavaFile.append("\treturn data;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse if (Double.isNaN(array[0])) {\n");
        currentJavaFile.append("\treturn new double[] {Double.NaN};\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\tint[] allIndexes = getDataLimits(ex[0]);\n");
        currentJavaFile.append("\tif ((allIndexes[0] >= 0) && (allIndexes[1] >= 0) && (!Double.isNaN(array[0]))){\n");
        currentJavaFile.append("\tdouble[] data = new double[allIndexes[1]-allIndexes[0]+1];\n");
        currentJavaFile.append("\tfor(int i=0; i<=allIndexes[1]-allIndexes[0]; i++){\n");
        currentJavaFile.append("\tdata[i] = getDataTable((int)array[0], i + allIndexes[0]);\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\treturn data;\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\telse {\n");
        currentJavaFile.append("\treturn new double[] {Double.NaN};\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        currentJavaFile.append("\t}\n");
        
        currentJavaFile.append("\tpublic void evaluateExpression(int idxExpr) {\n");
        currentJavaFile.append("\t\treturn;\n");
        currentJavaFile.append("\t}\n\n");
        
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
        currentJavaFile.append("\t\t\tSystem.err.println(\"GE result is DefaultValue.\");\n");
        currentJavaFile.append("\t\t\t\t\tresult = Double.POSITIVE_INFINITY;\n");
        currentJavaFile.append("\t\t\t}\n"); // End switch
        
        currentJavaFile.append("\t\t}\n"); // End try
        currentJavaFile.append("\t\tcatch (Exception ee) {\n");
        currentJavaFile.append("\t\t\tSystem.err.println(ee.fillInStackTrace());\n");
        currentJavaFile.append("\t\t\tSystem.err.println(\"Exception trying to calculate the GE result.\");\n");
        currentJavaFile.append("\t\t\tresult = Double.NaN;\n");
        
        currentJavaFile.append("\t\t}\n"); // End catch
        currentJavaFile.append("\t\tif(Double.isNaN(result)) {\n");
        //currentJavaFile.append("\t\t\tSystem.err.println(\"GE result is NaN.\");\n");
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
        
        // And now we evaluate all the solutions with the compiled file:
        evaluator = null;
        double cumulatedFitness = 0.0;
        
        try {
            evaluator = (AbstractPopEvaluator) (new MyLoader(compiler.getWorkDir())).loadClass("PopEvaluator" + threadId).newInstance();
        } catch (Exception ex) {
            logger.severe(ex.getLocalizedMessage());
        }
        
        for (int i = 0; i < solutions.size(); ++i) {
            Solution<Variable<Integer>> solution = solutions.get(i);
            
            // For every patient
            for (int p = 0; p < clinicalTable.size(); p++) {
                int fromP = -1;
                int toP = -1;
                
                //System.out.println("Patient: " + p);
                for (int ex=0; ex<patientsIdXs[p].length; ex++) {
                    if ((patientsIdXs[p][ex] >= 0) && (fromP < 0)) {
                        fromP = patientsIdXs[p][ex];
                    }
                    if ((patientsIdXs[p][ex] >= 0)) {
                        toP = patientsIdXs[p][ex];
                    }
                }
//System.out.println("From: " + fromP + ", to: " + toP);

                evaluator.setDataTable((ArrayList<double[]>) dataTable.getDataTable("training", fromP, toP));
                evaluator.setDataLimits(patientsIdXs[p]);
                
                cumulatedFitness = dataTable.evaluate(evaluator, solution, p, i);
                if (Double.isNaN(cumulatedFitness)) {
                    logger.info("I have a NaN number here");
                }
            }
            solution.getObjectives().set(0, cumulatedFitness);
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
        switch (properties.getProperty("Parallelization")) {
            case "yes":
                MasterWorkerThreads<Variable<Integer>> masterWorker = new MasterWorkerThreads<Variable<Integer>>(algorithm, problem, Integer.valueOf(properties.getProperty("NumCores")));
                masterWorker.execute();
            default:
                algorithm.initialize();
                algorithm.execute();
        }
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
