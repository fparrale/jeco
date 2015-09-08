package jeco.algorithm.ga;

import java.util.logging.Level;
import jeco.operator.comparator.SimpleDominance;
import jeco.operator.crossover.SBXCrossover;
import jeco.operator.mutation.PolynomialMutation;
import jeco.operator.selection.BinaryTournament;
import jeco.problem.Solution;
import jeco.problem.Solutions;
import jeco.problem.Variable;
import jeco.problems.Rastringin;
import jeco.util.logger.JecoLogger;

public class SimpleGeneticAlgorithm_example {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		JecoLogger.setup(Level.FINE);
		// First create the problem
		Rastringin problem = new Rastringin(4);
		// Second create the algorithm
		PolynomialMutation<Variable<Double>> mutationOp = new PolynomialMutation<Variable<Double>>(problem);
		SBXCrossover<Variable<Double>> crossoverOp = new SBXCrossover<Variable<Double>>(problem);
		SimpleDominance<Variable<Double>> comparator = new SimpleDominance<Variable<Double>>();
		BinaryTournament<Variable<Double>> selectionOp = new BinaryTournament<Variable<Double>>(comparator);
		SimpleGeneticAlgorithm<Variable<Double>> ga = new SimpleGeneticAlgorithm<Variable<Double>>(problem, 100, 5000, true, mutationOp, crossoverOp, selectionOp);
		ga.initialize();
		Solutions<Variable<Double>> solutions = ga.execute();
		for(Solution<Variable<Double>> solution : solutions) {
			System.out.println("Fitness = " + solution.getObjectives().get(0));
		}
		//System.out.println("solutions.size()="+ solutions.size());
		//System.out.println(solutions.toString());
		//System.out.println("solutions.size()="+ solutions.size());
	}
}
