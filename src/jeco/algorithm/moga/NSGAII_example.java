package jeco.algorithm.moga;

import java.util.logging.Logger;

import jeco.operator.crossover.SBXCrossover;
import jeco.operator.mutation.PolynomialMutation;
import jeco.operator.selection.BinaryTournamentNSGAII;
import jeco.problem.Solutions;
import jeco.problem.Variable;
import jeco.problems.dtlz.DTLZ1;
import jeco.util.logger.JecoLogger;

public class NSGAII_example {
	private static final Logger logger = Logger.getLogger(NSGAII_example.class.getName());
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		JecoLogger.setup();
		// First create the problem
		DTLZ1 problem = new DTLZ1(30);
		// Second create the algorithm
		NSGAII<Variable<Double>> algorithm = new NSGAII<Variable<Double>>(problem, 100, 250, new PolynomialMutation<Variable<Double>>(problem), new SBXCrossover<Variable<Double>>(problem), new BinaryTournamentNSGAII<Variable<Double>>());
		algorithm.initialize();
		Solutions<Variable<Double>> solutions = algorithm.execute();
		logger.info("solutions.size()="+ solutions.size());
		System.out.println(solutions.toString());
	}
}
