package jeco.algorithm.mopso;

import java.util.logging.Logger;


import jeco.problem.Solutions;
import jeco.problem.Variable;
import jeco.problems.zdt.ZDT1;
import jeco.util.logger.JecoLogger;

public class NSPSO_example {
	private static final Logger logger = Logger.getLogger(NSPSO_example.class.getName());
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		JecoLogger.setup();
		// First create the problem
		ZDT1 problem = new ZDT1(30);
		NSPSO<Variable<Double>> algorithm = new NSPSO<Variable<Double>>(problem, 100, 250, 0.4, 2.0, 2.0);
		algorithm.initialize();
		Solutions<Variable<Double>> solutions = algorithm.execute();
		logger.info("solutions.size()="+ solutions.size());
		System.out.println(solutions.toString());
	}
}
