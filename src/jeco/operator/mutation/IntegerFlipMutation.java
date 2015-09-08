package jeco.operator.mutation;

import jeco.problem.Problem;
import jeco.problem.Solution;
import jeco.problem.Variable;
import jeco.util.random.RandomGenerator;

//Solutions must be numeric
public class IntegerFlipMutation<T extends Variable<Integer>> extends MutationOperator<T> {
	protected Problem<T> problem;

	/**
	 * Constructor
	 * Creates a new IntegerFlipMutation mutation operator instance
	 */
	public IntegerFlipMutation(Problem<T> problem, double probability) {
		super(probability);
		this.problem = problem;
	} // IntegerFlipMutation

  @Override
	public Solution<T> execute(Solution<T> solution) {
		for (int i = 0; i < solution.getVariables().size(); i++) {
			if (RandomGenerator.nextDouble() < probability) {
				int lowerBound = (int)Math.round(problem.getLowerBound(i));
				int upperBound = (int)Math.round(problem.getUpperBound(i));
				solution.getVariables().get(i).setValue(RandomGenerator.nextInteger(lowerBound, upperBound));
			}
		}
		return solution;
	} // execute
} // IntegerFlipMutation

