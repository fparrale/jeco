package jeco.operator.mutation;

import jeco.problem.Solution;
import jeco.problem.Variable;
import jeco.util.random.RandomGenerator;

//Solutions must be numeric
public class BooleanMutation<T extends Variable<Boolean>> extends MutationOperator<T> {
	/**
	 * Constructor
	 * Creates a new IntegerFlipMutation mutation operator instance
	 */
	public BooleanMutation(double probability) {
		super(probability);
	} // IntegerFlipMutation

	public Solution<T> execute(Solution<T> solution) {
		for (int i = 0; i < solution.getVariables().size(); i++) {
			if (RandomGenerator.nextDouble() < probability) {
				solution.getVariables().get(i).setValue(!solution.getVariables().get(i).getValue());
			}
		}
		return solution;
	} // execute
}

