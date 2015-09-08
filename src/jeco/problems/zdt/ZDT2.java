package jeco.problems.zdt;

import java.util.ArrayList;

import jeco.operator.comparator.SolutionDominance;
import jeco.problem.Solution;
import jeco.problem.Solutions;
import jeco.problem.Variable;

public class ZDT2 extends ZDT {

    public ZDT2(Integer numberOfVariables) {
        super(numberOfVariables);
        for (int i = 0; i < numberOfVariables; i++) {
            lowerBound[i] = 0.0;
            upperBound[i] = 1.0;
        }
    } // ZDT1

    public ZDT2() {
        this(30);
    }

    public void evaluate(Solution<Variable<Double>> solution) {
        ArrayList<Variable<Double>> variables = solution.getVariables();
        double f1 = variables.get(0).getValue();
        double g = 0;
        for (int j = 1; j < numberOfVariables; ++j) {
            g += variables.get(j).getValue();
        }
        g *= 9.0;
        g /= numberOfVariables - 1;
        g += 1.0;
        double h = 1 - (f1 / g) * (f1 / g);
        solution.getObjectives().set(0, f1);
        solution.getObjectives().set(1, g * h);
    }

    public Solutions<Variable<Double>> computeParetoOptimalFront(int n) {
        Solutions<Variable<Double>> result = new Solutions<Variable<Double>>();

        double temp;
        for (int i = 0; i < n; ++i) {
            Solution<Variable<Double>> sol = new Solution<Variable<Double>>(numberOfObjectives);
            temp = 0.0 + (1.0 * i) / (n - 1);
            sol.getVariables().add(new Variable<Double>(temp));
            for (int j = 1; j < numberOfVariables; ++j) {
                sol.getVariables().add(new Variable<Double>(0.0));
            }

            evaluate(sol);
            result.add(sol);
        }

        result.reduceToNonDominated(new SolutionDominance<Variable<Double>>());
        return result;
    }
    public ZDT2 clone() {
    	ZDT2 clone = new ZDT2(this.numberOfVariables);
    	for(int i=0; i<numberOfVariables; ++i) {
    		clone.lowerBound[i] = lowerBound[i];
    		clone.upperBound[i] = upperBound[i];
    	}
    	return clone;
    }
}
