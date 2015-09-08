package jeco.operator.selection;

import java.util.Comparator;

import jeco.operator.comparator.SolutionDominance;
import jeco.problem.Solution;
import jeco.problem.Solutions;
import jeco.problem.Variable;
import jeco.util.random.RandomGenerator;

public class BinaryTournament<T extends Variable<?>> extends SelectionOperator<T> {

    protected Comparator<Solution<T>> comparator;

    public BinaryTournament(Comparator<Solution<T>> comparator) {
        this.comparator = comparator;
    } // BinaryTournament

    public BinaryTournament() {
        this(new SolutionDominance<T>());
    } // Constructor

    public Solutions<T> execute(Solutions<T> solutions) {
        Solutions<T> result = new Solutions<T>();
        Solution<T> s1, s2;
        s1 = solutions.get(RandomGenerator.nextInt(0, solutions.size()));
        s2 = solutions.get(RandomGenerator.nextInt(0, solutions.size()));

        int flag = comparator.compare(s1, s2);
        if (flag == -1) {
            result.add(s1);
        } else if (flag == 1) {
            result.add(s2);
        } else if (RandomGenerator.nextDouble() < 0.5) {
            result.add(s1);
        } else {
            result.add(s2);
        }
        return result;
    } // execute
} // BinaryTournament

