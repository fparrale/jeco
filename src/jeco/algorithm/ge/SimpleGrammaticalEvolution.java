package jeco.algorithm.ge;

import jeco.algorithm.Algorithm;
import jeco.algorithm.ga.SimpleGeneticAlgorithm;
import jeco.operator.comparator.SimpleDominance;
import jeco.operator.crossover.SinglePointCrossover;
import jeco.operator.mutation.IntegerFlipMutation;
import jeco.operator.selection.BinaryTournament;
import jeco.problem.Problem;
import jeco.problem.Solutions;
import jeco.problem.Variable;

/**
 * Grammatical evolution using just one objective.
 * 
 * @author J. M. Colmenar
 */
public class SimpleGrammaticalEvolution extends Algorithm<Variable<Integer>> {

    SimpleGeneticAlgorithm<Variable<Integer>> algorithm;
    
    public SimpleGrammaticalEvolution(Problem<Variable<Integer>> problem, int maxPopulationSize, int maxGenerations, double probMutation, double probCrossover) {
        super(problem);
        
        // Algorithm operators
        IntegerFlipMutation<Variable<Integer>> mutationOperator = new IntegerFlipMutation<Variable<Integer>>(problem, probMutation);
        SinglePointCrossover<Variable<Integer>> crossoverOperator = new SinglePointCrossover<Variable<Integer>>(problem, SinglePointCrossover.DEFAULT_FIXED_CROSSOVER_POINT, probCrossover, SinglePointCrossover.ALLOW_REPETITION);
        SimpleDominance<Variable<Integer>> comparator = new SimpleDominance<Variable<Integer>>();
        BinaryTournament<Variable<Integer>> selectionOp = new BinaryTournament<Variable<Integer>>(comparator);
        
        algorithm = new SimpleGeneticAlgorithm<Variable<Integer>>(problem, 
                maxPopulationSize, maxGenerations, true, mutationOperator, crossoverOperator, selectionOp);
        
    }

    @Override
    public void initialize() {
        algorithm.initialize();
    }

    @Override
    public void step() {
        algorithm.step();
    }

    @Override
    public Solutions<Variable<Integer>> execute() {
        return algorithm.execute();
    }

}
