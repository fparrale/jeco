package jeco.operator.crossover;

import jeco.problem.Solution;
import jeco.problem.Solutions;
import jeco.problem.Variable;
 /**
  * El cruce sí necesita plantillas porque se accede continuamente a las variables
  * @author jlrisco
  * @param <T>
  */

public abstract class CrossoverOperator<T extends Variable<?>> {
    abstract public Solutions<T> execute(Solution<T> parent1, Solution<T> parent2);
}
