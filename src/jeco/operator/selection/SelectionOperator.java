package jeco.operator.selection;

import jeco.problem.Solutions;
import jeco.problem.Variable;

public abstract class SelectionOperator<T extends Variable<?>> {
    abstract public Solutions<T> execute(Solutions<T> solutions);
}
