/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jeco.optimization.threads;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import jeco.problem.Problem;
import jeco.problem.Solution;
import jeco.problem.Solutions;
import jeco.problem.Variable;

/**
 *
 * @author jlrisco
 */
public class Worker<V extends Variable<?>> extends Thread {

    private static final Logger logger = Logger.getLogger(Worker.class.getName());

    protected Problem<V> problem;
    protected LinkedBlockingQueue<Solution<V>> sharedQueue = null;
    protected int numSolutions = 1;

    public Worker(Problem<V> problem, LinkedBlockingQueue<Solution<V>> sharedQueue, int numSolutions) {
        this.problem = problem;
        this.sharedQueue = sharedQueue;
        this.numSolutions = numSolutions;
    }

    @Override
    public void run() {
        Solutions<V> solutions = new Solutions<>();
        try {
            for (int i = 0; i < numSolutions; ++i) {
                Solution<V> solution = sharedQueue.poll(3, TimeUnit.SECONDS);
                if (solution != null) {
                    solutions.add(solution);
                }
            }
            problem.evaluate(solutions);
            solutions.clear();
        } catch (InterruptedException e) {
            logger.severe(e.getLocalizedMessage());
            logger.severe("Thread " + super.getId() + " has been interrupted. Shuting down ...");
        }

    }
}
