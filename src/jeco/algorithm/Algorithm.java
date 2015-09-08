/*
 * Copyright (C) 2010-2015 José Luis Risco Martín <jlrisco@ucm.es> and 
 * José Manuel Colmenar Verdugo <josemanuel.colmenar@urjc.es>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *  - José Luis Risco Martín
 *  - José Manuel Colmenar Verdugo
 */
package jeco.algorithm;

import java.util.Observable;
import jeco.problem.Problem;
import jeco.problem.Solutions;
import jeco.problem.Variable;

/**
 *
 * @author José L. Risco-Martín
 *
 */
public abstract class Algorithm<V extends Variable<?>> extends Observable {

  protected Problem<V> problem = null;
  // Attribute to stop execution of the algorithm.
  protected boolean stop = false;
  
  /**
   * Allows to stop execution after finishing the current generation; must be
   * taken into account in children classes.
   */
  public void stopExection() {
      stop = true;
  }

  public Algorithm(Problem<V> problem) {
    this.problem = problem;
  }

  public void setProblem(Problem<V> problem) {
    this.problem = problem;
  }

  public abstract void initialize();

  public abstract void step();

  public abstract Solutions<V> execute();
}
