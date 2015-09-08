package jeco.operator.comparator;

import java.util.Comparator;

import jeco.operator.assigner.CrowdingDistance;
import jeco.operator.assigner.FrontsExtractor;
import jeco.problem.Solution;
import jeco.problem.Variable;

public class ComparatorNSGAII<V extends Variable<?>> implements Comparator<Solution<V>> {

	public ComparatorNSGAII() {
	}

	public int compare(Solution<V> left, Solution<V> right) {
		Integer rankLeft = left.getProperties().get(FrontsExtractor.propertyRank).intValue();
		Integer rankRight = right.getProperties().get(FrontsExtractor.propertyRank).intValue();

		int comp = rankLeft.compareTo(rankRight);
		if (comp == -1) {
			return -1;
		} else if (comp == 0) {
			Double crowdedDistanceLeft = left.getProperties().get(CrowdingDistance.propertyCrowdingDistance).doubleValue();
			Double crowdedDistanceRight = right.getProperties().get(CrowdingDistance.propertyCrowdingDistance).doubleValue();
			comp = crowdedDistanceLeft.compareTo(crowdedDistanceRight);
			if (comp == 1) {
				return -1;
			} else if (comp == 0) {
				return 0;
			} else {
				return 1;
			}
		} else {
			return 1;
		}
	}
}
