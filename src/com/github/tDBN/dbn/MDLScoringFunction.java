package com.github.tDBN.dbn;

import java.util.List;

public class MDLScoringFunction extends LLScoringFunction {

	@Override
	public double evaluate(Observations observations, int transition, List<Integer> parentNodesPast,
			Integer parentNodePresent, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet) {
		
		LocalConfiguration c;
		if(observStatic == null && staticParentSet == null ) {
			c = new LocalConfiguration(observations.getAttributes(), observations.getMarkovLag(),
					parentNodesPast, parentNodePresent, childNode);
		} else {
			c = new LocalConfigurationWithStatic(observations.getAttributes(), observations.getMarkovLag(),
					parentNodesPast, parentNodePresent, childNode, observStatic.getAttributes(), staticParentSet);
		}
		
		double score = super.evaluate(observations, transition, parentNodesPast, parentNodePresent, childNode, observStatic, staticParentSet);

		// regularizer term
		score -= 0.5 * Math.log(observations.numObservations(transition)) * c.getNumParameters();

		return score;
	}

}
