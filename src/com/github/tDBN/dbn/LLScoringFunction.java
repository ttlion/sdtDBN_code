package com.github.tDBN.dbn;

import java.util.List;

public class LLScoringFunction implements ScoringFunction {

	@Override
	public double evaluate(Observations observations, int transition, List<Integer> parentNodesPast, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet) {
		return evaluate(observations, transition, parentNodesPast, null, childNode, observStatic, staticParentSet);
	}

	@Override
	public double evaluate(Observations observations, int transition, List<Integer> parentNodesPast,
			Integer parentNodePresent, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet) {
		
		LocalConfiguration c;
		if(observStatic == null && staticParentSet == null ) {
			c = new LocalConfiguration(observations.getAttributes(), observations.getMarkovLag(),
					parentNodesPast, parentNodePresent, childNode);
		} else {
			//System.out.println("Criando config com static!!");
			c = new LocalConfigurationWithStatic(observations.getAttributes(), observations.getMarkovLag(),
					parentNodesPast, parentNodePresent, childNode, observStatic.getAttributes(), staticParentSet);
			
			//System.out.println(c);
		}
		
		

		double score = 0;

		do {
			c.setConsiderChild(false);
			int Nij = observations.count(c, transition, observStatic);
			//System.out.println("Nij:" + Nij);
			c.setConsiderChild(true);
			do {
				int Nijk = observations.count(c, transition, observStatic);
				//System.out.println("\t Nijk:" + Nijk);
				if (Nijk != 0 && Nijk != Nij) {
					score += Nijk * (Math.log(Nijk) - Math.log(Nij));
				}
			} while (c.nextChild());
		} while (c.nextParents());

		return score;
	}

	@Override
	public double evaluate(Observations observations, List<Integer> parentNodesPast, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet) {
		return evaluate(observations, parentNodesPast, null, childNode, observStatic, staticParentSet);
	}

	@Override
	public double evaluate(Observations observations, List<Integer> parentNodesPast, Integer parentNodePresent,
			int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet) {
		return evaluate(observations, -1, parentNodesPast, parentNodePresent, childNode, observStatic, staticParentSet);
	}

}
