package com.github.tDBN.dbn;

import java.util.List;

/**
 * Interface for Scoring Fucntions
 * 
 * @author zlm
 * @author Tiago Leao
 * 
 */
public interface ScoringFunction {

	public abstract double evaluate(Observations observations, int transition, List<Integer> parentNodesPast,
			int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet);

	public abstract double evaluate(Observations observations, int transition, List<Integer> parentNodesPast,
			Integer parentNodePresent, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet);

	/**
	 * Calculate score when process is stationary.
	 */
	public abstract double evaluate(Observations observations, List<Integer> parentNodesPast,
			Integer parentNodePresent, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet);

	/**
	 * Calculate score when process is stationary without parents from present timestep.
	 */
	public abstract double evaluate(Observations observations, List<Integer> parentNodesPast, int childNode, ObservationsStatic observStatic, List<Integer> staticParentSet);

}
