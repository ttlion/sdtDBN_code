package com.github.tDBN.dbn;

import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

import com.github.tDBN.utils.Edge;
import com.github.tDBN.utils.Utils;

import au.com.bytecode.opencsv.CSVReader;

/**
 * Scores class.
 * It has all necessary methods to, from dynamic and, if given, static data, test all possible combinations
 * of nodes connections of a tDBN and create a tDBN that maximizes a certain scoring function given the
 * provided dynamic and static data.
 * 
 * @author zlm
 * @author Tiago Leao
 *
 */
public class Scores {

	private Observations observations;
	
	/**
	 * Static Observations may also be given
	 */
	private ObservationsStatic observStatic;

	/**
	 * scoresMatrix[t][i][j] is the score of the arc
	 * Xj[t+markovLag]->Xi[t+markovLag].
	 */
	private double[][][] scoresMatrix;

	/**
	 * parentNodes.get(t).get(i) is the list of optimal dynamic parents in
	 * {X[t],...,X[t+markovLag-1]} of Xi[t+markovLag] when there is no arc from
	 * X[t+markovLag] to X[t+markovLag].
	 */
	private List<List<List<Integer>>> parentNodesPast;
	
	/**
	 * parentStaticPast.get(t).get(i) is the list of optimal static parents in
	 * Xi[t+markovLag] when there is no arc from X[t+markovLag] to X[t+markovLag].
	 */
	private List<List<List<Integer>>> parentStaticPast;

	/**
	 * parentNodes.get(t).get(i).get(j) is the list of optimal dynamic parents in
	 * {X[t],...,X[t+markovLag-1]} of Xi[t+markovLag] when the arc
	 * Xj[t+markovLag]->Xi[t+markovLag] is present.
	 */
	private List<List<List<List<Integer>>>> parentNodes;
	
	/**
	 * parentStatic.get(t).get(i).get(j) is the list of optimal static parents 
	 * of Xi[t+markovLag] when the arc Xj[t+markovLag]->Xi[t+markovLag] is present.
	 */
	private List<List<List<List<Integer>>>> parentStatic;

	/**
	 * Upper limit on the number of dynamic parents from previous time slices.
	 */
	private int maxParents;
	
	/**
	 * Upper limit on the number of static parents from previous time slices.
	 */
	private int maxStaticParents;

	/**
	 * A list of all possible sets of dynamic parent nodes. Set cardinality lies within
	 * the range [1, maxParents].
	 */
	private List<List<List<List<Integer>>>> parentSets;
	
	/**
	 * A list of all possible sets of static parent nodes. Set cardinality lies within
	 * the range [1, maxStaticParents].
	 */
	private List<List<List<List<Integer>>>> staticSets;

	/**
	 * Each attribute has its lists of either dynamic (from previous timesteps or the same timestep)
	 * or static nodes that must be parents or that cannot be parents
	 */
	private List<List<List<Integer>>> forbiddenParentsPast; // This can simply be List<List<List<>>> because only used in the initial (no computationally expensive) part
	private List<List<List<Integer>>> mandatoryParentsPast;
	private List<List<Set<Integer>>> forbiddenParentsSameTimestep; // This must be List<List<Set<>>> because it will be used in each iteration of evaluate cycle
	private List<List<Set<Integer>>> mandatoryParentsSameTimestep;
	private List<List<List<Integer>>> forbiddenStaticParents;
	private List<List<List<Integer>>> mandatoryStaticParents;
	
	/*
	 *  Array with size n, where each position [i] is true if X[i] has at least 1 mandatory static parent, false otherwise
	 */
	private boolean[][] hasMandatoryStatic; 

	/**
	 * If true, evaluates only one score matrix for all transitions.
	 */
	private boolean stationaryProcess;

	private boolean evaluated = false;

	private boolean verbose;
	
	public Scores(Observations observations, int maxParents, boolean stationaryProcess, boolean verbose) {
		this(observations, maxParents, false, true, null, 0, null, null, null, null, null, null);
	}
	
	public Scores(Observations observations, int maxParents) {
		this(observations, maxParents, false, true);
	}
	
	/**
	 * Creates the Scores object, allocating space for every data structure necessary
	 */
	public Scores(Observations observations, int maxParents, boolean stationaryProcess, boolean verbose, ObservationsStatic observStatic, int maxStaticParents, String pathFileForbiddenDyn, String pathFileMandatoryDyn, String pathFileForbiddenStatic, String pathFileMandatoryStatic, String pathFileForbiddenDyn_sameTimestep, String pathFileMandatoryDyn_sameTimestep) {
		this.observations = observations; this.observStatic = observStatic;
		this.maxParents = maxParents; this.maxStaticParents = maxStaticParents;
		this.stationaryProcess = stationaryProcess;
		this.verbose = verbose;

		int n = this.observations.numAttributes(); int n_static = (observStatic != null) ? this.observStatic.numAttributes() : 0;
		int p = this.maxParents; int b = this.maxStaticParents;
		int markovLag = observations.getMarkovLag();
		
		// Check static parents error conditions
		if(observStatic != null) {
			if(maxStaticParents<0) {
				System.out.println("Maximum number of static parents must be >=0 !");
				System.exit(1);
			} else if (maxStaticParents > n_static) {
				System.out.println("Maximum number of static parents must be <= number static attributes!");
				System.exit(1);
			}	
		}
		
		boolean ret = fillForbiddenOrMandatoryLists(observations, observStatic, pathFileForbiddenDyn, pathFileMandatoryDyn, pathFileForbiddenStatic, pathFileMandatoryStatic);
		if(ret == false) {
			System.out.println("Error with forbidden or mandatory files!");
			System.exit(1);
		}
		
		ret = fillForbiddenOrMandatoryLists_sameTimestep(observations, pathFileForbiddenDyn_sameTimestep, pathFileMandatoryDyn_sameTimestep);
		if(ret == false) {
			System.out.println("Error with forbidden or mandatory files of the same timestep!");
			System.exit(1);
		}
		
		// calculate sum_i=1^k nCi
		int size = n * markovLag;
		for (int previous = n, i = 2; i <= p; i++) {
			int current = previous * (n - i + 1) / i;
			size += current;
			previous = current;
		}
		
		int numTransitions = stationaryProcess ? 1 : observations.numTransitions();
		
		// generate all possible dynamic parents sets
		parentSets = new ArrayList<List<List<List<Integer>>>>(numTransitions); // 1 set of lists per timestep
		
		for(int t=0 ; t<numTransitions; t++) { // For each transition
			List<List<List<Integer>>> listsInTransition = new ArrayList<List<List<Integer>>>(n); // 1 per node
			parentSets.add(listsInTransition);
			
			List<List<Integer>> forbiddenParentsPast_timestep = forbiddenParentsPast.get(t);
			List<List<Integer>> mandatoryParentsPast_timestep = mandatoryParentsPast.get(t);
			
			for(int i=0; i < n; i++) { // For each node
				List<List<Integer>> listsOfNodeInTimestep = new ArrayList<List<Integer>>(size);
				listsInTransition.add(listsOfNodeInTimestep);
				
				List<Integer> forbiddenParentsPast_node = forbiddenParentsPast_timestep.get(i);
				List<Integer> mandatoryParentsPast_node = mandatoryParentsPast_timestep.get(i);
				
				for (int k = 1; k <= p; k++) { // fill the created list
					generateCombinations(n * markovLag, k, listsOfNodeInTimestep, forbiddenParentsPast_node, mandatoryParentsPast_node);
				}
			}
		}
		
		if(observStatic != null) {
			
			int sizeStatic = n_static;
			for (int previous = n_static, i = 2; i <= b; i++) {
				int current = previous * (n_static - i + 1) / i;
				sizeStatic += current;
				previous = current;
			}
		
			hasMandatoryStatic = new boolean[numTransitions][n];
			
			// generate all possible static parents sets
			staticSets = new ArrayList<List<List<List<Integer>>>>(numTransitions); // 1 set of lists per timestep
			
			for(int t=0 ; t<numTransitions; t++) {
				List<List<List<Integer>>> listsInTransition_static = new ArrayList<List<List<Integer>>>(n); // 1 per node
				staticSets.add(listsInTransition_static);
				
				List<List<Integer>> forbiddenStaticParents_timestep = forbiddenStaticParents.get(t);
				List<List<Integer>> mandatoryStaticParents_timestep = mandatoryStaticParents.get(t);
				
				for(int i=0; i < n; i++) {
					List<List<Integer>> listsOfNodeInTimestep_static = new ArrayList<List<Integer>>(sizeStatic);
					listsInTransition_static.add(listsOfNodeInTimestep_static);
					
					List<Integer> forbiddenStaticParents_node = forbiddenStaticParents_timestep.get(i);
					List<Integer> mandatoryStaticParents_node = mandatoryStaticParents_timestep.get(i);
					
					hasMandatoryStatic[t][i] = (mandatoryStaticParents_node.size() != 0); // fill this array to be used in the evaluateWithStatic method
					
					for (int k = 1; k <= b; k++) { // fill the created list
						generateCombinations(n_static, k, listsOfNodeInTimestep_static, forbiddenStaticParents_node, mandatoryStaticParents_node);
					}
				}
			}
			
		}
		
		// ArrayList allocations
		parentNodesPast = new ArrayList<List<List<Integer>>>(numTransitions);
		parentNodes = new ArrayList<List<List<List<Integer>>>>(numTransitions);
		
		parentStaticPast = new ArrayList<List<List<Integer>>>(numTransitions);
		parentStatic = new ArrayList<List<List<List<Integer>>>>(numTransitions);

		for (int t = 0; t < numTransitions; t++) {

			parentNodesPast.add(new ArrayList<List<Integer>>(n));
			parentStaticPast.add(new ArrayList<List<Integer>>(n));
			
			// allocate all parentNodesPast and parentStaticPast positions also as ArrayLists
			List<List<Integer>> parentNodesPastTransition = parentNodesPast.get(t);
			List<List<Integer>> parentStaticNodesPastTransition = parentStaticPast.get(t);
			for (int i = 0; i < n; i++) {
				parentNodesPastTransition.add(new ArrayList<Integer>());
				parentStaticNodesPastTransition.add(new ArrayList<Integer>());
			}
			
			
			parentNodes.add(new ArrayList<List<List<Integer>>>(n));
			parentStatic.add(new ArrayList<List<List<Integer>>>(n));
			
			// allocate all parentNodes and parentStatic positions
			List<List<List<Integer>>> parentNodesTransition = parentNodes.get(t);
			List<List<List<Integer>>> parentStaticNodesTransition = parentStatic.get(t);
			for (int i = 0; i < n; i++) { // For each node 
				
				parentNodesTransition.add(new ArrayList<List<Integer>>(n));
				parentStaticNodesTransition.add(new ArrayList<List<Integer>>(n));
				
				List<List<Integer>> parentNodesTransitionHead = parentNodesTransition.get(i);
				List<List<Integer>> parentStaticNodesTransitionHead = parentStaticNodesTransition.get(i);
				for (int j = 0; j < n; j++) { // For each possible arc Xj->Xi
					parentNodesTransitionHead.add(new ArrayList<Integer>());
					parentStaticNodesTransitionHead.add(new ArrayList<Integer>());
				}
				
			}
		}

		// allocate scoresMatrix
		scoresMatrix = new double[numTransitions][n][n];

	}
	
	/**
	 * Evaluates all possible parents configurations, given the dynamic data, according to a certain scoring function.
	 * 
	 * @author zlmS
	 */
	public Scores evaluate(ScoringFunction sf) {
		
		if(observStatic!=null) { // If there is also static data to learn the tDBN, use evaluation with static data
			return evaluateWithStatic(sf);
		}
		
		int n = observations.numAttributes();
		int numTransitions = scoresMatrix.length;

		int[] numBestScoresPast = new int[n];
		int[][] numBestScores = new int[n][n];

		double bigNumber = 10E200;

		for (int t = 0; t < numTransitions; t++) {
			// System.out.println("evaluating score in transition " + t + "/" +
			// numTransitions);

			List<List<List<Integer>>> listsInTransition = parentSets.get(t);

			for (int i = 0; i < n; i++) {
				// System.out.println("evaluating node " + i + "/" + n);
				double bestScore = Double.NEGATIVE_INFINITY;
				for (List<Integer> parentSet : listsInTransition.get(i)) {
					double score = stationaryProcess ? sf.evaluate(observations, parentSet, i, null, null) : sf.evaluate(
							observations, t, parentSet, i, null, null);
					// System.out.println("Xi:" + i + " ps:" + parentSet +
					// " score:" + score);
					if (bestScore < score) {
						bestScore = score;
						parentNodesPast.get(t).set(i, parentSet);
						numBestScoresPast[i] = 1;
					} else if (bestScore == score)
						numBestScoresPast[i]++;
				}
				for (int j = 0; j < n; j++) {
					scoresMatrix[t][i][j] = -bestScore;
				}
			}

			List<Set<Integer>> forbiddenParents_timestep = forbiddenParentsSameTimestep.get(t);
			List<Set<Integer>> mandatoryParents_timestep = mandatoryParentsSameTimestep.get(t);
			
			for (int i = 0; i < n; i++) {
				Set<Integer> setForbiddenParents = forbiddenParents_timestep.get(i);
				Set<Integer> setMandatoryParents = mandatoryParents_timestep.get(i);
				for (int j = 0; j < n; j++) {
					if (i != j) {
						double bestScore = Double.NEGATIVE_INFINITY;
						for (List<Integer> parentSet : listsInTransition.get(i)) {
							double score = stationaryProcess ? sf.evaluate(observations, parentSet, j, i, null, null) : sf
									.evaluate(observations, t, parentSet, j, i, null, null);
							// System.out.println("Xi:" + i + " Xj:" + j +
							// " ps:" + parentSet + " score:" + score);
							if (bestScore < score) {
								bestScore = score;
								parentNodes.get(t).get(i).set(j, parentSet);
								numBestScores[i][j] = 1;
							} else if (bestScore == score)
								numBestScores[i][j]++;
						}

						if( setForbiddenParents.contains(j) ) { // bias Eij so that j->i does not appear in intra-slice
							scoresMatrix[t][i][j] += bestScore - bigNumber;
						} else if( setMandatoryParents.contains(j) ) { // bias Eij so that j->i appears in intra-slice
							scoresMatrix[t][i][j] += bestScore + bigNumber;
						} else { // do not bias if not neither forbidden nor mandatory
							scoresMatrix[t][i][j] += bestScore;
						}

					}
				}
			}

			if (verbose) {
				// System.out.println(Arrays.toString(numBestScoresPast));
				// System.out.println(Arrays.deepToString(numBestScores));
				long numSolutions = 1;
				for (int i = 0; i < n; i++)
					numSolutions *= numBestScoresPast[i];
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
						if (i != j)
							numSolutions *= numBestScores[i][j];
				System.out.println("Number of networks with max score: " + numSolutions);
			}

		}

		evaluated = true;

		return this;

	}

	// adapted from http://stackoverflow.com/a/7631893
	private void generateCombinations(int n, int k, List<List<Integer>> desiredList, List<Integer> listWithForbidden, List<Integer> listWithMandatory) {

		int[] comb = new int[k];
		for (int i = 0; i < comb.length; i++) {
			comb[i] = i;
		}

		boolean hasForbidden, hasMandatory;
		
		boolean done = false;
		while (!done) {

			hasForbidden = false;
			hasMandatory = true;

			List<Integer> intList = new ArrayList<Integer>(k);
			for (int i : comb) {
				intList.add(i);
				if(listWithForbidden.contains(i) == true) {
					hasForbidden = true;
					break;
				}
			}
			
			// Check if all mandatory parents are in the proposed intlist if this intlist is still valid (does not have forbidden)
			if(hasForbidden == false) {
				for(Integer elem : listWithMandatory) {
					if(intList.contains(elem) == false) {
						hasMandatory = false;
						break;
					}
				}
			}
			
			// intList may be added to desired list because it does not have any forbidden and has all mandatory
			if(hasForbidden == false && hasMandatory == true) {
				desiredList.add(intList);
			}

			int target = k - 1;
			comb[target]++;
			if (comb[target] > n - 1) {
				// carry the one
				while (comb[target] > ((n - 1 - (k - target)))) {
					target--;
					if (target < 0) {
						break;
					}
				}
				if (target < 0) {
					done = true;
				} else {
					comb[target]++;
					for (int i = target + 1; i < comb.length; i++) {
						comb[i] = comb[i - 1] + 1;
					}
				}
			}
		}
	}

	public double[][] getScoresMatrix(int transition) {
		return scoresMatrix[transition];
	}

	public DynamicBayesNet toDBN() {
		return toDBN(-1, false);
	}
	
	/**
	 * Creates a tDBN based on the evaluation performed of the several possible combinations
	 * of static and dynamic parents.
	 * To use this method, evaluate or evaluateWithStatic must first be called.
	 * 
	 * @author zlm
	 * @author Tiago Leao
	 */
	public DynamicBayesNet toDBN(int root, boolean spanning) {

		if (!evaluated)
			throw new IllegalStateException("Scores must be evaluated before being converted to DBN");

		int n = observations.numAttributes(); 
		int n_static = (observStatic!=null) ? observStatic.numAttributes() : 0;
		List<Attribute> staticAttributes = (observStatic != null) ? observStatic.getAttributes() : null;

		int numTransitions = scoresMatrix.length;

		List<BayesNet> transitionNets = new ArrayList<BayesNet>(numTransitions);

		for (int t = 0; t < numTransitions; t++) {

			// intra-slice relations are determined using Edmonds algorithm
			List<Edge> intraRelations = OptimumBranching.evaluate(scoresMatrix[t], root, spanning);

			if (verbose) {
				double score = 0;
				boolean[][] adj = Utils.adjacencyMatrix(intraRelations, n);

				for (int i = 0; i < n; i++) {
					boolean isRoot = true;
					for (int j = 0; j < n; j++) {
						if (adj[i][j]) {
							// score
							score += (scoresMatrix[t][i][j] - scoresMatrix[t][i][i]);
							isRoot = false;
						}
					}
					if (isRoot)
						// subtract since sign was inverted
						score -= scoresMatrix[t][i][i];
				}

				System.out.println("Network score: " + score);
			}

			List<Edge> interRelations = new ArrayList<Edge>(n * maxParents); // allocate space for dynamic parents
			
			List<Edge> staticParents = new ArrayList<Edge>(n_static * maxStaticParents); // allocate space for static parents too

			boolean[] hasParent = new boolean[n];

			for (Edge intra : intraRelations) {
				int tail = intra.getTail();
				int head = intra.getHead();
				List<List<List<Integer>>> parentNodesT = parentNodes.get(t);
				for (Integer nodePast : parentNodesT.get(head).get(tail)) {
					interRelations.add(new Edge(nodePast, head));
					hasParent[head] = true;
				}
				
				if(parentStatic != null) { // If static parents are given, also get the static parents
					List<List<List<Integer>>> staticparentNodesT = parentStatic.get(t);
					if(staticparentNodesT.get(head).get(tail)!=null) {
						for(Integer staticNodePast : staticparentNodesT.get(head).get(tail)) {
							staticParents.add(new Edge(staticNodePast, head));
						}
					}
				}
				
			}

			for (int i = 0; i < n; i++)
				if (!hasParent[i]) {
					List<List<Integer>> parentNodesPastT = parentNodesPast.get(t);
					for (int nodePast : parentNodesPastT.get(i)) {
						interRelations.add(new Edge(nodePast, i));
					}
					
					if(parentStaticPast != null) { // If static parents are given, also get the static parents
						List<List<Integer>> staticparentNodesPastT = parentStaticPast.get(t);
						if(staticparentNodesPastT.get(i)!=null) {
							for(Integer staticNodePast : staticparentNodesPastT.get(i)) {
								staticParents.add(new Edge(staticNodePast, i));
							}
						}
					}
				
				}

			BayesNet bt = new BayesNet(observations.getAttributes(), observations.getMarkovLag(), intraRelations,
					interRelations, staticAttributes, staticParents);

			transitionNets.add(bt);
		}

		return new DynamicBayesNet(observations.getAttributes(), transitionNets, staticAttributes);

	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		String ls = System.getProperty("line.separator");
		int n = scoresMatrix[0].length;
		DecimalFormat df = new DecimalFormat("0.00");

		int numTransitions = scoresMatrix.length;

		for (int t = 0; t < numTransitions; t++) {
			// sb.append("--- Transition " + t + " ---" + ls);
			// sb.append("Maximum number of parents in t: " + maxParents + ls);
			//
			// sb.append(ls);

			sb.append("Scores matrix:" + ls);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					sb.append(df.format(scoresMatrix[t][i][j]) + " ");
				}
				sb.append(ls);
			}

			// sb.append(ls);
			//
			// sb.append("Parents only in t:" + ls);
			// for (int i = 0; i < n; i++) {
			// sb.append(i + ": " + parentNodesPast.get(t).get(i) + ls);
			// }
			//
			// sb.append(ls);
			//
			// sb.append("Parents in t for each parent in t+1:" + ls);
			// sb.append("t+1:	");
			// for (int i = 0; i < n; i++) {
			// sb.append(i + "	");
			// }
			// sb.append(ls);
			// for (int i = 0; i < n; i++) {
			// sb.append(i + ":	");
			// for (int j = 0; j < n; j++) {
			// sb.append(parentNodes.get(t).get(i).get(j) + "	");
			// }
			// sb.append(ls);
			// }
			//
			// sb.append(ls);
		}

		return sb.toString();
	}
	
	
	/**
	 * Evaluates all possible parents configurations, given the dynamic and static data, according to a certain scoring function.
	 * 
	 * @author Tiago Leao
	 */
	public Scores evaluateWithStatic(ScoringFunction sf) {
		
		int n = observations.numAttributes();
		int numTransitions = scoresMatrix.length;

		int[] numBestScoresPast = new int[n];
		int[][] numBestScores = new int[n][n];

		double bigNumber = 10E200;

		for (int t = 0; t < numTransitions; t++) {
			//System.out.println("\n\n\nevaluating score in transition " + t + "/" + numTransitions);
			
			List<List<List<Integer>>> listsInTransition = parentSets.get(t);
			List<List<List<Integer>>> listsInTransition_static = staticSets.get(t);
			
			for (int i = 0; i < n; i++) {
				//System.out.println("\n\nevaluating node " + i + "/" + n);
				double bestScore = Double.NEGATIVE_INFINITY;
				for (List<Integer> parentSet : listsInTransition.get(i)) { // Try all possible configurations of dynamic parents
					
					double score;

					// First, try configuration without any static parents if possible
					if(hasMandatoryStatic[t][i] == false) { // In this case (no mandatory static parents), the configuration without static parents is possible
						score = stationaryProcess ? sf.evaluate(observations, parentSet, i, null, null) : sf.evaluate(
								observations, t, parentSet, i, null, null);
					
						if (bestScore < score) {
							bestScore = score;
							parentNodesPast.get(t).set(i, parentSet);
							parentStaticPast.get(t).set(i, null);
							numBestScoresPast[i] = 1;
						} else if (bestScore == score)
							numBestScoresPast[i]++;
					}					
					
					// Then try with all possible combinations of static parents
					for(List<Integer> staticParentSet : listsInTransition_static.get(i)) {
						
						score = stationaryProcess ? sf.evaluate(observations, parentSet, i, observStatic, staticParentSet) : sf.evaluate(
								observations, t, parentSet, i, observStatic, staticParentSet);
						
						if (bestScore < score) {
							bestScore = score;
							parentNodesPast.get(t).set(i, parentSet);
							parentStaticPast.get(t).set(i, staticParentSet);
							numBestScoresPast[i] = 1;
						} else if (bestScore == score)
							numBestScoresPast[i]++;
						
					}
					
				}
				for (int j = 0; j < n; j++) {
					scoresMatrix[t][i][j] = -bestScore;
				}
			}
			
			List<Set<Integer>> forbiddenParents_timestep = forbiddenParentsSameTimestep.get(t);
			List<Set<Integer>> mandatoryParents_timestep = mandatoryParentsSameTimestep.get(t);
			
			for (int i = 0; i < n; i++) {
				Set<Integer> setForbiddenParents = forbiddenParents_timestep.get(i);
				Set<Integer> setMandatoryParents = mandatoryParents_timestep.get(i);
				for (int j = 0; j < n; j++) {
					if (i != j) {
						double bestScore = Double.NEGATIVE_INFINITY;
						for (List<Integer> parentSet : listsInTransition.get(i)) { // Try all possible configurations of dynamic parents
							
							double score;
							
							// First, try configuration without any static parents if possible
							if(hasMandatoryStatic[t][i] == false) { // In this case (no mandatory static parents), the configuration without static parents is possible
								score = stationaryProcess ? sf.evaluate(observations, parentSet, j, i, null, null) : sf
										.evaluate(observations, t, parentSet, j, i, null, null);
							
								if (bestScore < score) {
									bestScore = score;
									parentNodes.get(t).get(i).set(j, parentSet);
									parentStatic.get(t).get(i).set(j, null);
									numBestScores[i][j] = 1;
								} else if (bestScore == score)
									numBestScores[i][j]++;
							}
							
							// Then try with all possible combinations of static parents
							for(List<Integer> staticParentSet : listsInTransition_static.get(i)) {

								score = stationaryProcess ? sf.evaluate(observations, parentSet, j, i, observStatic, staticParentSet) : sf
										.evaluate(observations, t, parentSet, j, i, observStatic, staticParentSet);
								
								if (bestScore < score) {
									bestScore = score;
									parentNodes.get(t).get(i).set(j, parentSet);
									parentStatic.get(t).get(i).set(j, staticParentSet);
									numBestScores[i][j] = 1;
								} else if (bestScore == score)
									numBestScores[i][j]++;
							}
							
						}

						if( setForbiddenParents.contains(j) ) { // bias Eij so that j->i does not appear in intra-slice
							scoresMatrix[t][i][j] += bestScore - bigNumber;
						} else if( setMandatoryParents.contains(j) ) { // bias Eij so that j->i appears in intra-slice
							scoresMatrix[t][i][j] += bestScore + bigNumber;
						} else { // do not bias if not neither forbidden nor mandatory
							scoresMatrix[t][i][j] += bestScore;
						}

					}
				}
			}

			if (verbose) {
				// System.out.println(Arrays.toString(numBestScoresPast));
				// System.out.println(Arrays.deepToString(numBestScores));
				long numSolutions = 1;
				for (int i = 0; i < n; i++)
					numSolutions *= numBestScoresPast[i];
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
						if (i != j)
							numSolutions *= numBestScores[i][j];
				System.out.println("Number of networks with max score: " + numSolutions);
			}

		}

		evaluated = true;

		return this;
	}


	
	public boolean fillForbiddenOrMandatoryLists(Observations observations, ObservationsStatic observationsStatic, String pathFileForbiddenDyn, String pathFileMandatoryDyn, String pathFileForbiddenStatic, String pathFileMandatoryStatic) {
		
		if(observations == null)
			return false;
		
		// Always init even if not filling, for the program not to crash
		int n = observations.numAttributes();
		int numTransitions = stationaryProcess ? 1 : observations.numTransitions();
		int markovLag = observations.getMarkovLag();
		
		forbiddenParentsPast = new ArrayList<List<List<Integer>>>(numTransitions);
		mandatoryParentsPast = new ArrayList<List<List<Integer>>>(numTransitions);
		
		for(int t=0; t<numTransitions; t++) {
			List<List<Integer>> forbiddenParentsPast_timestep = new ArrayList<List<Integer>>(n);
			List<List<Integer>> mandatoryParentsPast_timestep = new ArrayList<List<Integer>>(n);
			
			forbiddenParentsPast.add(forbiddenParentsPast_timestep);
			mandatoryParentsPast.add(mandatoryParentsPast_timestep);
			
			for(int i = 0; i < n; i++) {
				forbiddenParentsPast_timestep.add(new ArrayList<Integer>(n/3));
				mandatoryParentsPast_timestep.add(new ArrayList<Integer>(n/3));
			}
		}
		
		if(pathFileForbiddenDyn != null)
			fillForbiddenOrMandatoryLists(pathFileForbiddenDyn, observations, null, forbiddenParentsPast);
		if(pathFileMandatoryDyn != null)
			fillForbiddenOrMandatoryLists(pathFileMandatoryDyn, observations, null, mandatoryParentsPast);
		
		// Check error conditions
		for(int t=0; t<numTransitions; t++) {
			
			List<List<Integer>> forbidden_timestep = forbiddenParentsPast.get(t);
			List<List<Integer>> mandatory_timestep = mandatoryParentsPast.get(t);
			
			for(int i = 0; i < n; i++) {
				List<Integer> forbidden = forbidden_timestep.get(i);
				List<Integer> mandatory = mandatory_timestep.get(i);
				
				if(forbidden.size() >= n*markovLag) {
					System.err.println("Error: Cannot forbid all parents from past in att " + observations.getAttributes().get(i).getName());
					System.exit(1);
				}
				if(mandatory.size() > maxParents) {
					System.err.println("Error: Cannot make number of mandatory dynParentsFromPast > maxParentsFromPast in att " + observations.getAttributes().get(i).getName());
					System.exit(1);
				}
				for(Integer elem : forbidden) {
					if(mandatory.contains(elem)) {
						System.err.println("Error: Cannot have same dynamic node as mandatory and forbidden parent of att " + observations.getAttributes().get(i).getName());
						System.exit(1);
					}
				}
			}
			
		}
		
		
		if(observationsStatic != null) {
			// Always init even if not filling, for the program not to crash
			int n_static = observationsStatic.numAttributes();
			
			forbiddenStaticParents = new ArrayList<List<List<Integer>>>(numTransitions);
			mandatoryStaticParents = new ArrayList<List<List<Integer>>>(numTransitions);
			
			for(int t=0; t<numTransitions; t++) {
				List<List<Integer>> forbiddenStatic_timestep = new ArrayList<List<Integer>>(n);
				List<List<Integer>> mandatoryStatic_timestep = new ArrayList<List<Integer>>(n);
				
				forbiddenStaticParents.add(forbiddenStatic_timestep);
				mandatoryStaticParents.add(mandatoryStatic_timestep);
				
				for(int i = 0; i < n; i++) {
					forbiddenStatic_timestep.add(new ArrayList<Integer>(n_static/3));
					mandatoryStatic_timestep.add(new ArrayList<Integer>(n_static/3));
				}
			}
			
			if(pathFileForbiddenStatic != null)
				fillForbiddenOrMandatoryLists(pathFileForbiddenStatic, observations, observationsStatic, forbiddenStaticParents);
			if(pathFileMandatoryStatic != null)
				fillForbiddenOrMandatoryLists(pathFileMandatoryStatic, observations, observationsStatic, mandatoryStaticParents);
			
			// Check error conditions
			for(int t=0; t<numTransitions; t++) {
				
				List<List<Integer>> forbidden_timestep = forbiddenStaticParents.get(t);
				List<List<Integer>> mandatory_timestep = mandatoryStaticParents.get(t);
				
				for(int i = 0; i < n; i++) {
					List<Integer> forbidden = forbidden_timestep.get(i);
					List<Integer> mandatory = mandatory_timestep.get(i);
					
					if(mandatory.size() > maxStaticParents) {
						System.err.println("Error: Cannot make number of mandatory staticParents > maxStaticParents in att " + observations.getAttributes().get(i).getName());
						System.exit(1);
					}
					for(Integer elem : forbidden) {
						if(mandatory.contains(elem)) {
							System.err.println("Error: Cannot have same static node as mandatory and forbidden parent of att " + observations.getAttributes().get(i).getName());
							System.exit(1);
						}
					}
				}
			}
			
		}
		
		return true;
		
	}
	
	public void fillForbiddenOrMandatoryLists(String pathToFile, Observations observations, ObservationsStatic observationsStatic, List<List<List<Integer>>> listsToFill) {
		
		// Create map with keys the variable names and values their respective indexes in attributes array
		Map<String, Integer> nameToIndx = new HashMap<String, Integer>();
		int cnt=0;
		for (Attribute att : observations.getAttributes()) {
			nameToIndx.put(att.getName(), new Integer(cnt));
			cnt++;
		}
		
		// Create map for static attributes names if static exist
		Map<String, Integer> nameToIndxStatic = new HashMap<String, Integer>();
		if(observationsStatic != null) {
			cnt=0;
			for (Attribute att : observationsStatic.getAttributes()) {
				nameToIndxStatic.put(att.getName(), new Integer(cnt));
				cnt++;
			}
		}
		
		int n = observations.numAttributes();
		int numTransitions = stationaryProcess ? 1 : observations.numTransitions();
		int markovLag = observations.getMarkovLag();
		
		try {

			// open and parse the useful observations csv file
			CSVReader reader = new CSVReader(new FileReader(pathToFile));
			List<String[]> lines = reader.readAll();
			reader.close();

			ListIterator<String[]> li = lines.listIterator();
			
			String[] dataLine; // declare dataline
			
			while (li.hasNext()) {

				dataLine = li.next();

				// check for line sanity
				if(observationsStatic == null) {
					if (dataLine.length < 4) { // Dynamic file must have at least 4 elements in each line
						System.err.println("Error: File " + pathToFile + " badly formatted. Line has len < 4.");
						System.exit(1);
					}
					if ((dataLine.length)%2 != 0) { // Length of dynamic file must be even!
						System.err.println("Error: File " + pathToFile + " badly formatted. Line does not have an even number of elements.");
						System.exit(1);
					}
				} else { // Static file must have at least 3 elements in each line
					if (dataLine.length < 3) {
						System.err.println("Error: File " + pathToFile + " badly formatted. Line has len < 3.");
						System.exit(1);
					}
				}
				
				int timeStepInInput = -1;
				try {
					timeStepInInput = Integer.parseInt(dataLine[0]);
				} catch (NumberFormatException e) {
					System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[0] + " is not a valid timestep.");
					System.exit(1);
				}
				
				if( (timeStepInInput-markovLag) < 0 || (timeStepInInput-markovLag) >= numTransitions ) {
					System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[0] + " is not a valid timestep.");
					System.exit(1);
				}
				
				int t = timeStepInInput - markovLag;
				
				
				if(nameToIndx.containsKey(dataLine[1]) == false) {
					System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[1] + " is not a dynamic attribute.");
					System.exit(1);
				}
				
				int child = nameToIndx.get(dataLine[1]);
				
				if(observationsStatic == null) { // Case where putting the forbidden/mandatory dynamic parents
					
					for(int i=2; i<dataLine.length; i+=2) {
					
						if(nameToIndx.containsKey(dataLine[i]) == false) {
							System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[i] + " is not a valid attribute to consider as dynamic parent.");
							System.exit(1);
						}
						int parent = nameToIndx.get(dataLine[i]);
						
						int numberTimestepsBehind = -1000;
						try {
							numberTimestepsBehind = Integer.parseInt(dataLine[i+1]);
						} catch (NumberFormatException e) {
							System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[i+1] + " is not a valid timestep.");
							System.exit(1);
						}
						
						if(numberTimestepsBehind == 0) {
							System.err.println("Error: File " + pathToFile + " badly formatted. " + numberTimestepsBehind + " timesteps behind is not valid (must be at least 1 timestep behind).");
							System.exit(1);
						}
						
						if( (markovLag - Math.abs(numberTimestepsBehind)) < 0 ) {
							System.err.println("Error: File " + pathToFile + " badly formatted. " + Math.abs(numberTimestepsBehind) + " timesteps behind is not valid with markovLag " + markovLag);
							System.exit(1);
						}
						
						int parentId_withMarkovLagIncluded = parent + (markovLag - Math.abs(numberTimestepsBehind)) * n;
						
						// add forbidden/mandatory relation dataline[i]-->dataline[0]
						listsToFill.get(t).get(child).add(parentId_withMarkovLagIncluded);
					
					}
					
				} else { // Case where putting the forbidden/mandatory static parents
					
					for(int i=2; i<dataLine.length; i++) {
						if(nameToIndxStatic.containsKey(dataLine[i]) == false) {
							System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[i] + " is not a valid attribute to consider as static parent.");
							System.exit(1);
						}
						
						// add forbidden/mandatory relation dataline[i]-->dataline[0]
						listsToFill.get(t).get(child).add(nameToIndxStatic.get(dataLine[i]));
					}
					
				}
				
			}
			
		} catch (IOException e) {
			System.err.println("Error: File " + pathToFile + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}

		return;
	}
	
	
	public boolean fillForbiddenOrMandatoryLists_sameTimestep(Observations observations, String pathFileForbiddenDyn_sameTimestep, String pathFileMandatoryDyn_sameTimestep) {
		
		if(observations == null)
			return false;
		
		// Always init even if not filling, for the program not to crash
		int n = observations.numAttributes();
		int numTransitions = stationaryProcess ? 1 : observations.numTransitions();
		
		forbiddenParentsSameTimestep = new ArrayList<List<Set<Integer>>>(numTransitions);
		mandatoryParentsSameTimestep = new ArrayList<List<Set<Integer>>>(numTransitions);
		
		for(int t=0; t<numTransitions; t++) {
			List<Set<Integer>> forbiddenParentsSameTimestep_inTimestep = new ArrayList<Set<Integer>>(n);
			List<Set<Integer>> mandatoryParentsSameTimestep_inTimestep = new ArrayList<Set<Integer>>(n);
			
			forbiddenParentsSameTimestep.add(forbiddenParentsSameTimestep_inTimestep);
			mandatoryParentsSameTimestep.add(mandatoryParentsSameTimestep_inTimestep);
			
			for(int i = 0; i < n; i++) {
				forbiddenParentsSameTimestep_inTimestep.add(new HashSet<Integer>(n/3));
				mandatoryParentsSameTimestep_inTimestep.add(new HashSet<Integer>(n/3));
			}
		}
		
		if(pathFileForbiddenDyn_sameTimestep != null)
			fillForbiddenOrMandatoryLists_sameTimestep(pathFileForbiddenDyn_sameTimestep, observations, forbiddenParentsSameTimestep);
		if(pathFileMandatoryDyn_sameTimestep != null)
			fillForbiddenOrMandatoryLists_sameTimestep(pathFileMandatoryDyn_sameTimestep, observations, mandatoryParentsSameTimestep);
		
		// Check error conditions
		for(int t=0; t<numTransitions; t++) {
			List<Set<Integer>> forbidden_timestep = forbiddenParentsSameTimestep.get(t);
			List<Set<Integer>> mandatory_timestep = mandatoryParentsSameTimestep.get(t);
			
			for(int i = 0; i < n; i++) {
				Set<Integer> forbidden = forbidden_timestep.get(i);
				Set<Integer> mandatory = mandatory_timestep.get(i);
				
				if(forbidden.size() >= n-1) {
					if(forbidden.contains(i) == false || forbidden.size() == n ) {
						System.err.println("Error: Cannot forbid all parents from own timestep in att " + observations.getAttributes().get(i).getName() + ". Check proper parameter to define the root of the tree!");
						System.exit(1);
					}
				}
				
				if(mandatory.contains(i)) {
					System.err.println("Error: File " + pathFileMandatoryDyn_sameTimestep + " badly formatted. " + observations.getAttributes().get(i).getName() + " cannot be parent of itself in the same timestep.");
					System.exit(1);
				}
				
				for(Integer elem : forbidden) {
					if(mandatory.contains(elem)) {
						System.err.println("Error: Cannot have same dynamic node as mandatory and forbidden parent of att " + observations.getAttributes().get(i).getName());
						System.exit(1);
					}
				}
			}
		}
		
		
		return true;
	}
	
	public void fillForbiddenOrMandatoryLists_sameTimestep(String pathToFile, Observations observations, List<List<Set<Integer>>> setsToFill) {
		
		// Create map with keys the variable names and values their respective indexes in attributes array
		Map<String, Integer> nameToIndx = new HashMap<String, Integer>();
		int cnt=0;
		for (Attribute att : observations.getAttributes()) {
			nameToIndx.put(att.getName(), new Integer(cnt));
			cnt++;
		}
		
		int numTransitions = stationaryProcess ? 1 : observations.numTransitions();
		int markovLag = observations.getMarkovLag();
		
		try {

			// open and parse the useful observations csv file
			CSVReader reader = new CSVReader(new FileReader(pathToFile));
			List<String[]> lines = reader.readAll();
			reader.close();

			ListIterator<String[]> li = lines.listIterator();
			
			String[] dataLine; // declare dataline
			
			while (li.hasNext()) {

				dataLine = li.next();

				// check for line sanity
				if (dataLine.length < 3) {
					System.err.println("Error: File " + pathToFile + " badly formatted. Line has len < 3.");
					System.exit(1);
				}
				
				int timeStepInInput = -1;
				try {
					timeStepInInput = Integer.parseInt(dataLine[0]);
				} catch (NumberFormatException e) {
					System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[0] + " is not a valid timestep.");
					System.exit(1);
				}
				
				if( (timeStepInInput-markovLag) < 0 || (timeStepInInput-markovLag) >= numTransitions ) {
					System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[0] + " is not a valid timestep.");
					System.exit(1);
				}
				
				int t = timeStepInInput - markovLag;
				
				if(nameToIndx.containsKey(dataLine[1]) == false) {
					System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[1] + " is not a dynamic attribute.");
					System.exit(1);
				}
				
				int child = nameToIndx.get(dataLine[1]);
				
				for(int i=2; i<dataLine.length; i++) {
					
					if(nameToIndx.containsKey(dataLine[i]) == false) {
						System.err.println("Error: File " + pathToFile + " badly formatted. " + dataLine[i] + " is not a valid attribute to consider as dynamic parent.");
						System.exit(1);
					}
					
					// add forbidden/mandatory relation dataline[i]-->dataline[1]
					setsToFill.get(t).get(child).add(nameToIndx.get(dataLine[i]));
					
				}
				
			}
			
		} catch (IOException e) {
			System.err.println("Error: File " + pathToFile + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
		
		return;
	}
	

}




