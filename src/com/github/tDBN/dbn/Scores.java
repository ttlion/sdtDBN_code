package com.github.tDBN.dbn;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import com.github.tDBN.utils.Edge;
import com.github.tDBN.utils.Utils;

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
	private List<List<Integer>> parentSets;
	
	/**
	 * A list of all possible sets of static parent nodes. Set cardinality lies within
	 * the range [1, maxStaticParents].
	 */
	private List<List<Integer>> staticSets;

	/**
	 * If true, evaluates only one score matrix for all transitions.
	 */
	private boolean stationaryProcess;

	private boolean evaluated = false;

	private boolean verbose;
	
	public Scores(Observations observations, int maxParents, boolean stationaryProcess, boolean verbose) {
		this(observations, maxParents, false, true, null, 0);
	}
	
	public Scores(Observations observations, int maxParents) {
		this(observations, maxParents, false, true);
	}
	
	/**
	 * Creates the Scores object, allocating space for every data structure necessary
	 */
	public Scores(Observations observations, int maxParents, boolean stationaryProcess, boolean verbose, ObservationsStatic observStatic, int maxStaticParents) {
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
		
		// calculate sum_i=1^k nCi
		int size = n * markovLag;
		for (int previous = n, i = 2; i <= p; i++) {
			int current = previous * (n - i + 1) / i;
			size += current;
			previous = current;
		}
		
		// generate all possible dynamic parents sets
		parentSets = new ArrayList<List<Integer>>(size);
		for (int i = 1; i <= p; i++) {
			generateCombinations(n * markovLag, i, parentSets);
		}
		
		// Same thing but for static attributes
		
		int sizeStatic = n_static;
		for (int previous = n_static, i = 2; i <= b; i++) {
			int current = previous * (n_static - i + 1) / i;
			sizeStatic += current;
			previous = current;
		}
		
		// generate all possible static parents sets
		staticSets = new ArrayList<List<Integer>>(sizeStatic);
		for (int i = 1; i <= b; i++) {
			generateCombinations(n_static, i, staticSets);
		}
		
		// ArrayList allocations
		int numTransitions = stationaryProcess ? 1 : observations.numTransitions();
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

		for (int t = 0; t < numTransitions; t++) {
			// System.out.println("evaluating score in transition " + t + "/" +
			// numTransitions);
			for (int i = 0; i < n; i++) {
				// System.out.println("evaluating node " + i + "/" + n);
				double bestScore = Double.NEGATIVE_INFINITY;
				for (List<Integer> parentSet : parentSets) {
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

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i != j) {
						double bestScore = Double.NEGATIVE_INFINITY;
						for (List<Integer> parentSet : parentSets) {
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

						scoresMatrix[t][i][j] += bestScore;

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
	private void generateCombinations(int n, int k, List<List<Integer>> desiredList) {

		int[] comb = new int[k];
		for (int i = 0; i < comb.length; i++) {
			comb[i] = i;
		}

		boolean done = false;
		while (!done) {

			List<Integer> intList = new ArrayList<Integer>(k);
			for (int i : comb) {
				intList.add(i);
			}
			desiredList.add(intList);

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

		for (int t = 0; t < numTransitions; t++) {
			//System.out.println("\n\n\nevaluating score in transition " + t + "/" + numTransitions);
			for (int i = 0; i < n; i++) {
				//System.out.println("\n\nevaluating node " + i + "/" + n);
				double bestScore = Double.NEGATIVE_INFINITY;
				for (List<Integer> parentSet : parentSets) { // Try all possible configurations of dynamic parents
					
					// First, try configuration without any static parents
					double score = stationaryProcess ? sf.evaluate(observations, parentSet, i, null, null) : sf.evaluate(
							observations, t, parentSet, i, null, null);
					
					if (bestScore < score) {
						bestScore = score;
						parentNodesPast.get(t).set(i, parentSet);
						parentStaticPast.get(t).set(i, null);
						numBestScoresPast[i] = 1;
					} else if (bestScore == score)
						numBestScoresPast[i]++;
					
					// Then try with all possible combinations of static parents
					for(List<Integer> staticParentSet : staticSets) {
						
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
			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i != j) {
						double bestScore = Double.NEGATIVE_INFINITY;
						for (List<Integer> parentSet : parentSets) { // Try all possible configurations of dynamic parents
							
							// First, try configuration without any static parents
							double score = stationaryProcess ? sf.evaluate(observations, parentSet, j, i, null, null) : sf
									.evaluate(observations, t, parentSet, j, i, null, null);
							
							if (bestScore < score) {
								bestScore = score;
								parentNodes.get(t).get(i).set(j, parentSet);
								parentStatic.get(t).get(i).set(j, null);
								numBestScores[i][j] = 1;
							} else if (bestScore == score)
								numBestScores[i][j]++;
							
							// Then try with all possible combinations of static parents
							for(List<Integer> staticParentSet : staticSets) {
								
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

						scoresMatrix[t][i][j] += bestScore;

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

}
