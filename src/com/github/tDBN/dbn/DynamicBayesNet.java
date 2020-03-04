package com.github.tDBN.dbn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;

import com.github.tDBN.utils.Edge;

public class DynamicBayesNet {

	private List<Attribute> attributes;
	
	private List<Attribute> staticAttributes;

	private int markovLag;

	private BayesNet initialNet;

	private List<BayesNet> transitionNets;

	public DynamicBayesNet(List<Attribute> attributes, BayesNet initialNet, List<BayesNet> transitionNets, List<Attribute> staticAttributes) {
		this.attributes = attributes; this.staticAttributes = staticAttributes;
		this.initialNet = initialNet;
		this.transitionNets = transitionNets;
		this.markovLag = transitionNets.get(0).getMarkovLag();
	}

	public DynamicBayesNet(List<Attribute> attributes, List<BayesNet> transitionNets) {
		this(attributes, null, transitionNets, null);
	}

	public DynamicBayesNet(List<Attribute> attributes, List<BayesNet> transitionNets,List<Attribute> staticAttributes) {
		this(attributes, null, transitionNets, staticAttributes);
	}

	public DynamicBayesNet generateParameters() {
		initialNet.generateParameters();
		for (BayesNet transitionNet : transitionNets)
			transitionNet.generateParameters();
		return this;
	}

	public void learnParameters(Observations o) {
		learnParameters(o, false);
	}

	public void learnParameters(Observations o, ObservationsStatic observStatic) {
		learnParameters(o, false, observStatic);
	}
	
	public String learnParameters(Observations o, boolean stationaryProcess) {
		return learnParameters(o, stationaryProcess, null);
	}

	public String learnParameters(Observations o, boolean stationaryProcess, ObservationsStatic observStatic) {

		if (stationaryProcess) {
			// assert there is only one transition network
			if (transitionNets.size() > 1)
				throw new IllegalArgumentException("DBN has more than one transition network, cannot "
						+ "learn parameters considering a stationary process");

			return transitionNets.get(0).learnParameters(o, -1, observStatic);

		} else {
			int T = transitionNets.size();
			for (int t = 0; t < T; t++) {
				transitionNets.get(t).learnParameters(o, t, observStatic);
			}
		}
		return null;
	}

	public Observations generateObservations(int numIndividuals) {
		return generateObservations(numIndividuals, transitionNets.size(), false);
	}

	public Observations generateObservations(int numIndividuals, int numTransitions, boolean stationaryProcess) {
		int[][][] obsMatrix = generateObservationsMatrix(null, numIndividuals, numTransitions, stationaryProcess, false);
		return new Observations(attributes, obsMatrix);
	}

	public Observations forecast(Observations originalObservations, int numTransitions, boolean stationaryProcess,
			boolean mostProbable) {
		if (stationaryProcess) {
			// assert there is only one transition network
			if (transitionNets.size() > 1)
				throw new IllegalArgumentException("DBN has more than one transition network, cannot "
						+ "learn parameters considering a stationary process");

		}
		List<int[]> initialObservations = originalObservations.getFirst();
		int[][][] obsMatrix = generateObservationsMatrix(initialObservations, initialObservations.size(),
				numTransitions, stationaryProcess, mostProbable);
		return new Observations(originalObservations, obsMatrix);
	}

	public Observations forecast(Observations originalObservations) {
		return forecast(originalObservations, transitionNets.size(), false, false);
	}

	/**
	 * 
	 * @param initialObservations
	 * @param numTransitions
	 * @param numIndividuals
	 * @return
	 */
	private int[][][] generateObservationsMatrix(List<int[]> initialObservations, int numIndividuals,
			int numTransitions, boolean stationaryProcess, boolean mostProbable) {
		// System.out.println("generating observations");

		if (!stationaryProcess && numTransitions > transitionNets.size())
			throw new IllegalArgumentException("DBN only has " + transitionNets.size() + " "
					+ "transitions defined, cannot generate " + numTransitions + ".");

		int n = attributes.size();
		//
		int[][][] obsMatrix = new int[numTransitions][numIndividuals][(markovLag + 1) * n];

		for (int subject = 0; subject < numIndividuals; subject++) {
			int[] observation0 = initialObservations != null ? initialObservations.get(subject) : initialNet
					.nextObservation(null, mostProbable);
			int[] observationT = observation0;
			for (int transition = 0; transition < numTransitions; transition++) {
				System.arraycopy(observationT, 0, obsMatrix[transition][subject], 0, n * markovLag);

				int[] observationTplus1 = stationaryProcess ? transitionNets.get(0).nextObservation(observationT,
						mostProbable) : transitionNets.get(transition).nextObservation(observationT, mostProbable);
				System.arraycopy(observationTplus1, 0, obsMatrix[transition][subject], n * markovLag, n);

				observationT = Arrays.copyOfRange(obsMatrix[transition][subject], n, (markovLag + 1) * n);
			}
		}

		return obsMatrix;
	}

	public static List<int[]> compare(DynamicBayesNet original, DynamicBayesNet recovered) {
		return compare(original, recovered, false);
	}

	public static List<int[]> compare(DynamicBayesNet original, DynamicBayesNet recovered, boolean verbose) {
		assert (original.transitionNets.size() == recovered.transitionNets.size());
		int numTransitions = original.transitionNets.size();
		List<int[]> counts = new ArrayList<int[]>(numTransitions);
		for (int t = 0; t < numTransitions; t++) {
			counts.add(BayesNet.compare(original.transitionNets.get(t), recovered.transitionNets.get(t), verbose));
		}
		return counts;
	}

	public String toDot(boolean compactFormat) {
		StringBuilder sb = new StringBuilder();
		String ls = System.getProperty("line.separator");
		String dl = ls + ls;
		int n = attributes.size(); int n_static = (staticAttributes!=null) ? staticAttributes.size():0;
		int T = transitionNets.size();

		if (compactFormat && (T != 1 || markovLag != 1))
			throw new IllegalStateException(
					"More than one transition network or Markov lag larger than 1, cannot create compact graph.");

		boolean hasStaticArrows = false;
		for (int t = 0; t < T; t++) {
			if(transitionNets.get(t).hasStaticArrows() == true)
				hasStaticArrows = true;
		}
		
		// digraph init
		sb.append("digraph dbn{" + dl);

		sb.append("rankdir=LR;"+ls); // To have horizontal image

		if (compactFormat) {
			for (int i = 0; i < n; i++) {
				sb.append("X" + i);
				String attributeName = attributes.get(i).getName();
				if (attributeName != null)
					sb.append("[label=\"" + attributeName);
				else
					sb.append("[label=\"X" + i);
				sb.append("\"];" + ls);
			}
			sb.append(ls);
		} else {
			for (int t = 0; t < T + markovLag; t++) {
				// slice t attributes
				sb.append("{" + ls + "rank=same;"+ls);
				
				sb.append("rank" + t +  " [style=invisible];" + ls);
				
				for (int i = 0; i < n; i++) {
					sb.append("X" + i + "_" + t);
					String attributeName = attributes.get(i).getName();
					if (attributeName != null)
						sb.append("[label=\"" + attributeName);
					else
						sb.append("[label=\"X" + i);
					sb.append("[" + t + "]\"];" + ls);
				}
				sb.append("};"+ ls);
			}
			
			// static attributes
			if(hasStaticArrows == true) {
				sb.append("{" + ls + "rank=source;"+ls);
			} else {
				sb.append("{" + ls + "rank=min;"+ls);
			}
			
			// just some graphical definitions for the static attributes
			for (int i = 0; i < n_static; i++) {
				sb.append(staticAttributes.get(i).getName());
				sb.append("[shape=polygon, sides=5];" + ls);
			}
			sb.append("};"+ ls);
			
			for (int t = 0; t < T + markovLag - 1; t++) {
				sb.append("rank" + t + " -> rank" + (t+1) + " [color=white];" + ls);
			}
			
			for (int t = 0; t < T + markovLag; t++) {
				
				sb.append("rank" + t);
				
				for (int i = 0; i < n; i++) {
					sb.append(" -> X" + i + "_" + t);
				}
				
				sb.append(" [style=invis]" + ls);
				
			}
			
			
		}
		sb.append(ls);
		
		
		
		
		// transition and intra-slice (t>0) edges
		for (int t = 0; t < T; t++)
			sb.append(transitionNets.get(t).toDot(t, compactFormat));

		sb.append(ls + "}" + ls);

		return sb.toString();
	}

	public String toString(boolean printParameters) {
		StringBuilder sb = new StringBuilder();
		String ls = System.getProperty("line.separator");

		if (initialNet != null)
			sb.append(initialNet.toString(-1, false));

		int i = 0;
		for (Iterator<BayesNet> iter = transitionNets.iterator(); iter.hasNext();) {
			sb.append(iter.next().toString(i, printParameters));
			i++;
			if (iter.hasNext())
				sb.append("-----------------" + ls + ls);
		}

		return sb.toString();
	}

	public String toString() {
		return toString(false);
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) {
		// a really simple network

		Attribute a1 = new NominalAttribute();
		a1.setName("me");
		a1.add("yes");
		a1.add("no");

		Attribute a2 = new NumericAttribute();
		a2.add("10");
		a2.add("20");

		Attribute a3 = new NumericAttribute();
		a3.add("0");
		a3.add("1");

		Attribute a4 = new NumericAttribute();
		a4.add("0");
		a4.add("1");

		Attribute a5 = new NumericAttribute();
		a5.add("0");
		a5.add("1");

		Attribute a6 = new NumericAttribute();
		a6.add("0");
		a6.add("1");

		List<Attribute> a = Arrays.asList(a1, a2, a3, a4, a5, a6);

		Edge e1 = new Edge(0, 0);
		Edge e2 = new Edge(0, 1);
		Edge e3 = new Edge(0, 2);
		Edge e4 = new Edge(0, 3);
		Edge e5 = new Edge(0, 4);
		Edge e6 = new Edge(0, 5);
		Edge e7 = new Edge(1, 0);
		Edge e8 = new Edge(1, 1);
		Edge e9 = new Edge(1, 2);
		Edge e10 = new Edge(1, 3);
		Edge e11 = new Edge(1, 4);
		Edge e12 = new Edge(1, 5);
		Edge e13 = new Edge(2, 0);
		Edge e14 = new Edge(2, 1);
		Edge e15 = new Edge(2, 2);
		Edge e16 = new Edge(2, 3);
		Edge e17 = new Edge(2, 4);
		Edge e18 = new Edge(2, 5);
		Edge e19 = new Edge(3, 0);
		Edge e20 = new Edge(3, 1);
		Edge e21 = new Edge(3, 2);
		Edge e22 = new Edge(3, 3);
		Edge e23 = new Edge(3, 4);
		Edge e24 = new Edge(3, 5);
		Edge e25 = new Edge(4, 0);
		Edge e26 = new Edge(4, 1);
		Edge e27 = new Edge(4, 2);
		Edge e28 = new Edge(4, 3);
		Edge e29 = new Edge(4, 4);
		Edge e30 = new Edge(4, 5);
		Edge e31 = new Edge(5, 0);
		Edge e32 = new Edge(5, 1);
		Edge e33 = new Edge(5, 2);
		Edge e34 = new Edge(5, 3);
		Edge e35 = new Edge(5, 4);
		Edge e36 = new Edge(5, 5);

		// 0 -> 1 -> 2 -> 3 -> 4 -> 5
		List<Edge> intra = Arrays.asList(e2, e9, e16, e23, e30);

		//
		List<Edge> inter = Arrays.asList(e1, e4, e8, e15, e16, e20, e22, e26, e28, e31, e32, e34, e35);

		BayesNet b0 = new BayesNet(a, intra);
		b0.generateParameters();
		BayesNet bt = new BayesNet(a, intra, inter);
		bt.generateParameters();

		DynamicBayesNet dbn1 = new DynamicBayesNet(a, b0, Arrays.asList(bt), null);

		Observations o = dbn1.generateObservations(100);
		System.out.println(o);

		Scores s = new Scores(o, 2);
		s.evaluate(new LLScoringFunction());

		// System.out.println(s);

		// System.out.println(OptimumBranching.evaluate(s.getScoresMatrix()));

		DynamicBayesNet dbn2 = s.toDBN();

		System.out.println(dbn1.toDot(false));

	}
	
	/**
	 * Determine the most probable value for a certain attribute in a certain timestep.
	 * If some parents are needed to be determined they are, using a STACK aproach: a node
	 * is put in the stack if it needs to be determined. If a node needs a parent to be determined, then
	 * put the parent into the STACK for the parent to be determined before.
	 * 
	 * @param transitionNetID
	 * 		The transitionNetID (timestep-markovLag)
	 * @param attributeID
	 * 		The attribute ID in the list of attributes
	 * @param stationary
	 * 		Whether the dbn is stationary or not 
	 * @param observations
	 * 		Three dimensional array with dynamic observations
	 * @param subject
	 * 		The desired subject
	 * @param staticObservations
	 * 		Two dimensional arrat with static observations
	 * 
	 * @param infMode
	 * 		If 1, estimates values with most probable, otherwise estimates with random sampling according to probability distribution
	 * 
	 * @author Tiago Leao
	 * 
	 */
	public void getMostProbable(int transitionNetID, int attributeID, boolean stationary, int[][][] observations, int subject, int[][] staticObservations, int infMode) {
		
		// if timestep < markovLag, it will not be possible to determine the desired value
		if(transitionNetID<0) return;

		int n = attributes.size(), i=0;
		int config[] = new int[(markovLag + 1) * n]; int staticConfig[] = (staticAttributes!=null) ? new int[staticAttributes.size()] : null;
		Configuration desiredConfig;
		List<Double> distribution;
		Map<Configuration, List<Double>> cpt;
		double sum;
		boolean allParentsSpecified;
		
		Random randGen = new Random(); // for generating random numbers
		double rand, cummulative;
		int indxSelected, currIndx;
		
		// Get all static and dynamic parents of the attribute
		List<Integer> parents, staticParents;
		if (stationary == true) {
			parents = this.getTransitionNet(0).getParents(attributeID);
			staticParents = this.getTransitionNet(0).getStaticParents(attributeID);
		} else {
			parents = this.getTransitionNet(transitionNetID).getParents(attributeID);
			staticParents = this.getTransitionNet(transitionNetID).getStaticParents(attributeID);
		}
		
		// Get a stack for the attributes and a stack for the timesteps (these 2 stacks work as 1 stack with tuples (attribute, timestep), 
		// as every action on one stack is always also performed on the other stack)
		Stack<Integer> stackAttID = new Stack<Integer>();
		Stack<Integer> stackTimesteps = new Stack<Integer>();
		
		// Put the attribute specified in the input in the stack
		if(observations[transitionNetID][subject][attributeID + n*markovLag] < 0 ) {
			stackAttID.push(attributeID);
			stackTimesteps.push(transitionNetID);
		}
		
		// If some static parent is not specified in the observation file, inference is not possible
		for(Integer staticParent : staticParents) { 
			
			if(staticObservations == null) return; // There are static parents but no static observations were given
			
			if(staticObservations[subject][staticParent] < 0) { // There are static parents but no static observations were given for the static parent needed
				return;
			}
		}
		
		// Check all parents of the desired attribute. Put in stack the ones that need to be estimated
		for(Integer parent : parents) {
			if(observations[transitionNetID][subject][parent] < 0) {
				
				if(transitionNetID + parent/n - markovLag < 0) return; // If it is not possible to determine the parent, just return
					
				stackAttID.push(parent%n);
				stackTimesteps.push(transitionNetID + parent/n - markovLag);
			}
		}
		
		// Use the stack approach explained in the beginning of method
		while(stackAttID.isEmpty() == false) {
			
			// Get attribute and transitionNetID, only peeking from stack
			attributeID = stackAttID.peek();
			transitionNetID = stackTimesteps.peek();
			
			// Check if the attribute was already determined. If yes remove it from stack
			if ( observations[transitionNetID][subject][n * markovLag + attributeID] >= 0) {
				attributeID = stackAttID.pop();
				transitionNetID = stackTimesteps.pop();
				continue;
			}
			
			// Get all static and dynamic parents of the attribute
			if (stationary == true) {
				parents = this.getTransitionNet(0).getParents(attributeID);
				staticParents = this.getTransitionNet(0).getStaticParents(attributeID);
			} else {
				parents = this.getTransitionNet(transitionNetID).getParents(attributeID);
				staticParents = this.getTransitionNet(transitionNetID).getStaticParents(attributeID);
			}
			
			// If some static parent is not specified in the observation file, inference is not possible
			for(Integer staticParent : staticParents) {
				
				if(staticObservations == null) return; // There are static parents but no static observations were given
				
				if(staticObservations[subject][staticParent] < 0) { // There are static parents but no static observations were given for the static parent needed
					return;
				}
			}
			
			allParentsSpecified = true;
			
			// Check all parents of the desired attribute. Put in stack the ones that need to be estimated
			for(Integer parent : parents) {
				if(observations[transitionNetID][subject][parent] < 0) {
					
					if(transitionNetID + parent/n - markovLag < 0) return; // If it is not possible to determine the parent, just return
					
					stackAttID.push(parent%n);
					stackTimesteps.push(transitionNetID + parent/n - markovLag);
					
					allParentsSpecified = false; // Mark that at least 1 parent still needs to be estimated
				}
			}
			
			if (allParentsSpecified == false) { 
				continue;
				
			} else { // If all parents are already specified, attribute being analysed can be determined
				
				// Remove attribute and transitionNetID from stack
				attributeID = stackAttID.pop();
				transitionNetID = stackTimesteps.pop();
				
				// Initialize the dynamic attributes values (-1 if not parent, 0 if the child attribute)
				for(i=0; i<config.length; i++)
					config[i] = -1; // By default store an attribute as not being a parent
				
				config[n * markovLag + attributeID] = 0; // The child must have 0 in its array position
				
				// If needed (if there are static parents), initialize the static attributes values (-1 if not parent)
				if(staticParents.isEmpty() == false) {
					for(i=0; i<staticConfig.length; i++)
						staticConfig[i] = -1; // By default store an attribute as not being a parent
				}
				
				// Get the proper values for each dynamic and static parents
				for(Integer parent : parents) 
					config[parent] = observations[transitionNetID][subject][parent];
				
				for(Integer staticParent : staticParents)
					staticConfig[staticParent] = staticObservations[subject][staticParent];
				
				
				// Create configuration of parents
				if(staticParents.isEmpty() == true) {
					desiredConfig = new LocalConfiguration(attributes, config); // Only dynamic parents
				} else {
					desiredConfig = new LocalConfigurationWithStatic(attributes, config, staticAttributes, staticConfig); // Dynamic and static parents
				}
				
				// Get the proper CPT
				if (stationary == true) {
					cpt = this.getTransitionNet(0).getCPT(attributeID);
				} else {
					cpt = this.getTransitionNet(transitionNetID).getCPT(attributeID);
				}
				
				// Get the distribution
				distribution = cpt.get(desiredConfig);
				
				System.out.println("---------");
				System.out.println("Estimating: " + attributes.get(attributeID).getName() + "[" +  (transitionNetID+markovLag) + "]");
				System.out.println(distribution);
				
				if(infMode == 1) { // Estimate with most probable value
					
					// Get the index of the most probable value
					
					sum = 0; for(double d : distribution) sum += d;
					
					distribution.add(1.0-sum);
					
					indxSelected = distribution.indexOf(Collections.max(distribution)); // Get the index
					
					distribution.remove(distribution.size()-1);
					
					System.out.println("Indx escolhido pq max: " + indxSelected);
					
				} else {
					// Randomly estimate a value, using the probabilities of distribution
					rand = randGen.nextDouble();
					cummulative = 0.0;
					indxSelected = -1;
					currIndx = 0;
					
					for(double d : distribution) {
						cummulative += d;
						if(rand < cummulative) {
							indxSelected = currIndx;
							break;
						}
						currIndx++;
					}
					
					if(indxSelected == -1) 
						indxSelected = currIndx;
					
					System.out.println("Random nb: " + rand + "|| Indx escolhido: " + indxSelected);
					
				}
				
				System.out.println("---------");
				
				observations[transitionNetID][subject][n * markovLag + attributeID]  = indxSelected; // Put that value in the observations matrix
				
				// Propagate determined value to other observations matrices, according to markovLag
				for(i=1; i < (markovLag+1); i++) {
					if (transitionNetID+i >= observations.length) break;
					observations[transitionNetID+i][subject][n * (markovLag-i) + attributeID] = indxSelected;
				}
				
			}
		}
	}

	public BayesNet getTransitionNet(int t) {
		return transitionNets.get(t);
	}
	
	public int getNumberTransitionNets() {
		return transitionNets.size();
	}
	
}
