package com.github.tDBN.dbn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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

	public String learnParameters(Observations o, boolean stationaryProcess) {

		if (stationaryProcess) {
			// assert there is only one transition network
			if (transitionNets.size() > 1)
				throw new IllegalArgumentException("DBN has more than one transition network, cannot "
						+ "learn parameters considering a stationary process");

			return transitionNets.get(0).learnParameters(o, -1);

		} else {
			int T = transitionNets.size();
			for (int t = 0; t < T; t++) {
				transitionNets.get(t).learnParameters(o, t);
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
		int n = attributes.size(); int n_static = staticAttributes.size();
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

		sb.append("rankdir=LR;"+ls); // Para se ter imagem na horizontal

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
			
			for (int i = 0; i < n_static; i++) {
				sb.append(staticAttributes.get(i).getName());
				sb.append("[shape=polygon, sides=5];" + ls);
			}
			sb.append("};"+ ls);
			
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

	public void getMostProbable(int transitionNetID, int attributeID, boolean stationary, int[][][] observations, int subject) {
		
		if(transitionNetID<0) return; // Neste caso nao se conseguira prever nada porque timestep do no e menor do que markovLag

		int n = attributes.size(), i=0, value;
		int config[] = new int[(markovLag + 1) * n];
		Configuration desiredConfig;
		List<Double> distribution;
		Map<Configuration, List<Double>> cpt;
		double sum;
		boolean allParentsSpecified;
		
		// Obter todos os pais
		List<Integer> parents;
		if (stationary == true) {
			parents = this.getTransitionNet(0).getParents(attributeID);
		} else {
			parents = this.getTransitionNet(transitionNetID).getParents(attributeID);
		}
		
		Stack<Integer> stackAttID = new Stack<Integer>();
		Stack<Integer> stackTimesteps = new Stack<Integer>();
		
		if(observations[transitionNetID][subject][attributeID + n*markovLag] < 0 ) {
			//System.out.println("Pondo na STACKaaaaaa: " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]");
			stackAttID.push(attributeID);
			stackTimesteps.push(transitionNetID);
		}

		for(Integer parent : parents) {
			if(observations[transitionNetID][subject][parent] < 0) {
				
				if(transitionNetID + parent/n - markovLag < 0) return; // Nao e possivel prever este pai, entao abortar
					
				//System.out.println("Pondo na STACK: " + attributes.get(parent%n).getName() + "[" + (transitionNetID + parent/n) + "]");
				stackAttID.push(parent%n);
				stackTimesteps.push(transitionNetID + parent/n - markovLag);
			}
		}
		
		while(stackAttID.isEmpty() == false) {
			
			attributeID = stackAttID.peek();
			transitionNetID = stackTimesteps.peek();
			//System.out.println("Olhando na STACK: " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]");
			
			if ( observations[transitionNetID][subject][n * markovLag + attributeID] >= 0) {
				//System.out.println("Ja estava preenchido! Tirando da STACK: " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]");
				attributeID = stackAttID.pop();
				transitionNetID = stackTimesteps.pop();
				continue;
			}
			
			if (stationary == true) {
				parents = this.getTransitionNet(0).getParents(attributeID);
			} else {
				parents = this.getTransitionNet(transitionNetID).getParents(attributeID);
			}
			
			allParentsSpecified = true;
			for(Integer parent : parents) {
				if(observations[transitionNetID][subject][parent] < 0) {
					
					if(transitionNetID + parent/n - markovLag < 0) return; // Nao e possivel prever este pai, entao abortar
					
					//System.out.println("Pondo na STACK: " + attributes.get(parent%n).getName() + "[" + (transitionNetID + parent/n) + "]");
					stackAttID.push(parent%n);
					stackTimesteps.push(transitionNetID + parent/n - markovLag);
					allParentsSpecified = false;
				}
			}
			
			if (allParentsSpecified == false) {
				continue;
			} else {
				attributeID = stackAttID.pop();
				transitionNetID = stackTimesteps.pop();
				
				//System.out.println("Tirando da STACK: " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]");
				
				// Inicializar os valores dos atributos (-1 se nao e pai, 0 se proprio attributo)
				for(i=0; i<config.length; i++)
					config[i] = -1; // Por default meter que atributo nao e pai, os pais serao postos num ciclo a seguir
				
				config[n * markovLag + attributeID] = 0; // O proprio no child tem que ter 0 na sua posicao
				
				for(Integer parent : parents) 
					config[parent] = observations[transitionNetID][subject][parent];
				
				desiredConfig = new Configuration(attributes, config);
				
				// Obter a CPT correta relativa ao atributo em questao na transicao em questao
				if (stationary == true) {
					cpt = this.getTransitionNet(0).getCPT(attributeID);
				} else {
					cpt = this.getTransitionNet(transitionNetID).getCPT(attributeID);
				}
				
				distribution = cpt.get(desiredConfig);
				
				sum = 0; for(double d : distribution) sum += d;
				
				distribution.add(1.0-sum);
				
				value = distribution.indexOf(Collections.max(distribution));
				
				distribution.remove(distribution.size()-1);
				
				observations[transitionNetID][subject][n * markovLag + attributeID]  = value; // Meter valor no timestep que se esta a analisar
				
				//System.out.println("Valor escolhido: " + attributes.get(attributeID).get(value));
				
				// Propagar valor para as outras matrizes de observacoes
				for(i=1; i < (markovLag+1); i++) {
					if (transitionNetID+i >= observations.length) break;
					observations[transitionNetID+i][subject][n * (markovLag-i) + attributeID] = value;
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
