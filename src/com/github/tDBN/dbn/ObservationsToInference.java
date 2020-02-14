package com.github.tDBN.dbn;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Map.Entry;

import com.github.tDBN.utils.Utils;

import au.com.bytecode.opencsv.CSVReader;

public class ObservationsToInference {

	/**
	 * Three-dimensional matrix of coded observation data which will be used for
	 * learning a dynamic Bayesian network.
	 * <ul>
	 * <li>the 1st index refers to the transition {t - markovLag + 1, ...
	 * ,t}->t+1;
	 * <li>the 2nd index refers to the the subject (set of observed attributes);
	 * <li>the 3rd index refers to the attribute and lies within the range [0,
	 * (1 + markovLag)*n[, where [0, markovLag*n[ refers to attributes in the
	 * past andand [markovLag*n, (1 + markovLag)*n[ refers to attributes in time
	 * t+1.
	 * </ul>
	 */
	private int[][][] usefulObservations;

	private int[][] staticUsefulObservations;

	/**
	 * Each column of the useful observation data refers to an attribute.
	 */
	private List<Attribute> attributes;

	private List<Attribute> staticAttributes;

	/**
	 * Number of subjects per transition. Only those who have complete data for
	 * a transition are stored.
	 */
	private int[] numSubjects;

	/**
	 * File that contains observations that will be converted to attributes and
	 * from which one can learn a DBN.
	 */
	private String observationsFileName;

	private String staticObservationsFileName;

	/**
	 * Header line of input useful observations CSV file.
	 */
	private String[] observationsHeader;

	private String[] staticObservationsHeader;

	/**
	 * Order of the Markov process, which is the number of previous time slices
	 * that influence the values in the following slice. Default is first-order
	 * Markov.
	 */
	private int markovLag = 1;
	
	private int inferenceAtts[][];
	private String subjects[];
	private int matrixWithAllPredictions[][][];
	private String mostProbable[][];
	
	// AQUI COMO CADA SUBJECT ESTA SEMPRE NA MESMA LINHA, SO E PRECISO MAPA  
	// DE STRING PRA INTEGER (EM VEZ DE PARA ARRAY COMO E FEITO EM Observations.java)
	private Map<String, Integer> subjectLineInMatrices; 
	
	public ObservationsToInference(String obsFileName, Integer markovLag, List<Attribute> inputAttributes, String staticObsFileName, List<Attribute> inputStaticAttributes) {
		this.observationsFileName = obsFileName; this.staticObservationsFileName = staticObsFileName;
		this.markovLag = markovLag != null ? markovLag : 1;
		
		if (inputAttributes == null) {
			System.out.println("Error!");
			System.exit(1);
		} else {
			this.attributes = inputAttributes;
			this.staticAttributes = inputStaticAttributes;
		}
		readFromFiles();
	}
	
	
	private void readFromFiles() {

		try {

			// open and parse the observations csv file
			CSVReader reader = new CSVReader(new FileReader(observationsFileName));
			List<String[]> lines = reader.readAll();
			reader.close();

			ListIterator<String[]> li = lines.listIterator();

			// get first line
			String[] header = li.next();

			int numTimeSlices = parseNumTimeSlices(header);
			int numTransitions = numTimeSlices - markovLag;

			int numAttributes = (header.length - 1) / numTimeSlices;

			observationsHeader = processHeader(header, numAttributes);
			
			if ( checkHeaderFormat(observationsHeader, attributes) == false ) {
				System.err.println("Format (order of variables) of observation file to make inference is not the same of learning file format");
				System.exit(1);
			}

			// allocate observations matrix
			int totalNumSubjects = lines.size() - 1;
			usefulObservations = new int[numTransitions][totalNumSubjects][(markovLag + 1) * numAttributes];
			numSubjects = new int[numTransitions];
			subjects = new String[totalNumSubjects];

			subjectLineInMatrices = new LinkedHashMap<String, Integer>(totalNumSubjects);

			String[] dataLine;
			int i=0;
			while (li.hasNext()) {

				dataLine = li.next();

				// check for line sanity
				if (dataLine.length != numTimeSlices * numAttributes + 1) {
					System.err.println(Arrays.deepToString(dataLine));
					System.err
							.println("Observations file: input data line does not have the correct number of columns.");
					System.err.println("Line length: " + dataLine.length);
					System.err.println("Number of time slices: " + numTimeSlices);
					System.err.println("Number of attributes: " + numAttributes);
					System.exit(1);
				}
				
				subjects[i++] = dataLine[0]; // Store subjects IDs

				subjectLineInMatrices.put(dataLine[0], new Integer(numSubjects[0]));

				for (int t = 0; t < numTransitions; t++) {

					String[] transition = Arrays.copyOfRange(dataLine, 1 + t * numAttributes, 1 + (t + markovLag + 1) * numAttributes);
					for (int j = 0; j < (markovLag + 1) * numAttributes; j++) {
						String value = transition[j];
						if (value.length() == 0 || value.equals("?")) { // Missing value
							usefulObservations[t][numSubjects[t]][j] = -1;
							continue;
						}
						
						int attributeId = j % numAttributes;
						Attribute attribute = attributes.get(attributeId);
						
						if (attribute.hasValue(value) == false) { // Value not in training data
							usefulObservations[t][numSubjects[t]][j] = -2;
							continue;
						}
						
						usefulObservations[t][numSubjects[t]][j] = attribute.getIndex(value);
					}
					numSubjects[t]++;
				}
			}
			
//			System.out.println("matrix Dynamic criada!!:");
//			for(int[][] matriz: usefulObservations) {
//				for(int[] linha : matriz) {
//					for(int valor : linha){
//						System.out.print(valor + " ");
//					}
//					System.out.println();
//				}
//				System.out.println();
//			}

		} catch (IOException e) {
			System.err.println("File " + observationsFileName + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
		
		// Return if there are no static observations
		if (staticObservationsFileName == null) 
			return;
		
		// Parse static observations file
		try {
			
			// open and parse the static observations csv file
			CSVReader reader = new CSVReader(new FileReader(staticObservationsFileName));
			List<String[]> lines = reader.readAll();
			reader.close();

			ListIterator<String[]> li = lines.listIterator();

			// get first line
			String[] header = li.next();

			int numStaticAttributes = (header.length - 1);
			
			// Parse static header	
			staticObservationsHeader = new String[numStaticAttributes];
			String auxiliarStaticHeder[] = Arrays.copyOfRange(header, 1, numStaticAttributes + 1);
			int i = 0;
			for (String attName : auxiliarStaticHeder)
				staticObservationsHeader[i++] = attName;
			
			// Check if same format of learning static file
			if ( checkHeaderFormat(staticObservationsHeader, staticAttributes) == false ) {
				System.err.println("Format (order of variables) of static observation file is not the same of static learning file format");
				System.exit(1);
			}
			
			// allocate static observations matrix
			staticUsefulObservations = new int[subjects.length][numStaticAttributes];
			
			String[] dataLine;
			while (li.hasNext()) {

				dataLine = li.next();

				// check for line sanity
				if (dataLine.length != numStaticAttributes + 1) {
					System.err.println(Arrays.deepToString(dataLine));
					System.err.println("Static Observations file: input data line does not have the correct number of columns.");
					System.err.println("Line length: " + dataLine.length);
					System.err.println("Number of attributes: " + numStaticAttributes);
					System.exit(1);
				}
				
				String subject = dataLine[0];
				
				Integer subjectLine = subjectLineInMatrices.get(subject);
				
				if( subjectLine == null) {
					System.err.println("Subject "+ subject + " is not in dynamic observations file! Aborting!");
					System.exit(1);
				}
				
				for (int j = 1; j < (numStaticAttributes + 1); j++) {
					
					String value = dataLine[j];
					int staticAttributeId = j-1;
					
					if (value.length() == 0 || value.equals("?")) { // Missing value
						staticUsefulObservations[subjectLine][staticAttributeId] = -1;
						continue;
					}
					
					Attribute staticAttribute = staticAttributes.get(staticAttributeId);
					
					if (staticAttribute.hasValue(value) == false) { // Value not in static training data
						staticUsefulObservations[subjectLine][staticAttributeId] = -2;
						continue;
					}
					
					staticUsefulObservations[subjectLine][staticAttributeId] = staticAttribute.getIndex(value);
				}
				
			}
			
//			System.out.println("matrix Static criada!!:");
//			for(int[] linha : staticUsefulObservations) {
//				for(int valor : linha){
//					System.out.print(valor + " ");
//				}
//				System.out.println();
//			}
			
		} catch (IOException e) {
			System.err.println("File " + staticObservationsFileName + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	/**
	 * Reads the second and last column of the header, parses the integer time
	 * value and returns the difference between the two, plus one. If parsing is
	 * not possible, exits. Also performs error checking on the number of
	 * columns.
	 * 
	 * @return the number of time slices in input file
	 */
	private static int parseNumTimeSlices(String[] header) {

		int timeFirstColumn = 0, timeLastColumn = 0;

		try {
			// get first and last column time identifier
			timeFirstColumn = Integer.parseInt(header[1].split("__")[1]);
			timeLastColumn = Integer.parseInt(header[header.length - 1].split("__")[1]);

		} catch (ArrayIndexOutOfBoundsException e) {
			System.err.println(Arrays.deepToString(header));
			System.err.println("Input file header does not comply to the 'attribute__t' format.");
			System.exit(1);
		} catch (NumberFormatException e) {
			System.err.println(Arrays.deepToString(header));
			System.err.println("Input file header does not comply to the 'attribute__t' format.");
			System.exit(1);
		}

		int numTimeSlices = timeLastColumn - timeFirstColumn + 1;

		// the number of columns per time slice must be constant
		// header contains an extra column with subject id
		if ((header.length - 1) % numTimeSlices != 0) {
			System.err.println(Arrays.deepToString(header));
			System.err.println("Input file header does not have a number of columns"
					+ " compatible with the number of time slices.");
			System.err.println("Header length: " + header.length);
			System.err.println("Number of time slices: " + numTimeSlices);
			System.exit(1);
		}

		return numTimeSlices;
	}
	
	/**
	 * Gets the name of the attributes from an input header line and the number
	 * of attributes.
	 */
	private String[] processHeader(String[] header, int numAttributes) {
		String[] newHeader = new String[numAttributes];
		String stripFirstHeader[] = Arrays.copyOfRange(header, 1, numAttributes + 1);
		int i = 0;
		for (String column : stripFirstHeader) {
			String[] columnParts = column.split("__");
			newHeader[i++] = columnParts[0];
		}
		return newHeader;
	}
	
	private boolean checkHeaderFormat(String[] inputHeader, List<Attribute> attributesToCompare) {
		int i=0;
		for(Attribute att : attributesToCompare) {
			if(att.getName().equals(inputHeader[i]) == false ) return false;
			i++;
		}
		return true;
	}
	
	public void parseAttributes(String attFileName) {
		
		int newMatrixDimension = -1;
		
		try {
			// open and parse the observations csv file
			CSVReader reader = new CSVReader(new FileReader(attFileName));
			List<String[]> lines = reader.readAll();
			reader.close();
			
			inferenceAtts = new int[lines.size()][2]; // [0] tem numero do atributo de acordo com a ordem de "List<Attribute> attributes" , [1] tem timestep
			
			int i=-1, desiredtimestep, found = 0;
			for(String[] line : lines) {
				i++;
				
				if (line.length != 2) {
					System.out.println("File with variables to make inference is not in proper format");
					System.exit(1);
				}
				
				found = 0;
				for ( Attribute att : attributes ) {
					if(att.getName().equals(line[0])) { // Found the proper attribute by its name!
						found = 1;
						
						try {
							desiredtimestep = Integer.parseInt(line[1]);
							found = 2;
						} catch (NumberFormatException e) {
							break;
						}
						
						if(desiredtimestep < 0) break;

						if( (desiredtimestep - markovLag + 1) > usefulObservations.length) {
								if((desiredtimestep - markovLag + 1) > newMatrixDimension) {
									newMatrixDimension = desiredtimestep - markovLag + 1;
									//System.out.println("Nova Dimensao!!: " + newMatrixDimension);
								}
						}

						found=3; // 3 = success
						inferenceAtts[i][0] = attributes.indexOf(att);
						inferenceAtts[i][1] = desiredtimestep;
						
					}
				}
				
				if(found == 0) {
					System.out.println("File with variables to make inference is not in proper format, specified attribute does not exist!");
					System.exit(1);
				} else if (found == 1) {
					System.out.println("File with variables to make inference is not in proper format, timestep not an integer!");
					System.exit(1);
				} else if (found == 2) {
					System.out.println("File with variables to make inference is not in proper format, negative timesteps not valid!");
					System.exit(1);
				}
				
			}
			
			if(newMatrixDimension!=-1) // E necessario aumentar a matriz das observacoes para ter os novos timesteps desejados
				usefulObservations = increaseObservationMatrix(usefulObservations, newMatrixDimension);
			
			mostProbable = new String[subjects.length][inferenceAtts.length];
			
		} catch (IOException e) {
			System.err.println("File " + attFileName + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	
	public void makeInference(boolean stationary, DynamicBayesNet dbn, String fileToWrite, boolean presentDistr) {
		
		StringBuilder sb = new StringBuilder();
		
		int attributeID = 0, i=0, transitionNetID = 0, attToPrintDim = 0, currLineIndx=0, maxProbIndx;
		
		int configs[][] = new int[usefulObservations[0].length][(markovLag + 1) * attributes.size()]; // Tem que ter todos os att do mesmo instante e os de tras necessarios consoante markovLag
		int staticConfigs[][] = ( staticUsefulObservations != null ) ? new int[usefulObservations[0].length][staticAttributes.size()] : null ;
		
		Configuration desiredConfigs[] = new Configuration[usefulObservations[0].length];
		
		List<Integer> parents, staticParents;
		Map<Configuration, List<Double>> cpt;
		List<Double> distribution;
		double count, maxProb;
		int canMakeInference[] = new int[usefulObservations[0].length];
		Attribute attToPrint;
		String valueToStore = null;
		
		currLineIndx=-1;
		for(int[] infLine : inferenceAtts) {
			currLineIndx++;
			
			attributeID = infLine[0];
			transitionNetID = infLine[1] - markovLag;
			
			//System.out.println("\n\nDetermining " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]" );
			
			if(transitionNetID < 0) {
				sb.append("Distributions " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]:\n");
				sb.append("Not possible to determine because timestep " + (transitionNetID+markovLag) + " < markovLag\n\n");
				continue;
			}
			
			if(stationary == false && transitionNetID >= dbn.getNumberTransitionNets()) { // Not possible to make inference in this case
				sb.append("Distributions " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]:\n");
				sb.append("Not possible to determine because non-stationary DBN only until timestep " + (dbn.getNumberTransitionNets()+markovLag-1) + "\n\n");
				continue;
			}

			// Obter a CPT correta relativa ao atributo em questao na transicao em questao
			if (stationary == true) {
				cpt = dbn.getTransitionNet(0).getCPT(attributeID);
			} else {
				cpt = dbn.getTransitionNet(transitionNetID).getCPT(attributeID);
			}
			
			// Obter todos os pais
			if (stationary == true) {
				parents = dbn.getTransitionNet(0).getParents(attributeID);
				staticParents = dbn.getTransitionNet(0).getStaticParents(attributeID);
			} else {
				parents = dbn.getTransitionNet(transitionNetID).getParents(attributeID);
				staticParents = dbn.getTransitionNet(transitionNetID).getStaticParents(attributeID);
			}
			
			// Inicializar os valores dos atributos (-1 se nao e pai, 0 se proprio attributo)
			for(int[] line : configs) {
				for(i=0; i<line.length; i++)
					line[i] = -1; // Por default meter que atributo nao e pai, os pais serao postos num ciclo a seguir
				
				line[attributes.size() * markovLag + attributeID] = 0; // O proprio no child tem que ter 0 na sua posicao
			}
			
			// Inicializar os valores dos atributos estaticos (todos -1) se for necessario
			if( staticParents.isEmpty() == false && staticConfigs != null ) {
				for(int[] line : staticConfigs) {
					for(i=0; i<line.length; i++)
						line[i] = -1;
				}
			}
				
			i=-1;
			for(int[] observation : usefulObservations[transitionNetID]) {
				i++;
				
				canMakeInference[i] = 0;
				//System.out.println("Values for subject " + i + ": ");
				
				// Check if all necessary staticParents are given in the observation being analyzed, if not, inference is not possible
				for(Integer staticParent : staticParents) {
					
					if(staticUsefulObservations == null) { // Has static parents, but no static observations were given, inference not possible
						canMakeInference[i] = -1;
						break;
					}
					
					if(staticUsefulObservations[i][staticParent] < 0) { // The necessary static observation was not given, inference not possible
						canMakeInference[i] = -1;
						break;
					}
					staticConfigs[i][staticParent] = staticUsefulObservations[i][staticParent];
				}
				
				if(canMakeInference[i] == 0) {
					for(Integer parent : parents) {
						
						if(observation[parent]<0) {
							//System.out.println("CALLING MOST PROBABLE WITH (netID=" + (transitionNetID + parent/(attributes.size()) - markovLag) + ", attID=" + parent%(attributes.size()) + ", " + stationary + ", tag= " + i + ")");
							dbn.getMostProbable(transitionNetID + parent/(attributes.size()) - markovLag, parent%(attributes.size()), stationary, usefulObservations, i, staticUsefulObservations);
							//System.out.println("\tEstimado " + attributes.get(parent%(attributes.size())).getName() + "[" + (transitionNetID + parent/(attributes.size())) + "]" + " -- Value: " + observation[parent]);
						}
						
						configs[i][parent] = observation[parent];
						//System.out.println("\tParent " + attributes.get(parent%(attributes.size())).getName() + "[" + (transitionNetID + parent/(attributes.size())) + "]" + " -- Value: " + observation[parent]);
						
						if(configs[i][parent] < 0) { // Nao foi possivel prever porque nao havia rede
							canMakeInference[i] = -1;
							break;
						}
					}
				}
				
				// Configuration e feita com a configs[i], que representa a configuracao dos pais no caso que se esta a testar
				if(canMakeInference[i] == 0) {
					
					if(staticParents.isEmpty() == false) { // Criar configuracao usando static attributes
						desiredConfigs[i] = new LocalConfigurationWithStatic(attributes, configs[i], staticAttributes, staticConfigs[i]);
						//System.out.println("CRiada config STATIC para no "+subjects[i]+":" + desiredConfigs[i]);
					} else { // Criar configuracao nao usando static attributes
						desiredConfigs[i] = new LocalConfiguration(attributes, configs[i]);
						//System.out.println("Criada config NAO STATIC para no "+subjects[i]+":" + desiredConfigs[i]);
					}
					
				}
					
			}
			
			// Imprimir atributo e instante temporal que se esta a fazer inferencia e os seus valores possiveis
			sb.append("Distributions " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]:\n");
			sb.append("id");
			
			attToPrint = attributes.get(attributeID);
			for(attToPrintDim = 0 ; attToPrintDim < attToPrint.size(); attToPrintDim++) {
				sb.append("," + attToPrint.get(attToPrintDim));
			}
			sb.append("\n");
			
			// Para cada caso, obter a linha correta da CPT, de acordo com a configuracao dos pais, que esta em desiredConfig
			i=-1;
			for(Configuration desiredConfig : desiredConfigs) {
				
				i++;
				sb.append(subjects[i]);
				if(canMakeInference[i] == -1 ) {
					for(int k=0; k<attToPrintDim; k++) sb.append(",-1");
					
					sb.append("\n");
					
					mostProbable[i][currLineIndx] = null;
					
					continue;
				}
				
				//System.out.println("Config subject " + subjects[i]  + ":" + desiredConfig);
				distribution = cpt.get(desiredConfig);
				
				count = 0.0; 
				for(double d : distribution) {
					sb.append(String.format(",%.3f", d));
					count += d;
				}
				sb.append(String.format(",%.3f\n", 1.0-count));
				
				maxProb = Collections.max(distribution);
				maxProbIndx = distribution.indexOf(maxProb);
				
				if(maxProb < (1.0-count))
					maxProbIndx = attToPrintDim-1;
				
				valueToStore = attributes.get(attributeID).get(maxProbIndx);
				
				mostProbable[i][currLineIndx] = valueToStore;
			}
			sb.append("\n");
		}
		
		if (presentDistr == false )
			sb = MostProbableValuesFromInference();
		
		
		if(fileToWrite != null) {
			try {
				Utils.writeToFile(fileToWrite, sb.toString());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		} else {
			System.out.println(sb.toString());
		}
	}
	
	private int[][][] increaseObservationMatrix(int[][][] currentMatrix, int newMatrixDimension) {
		
		int transitionNetID, subject, att, n = attributes.size();
		
		int [][][] auxMatrix = new int[newMatrixDimension][subjects.length][(markovLag + 1) * n];
		
		for(transitionNetID=0; transitionNetID<currentMatrix.length; transitionNetID++)
			for(subject=0; subject<subjects.length; subject++)
				System.arraycopy(currentMatrix[transitionNetID][subject], 0, auxMatrix[transitionNetID][subject], 0, (markovLag + 1) * n);
		
		for(transitionNetID=currentMatrix.length; transitionNetID < newMatrixDimension; transitionNetID++) // Fill with -1 the new positions
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * n; att++)
					auxMatrix[transitionNetID][subject][att] = -1;
		
		int maxNetIDtoProp = Math.min((currentMatrix.length + markovLag), newMatrixDimension);
		
		// Propagar valores passados para as outras matrizes de observacoes
		for(transitionNetID = currentMatrix.length; transitionNetID < maxNetIDtoProp; transitionNetID++) {
			for(subject=0; subject < subjects.length; subject++) {
				for(att=0; att < markovLag * n; att++) { // So vai ate aos que nao sao deste timestep
					auxMatrix[transitionNetID][subject][att] = auxMatrix[transitionNetID - 1][subject][att + n];
				}
			}
		}
		
//		System.out.println("Nova Matriz aumentada!!:");
//		for(int[][] matriz: auxMatrix) {
//			for(int[] linha : matriz) {
//				for(int valor : linha){
//					System.out.print(valor + " ");
//				}
//				System.out.println();
//			}
//			System.out.println();
//		}
		
		return auxMatrix;
	}
	
	public void getMostProbableTrajectory(int timestepMax, DynamicBayesNet dbn, boolean stationary) {
		
		int transitionNetID, subject, att;
		
		int maxNumTransNet = timestepMax - markovLag + 1;
		
		if(stationary == false) {
			if(maxNumTransNet > dbn.getNumberTransitionNets()) {
				System.err.println("Cannot predict until timestep " + timestepMax + " because non-stationary DBN learnt with markovLag=" + markovLag + " only has " + dbn.getNumberTransitionNets() +" transition networks!");
				System.exit(1);
			}
		}
		
		if(maxNumTransNet > usefulObservations.length) {
			matrixWithAllPredictions = increaseObservationMatrix(usefulObservations, maxNumTransNet);
		} else {
			matrixWithAllPredictions = usefulObservations;
		}
		
		// Predict all nodes where the are still observations
		for(transitionNetID=0; transitionNetID < usefulObservations.length; transitionNetID++)
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * attributes.size(); att++)
					if(matrixWithAllPredictions[transitionNetID][subject][att] < 0) // Case where value not given in original data
						dbn.getMostProbable(transitionNetID + att/(attributes.size()) - markovLag, att%(attributes.size()), stationary, matrixWithAllPredictions, subject, staticUsefulObservations);
		
		if(maxNumTransNet <= usefulObservations.length) return;
		
		// Predict all values of remaining timesteps
		for(transitionNetID=usefulObservations.length; transitionNetID < maxNumTransNet; transitionNetID++)
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * attributes.size(); att++)
					dbn.getMostProbable(transitionNetID + att/(attributes.size()) - markovLag, att%(attributes.size()), stationary, matrixWithAllPredictions, subject, staticUsefulObservations);
		
	}
	
	
	public void printMostProbableTrajectory(boolean isToFile, String fileName) {
		
		String output = toTimeSeriesHorizontal(matrixWithAllPredictions);
		
		if(isToFile == true) {
			try {
				Utils.writeToFile(fileName, output);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		} else {
			System.out.println(output);
		}
		
		return;
	}
	
	
	public String toTimeSeriesHorizontal(int[][][] desiredMatrix) {
		
		StringBuilder sb = new StringBuilder();
		String ls = System.getProperty("line.separator");
		int numTransitions = desiredMatrix.length + markovLag;
		int numAttributes = attributes.size();
		int numSubjects = subjects.length;
		int transition, att, subject, transitionNetID;
		int startCurrentAtt = markovLag*attributes.size();
		
		// Create header line
		sb.append("id");
		for (transition = 0; transition < numTransitions; transition++) {
			for (att=0; att < numAttributes; att++) {
				sb.append(","+ attributes.get(att).getName() +"__"+transition);
			}
		}
		sb.append(ls);
		
		// Create every line
		for(subject=0; subject<numSubjects;  subject++ ) {
			sb.append(subjects[subject]);
			
			// Create 1st markovLag transitions
			for (att=0; att < numAttributes * markovLag; att++) {
				if(desiredMatrix[0][subject][att]<0) {
					sb.append(",");
				} else {
					sb.append(","+   attributes.get(att%numAttributes).get(desiredMatrix[0][subject][att]));
				}
			}
			
			// Create remaining transitions
			for (transitionNetID = 0; transitionNetID < desiredMatrix.length; transitionNetID++) {
				for (att = startCurrentAtt; att < numAttributes*(markovLag+1); att++) {
					if(desiredMatrix[transitionNetID][subject][att]<0){
						sb.append(",");
					} else {
						sb.append(","+   attributes.get(att%numAttributes).get(desiredMatrix[transitionNetID][subject][att]));
					}
				}
			}
			sb.append(ls);
		}
		
		return sb.toString();
	}
	
	private StringBuilder MostProbableValuesFromInference() {
		
		if (mostProbable==null || inferenceAtts==null) {
			System.out.println("Must call before the proper methods to create and fill the most probable values for the desired attributes!");
			return null;
		}
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("id");
		for(int infAtt[] : inferenceAtts)
			sb.append(","+ attributes.get(infAtt[0]).getName() + "[" + infAtt[1] + "]");
		sb.append("\n");
		
		int i=0;
		for(String[] linha : mostProbable) {
			sb.append(subjects[i++]);
			for(String valor : linha){
				sb.append("," + valor);
			}
			sb.append("\n");
		}
		
		return sb;
	}
	
}
