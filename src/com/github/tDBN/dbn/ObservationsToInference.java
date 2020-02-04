package com.github.tDBN.dbn;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

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

	/**
	 * Each column of the useful observation data refers to an attribute.
	 */
	private List<Attribute> attributes;

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

	/**
	 * Header line of input useful observations CSV file.
	 */
	private String[] observationsHeader;

	/**
	 * Order of the Markov process, which is the number of previous time slices
	 * that influence the values in the following slice. Default is first-order
	 * Markov.
	 */
	private int markovLag = 1;
	
	private int inferenceAtts[][];
	private String subjects[];
	private int matrixWithAllPredictions[][][];
	
	
	public ObservationsToInference(String obsFileName, Integer markovLag, List<Attribute> inputAttributes) {
		this.observationsFileName = obsFileName;
		this.markovLag = markovLag != null ? markovLag : 1;
		
		if (inputAttributes == null) {
			System.out.println("Error!");
			System.exit(1);
		} else {
			this.attributes = inputAttributes;
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
			
			if ( checkHeaderFormat(observationsHeader) == false ) {
				System.err.println("Format (order of variables) of observation file to make inference is not the same of learning file format");
				System.exit(1);
			}

			// allocate observations matrix
			int totalNumSubjects = lines.size() - 1;
			usefulObservations = new int[numTransitions][totalNumSubjects][(markovLag + 1) * numAttributes];
			numSubjects = new int[numTransitions];
			subjects = new String[totalNumSubjects];

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

		} catch (IOException e) {
			System.err.println("File " + observationsFileName + " could not be opened.");
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
	
	private boolean checkHeaderFormat(String[] inputHeader) {
		int i=0;
		for(Attribute att : attributes) {
			if(att.getName().equals(inputHeader[i]) == false ) return false;
			i++;
		}
		return true;
	}
	
	public int[][][] getObservationsMatrix() {
		return usefulObservations;
	}
	
	public void parseAttributes(String attFileName) {
		
		int timestepMaxDesired = -1;
		
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
						
						if(desiredtimestep > usefulObservations.length) {
							if(desiredtimestep > timestepMaxDesired)
								timestepMaxDesired = desiredtimestep;
						}
						
						if(desiredtimestep<1 ) break;
						
						found=3; // 3 = success
						inferenceAtts[i][0] = attributes.indexOf(att);
						inferenceAtts[i][1] = desiredtimestep - 1;
						
					}
				}
				
				if(found == 0) {
					System.out.println("File with variables to make inference is not in proper format, specified attribute does not exist!");
					System.exit(1);
				} else if (found == 1) {
					System.out.println("File with variables to make inference is not in proper format, timestep not an integer!");
					System.exit(1);
				} else if (found == 2) {
					System.out.println("File with variables to make inference is not in proper format, timestep not valid!");
					System.exit(1);
				}
				
			}
			
			if(timestepMaxDesired!=-1) // E necessario aumentar a matriz das observacoes para ter os novos timesteps desejados
				usefulObservations = increaseObservationMatrix(usefulObservations, timestepMaxDesired);
			
		} catch (IOException e) {
			System.err.println("File " + attFileName + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	
	public void makeInference(boolean stationary, DynamicBayesNet dbn, String fileToWrite) {
		
		StringBuilder sb = new StringBuilder();
		
		int attributeID = 0, timestep = 0, i=0;
		
		int configs[][] = new int[usefulObservations[0].length][(markovLag + 1) * attributes.size()]; // Tem que ter todos os att do mesmo instante e os de tras necessarios consoante markovLag
		Configuration desiredConfigs[] = new Configuration[usefulObservations[0].length];
		
		List<Integer> parents;
		Map<Configuration, List<Double>> cpt;
		List<Double> distribution;
		double count;
		int canMakeInference[] = new int[usefulObservations[0].length];
		
		for(int[] infLine : inferenceAtts) {
			attributeID = infLine[0];
			timestep = infLine[1];
		
			// Obter a CPT correta relativa ao atributo em questao na transicao em questao
			if (stationary == true) {
				cpt = dbn.getTransitionNet(0).getCPT(attributeID);
			} else {
				cpt = dbn.getTransitionNet(timestep).getCPT(attributeID);
			}
			
			// Obter todos os pais
			if (stationary == true) {
				parents = dbn.getTransitionNet(0).getParents(attributeID);
			} else {
				parents = dbn.getTransitionNet(timestep).getParents(attributeID);
			}
			
			// Inicializar os valores dos atributos (-1 se nao e pai, 0 se proprio attributo)
			for(int[] line : configs) {
				for(i=0; i<line.length; i++)
					line[i] = -1; // Por default meter que atributo nao � pai, os pais serao postos num ciclo a seguir
				
				line[attributes.size() * markovLag + attributeID] = 0; // O proprio no child tem que ter 0 na sua posicao
			}
			
			i=-1;
			for(int[] observation : usefulObservations[timestep]) {
				i++;
				
				canMakeInference[i] = 0;
				System.out.println("Values for subject " + i + ": ");
				for(Integer parent : parents) {
					
					if(observation[parent]<0) {
						dbn.getMostProbable(timestep - markovLag + parent/(attributes.size()), parent%(attributes.size()), stationary, usefulObservations, i);
						System.out.println("estimado " + attributes.get(parent%(attributes.size())).getName() + "[" + (timestep - markovLag + parent/(attributes.size())+1) + "]");
						System.out.println("parent: " + parent);
					}
					
					configs[i][parent] = observation[parent];
					System.out.println("\tParent " + parent + " -- Value: " + observation[parent]);
				}
				
				// Configuration e feita com a configs[i], que representa a configuracao dos pais no caso que se esta a testar
				if(canMakeInference[i] == 0)
					desiredConfigs[i] = new Configuration(attributes, configs[i]);
			}
			
			// Imprimir atributo e instante temporal que se esta a fazer inferencia e os seus valores possiveis
			sb.append("\nDistributions of node " + attributes.get(attributeID).getName() + "[" + (timestep+1) + "]:\n");
			sb.append("Posssible values: " + attributes.get(attributeID)+ "\n");
			
			// Para cada caso, obter a linha correta da CPT, de acordo com a configuracao dos pais, que esta em desiredConfig
			i=-1;
			for(Configuration desiredConfig : desiredConfigs) {
				
				i++;
				sb.append("\tDistribution subject " + subjects[i] + ": [");
				if(canMakeInference[i] == -1 ) {
					sb.append("Cannot make inference because one of the parents was not specified]\n");
					continue;
				} else if (canMakeInference[i] == -2) {
					sb.append("All probabilities 0 because one of parents has a value that never appeared in training]\n");
					continue;
				}
				
				distribution = cpt.get(desiredConfig);
				
				count = 0.0; 
				for(double d : distribution) {
					sb.append(String.format("%.3f, ", d));
					count += d;
				}
				sb.append(String.format("%.3f]\n", 1.0-count));
			}
		}
		
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
	
	private int[][][] increaseObservationMatrix(int[][][] currentMatrix, int newMaxTimesteps) {
		
		int i, j;
		int timestep, subject, att;
		
		int [][][] auxMatrix = new int[newMaxTimesteps][subjects.length][(markovLag + 1) * attributes.size()];
		
		for(i=0; i<currentMatrix.length; i++)
			for(j=0; j<subjects.length; j++)
				System.arraycopy(currentMatrix[i][j], 0, auxMatrix[i][j], 0, (markovLag + 1) * attributes.size());
		
		for(timestep=currentMatrix.length; timestep < newMaxTimesteps; timestep++) // Fill with -1 the new positions
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * attributes.size(); att++)
					auxMatrix[timestep][subject][att] = -1;
		
		// Propagar valores passados para as outras matrizes de observacoes
		for(timestep = currentMatrix.length; timestep < newMaxTimesteps; timestep++) {
			for(subject=0; subject < subjects.length; subject++) {
				for(att=0; att < markovLag * attributes.size(); att++) { // So vai ate aos que nao sao deste timestep
					auxMatrix[timestep][subject][att] = auxMatrix[timestep - markovLag + att/(attributes.size()) ][subject][attributes.size() * markovLag + att%(attributes.size())];
				}
			}
		}
		
		return auxMatrix;
	}
	
	public void getMostProbableTrajectory(int timestepMax, DynamicBayesNet dbn, boolean stationary) {
		
		int timestep, subject, att;
		
		if(stationary == false) {
			if(timestepMax > dbn.getNumberTransitionNets()) {
				System.err.println("Cannot predict until timestep " + timestepMax + " because DBN learnt only has " + dbn.getNumberTransitionNets() +" timesteps!");
				System.exit(1);
			}
		}
		
		matrixWithAllPredictions = increaseObservationMatrix(usefulObservations, timestepMax);
		
		// Predict all nodes where the are still observations
		for(timestep=0; timestep < usefulObservations.length; timestep++)
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * attributes.size(); att++)
					if(matrixWithAllPredictions[timestep][subject][att] < 0) // Case where value not given in original data
						dbn.getMostProbable(timestep - markovLag + att/(attributes.size()), att%(attributes.size()), stationary, matrixWithAllPredictions, subject);
		
		// Predict all values of remaining timesteps
		for(timestep=usefulObservations.length; timestep < timestepMax; timestep++)
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * attributes.size(); att++)
					dbn.getMostProbable(timestep - markovLag + att/(attributes.size()), att%(attributes.size()), stationary, matrixWithAllPredictions, subject);
						

		System.out.println("Matrix obtidaaaaa:");
		for(int[][] matriz: matrixWithAllPredictions) {
			for(int[] linha : matriz) {
				for(int valor : linha){
					System.out.print(valor + " ");
				}
				System.out.println();
			}
			System.out.println();
			System.out.println();
		}
		
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
		int numTransitions = desiredMatrix.length;
		int numAttributes = attributes.size();
		int numSubjects = subjects.length;
		int transition, att, subject;
		int startCurrentAtt = markovLag*attributes.size();
		
		// Create header line
		sb.append("id");
		for (transition = 0; transition <= numTransitions; transition++) {
			for (att=0; att < numAttributes; att++) {
				sb.append(","+ attributes.get(att).getName() +"__"+transition);
			}
		}
		sb.append(ls);
		
		// Create every line
		for(subject=0; subject<numSubjects;  subject++ ) {
			sb.append(subjects[subject]);
			
			// Create 1st transition
			for (att=0; att < numAttributes; att++) {
				if(desiredMatrix[0][subject][att]<0) {
					sb.append(",");
				} else {
					sb.append(","+   attributes.get(att).get(desiredMatrix[0][subject][att]));
				}
			}
			
			// Create remaining transitions
			for (transition = 0; transition < numTransitions; transition++) {
				for (att=0; att < numAttributes; att++) {
					if(desiredMatrix[transition][subject][startCurrentAtt+att]<0){
						sb.append(",");
					} else {
						sb.append(","+   attributes.get(att).get(desiredMatrix[transition][subject][startCurrentAtt+att]));
					}
				}
			}
			sb.append(ls);
		}
		
		return sb.toString();
	}
	
}