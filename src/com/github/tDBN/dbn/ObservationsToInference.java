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
import java.util.Random;
import java.util.Map.Entry;

import com.github.tDBN.utils.Utils;

import au.com.bytecode.opencsv.CSVReader;

/**
 * 
 * Class to store the static and dynamic observations useful for inference in a created tDBN.
 * Stores the desired variables to perform inference.
 * It also has the methods to perform inference.
 * 
 * @author Tiago Leao
 * 
 */
public class ObservationsToInference {

	/**
	 * Three-dimensional matrix of coded dynamic observation data which will be used for
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
	 * Three-dimensional matrix of coded static observation data which will be used for
	 * learning a dynamic Bayesian network.
	 * <ul>
	 * <li>the 1st index refers to the the subject (set of observed attributes);
	 * <li>the 2nd index refers to the attribute
	 * </ul>
	 */
	private int[][] staticUsefulObservations;

	/**
	 * Each column of the useful dynamic observation data refers to a dynamic attribute.
	 */
	private List<Attribute> attributes;
	
	/**
	 * Each column of the useful static observation data refers to a static attribute.
	 */
	private List<Attribute> staticAttributes;

	/**
	 * Number of dynamic subjects per transition. All subjects are stored, even if having missing data.
	 */
	private int[] numSubjects;

	/**
	 * File that contains dynamic observations that will be used for inference.
	 */
	private String observationsFileName;
	
	/**
	 * File that contains static observations that will be used for inference.
	 */
	private String staticObservationsFileName;

	/**
	 * Header line of input useful dynamic observations CSV file.
	 */
	private String[] observationsHeader;
	
	/**
	 * Header line of input useful static observations CSV file.
	 */
	private String[] staticObservationsHeader;

	/**
	 * Order of the Markov process, which is the number of previous time slices
	 * that influence the values in the following slice. Default is first-order
	 * Markov.
	 */
	private int markovLag = 1;
	
	/**
	 * Matrix to store the attributes in which user specified to make inference.
	 * <ul>
	 * <li>the 1st index refers to attribute;
	 * <li>the 2nd index refers to the timestep
	 * </ul>
	 * 
	 */
	private int inferenceAtts[][];
	
	/**
	 * Store all subjects
	 */
	private String subjects[];
	
	/**
	 * Three-dimensional matrix with the data of all predictions made.
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
	private int matrixWithAllPredictions[][][];
	
	/**
	 * Two-dimensional matrix with the data of all predictions made in the attributes asked by user.
	 * <ul>
	 * <li>the 1st index refers to attribute;
	 * <li>the 2nd index refers to the timestep
	 * </ul>
	 * 
	 */
	private String predictions[][];
	
	/**
	 * Store line of each subject. As missing values are also stored, this line is the same in all matrices.
	 * <ul>
	 * <li> The key is the subject;
	 * <li> The value is the respective line.
	 * </ul>
	 */
	private Map<String, Integer> subjectLineInMatrices;
	
	
	/**
	 * Creates ObservationsToInference object from static and dynamic observations files and static and dynamic observations attributes already taken from the learned tdBN.
	 */
	public ObservationsToInference(String obsFileName, Integer markovLag, List<Attribute> inputAttributes, String staticObsFileName, List<Attribute> inputStaticAttributes) {
		this.observationsFileName = obsFileName; this.staticObservationsFileName = staticObsFileName;
		this.markovLag = markovLag != null ? markovLag : 1;
		
		if (inputAttributes == null) { // Attributes must be given
			System.out.println("Error!");
			System.exit(1);
		} else {
			this.attributes = inputAttributes;
			this.staticAttributes = inputStaticAttributes;
		}
		readFromFiles();
	}
	
	private void readFromFiles() {
		
		// Parse dynamic observations file
		try {

			// open and parse the dynamic observations csv file
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
			
			// Store the line of each subject in all matrices
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

				subjectLineInMatrices.put(dataLine[0], new Integer(numSubjects[0])); // Store the line of the subject in all matrices

				for (int t = 0; t < numTransitions; t++) {

					String[] transition = Arrays.copyOfRange(dataLine, 1 + t * numAttributes, 1 + (t + markovLag + 1) * numAttributes);
					for (int j = 0; j < (markovLag + 1) * numAttributes; j++) {
						String value = transition[j];
						if (value.length() == 0 || value.equals("?")) { // Missing value, put -1
							usefulObservations[t][numSubjects[t]][j] = -1;
							continue;
						}
						
						int attributeId = j % numAttributes;
						Attribute attribute = attributes.get(attributeId);
						
						if (attribute.hasValue(value) == false) { // Value not in training data, put -2
							usefulObservations[t][numSubjects[t]][j] = -2;
							continue;
						}
						
						usefulObservations[t][numSubjects[t]][j] = attribute.getIndex(value); // Value in training data, put proper value
					}
					numSubjects[t]++;
				}
			}

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
				
				Integer subjectLine = subjectLineInMatrices.get(subject); // Get the subject line (dynamic and static matrices must be coherent!)
				
				if( subjectLine == null) {
					System.err.println("Subject "+ subject + " is not in dynamic observations file! Aborting!");
					System.exit(1);
				}
				
				for (int j = 1; j < (numStaticAttributes + 1); j++) {
					
					String value = dataLine[j];
					int staticAttributeId = j-1;
					
					if (value.length() == 0 || value.equals("?")) { // Missing value, put -1
						staticUsefulObservations[subjectLine][staticAttributeId] = -1;
						continue;
					}
					
					Attribute staticAttribute = staticAttributes.get(staticAttributeId);
					
					if (staticAttribute.hasValue(value) == false) { // Value not in static training data, put -2
						staticUsefulObservations[subjectLine][staticAttributeId] = -2;
						continue;
					}
					
					staticUsefulObservations[subjectLine][staticAttributeId] = staticAttribute.getIndex(value); // Value in training data, put proper value
				}
				
			}
			
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
	
	/**
	 * Check if header line matches a specified the attribute order
	 */
	private boolean checkHeaderFormat(String[] inputHeader, List<Attribute> attributesToCompare) {
		int i=0;
		for(Attribute att : attributesToCompare) {
			if(att.getName().equals(inputHeader[i]) == false ) return false;
			i++;
		}
		return true;
	}
	
	/**
	 * Store the attributes and respective timesteps, on which the user specified inference to be made.
	 * 
	 * @param attfileName 
	 * 		file given by user, with the attributes and respective timesteps
	 */
	public void parseAttributes(String attFileName) {
		
		int newMatrixDimension = -1;
		
		try {
			// open and parse the observations csv file
			CSVReader reader = new CSVReader(new FileReader(attFileName));
			List<String[]> lines = reader.readAll();
			reader.close();
			
			inferenceAtts = new int[lines.size()][2];
			
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
						
						try { // Check if second element of line is an integer
							desiredtimestep = Integer.parseInt(line[1]);
							found = 2;
						} catch (NumberFormatException e) {
							break;
						}
						
						if(desiredtimestep < 0) break; // timestep must be positive
						
						// Check if observation matrix must be increased
						if( (desiredtimestep - markovLag + 1) > usefulObservations.length) {
								if( (desiredtimestep - markovLag + 1) > newMatrixDimension) {
									newMatrixDimension = desiredtimestep - markovLag + 1;
								}
						}

						found=3; // 3 = success
						inferenceAtts[i][0] = attributes.indexOf(att);
						inferenceAtts[i][1] = desiredtimestep;
						
					}
				}
				
				// Check for error conditions
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
			
			if(newMatrixDimension!=-1) // Check if necessary to increase observations matrix size
				usefulObservations = increaseObservationMatrix(usefulObservations, newMatrixDimension);
			
			predictions = new String[subjects.length][inferenceAtts.length];
			
		} catch (IOException e) {
			System.err.println("File " + attFileName + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	/**
	 * Make inference on the specified set of attributes, using a tDBN learned
	 * 
	 * @param stationary 
	 * 		whether the tDBN is stationary or not
	 * @param dbn 
	 * 		the learned tDBN
	 * @param fileToWrite 
	 * 		a file path where to write the inference perfomed. If null, output written to terminal
	 * @param presentDistr 
	 * 		If true outputs all the distribution, if false only outputs most probable value for each attribute
	 * 
	 * @param infFormat
	 * 		If 1 estimates with most probable, otherwise estimates with random sampling according to distribution
	 */
	public void makeInference(boolean stationary, DynamicBayesNet dbn, String fileToWrite, int infFormat) {
		
		StringBuilder sb = new StringBuilder();
		
		int attributeID = 0, i=0, transitionNetID = 0, attToPrintDim = 0, currLineIndx=0;
		
		// Allocate matrices for each static and dynamic configuration of each subject
		int configs[][] = new int[usefulObservations[0].length][(markovLag + 1) * attributes.size()];
		int staticConfigs[][] = ( staticUsefulObservations != null ) ? new int[usefulObservations[0].length][staticAttributes.size()] : null;
		
		// Alocate array for each Configuration object associated to each subject
		Configuration desiredConfigs[] = new Configuration[usefulObservations[0].length];
		
		List<Integer> parents, staticParents;
		Map<Configuration, List<Double>> cpt;
		List<Double> distribution;
		double count, maxProb;
		int canMakeInference[] = new int[usefulObservations[0].length];
		Attribute attToPrint;
		String valueToStore = null;
		
		Random randGen = new Random(); // for generating random numbers
		double rand, cummulative;
		int indxSelected, currIndx;
		
		currLineIndx=-1;
		for(int[] infLine : inferenceAtts) { // For each attribute inference is desired
			currLineIndx++;
			
			attributeID = infLine[0];
			transitionNetID = infLine[1] - markovLag;
			
			if(transitionNetID < 0) { // Error condition
				sb.append("Distributions " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]:\n");
				sb.append("Not possible to determine because timestep " + (transitionNetID+markovLag) + " < markovLag\n\n");
				continue;
			}
			
			if(stationary == false && transitionNetID >= dbn.getNumberTransitionNets()) { // Error condition
				sb.append("Distributions " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]:\n");
				sb.append("Not possible to determine because non-stationary DBN only until timestep " + (dbn.getNumberTransitionNets()+markovLag-1) + "\n\n");
				continue;
			}

			// Get proper cpt
			if (stationary == true) {
				cpt = dbn.getTransitionNet(0).getCPT(attributeID);
			} else {
				cpt = dbn.getTransitionNet(transitionNetID).getCPT(attributeID);
			}
			
			// Get proper dynamic and static parents
			if (stationary == true) {
				parents = dbn.getTransitionNet(0).getParents(attributeID);
				staticParents = dbn.getTransitionNet(0).getStaticParents(attributeID);
			} else {
				parents = dbn.getTransitionNet(transitionNetID).getParents(attributeID);
				staticParents = dbn.getTransitionNet(transitionNetID).getStaticParents(attributeID);
			}
			
			// Initialize the dynamic attributes values (-1 if not parent, 0 if the child attribute)
			for(int[] line : configs) {
				for(i=0; i<line.length; i++)
					line[i] = -1; // By default store an attribute as not being a parent
				
				line[attributes.size() * markovLag + attributeID] = 0; // The child must have 0 in its array position
			}
			
			// If needed (if there are static parents), initialize the static attributes values (-1 if not parent)
			if( staticParents.isEmpty() == false && staticConfigs != null ) {
				for(int[] line : staticConfigs) {
					for(i=0; i<line.length; i++)
						line[i] = -1; // By default store an attribute as not being a parent
				}
			}
				
			i=-1;
			for(int[] observation : usefulObservations[transitionNetID]) {
				i++;
				
				canMakeInference[i] = 0;
				
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
					
					// Check the dynamic parents in the observations file
					for(Integer parent : parents) {
						
						// If dynamic parent not in observations file, must be estimated according to tDBN learned
						if(observation[parent]<0) {
							dbn.getMostProbable(transitionNetID + parent/(attributes.size()) - markovLag, parent%(attributes.size()), stationary, usefulObservations, i, staticUsefulObservations, infFormat);
						}
						
						configs[i][parent] = observation[parent];
						
						if(configs[i][parent] < 0) { // If, accoring to the tDBN, it was not possible to estimate parent value, inference is not possible
							canMakeInference[i] = -1;
							break;
						}
					}
				}
				
				// If inference if possible for this attribute in this timestep, create configuration of parents
				if(canMakeInference[i] == 0) {
					if(staticParents.isEmpty() == false) { // Create configuration using static attributes (node has static parents)
						desiredConfigs[i] = new LocalConfigurationWithStatic(attributes, configs[i], staticAttributes, staticConfigs[i]);
						
					} else { // Create configuration not using static attributes (node does not have static parents)
						desiredConfigs[i] = new LocalConfiguration(attributes, configs[i]);
					}
				}
					
			}
			
			// Get attribute and proper timestep in which inference is being made
			sb.append("Distributions " + attributes.get(attributeID).getName() + "[" + (transitionNetID+markovLag) + "]:\n");
			sb.append("id");
			
			// Get possible values of the specified attribute
			attToPrint = attributes.get(attributeID);
			for(attToPrintDim = 0 ; attToPrintDim < attToPrint.size(); attToPrintDim++) {
				sb.append("," + attToPrint.get(attToPrintDim));
			}
			sb.append("\n");
			
			// For each case given in observation files, get the proper CPT line of the DBN and perform inference using it 
			i=-1;
			for(Configuration desiredConfig : desiredConfigs) {
				
				i++;
				sb.append(subjects[i]);
				if(canMakeInference[i] == -1 ) { // If inference was found to be impossible, just continue
					for(int k=0; k<attToPrintDim; k++) sb.append(",-1");
					
					sb.append("\n");
					
					predictions[i][currLineIndx] = null;
					
					continue;
				}
				
				distribution = cpt.get(desiredConfig);
				
				count = 0.0; 
				for(double d : distribution) {
					sb.append(String.format(",%.3f", d));
					count += d;
				}
				sb.append(String.format(",%.3f\n", 1.0-count));
				
				System.out.println("---------");
				System.out.println("Estimating: " + attributes.get(attributeID).getName() + "[" +  (transitionNetID+markovLag) + "]");
				System.out.println(distribution);
				
				if(infFormat == 1) { // Select most probable
					maxProb = Collections.max(distribution);
					indxSelected = distribution.indexOf(maxProb);
					
					if(maxProb < (1.0-count))
						indxSelected = attToPrintDim-1;
					
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
				System.out.println("---------");
				System.out.println("---------");
				System.out.println("---------");
				
				valueToStore = attributes.get(attributeID).get(indxSelected);
				
				predictions[i][currLineIndx] = valueToStore;
			}
			sb.append("\n");
		}
		
		// If user specified to only print the most probable or random sampled value instead of all distribution, create proper information to be provided to user
		if (infFormat != 2 )
			sb = PredictedValuesFromInference();
		
		// Write information to user
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
	
	/**
	 * Write in csv tabular form the predicted value for each attribute inference was made for each subject
	 * 
	 */
	private StringBuilder PredictedValuesFromInference() {
		
		if (predictions==null || inferenceAtts==null) { // Error condition
			System.out.println("Must call before the proper methods to create and fill the most probable values for the desired attributes!");
			return null;
		}
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("id"); // Write each attribute name and timestep
		for(int infAtt[] : inferenceAtts)
			sb.append(","+ attributes.get(infAtt[0]).getName() + "[" + infAtt[1] + "]");
		sb.append("\n");
		
		// Write the determined value for attribute for each subject
		int i=0;
		for(String[] line : predictions) {
			sb.append(subjects[i++]);
			for(String value : line){
				sb.append("," + value);
			}
			sb.append("\n");
		}
		
		return sb;
	}
	
	/**
	 * Increase observation matrix.
	 * Some values must be propagated to the newly created extra positions, according to the markovLag, remaining values
	 * are initialized with -1.
	 * 
	 */
	private int[][][] increaseObservationMatrix(int[][][] currentMatrix, int newMatrixDimension) {
		
		int transitionNetID, subject, att, n = attributes.size();
		
		// Create new matrix
		int [][][] auxMatrix = new int[newMatrixDimension][subjects.length][(markovLag + 1) * n];
		
		// Copy values
		for(transitionNetID=0; transitionNetID<currentMatrix.length; transitionNetID++)
			for(subject=0; subject<subjects.length; subject++)
				System.arraycopy(currentMatrix[transitionNetID][subject], 0, auxMatrix[transitionNetID][subject], 0, (markovLag + 1) * n);
		
		// Fill with -1 the new positions
		for(transitionNetID=currentMatrix.length; transitionNetID < newMatrixDimension; transitionNetID++)
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * n; att++)
					auxMatrix[transitionNetID][subject][att] = -1;
		
		// Propagate values to the newly created positions, according to markovLag
		
		int maxNetIDtoProp = Math.min((currentMatrix.length + markovLag), newMatrixDimension);
		
		for(transitionNetID = currentMatrix.length; transitionNetID < maxNetIDtoProp; transitionNetID++) {
			for(subject=0; subject < subjects.length; subject++) {
				for(att=0; att < markovLag * n; att++) { // So vai ate aos que nao sao deste timestep
					auxMatrix[transitionNetID][subject][att] = auxMatrix[transitionNetID - 1][subject][att + n];
				}
			}
		}
		
		return auxMatrix;
	}
	
	
	/**
	 * Get most probable trajectory of all attributes until a certain timestep
	 *
	 * @param timestepMax
	 * 		the timestep until which to estimate the attributes values
	 * 
	 * @param dbn
	 * 		the DBN used to predict the values
	 * 
	 * @param stationary
	 * 		whether the dbn is stationary or not
	 * 
	 * @param infFormat
	 * 		If 1 estimates with most probable, otherwise estimates with random sampling according to distribution
	 * 
	 */
	public void getMostProbableTrajectory(int timestepMax, DynamicBayesNet dbn, boolean stationary, int infFormat) {
		
		int transitionNetID, subject, att;
		
		int maxNumTransNet = timestepMax - markovLag + 1;
		
		if(stationary == false) { // Check for error condition
			if(maxNumTransNet > dbn.getNumberTransitionNets()) {
				System.err.println("Cannot predict until timestep " + timestepMax + " because non-stationary DBN learnt with markovLag=" + markovLag + " only has " + dbn.getNumberTransitionNets() +" transition networks!");
				System.exit(1);
			}
		}
		
		// Check if needed to increase size of observation matrix
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
						dbn.getMostProbable(transitionNetID + att/(attributes.size()) - markovLag, att%(attributes.size()), stationary, matrixWithAllPredictions, subject, staticUsefulObservations, infFormat);
		
		// If all values already predicted, return
		if(maxNumTransNet <= usefulObservations.length) return;
		
		// Predict all values of remaining timesteps
		for(transitionNetID=usefulObservations.length; transitionNetID < maxNumTransNet; transitionNetID++)
			for(subject=0; subject < subjects.length; subject++)
				for(att=0; att < (markovLag + 1) * attributes.size(); att++)
					dbn.getMostProbable(transitionNetID + att/(attributes.size()) - markovLag, att%(attributes.size()), stationary, matrixWithAllPredictions, subject, staticUsefulObservations, infFormat);
		
		return;
	}
	
	/**
	 * Print to user, in csv tabular form, all values determined for each variable in the most probable trajectory
	 */
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
	
	/**
	 * Write observation three-dimensional matrix in tabular csv form.
	 */
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
	
}
