package com.github.tDBN.dbn;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;

import com.github.tDBN.utils.Utils;

public class ObservationsStatic {

	/**
	 * Three-dimensional matrix of coded observation data which will be used for
	 * learning a dynamic Bayesian network.
	 * <ul>
	 * <li>the 1st index refers to the transition t+MarkovLag;
	 * <li>the 2nd index refers to the the subject;
	 * <li>the 3rd index refers to the static attribute
	 * </ul>
	 */
	private int[][][] usefulObservations;

	/**
	 * Each column of the useful observation data refers to an attribute.
	 */
	private List<Attribute> attributes;

	/**
	 * File that contains the static observations that will be converted to attributes and
	 * from which one can learn a DBN.
	 */
	private String usefulObservationsFileName;

	/**
	 * Header line of input useful observations CSV file.
	 */
	private String[] usefulObservationsHeader;
	
	private int[] numSubjects;

	public ObservationsStatic(String usefulObsFileName, Map<String, boolean[]> subjectIsPresent, int numTransInTemporal) {
		this.usefulObservationsFileName = usefulObsFileName;
		
		try {

			// open and parse the useful observations csv file
			CSVReader reader = new CSVReader(new FileReader(usefulObservationsFileName));
			List<String[]> lines = reader.readAll();
			reader.close();

			ListIterator<String[]> li = lines.listIterator();

			// get first line
			String[] header = li.next();

			int numAttributes = (header.length - 1);
			attributes = new ArrayList<Attribute>(numAttributes);

			usefulObservationsHeader = processHeader(header, numAttributes);

			// allocate observations matrix
			int totalNumSubjects = lines.size() - 1;
			usefulObservations = new int[numTransInTemporal][totalNumSubjects][numAttributes];
			numSubjects = new int[numTransInTemporal];

			String[] dataLine = li.next();

			// fill attributes from first subject (it must not have missing values)
			String[] firstObservation = Arrays.copyOfRange(dataLine, 1, numAttributes + 1);
			if (countMissingValues(firstObservation) > 0) {
				System.err.println(firstObservation);
				System.err.println("First subject contains missing values.");
				System.exit(1);
			}
			int i = 0, numsubjects = 0;
			for (String value : firstObservation) {
				Attribute attribute;
				// numeric attribute
				if (Utils.isNumeric(value))
					attribute = new NumericAttribute();
				// nominal attribute
				else
					attribute = new NominalAttribute();
				attribute.setName(usefulObservationsHeader[i++]);
				attributes.add(attribute);
			}

			// rewind one line
			li.previous();

			while (li.hasNext()) {

				dataLine = li.next();

				// check for line sanity
				if (dataLine.length != numAttributes + 1) {
					System.err.println(Arrays.deepToString(dataLine));
					System.err
							.println("Observations file: input data line does not have the correct number of columns.");
					System.err.println("Line length: " + dataLine.length);
					System.err.println("Number of attributes: " + numAttributes);
					System.exit(1);
				}

				// record subject id
				String subject = dataLine[0];
				
				boolean existsInTemporal[] = subjectIsPresent.get(subject);
				
				if( existsInTemporal == null) {
					System.err.println("Subject "+ subject + " is not in temporal file! Aborting!");
					System.exit(1);
				}
				
				for (int t = 0; t < numTransInTemporal; t++) {

					boolean observationsOk = existsInTemporal[t];
					
					// if obsservationsOk==false, do not put in the variables associated to this 
					// timestep (cannot use static att in this timestep because there where no dynamic att used)
					if(observationsOk == false) {
						continue;
					}
					
					for (int j = 0; j < numAttributes; j++) {
						String value = dataLine[j+1];
						
						if( value.equals("?") == true ) { // Vamos ver se consigo pondo aqui -1 nos missings consigo depois calcular as probabilidades bem
							usefulObservations[t][numSubjects[t]][j] = -1;
							continue;
						}
						
						int attributeId = j % numAttributes;
						Attribute attribute = attributes.get(attributeId);
						attribute.add(value);
						usefulObservations[t][numSubjects[t]][j] = attribute.getIndex(value);
					}
					numSubjects[t]++;
				}
				
			}
			
			i=0;
			System.out.println("Matriz estatica:");
			for(int[][] matriz : usefulObservations) {
				for(int[] linha : matriz) {
					for(int value : linha) {
						if(value != -1)
							System.out.print(attributes.get(i).get(value) + "  ");
						else
							System.out.print(-1);
						i++;
					}
					i=0;
					System.out.println("");
				}
				System.out.println("");
			}

		} catch (IOException e) {
			System.err.println("File " + usefulObservationsFileName + " could not be opened.");
			e.printStackTrace();
			System.exit(1);
		}
	}

	/**
	 * Gets the name of the attributes from an input header line and the number
	 * of attributes.
	 */
	private String[] processHeader(String[] header, int numAttributes) {
		String[] newHeader = new String[numAttributes];
		String stripFirstHeader[] = Arrays.copyOfRange(header, 1, numAttributes + 1);
		int i = 0;
		for (String name : stripFirstHeader) {
			newHeader[i++] = new String(name);
		}
		return newHeader;
	}
	
	private static int countMissingValues(String[] dataLine) {

		int missing = 0;

		for (String value : dataLine)
			if (value.length() == 0 || value.equals("?"))
				missing++;

		return missing;
	}
	
	private boolean checkIfInTemporal(String subject, String subjectsInTemporal[]) {
		for(String subjInTemp : subjectsInTemporal )
			if(subjInTemp.equals(subject) ) 
				return true;
		
		return false;
	}

//	public int numObservations(int transition) {
//
//		// stationary process
//		if (transition < 0) {
//			int numObs = 0;
//			int T = numTransitions();
//			for (int t = 0; t < T; t++)
//				numObs += numSubjects[t];
//			return numObs;
//		}
//
//		// time-varying process
//		return numSubjects[transition];
//	}

	public int numAttributes() {
		return attributes.size();
	}

	public List<Attribute> getAttributes() {
		return attributes;
	}

	public int[][][] getObservationsMatrix() {
		return usefulObservations;
	}

//	/**
//	 * Given a network configuration (parents and child values), counts all
//	 * observations in some transition that are compatible with it. If
//	 * transition is negative, counts matches in all transitions.
//	 */
//	public int count(LocalConfiguration c, int transition) {
//
//		// stationary process
//		if (transition < 0) {
//			int allMatches = 0;
//			int T = numTransitions();
//			for (int t = 0; t < T; t++)
//				allMatches += count(c, t);
//			return allMatches;
//		}
//
//		// time-varying process
//		int matches = 0;
//		int N = numObservations(transition);
//		for (int i = 0; i < N; i++)
//			if (c.matches(usefulObservations[transition][i]))
//				matches++;
//		return matches;
//	}

}
