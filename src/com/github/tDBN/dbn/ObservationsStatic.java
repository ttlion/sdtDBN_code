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

/**
 * 
 * Class to store the static observations useful for the tDBN learning and manage them.
 * 
 * @author Tiago Leao
 * 
 */
public class ObservationsStatic {

	/**
	 * Three-dimensional matrix of coded static observation data which will be used for
	 * learning a dynamic Bayesian network.
	 * <ul>
	 * <li>the 1st index refers to the transition t+MarkovLag;
	 * <li>the 2nd index refers to the the subject;
	 * <li>the 3rd index refers to the static attribute
	 * </ul>
	 */
	private int[][][] usefulObservations;

	/**
	 * Each column of the useful observation data refers to a static attribute.
	 */
	private List<Attribute> attributes;

	/**
	 * File that contains the static observations that will be converted to attributes and from which one can learn a DBN.
	 * 
	 */
	private String usefulObservationsFileName;

	/**
	 * Header line of input useful observations CSV file.
	 */
	private String[] usefulObservationsHeader;
	
	private int[] numSubjects;

	public ObservationsStatic(String usefulObsFileName, Map<String, int[]> subjectLinePerTemporalMatrix, int numTransInTemporal, int numTempSubjects) {
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
			int totalNumSubjects = numTempSubjects;
			
			if(totalNumSubjects==0) {
				System.out.println("There are 0 temporal subjects! Aborting");
				System.exit(1);
			}

			usefulObservations = new int[numTransInTemporal][totalNumSubjects][numAttributes];
			numSubjects = new int[numTransInTemporal];
			
			for(int i=0; i< usefulObservations.length; i++) { // Init usefulObservations to -1 in all positions
				for(int j=0; j<usefulObservations[i].length; j++) {
					for(int k=0; k<usefulObservations[i][j].length; k++) {
						usefulObservations[i][j][k] = -1;
					}
				}
			}

			String[] dataLine = li.next();

			// fill attributes from first subject (it must not have missing values)
			String[] firstObservation = Arrays.copyOfRange(dataLine, 1, numAttributes + 1);
			if (countMissingValues(firstObservation) > 0) {
				System.err.println(firstObservation);
				System.err.println("First subject contains missing values.");
				System.exit(1);
			}
			int i = 0;
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
				
				// Get array with the lines of each temporal matrix where the subject is
				int subjectLinesPerSlice[] = subjectLinePerTemporalMatrix.get(subject);
				
				if( subjectLinesPerSlice == null) {
					System.err.println("Subject "+ subject + " is not in temporal file! Aborting!");
					System.exit(1);
				}
				
				for (int t = 0; t < numTransInTemporal; t++) {

					int subjectLine = subjectLinesPerSlice[t]; // Get the proper line of the subject in this temporal transition
					
					// if subjectLine == -1, do not put in the variables associated to this 
					// timestep (cannot use static att in this timestep because there where no dynamic att used)
					if(subjectLine == -1) {
						continue;
					}
					
					for (int j = 0; j < numAttributes; j++) {
						String value = dataLine[j+1];
						
						if( value.length() == 0 || value.equals("?") == true ) { // If there is no static attribute value given (missing value), put -1
							usefulObservations[t][subjectLine][j] = -1;
							continue;
						}
						
						int attributeId = j % numAttributes;
						Attribute attribute = attributes.get(attributeId);
						attribute.add(value);
						usefulObservations[t][subjectLine][j] = attribute.getIndex(value); // Store the index of the value of the static attribute
					}
					numSubjects[t]++;
				}
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

	public int numAttributes() {
		return attributes.size();
	}

	public List<Attribute> getAttributes() {
		return attributes;
	}

	public int[][][] getObservationsMatrix() {
		return usefulObservations;
	}

}
