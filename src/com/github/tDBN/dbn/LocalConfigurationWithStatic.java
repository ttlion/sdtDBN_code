package com.github.tDBN.dbn;

import java.util.Arrays;
import java.util.List;

/**
 * 
 * Extends the LocalConfiguration to represent a configuration also considering static parents
 * 
 * @author Tiago Leao
 * 
 */
public class LocalConfigurationWithStatic extends LocalConfiguration {
	
	private List<Attribute> staticAttributes;

	private int[] staticConfiguration;
	
	private int[] staticParentIndices;
	
	/**
	 * For easily cloning a LocalConfigurationWithStatic object.
	 * 
	 */
	public LocalConfigurationWithStatic(LocalConfigurationWithStatic original) {
		super(original);
		this.staticAttributes = original.staticAttributes;
		this.staticConfiguration = original.staticConfiguration.clone();
		this.staticParentIndices = original.staticParentIndices.clone();		
	}
	
	/**
	 * Allocates LocalConfigurationWithStatic only based on lists of dynamic and static attributes and the configuration arrays of static and dynamic attributes.
	 * 
	 * @param attributes
	 *            list of dynamic  attributes characterizing the nodes
	 * @param configuration
	 *            array of integers, each integer representing a value of the proper dynamic attribute. -1 means the attribute is not in configuration.
	 * @param staticAttributes
	 *            list of static  attributes characterizing the nodes
	 * @param staticConfiguration
	 *            array of integers, each integer representing a value of the proper static attribute. -1 means the attribute is not in configuration.
	 */
	public LocalConfigurationWithStatic(List<Attribute> attributes, int[] configuration, List<Attribute> staticAttributes, int[] staticConfiguration) {
		super(attributes,configuration);
		this.staticAttributes = staticAttributes;
		this.staticConfiguration = staticConfiguration;
	}
	
	public LocalConfigurationWithStatic(List<Attribute> attributes, int markovLag, List<Integer> parentNodesPast,
			int childNode, List<Attribute> staticAttributes, List<Integer> staticParentSet) {
		this(attributes, markovLag, parentNodesPast, null, childNode, staticAttributes, staticParentSet);
	}
	
	/**
	 * Calls LocalConfiguration constructor to allocate the dynamic parents configuration array.
	 * Allocates static parents configuration array, setting all parents to their first value.
	 * 
	 */
	public LocalConfigurationWithStatic(List<Attribute> attributes, int markovLag, List<Integer> parentNodesPast,
			Integer parentNodesPresent, int childNode, List<Attribute> staticAttributes, List<Integer> staticParentSet) {
		
		super(attributes, markovLag, parentNodesPast, parentNodesPresent, childNode);
		
		this.staticAttributes = staticAttributes;
		
		int numStaticParents = (staticParentSet != null) ? staticParentSet.size() : 0;

		staticParentIndices = new int[numStaticParents];
		staticConfiguration = new int[staticAttributes.size()];
		
		for(int j=0; j < staticAttributes.size(); j++) { // Reset staticConfiguration array
			staticConfiguration[j]=-1; 
		}
		
		// Store indices of static parents in proper array (just as it is done with dynamic parents)
		int i = 0;
		if (staticParentSet != null)
			for (Integer parentNode : staticParentSet) {
				staticParentIndices[i++] = parentNode;
			}
		
		// Only needed to reset static parents (LocalConfiguration constructor reseted dynamic parents)
		resetStaticParents();
		
	}
	
	/**
	 * Method overridden to also reset the static parents indices.
	 */
	@Override
	public void resetParents() {
		
		for (int i = 0; i < parentIndices.length; i++) {
			configuration[parentIndices[i]] = 0;
		}
		
		if(staticConfiguration!=null && staticParentIndices!= null)
			resetStaticParents();
	}
	
	/**
	 * Method to only reset static parents.
	 */
	public void resetStaticParents() {
		for (int i = 0; i < staticParentIndices.length; i++) {
			staticConfiguration[staticParentIndices[i]] = 0;
		}
	}
	
	/**
	 * Method overridden to also take into account the static parents possible values.
	 */
	@Override
	public int getParentsRange() {
		if (parentIndices.length == 0 && staticParentIndices.length == 0 ) {
			return 0;
		}
		
		int n = attributes.size();
		int n_static = staticAttributes.size();
		
		int result = 1;
		for (int i = 0; i < parentIndices.length; i++) {
			result *= attributes.get(parentIndices[i] % n).size();
		}
		for (int i = 0; i < staticParentIndices.length; i++) {
			result *= staticAttributes.get(staticParentIndices[i] % n_static).size();
		}
		
		return result;
	}
	
	/**
	 * Method overridden to also take into account the static parents.
	 */
	@Override
	public int getNumParameters() {
		return getParentsRange() * (getChildRange() - 1);
	}
	
	/**
	 * Method overridden to also take into account the static parents.
	 */
	@Override
	public boolean nextParents() { 
		
		int n = attributes.size();
		int n_static = staticAttributes.size();
		
		// Try to change the value of a dynamic parent
		for (int i = 0; i < parentIndices.length; i++) {
			if (++configuration[parentIndices[i]] < attributes.get(parentIndices[i] % n).size()) {
				return true;
			} else {
				configuration[parentIndices[i]] = 0;
			}
		}
		
		// If all dynamic parents values already tried, try to change the value of a static parent
		for (int i = 0; i < staticParentIndices.length; i++) {
			if (++staticConfiguration[staticParentIndices[i]] < staticAttributes.get(staticParentIndices[i] % n_static).size()) {
				break;
			} else {
				staticConfiguration[staticParentIndices[i]] = 0;
				if (i == staticParentIndices.length - 1) { // Case where all combinations of static and dynamic parents were checked
					resetParents();	
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Method overridden to also take into account the static parents.
	 */
	@Override
	public boolean matches(int[] observationDyn, int[] observationStatic) {
		
		
		int n = attributes.size();
		
		// Check dynamic parents values
		for (int i = 0; i < configuration.length; i++) {
			if (configuration[i] > -1) {
				if (observationDyn[i] != configuration[i]) {
					if (considerChild || i != childNode + n * markovLag) {
						return false;
					}
				}
			}
		}
		
		// Check static parents values
		for (int i = 0; i < staticConfiguration.length; i++) {
			if (staticConfiguration[i] > -1) {
				if (observationStatic[i] != staticConfiguration[i]) {
						return false;
				}
			}
		}
		
		return true;
	}
	
	/**
	 * Method overridden to also take into account the static parents.
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();	
		
		sb.append("[");
		
		// Print static config distribution
		for (int i = 0; i < staticConfiguration.length; i++) {
			if (staticConfiguration[i] != -1) {
				sb.append(staticAttributes.get(i).getName() + "=" + staticAttributes.get(i).get(staticConfiguration[i]));
				sb.append(", ");
			}
		}
		
		// Print dynamic config distribution
		int n = attributes.size();
		for (int i = 0; i < configuration.length; i++) {
			if (configuration[i] != -1 && i != n * markovLag + childNode) {
				int lag = i / n;
				int id = i % n;
				sb.append(attributes.get(id).getName() + "[" + lag + "]=" + attributes.get(id).get(configuration[i]));
				sb.append(", ");
			}
		}
		
		// Readable version
		if (sb.length() > 0) {
			sb.setLength(sb.length() - 2);
		}
		
		sb.append("]");
		return sb.toString();
		
	}
	
	/**
	 * Method overridden to also take into account the static configuration parents.
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + Arrays.hashCode(staticConfiguration);
		return result;
	}
	
	/**
	 * Checks the dynamic configuration and then checks the static configuration also.
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		LocalConfigurationWithStatic other = (LocalConfigurationWithStatic) obj;
		if (!Arrays.equals(staticConfiguration, other.staticConfiguration))
			return false;
		return true;
	}
	
}
