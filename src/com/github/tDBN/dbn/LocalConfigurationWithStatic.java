package com.github.tDBN.dbn;

import java.util.Arrays;
import java.util.List;

public class LocalConfigurationWithStatic extends LocalConfiguration {
	
	// Meter variaveis para as configuracoes dos pais static
	
	private List<Attribute> staticAttributes;

	private int[] staticConfiguration;
	
	private int[] staticParentIndices;
	
	public LocalConfigurationWithStatic(LocalConfigurationWithStatic original) {
		super(original);
		this.staticAttributes = original.staticAttributes;
		this.staticConfiguration = original.staticConfiguration.clone();
		this.staticParentIndices = original.staticParentIndices.clone();		
	}
	
	public LocalConfigurationWithStatic(List<Attribute> attributes, int[] configuration, List<Attribute> staticAttributes, int[] staticConfiguration) {
		super(attributes,configuration);
		this.staticAttributes = staticAttributes;
		this.staticConfiguration = staticConfiguration;
	}
	
	public LocalConfigurationWithStatic(List<Attribute> attributes, int markovLag, List<Integer> parentNodesPast,
			int childNode, List<Attribute> staticAttributes, List<Integer> staticParentSet) {
		this(attributes, markovLag, parentNodesPast, null, childNode, staticAttributes, staticParentSet);
	}
	
	public LocalConfigurationWithStatic(List<Attribute> attributes, int markovLag, List<Integer> parentNodesPast,
			Integer parentNodesPresent, int childNode, List<Attribute> staticAttributes, List<Integer> staticParentSet) {
		
		super(attributes, markovLag, parentNodesPast, parentNodesPresent, childNode);
		
		this.staticAttributes = staticAttributes;
		
		int numStaticParents = (staticParentSet != null) ? staticParentSet.size() : 0;

		staticParentIndices = new int[numStaticParents];
		staticConfiguration = new int[staticAttributes.size()];
		
		for(int j=0; j < staticAttributes.size(); j++) { // Reset staticConfiguration
			staticConfiguration[j]=-1; 
		}
		
		int i = 0;

		if (staticParentSet != null)
			for (Integer parentNode : staticParentSet) { // ATENCAO!! E AQUI QUE TOU A POR O INDICE DO PAI ESTATICO!!
				staticParentIndices[i++] = parentNode;
			}
		
		resetStaticParents(); // So preciso de reset aos static porque os nao static ja tao reset (constructor de LocalConfiguration foi chamado)
		
	}
	
	@Override
	public void resetParents() {
		
		for (int i = 0; i < parentIndices.length; i++) {
			configuration[parentIndices[i]] = 0;
		}
		
		if(staticConfiguration!=null && staticParentIndices!= null)
			resetStaticParents();
	}
	
	public void resetStaticParents() {
		for (int i = 0; i < staticParentIndices.length; i++) {
			staticConfiguration[staticParentIndices[i]] = 0;
		}
	}
	
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
	
	@Override
	public int getNumParameters() {
		return getParentsRange() * (getChildRange() - 1);
	}
	
	@Override
	public boolean nextParents() { 
	
		//TODO CONFIRMAR QUE ESTA BEM
//		System.out.print("Pais antes:");
//		for(int value : staticConfiguration)
//			System.out.print(value + ", ");
//		for(int value : configuration)
//			System.out.print(value + ", ");
//		System.out.println("");
		
		
		int n = attributes.size();
		int n_static = staticAttributes.size();

		for (int i = 0; i < parentIndices.length; i++) {
			if (++configuration[parentIndices[i]] < attributes.get(parentIndices[i] % n).size()) {
//				System.out.print("Pais dps 1oFOR:");
//				for(int value : staticConfiguration)
//					System.out.print(value + ", ");
//				for(int value : configuration)
//					System.out.print(value + ", ");
//				System.out.println("");
				return true;
			} else {
				configuration[parentIndices[i]] = 0;
			}
		}
		
		for (int i = 0; i < staticParentIndices.length; i++) {
			if (++staticConfiguration[staticParentIndices[i]] < staticAttributes.get(staticParentIndices[i] % n_static).size()) {
				break;
			} else {
				staticConfiguration[staticParentIndices[i]] = 0;
				if (i == staticParentIndices.length - 1) {
					resetParents(); // SO AQUI DEPOIS DE SE VER OS STATIC TB E QUE SE FAZ RESET DOS PARENTS!!
					
//					System.out.print("Pais dps 2oFOR:");
//					for(int value : staticConfiguration)
//						System.out.print(value + ", ");
//					for(int value : configuration)
//						System.out.print(value + ", ");
//					System.out.println("");
					
					return false;
				}
			}
		}
		
//		System.out.print("Pais dps ENDs:");
//		for(int value : staticConfiguration)
//			System.out.print(value + ", ");
//		for(int value : configuration)
//			System.out.print(value + ", ");
//		System.out.println("");
		
		return true;
	}
	
	@Override
	public boolean matches(int[] observationDyn, int[] observationStatic) {
		
		//TODO ESTARA TUDO?
		
		int n = attributes.size();
//		int n_static = staticAttributes.size();

		for (int i = 0; i < configuration.length; i++) {
			if (configuration[i] > -1) {
				if (observationDyn[i] != configuration[i]) {
					if (considerChild || i != childNode + n * markovLag) {
						return false;
					}
				}
			}
		}
		
		for (int i = 0; i < staticConfiguration.length; i++) {
			if (staticConfiguration[i] > -1) {
				if (observationStatic[i] != staticConfiguration[i]) {
						return false;
				}
			}
		}
		
		return true;
		
	}
	
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

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + Arrays.hashCode(staticConfiguration);
		return result;
	}

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
