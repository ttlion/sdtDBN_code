package com.github.tDBN.dbn;

import java.io.Serializable;

import com.github.tDBN.utils.BidirectionalArray;

public class NominalAttribute implements Attribute, Serializable {
	
	static final long serialVersionUID = 42L;

	private String name;

	private BidirectionalArray<String> values = new BidirectionalArray<String>();

	@Override
	public boolean isNumeric() {
		return false;
	}

	@Override
	public boolean isNominal() {
		return true;
	}

	@Override
	public int size() {
		return values.size();
	}

	@Override
	public boolean add(String value) {
		return values.add(value);
	}

	@Override
	public String toString() {
		return "" + values;
	}

	@Override
	public int getIndex(String value) {
		return values.getIndex(value);
	}

	@Override
	public String get(int index) {
		return values.get(index);
	}

	public void setName(String name) {
		this.name = name;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public boolean hasValue(String value) {
		return values.containsValue(value);
	}

}
