package org.broadinstitute.sting.indels;

public class Matrix<T> {
	private int size;
	private Object [][] data;
	
	public Matrix(int n) {
		size = n;
		data = new Object[n][n];
	}
	
	@SuppressWarnings("unchecked")
	public T get(int i, int j) {
		assert (i<size) && (j < size) : "Matrix index is out of bounds"; 
		return (T) data[i][j]; 
	}

	public void set(int i, int j, T value) {
		assert (i<size) && (j < size) : "Matrix index is out of bounds"; 
		data[i][j] = value;
	}
}
