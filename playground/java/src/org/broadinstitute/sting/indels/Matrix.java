package org.broadinstitute.sting.indels;

public class Matrix<T> {
	private int nRows;
    private int nCols;
	private Object [][] data;

    /** Instantiates a generic matrix of objects with n rows and m cols.
     *
     * @param n number of rows
     * @param m number of columns
     */
	public Matrix(int n, int m) {
		nRows = n;
        nCols = m;
		data = new Object[n][m];
	}

    /** Instantiates a square n x n matrix of objects.
     *
      * @param n size of the matrix
     */
    public Matrix (int n) {
        this(n,n);
    }

	@SuppressWarnings("unchecked")
	public T get(int i, int j) {
		assert (i < nRows ) && (j < nCols) : "Matrix index is out of bounds";
		return (T) data[i][j]; 
	}

	public void set(int i, int j, T value) {
		assert ( i < nRows ) && ( j < nCols ) : "Matrix index is out of bounds"; 
		data[i][j] = value;
	}
}
