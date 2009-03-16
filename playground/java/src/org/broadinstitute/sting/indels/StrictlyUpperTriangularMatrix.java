package org.broadinstitute.sting.indels;

public class StrictlyUpperTriangularMatrix {
	private double [] data;
	private int size;
	
	public StrictlyUpperTriangularMatrix(int dimension) {
		assert dimension >=2 : "Distance matrix can not be smaller than 2x2";
		if ( dimension % 2 == 0 ) {
				int k = dimension >> 1; // dimension/2
				data = new double[k*(dimension-1)];
		} else {
				int k = ( dimension -1 ) >> 1; // (dimension -1)/2
				data = new double[k*dimension];
		}
		size = dimension;
	}
	
	public double get(int i, int j) {
		assert (i < size) && (j<size) : "Out of bound index into distance matrix detected";
		if ( i >= j ) return 0.0;
		
		// we are guaranteed now that i < j; now
		// translate i,j into the linear offset into our internal data array and return the value:
		return data[ linearOffset(i,j) ];
	}

	public void set(int i, int j, double value) {
		assert (i < size) && (j<size) : "Out of bound index into distance matrix detected";
		assert  i < j : "Only i < j elements can be set in strictly upper diagonal matrix" ;
		
		// we are guaranteed now that i < j; now
		// translate i,j into the linear offset into our internal data array and set the value
		data[ linearOffset(i,j) ] = value;
	}
	
	public int size() { return size; }
	
	/** Returns ready-to-print string representing the full matrix (don't use for 1000x1000 matrices!!);
	 *  each element is formatted according to a default format. 
	 * @return
	 * @see format(String f)
	 */
	public String format() {
		return format("%6.3f ");
	}
	
	/** Returns ready-to-print string representing the full matrix (don't use for 1000x1000 matrices!!);
	 *  each element is formatted according to a specified format string (note: format string must include all desired
	 *  whitespaces before and/or after an element, as this method itself does not add any spaces between the formatted elements). 
	 * @return
	 * @see format()
	 */
	public String format(String f) {
		StringBuilder b = new StringBuilder();
		java.util.Formatter frm = new java.util.Formatter(b);
		for ( int i = 0 ; i < size ; i++ ) {
			for ( int j = 0 ; j < size ; j++ ) {
				if ( i < j ) frm.format(f, get(i,j));
				else frm.format(f, 0.0);
			}
			b.append('\n');
		}
		return b.toString();
	}
	
	/** Computes linear offset into the array that keeps actual data given "row" and "column" indices into the matrix
	 * this class represents; this method is unchecked, but it expects i < j otherwise the result will be incorrect.
	 * @param i row index
	 * @param j column index
	 * @return linear offset into the data[] member of this class
	 */
	private int linearOffset(int i, int j) {
		int k = (( size << 1 ) - i - 1)*i; // [ 2*d - (i+1) ] * i 
		k >>= 1; // k/=2
		// now k is the offset of the first stored element in row i
		return ( k + (j - i - 1));
	}

	private static void testMe() {
		StrictlyUpperTriangularMatrix m = new StrictlyUpperTriangularMatrix(3);
		
		m.set(0,1,0.54321);
		m.set(0,2,0.43215);
		m.set(1,2,0.321);

		System.out.println( m.format());
	}
	
	public static void main(String[] argv) {
		testMe();
	}
}
