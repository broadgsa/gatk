package org.broadinstitute.sting.indels;

public class StrictlyUpperTriangularMatrix extends SymmetricMatrix {

	public StrictlyUpperTriangularMatrix(int dimension) {
        super(dimension);
		assert dimension >=2 : "Distance matrix can not be smaller than 2x2";
	}
	
	public double get(int i, int j) {
		if ( i >= j ) return 0.0;
		return super.get(i,j);
	}

	public void set(int i, int j, double value) {
		assert  i < j : "Only i < j elements can be set in strictly upper diagonal matrix" ;
        super.set(i,j,value);
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
