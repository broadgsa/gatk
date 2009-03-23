package org.broadinstitute.sting.indels;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 22, 2009
 * Time: 3:43:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class SymmetricMatrix {

    protected double [] data;
    protected int size;

    public SymmetricMatrix(int dimension) {
        assert dimension >= 0 : "Matrix size can not be negative";
        if ( dimension % 2 == 0 ) {
                int k = dimension >> 1; // dimension/2
                data = new double[k*(dimension+1)];
        } else {
                int k = ( dimension + 1 ) >> 1; // (dimension + 1)/2
                data = new double[k*dimension];
        }
        size = dimension;
    }

    public double get(int i, int j) {
        assert (i < size) && ( j < size) : "Out of bound index into matrix detected";
        if ( i >= j ) return data[linearOffset(j,i)]; // we store only the upper triangle in memory		
        return data[ linearOffset(i,j) ];
    }

    public void set(int i, int j, double value) {
        assert (i < size) && (j < size) : "Out of bound index into matrix detected";

        if ( i >= j ) data[ linearOffset(j,i) ] = value;
        else data[ linearOffset(i,j) ] = value;
    }

    public int size() { return size; }

    public int nRows() { return size; }
    public int nCols() { return size; }

    /** Returns ready-to-print string representing the full matrix (don't use for 1000x1000 matrices!!);
     *  each element is formatted according to a default format.
     * @return
     * @see #format(String f)
     */
    public String format() {
        return format("%6.3f ");
    }


    /** Returns ready-to-print string representing the full matrix (don't use for 1000x1000 matrices!!);
     *  each element is formatted according to a specified format string (note: format string must include all desired
     *  whitespaces before and/or after an element, as this method itself does not add any spaces between the formatted elements).
     * @return
     * @see #format()
     */
    public String format(String f) {
        StringBuilder b = new StringBuilder();
        java.util.Formatter frm = new java.util.Formatter(b);
        for ( int i = 0 ; i < size ; i++ ) {
            for ( int j = 0 ; j < size ; j++ ) {
                frm.format(f, get(i,j));
            }
            b.append('\n');
        }
        return b.toString();
    }


    /** Computes linear offset into the internal array that keeps actual data, given "row" and "column" indices
     * into the matrix; this method is unchecked, but it expects i <= j otherwise the result is unspecified.
     * @param i row index
     * @param j column index
     * @return linear offset into the data[] member of this class
     */
    private int linearOffset(int i, int j) {
        int k = (( size << 1 ) - i + 1)*i; // [ 2*d - (i+1) ] * i
        k >>= 1; // k/=2
        // now k is the offset of the first stored element in row i
        return ( k + j - i );
    }

}
