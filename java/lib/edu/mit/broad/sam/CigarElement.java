/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

/**
 * One component of a cigar string.  The component comprises the operator, and the number of bases to which
 * the  operator applies.
 */
public class CigarElement {
    private final int length;
    private final CigarOperator operator;

    public CigarElement(final int length, final CigarOperator operator) {
        this.length = length;
        this.operator = operator;
    }

    public int getLength() {
        return length;
    }

    public CigarOperator getOperator() {
        return operator;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof CigarElement)) return false;

        final CigarElement that = (CigarElement) o;

        if (length != that.length) return false;
        if (operator != that.operator) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = length;
        result = 31 * result + (operator != null ? operator.hashCode() : 0);
        return result;
    }
}
