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

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * A list of CigarElements, which describes how a read aligns with the reference.
 * E.g. the Cigar string 10M1D25M means
 * * match or mismatch for 10 bases
 * * deletion of 1 base
 * * match or mismatch for 25 bases
 */
public class Cigar {
    private final List<CigarElement> cigarElements = new ArrayList<CigarElement>();

    public Cigar() {
    }

    public Cigar(final List<CigarElement> cigarElements) {
        this.cigarElements.addAll(cigarElements);
    }

    public List<CigarElement> getCigarElements() {
        return Collections.unmodifiableList(cigarElements);
    }

    public CigarElement getCigarElement(final int i) {
        return cigarElements.get(i);
    }

    public void add(final CigarElement cigarElement) {
        cigarElements.add(cigarElement);
    }

    public int numCigarElements() {
        return cigarElements.size();
    }

    public int getReferenceLength() {
        int length = 0;
        for (CigarElement element : cigarElements) {
            switch (element.getOperator()) {
                case M:
                case D:
                case N:
                    length += element.getLength();
            }
        }
        return length;
    }

    public int getPaddedReferenceLength() {
        int length = 0;
        for (CigarElement element : cigarElements) {
            switch (element.getOperator()) {
                case M:
                case D:
                case N:
                case P:
                    length += element.getLength();
            }
        }
        return length;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof Cigar)) return false;

        final Cigar cigar = (Cigar) o;

        if (cigarElements != null ? !cigarElements.equals(cigar.cigarElements) : cigar.cigarElements != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        return cigarElements != null ? cigarElements.hashCode() : 0;
    }
}
