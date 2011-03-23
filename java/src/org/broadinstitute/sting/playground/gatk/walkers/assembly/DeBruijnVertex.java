package org.broadinstitute.sting.playground.gatk.walkers.assembly;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */
// simple node class for storing kmer sequences
public class DeBruijnVertex {

    // used for equals()
    protected byte[] actualSequence;

    // used for printing and traversing graphs
    protected byte[] printableSequence;

    public DeBruijnVertex(byte[] sequence) {
        actualSequence = sequence;
        printableSequence = new byte[sequence.length];
        System.arraycopy(sequence, 0, printableSequence, 0, sequence.length);
    }

    public boolean equals(DeBruijnVertex v) {
        return Arrays.equals(actualSequence, v.actualSequence);
    }

    public String toString() {
        return new String(printableSequence);
    }    

    public void addPrefix(byte[] prefix, boolean justPrintableSequence) {
        printableSequence = addPrefix(printableSequence, prefix);
        if ( !justPrintableSequence )
            actualSequence = addPrefix(actualSequence, prefix);
    }

    private static byte[] addPrefix(byte[] sequence, byte[] prefix) {
        byte[] newSequence = new byte[sequence.length + prefix.length];
        System.arraycopy(prefix, 0, newSequence, 0, prefix.length);
        System.arraycopy(sequence, 0, newSequence, prefix.length, sequence.length);
        return newSequence;
    }

    public void removePrefix(int prefixLength, boolean justPrintableSequence) {
        printableSequence = removePrefix(printableSequence, prefixLength);
        if ( !justPrintableSequence )
            actualSequence = removePrefix(actualSequence, prefixLength);
    }

    private static byte[] removePrefix(byte[] sequence, int prefixLength) {
        int newLength = sequence.length - prefixLength;
        byte[] newSequence = new byte[newLength];
        System.arraycopy(sequence, prefixLength, newSequence, 0, newLength);
        return newSequence;
    }

    public void removeSuffix(int suffixLength, boolean justPrintableSequence) {
        printableSequence = removeSuffix(printableSequence, suffixLength);
        if ( !justPrintableSequence )
            actualSequence = removeSuffix(actualSequence, suffixLength);
    }

    private static byte[] removeSuffix(byte[] sequence, int suffixLength) {
        int newLength = sequence.length - suffixLength;
        byte[] newSequence = new byte[newLength];
        System.arraycopy(sequence, 0, newSequence, 0, newLength);
        return newSequence;
    }
}
