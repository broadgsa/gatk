package org.broadinstitute.sting.playground.fourbasecaller;

import java.lang.String;

public enum Nucleotide {
    A, C, G, T;

    Nucleotide complement() {
        switch (this) {
            case A: return T;
            case C: return G;
            case G: return C;
            case T: return A;
        }

        throw new IllegalArgumentException("Nucleotide must be A, C, G, or T.");
    }

    public char asChar() {
        switch(this) {
            case A: return 'A';
            case C: return 'C';
            case G: return 'G';
            case T: return 'T';
        }

        throw new IllegalArgumentException("Nucleotide must be A, C, G, or T.");
    }

    public String toString() {
        char[] value = new char[1];
        value[0] = this.asChar();

        return new String(value);
    }
}
