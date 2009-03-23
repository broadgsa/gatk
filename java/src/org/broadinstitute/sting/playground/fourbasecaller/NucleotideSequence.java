package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.playground.fourbasecaller.Nucleotide;

import java.lang.String;
import java.util.Vector;

class NucleotideSequence extends java.util.Vector<Nucleotide> {
    public NucleotideSequence(int initialCapacity) {
        this.setSize(initialCapacity);
    }

    public NucleotideSequence(NucleotideSequence ns) {
        this.setSize(ns.size());

        for (int locus = 0; locus < ns.size(); locus++) {
            switch (ns.get(locus)) {
                case A: this.set(locus, Nucleotide.A); break;
                case C: this.set(locus, Nucleotide.C); break;
                case G: this.set(locus, Nucleotide.G); break;
                case T: this.set(locus, Nucleotide.T); break;
            }
        }
    }

    public NucleotideSequence reverseComplement() {
        if (this.size() % 2 != 0) {
            this.get(this.size()/2 + 1).complement();
        }

        for (int locus = 0; locus < this.size()/2; locus++) {
            Nucleotide temp = this.get(locus).complement();
            this.set(locus, this.get(this.size() - locus - 1).complement());
            this.set(this.size() - locus - 1, temp);
        }

        return this;
    }

    public int editDistance(NucleotideSequence ns) { return editDistance(ns, this.size() < ns.size() ? this.size() : ns.size()); }

    public int editDistance(NucleotideSequence ns, int searchLength) {
        int editDistance = 0;
        for (int locus = 0; locus < searchLength; locus++) {
            editDistance += (this.get(locus) != ns.get(locus)) ? 1 : 0;
        }

        return editDistance;
    }

    public String toString() {
        char[] charbuffer = new char[this.size()];
        for (int locus = 0; locus < this.size(); locus++) {
            charbuffer[locus] = this.get(locus).asChar();
        }

        return new String(charbuffer);
    }
}
