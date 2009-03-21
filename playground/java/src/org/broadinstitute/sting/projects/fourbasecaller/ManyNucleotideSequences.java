package org.broadinstitute.sting.projects.fourbasecaller;

import java.io.File;
import java.util.Vector;

import org.broadinstitute.sting.projects.fourbasecaller.NucleotideSequence;

import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;

class ManyNucleotideSequences extends java.util.Vector<NucleotideSequence> {
    public ManyNucleotideSequences() {}

    public ManyNucleotideSequences(int initialCapacity) {
        this.setSize(initialCapacity);
    }

    public ManyNucleotideSequences(String seqpath) {
        ReferenceSequenceFile seqfile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(seqpath));

        ReferenceSequence rs;
        while ((rs = seqfile.nextSequence()) != null) {
            byte[] bases = rs.getBases();

            NucleotideSequence ns = new NucleotideSequence(bases.length);
            for (int locus = 0; locus < bases.length; locus++) {
                Nucleotide nt;
                switch(bases[locus]) {
                    case 'A': nt = Nucleotide.A; break;
                    case 'C': nt = Nucleotide.C; break;
                    case 'G': nt = Nucleotide.G; break;
                    case 'T': nt = Nucleotide.T; break;
                    default: nt = Nucleotide.T; break;
                }

                ns.set(locus, nt);
            }

            this.add(ns);
        }
    }

    public ManyNucleotideSequences reverseComplement() {
        for (int contig = 0; contig < this.size(); contig++) {
            this.get(contig).reverseComplement();
        }

        return this;
    }
}
