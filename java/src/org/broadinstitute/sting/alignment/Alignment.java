package org.broadinstitute.sting.alignment;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.CigarElement;

/**
 * Represents an alignment of a read to a site in the reference genome.
 *
 * @author mhanna
 * @version 0.1
 */
public class Alignment {
    protected int contigIndex;
    protected long alignmentStart;
    protected boolean negativeStrand;
    protected int mappingQuality;

    private char[] cigarOperators;
    private int[] cigarLengths;

    /**
     * Gets the index of the given contig.
     * @return the inde
     */
    public int getContigIndex() { return contigIndex; }

    /**
     * Gets the starting position for the given alignment.
     * @return Starting position.
     */
    public long getAlignmentStart() { return alignmentStart; }

    /**
     * Is the given alignment on the reverse strand?
     * @return True if the alignment is on the reverse strand.
     */
    public boolean isNegativeStrand() { return negativeStrand; }

    /**
     * Gets the score of this alignment.
     * @return The score.
     */
    public int getMappingQuality() { return mappingQuality; }

    /**
     * Gets the cigar for this alignment.
     * @return sam-jdk formatted alignment.
     */
    public Cigar getCigar() {
        Cigar cigar = new Cigar();
        for(int i = 0; i < cigarOperators.length; i++) {
            CigarOperator operator = CigarOperator.characterToEnum(cigarOperators[i]);
            cigar.add(new CigarElement(cigarLengths[i],operator));
        }
        return cigar;
    }

    /**
     * Temporarily implement getCigarString() for debugging; the TextCigarCodec is unfortunately
     * package-protected.
     * @return
     */
    public String getCigarString() {
        Cigar cigar = getCigar();
        if(cigar.isEmpty()) return "*";

        StringBuilder cigarString = new StringBuilder();
        for(CigarElement element: cigar.getCigarElements()) {
            cigarString.append(element.getLength());
            cigarString.append(element.getOperator());
        }
        return cigarString.toString();
    }

    /**
     * Stub for inheritance.
     */
    public Alignment() {}    

    /**
     * Create a new alignment object.
     * @param contigIndex The contig to which this read aligned.
     * @param alignmentStart The point within the contig to which this read aligned.
     * @param negativeStrand Forward or reverse alignment of the given read.
     * @param mappingQuality How good does BWA think this mapping is?
     * @param cigarOperators The ordered operators in the cigar string.
     * @param cigarLengths The lengths to which each operator applies.
     */
    public Alignment(int contigIndex, int alignmentStart, boolean negativeStrand, int mappingQuality,
                     char[] cigarOperators, int[] cigarLengths) {
        this.contigIndex = contigIndex;
        this.alignmentStart = alignmentStart;
        this.negativeStrand = negativeStrand;
        this.mappingQuality = mappingQuality;
        this.cigarOperators = cigarOperators;
        this.cigarLengths = cigarLengths;
    }
}
