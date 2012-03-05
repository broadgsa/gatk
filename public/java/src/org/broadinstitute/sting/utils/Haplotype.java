/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

public class Haplotype {
    protected final byte[] bases;
    protected final double[] quals;
    private GenomeLoc genomeLocation = null;
    private boolean isReference = false;
 
    /**
     * Create a simple consensus sequence with provided bases and a uniform quality over all bases of qual
     *
     * @param bases bases
     * @param qual  qual
     */
    public Haplotype(byte[] bases, int qual) {
        this.bases = bases;
        quals = new double[bases.length];
        Arrays.fill(quals, (double)qual);
    }

    public Haplotype(byte[] bases, double[] quals) {
        this.bases = bases;
        this.quals = quals;
    }

    public Haplotype(byte[] bases) {
        this(bases, 0);
    }

    public Haplotype(byte[] bases, GenomeLoc loc) {
        this(bases);
        this.genomeLocation = loc;
    }

    public Haplotype(byte[] bases, GenomeLoc loc, boolean isRef) {
        this(bases, loc);
        this.isReference = isRef;
    }

    @Override
    public boolean equals( Object h ) {
        return h instanceof Haplotype && Arrays.equals(bases, ((Haplotype) h).bases);
    }

    public double getQualitySum() {
        double s = 0;
        for (int k=0; k < bases.length; k++) {
            s += quals[k];
        }
        return s;
    }

    public String toString() {
        String returnString = "";
        for(int iii = 0; iii < bases.length; iii++) {
            returnString += (char) bases[iii];
        }
        return returnString;
    }

    public double[] getQuals() {
        return quals;
    }
    public byte[] getBases() {
        return bases;
    }

    public long getStartPosition() {
        return genomeLocation.getStart();
    }

    public long getStopPosition() {
        return genomeLocation.getStop();
    }

    public boolean isReference() {
        return isReference;
    }

    public byte[] insertAllele( final Allele refAllele, final Allele altAllele, int refInsertLocation, final int hapStart, final Cigar haplotypeCigar ) {
        
        if( refAllele.length() != altAllele.length() ) { refInsertLocation++; }
        int haplotypeInsertLocation = getHaplotypeCoordinateForReferenceCoordinate(hapStart, haplotypeCigar, refInsertLocation);
        if( haplotypeInsertLocation == -1 ) { // desired change falls inside deletion so don't bother creating a new haplotype
            return bases.clone();
        }
        byte[] newHaplotype;

        try {
            if( refAllele.length() == altAllele.length() ) { // SNP or MNP
                newHaplotype = bases.clone();
                for( int iii = 0; iii < altAllele.length(); iii++ ) {
                    newHaplotype[haplotypeInsertLocation+iii] = altAllele.getBases()[iii];
                }
            } else if( refAllele.length() < altAllele.length() ) { // insertion                
                final int altAlleleLength = altAllele.length();
                newHaplotype = new byte[bases.length + altAlleleLength];
                for( int iii = 0; iii < bases.length; iii++ ) {
                    newHaplotype[iii] = bases[iii];
                }
                for( int iii = newHaplotype.length - 1; iii > haplotypeInsertLocation + altAlleleLength - 1; iii-- ) {
                    newHaplotype[iii] = newHaplotype[iii-altAlleleLength];
                }
                for( int iii = 0; iii < altAlleleLength; iii++ ) {
                    newHaplotype[haplotypeInsertLocation+iii] = altAllele.getBases()[iii];
                }
            } else { // deletion
                final int shift = refAllele.length() - altAllele.length();
                newHaplotype = new byte[bases.length - shift];
                for( int iii = 0; iii < haplotypeInsertLocation + altAllele.length(); iii++ ) {
                    newHaplotype[iii] = bases[iii];
                }
                for( int iii = haplotypeInsertLocation + altAllele.length(); iii < newHaplotype.length; iii++ ) {
                    newHaplotype[iii] = bases[iii+shift];
                }
            }
        } catch (Exception e) { // event already on haplotype is too large/complex to insert another allele, most likely because of not enough reference padding
            return bases.clone();
        }
        
        return newHaplotype;
    }

    public static LinkedHashMap<Allele,Haplotype> makeHaplotypeListFromAlleles(List<Allele> alleleList, int startPos, ReferenceContext ref,
                                                               final int haplotypeSize, final int numPrefBases) {


        LinkedHashMap<Allele,Haplotype> haplotypeMap = new LinkedHashMap<Allele,Haplotype>();

        Allele refAllele = null;

        for (Allele a:alleleList) {
            if (a.isReference()) {
                refAllele = a;
                break;
            }

        }

        if (refAllele == null)
            throw new ReviewedStingException("BUG: no ref alleles in input to makeHaplotypeListfrom Alleles at loc: "+ startPos);

        byte[] refBases = ref.getBases();


        int startIdxInReference = (int)(1+startPos-numPrefBases-ref.getWindow().getStart());
        //int numPrefBases = (int)(vc.getStart()-ref.getWindow().getStart()+1); // indel vc starts one before event


        byte[] basesBeforeVariant = Arrays.copyOfRange(refBases,startIdxInReference,startIdxInReference+numPrefBases);
        int startAfter = startIdxInReference+numPrefBases+ refAllele.getBases().length;
        // protect against long events that overrun available reference context
        if (startAfter > refBases.length)
            startAfter = refBases.length;
        byte[] basesAfterVariant = Arrays.copyOfRange(refBases,
                startAfter, refBases.length);


        // Create location for all haplotypes
        int startLoc = ref.getWindow().getStart() + startIdxInReference;
        int stopLoc = startLoc + haplotypeSize-1;

        GenomeLoc locus = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),startLoc,stopLoc);


        for (Allele a : alleleList) {

            byte[] alleleBases = a.getBases();
            // use string concatenation
            String haplotypeString = new String(basesBeforeVariant) + new String(alleleBases) + new String(basesAfterVariant);
            haplotypeString = haplotypeString.substring(0,haplotypeSize);

           haplotypeMap.put(a,new Haplotype(haplotypeString.getBytes(), locus, a.isReference()));

        }

        return haplotypeMap;
    }

    private static Integer getHaplotypeCoordinateForReferenceCoordinate( final int haplotypeStart, final Cigar haplotypeCigar, final int refCoord ) {
        int readBases = 0;
        int refBases = 0;
        boolean fallsInsideDeletion = false;

        int goal = refCoord - haplotypeStart;  // The goal is to move this many reference bases
        boolean goalReached = refBases == goal;

        Iterator<CigarElement> cigarElementIterator = haplotypeCigar.getCigarElements().iterator();
        while (!goalReached && cigarElementIterator.hasNext()) {
            CigarElement cigarElement = cigarElementIterator.next();
            int shift = 0;

            if (cigarElement.getOperator().consumesReferenceBases() || cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
                if (refBases + cigarElement.getLength() < goal)
                    shift = cigarElement.getLength();
                else
                    shift = goal - refBases;

                refBases += shift;
            }
            goalReached = refBases == goal;

            if (!goalReached && cigarElement.getOperator().consumesReadBases())
                readBases += cigarElement.getLength();

            if (goalReached) {
                // Is this base's reference position within this cigar element? Or did we use it all?
                boolean endsWithinCigar = shift < cigarElement.getLength();

                // If it isn't, we need to check the next one. There should *ALWAYS* be a next one
                // since we checked if the goal coordinate is within the read length, so this is just a sanity check.
                if (!endsWithinCigar && !cigarElementIterator.hasNext())
                    return -1;

                CigarElement nextCigarElement;

                // if we end inside the current cigar element, we just have to check if it is a deletion
                if (endsWithinCigar)
                    fallsInsideDeletion = cigarElement.getOperator() == CigarOperator.DELETION;

                    // if we end outside the current cigar element, we need to check if the next element is an insertion or deletion.
                else {
                    nextCigarElement = cigarElementIterator.next();

                    // if it's an insertion, we need to clip the whole insertion before looking at the next element
                    if (nextCigarElement.getOperator() == CigarOperator.INSERTION) {
                        readBases += nextCigarElement.getLength();
                        if (!cigarElementIterator.hasNext())
                            return -1;

                        nextCigarElement = cigarElementIterator.next();
                    }

                    // if it's a deletion, we will pass the information on to be handled downstream.
                    fallsInsideDeletion = nextCigarElement.getOperator() == CigarOperator.DELETION;
                }

                // If we reached our goal outside a deletion, add the shift
                if (!fallsInsideDeletion && cigarElement.getOperator().consumesReadBases())
                    readBases += shift;

                    // If we reached our goal inside a deletion, but the deletion is the next cigar element then we need
                    // to add the shift of the current cigar element but go back to it's last element to return the last
                    // base before the deletion (see warning in function contracts)
                else if (fallsInsideDeletion && !endsWithinCigar)
                    readBases += shift - 1;

                    // If we reached our goal inside a deletion then we must backtrack to the last base before the deletion
                else if (fallsInsideDeletion && endsWithinCigar)
                    readBases--;
            }
        }

        if (!goalReached)
            return -1;

        return (fallsInsideDeletion ? -1 : readBases);
    }

}
