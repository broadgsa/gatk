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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public class Haplotype {
    protected final byte[] bases;
    protected final double[] quals;
    private GenomeLoc genomeLocation = null;
    private HashMap<String, double[]> readLikelihoodsPerSample = null;
    private HashMap<Integer, VariantContext> eventMap = null;
    private boolean isRef = false;
    private Cigar cigar;
    private int alignmentStartHapwrtRef;
 
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

    @Override
    public boolean equals( Object h ) {
        return h instanceof Haplotype && Arrays.equals(bases, ((Haplotype) h).bases);
    }
    
    @Override
    public int hashCode() {
        return Arrays.hashCode(bases);
    }

    public void addReadLikelihoods( final String sample, final double[] readLikelihoods ) {
        if( readLikelihoodsPerSample == null ) {
            readLikelihoodsPerSample = new HashMap<String, double[]>();
        }
        readLikelihoodsPerSample.put(sample, readLikelihoods);
    }

    @Ensures({"result != null"})
    public double[] getReadLikelihoods( final String sample ) {
        return readLikelihoodsPerSample.get(sample);
    }
    
    public Set<String> getSampleKeySet() {
        return readLikelihoodsPerSample.keySet();
    }

    public HashMap<Integer, VariantContext> getEventMap() {
        return eventMap;
    }

    public void setEventMap( final HashMap<Integer, VariantContext> eventMap ) {
        this.eventMap = eventMap;
    }

    public boolean isReference() {
        return isRef;
    }

    public void setIsReference( boolean isRef ) {
        this.isRef = isRef;
    }

    public double getQualitySum() {
        double s = 0;
        for (int k=0; k < bases.length; k++) {
            s += quals[k];
        }
        return s;
    }

    @Override
    public String toString() {
        return new String(bases);
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

    public int getAlignmentStartHapwrtRef() {
        return alignmentStartHapwrtRef;
    }

    public void setAlignmentStartHapwrtRef( final int alignmentStartHapwrtRef ) {
        this.alignmentStartHapwrtRef = alignmentStartHapwrtRef;
    }

    public Cigar getCigar() {
        return cigar;
    }

    public void setCigar( final Cigar cigar ) {
        this.cigar = cigar;
    }

    @Requires({"refInsertLocation >= 0"})
    public Haplotype insertAllele( final Allele refAllele, final Allele altAllele, int refInsertLocation ) {
        
        if( refAllele.length() != altAllele.length() ) { refInsertLocation++; }
        int haplotypeInsertLocation = ReadUtils.getReadCoordinateForReferenceCoordinate(alignmentStartHapwrtRef, cigar, refInsertLocation, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if( haplotypeInsertLocation == -1 ) { // desired change falls inside deletion so don't bother creating a new haplotype
            return new Haplotype(bases.clone());
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
            return new Haplotype(bases.clone());
        }
        
        return new Haplotype(newHaplotype);
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
        final int startLoc = ref.getWindow().getStart() + startIdxInReference;
        final int stopLoc = startLoc + haplotypeSize-1;

        final GenomeLoc locus = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),startLoc,stopLoc);


        for (final Allele a : alleleList) {

            byte[] alleleBases = a.getBases();
            // use string concatenation
            String haplotypeString = new String(basesBeforeVariant) + new String(alleleBases) + new String(basesAfterVariant);
            haplotypeString = haplotypeString.substring(0,haplotypeSize);

           haplotypeMap.put(a,new Haplotype(haplotypeString.getBytes(), locus));

        }

        return haplotypeMap;
    }
}
