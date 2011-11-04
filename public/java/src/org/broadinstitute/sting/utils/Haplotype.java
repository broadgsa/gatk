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

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

public class Haplotype {
    protected byte[] bases = null;
    protected double[] quals = null;
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
    public byte[] getBasesAsBytes() {
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

}
