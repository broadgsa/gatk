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

package org.broadinstitute.sting.utils.genotype;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

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

    public Haplotype(String bases, double[] quals) {
        this.bases = bases.getBytes();
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


    public String toString() { return new String(this.bases); }

    public double getQualitySum() {
        double s = 0;
        for (int k=0; k < bases.length; k++) {
            s += quals[k];
        }
        return s;
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

    public boolean isReference() {
        return isReference;
    }


    public static List<Haplotype> makeHaplotypeListFromVariantContextAlleles(VariantContext vc, ReferenceContext ref, final int haplotypeSize) {


        List<Haplotype> haplotypeList = new ArrayList<Haplotype>();

        byte[] refBases = ref.getBases();

        int numPrefBases = haplotypeSize/2;

        int startIdxInReference = (int)(1+vc.getStart()-numPrefBases-ref.getWindow().getStart());
        //int numPrefBases = (int)(vc.getStart()-ref.getWindow().getStart()+1); // indel vc starts one before event

 
        byte[] basesBeforeVariant = Arrays.copyOfRange(refBases,startIdxInReference,startIdxInReference+numPrefBases);
        byte[] basesAfterVariant = Arrays.copyOfRange(refBases,
                startIdxInReference+numPrefBases+vc.getReference().getBases().length, refBases.length);


        // Create location for all haplotypes
        long startLoc = ref.getWindow().getStart() + startIdxInReference;
        long stopLoc = startLoc + haplotypeSize-1;

        GenomeLoc locus = GenomeLocParser.createGenomeLoc(ref.getLocus().getContigIndex(),startLoc,
                 stopLoc);


        for (Allele a : vc.getAlleles()) {

            byte[] alleleBases = a.getBases();
            // use string concatenation
            String haplotypeString = new String(basesBeforeVariant) + new String(alleleBases) + new String(basesAfterVariant);
            haplotypeString = haplotypeString.substring(0,haplotypeSize);

           haplotypeList.add(new Haplotype(haplotypeString.getBytes(), locus, a.isReference()));

        }

        return haplotypeList;
    }

}
