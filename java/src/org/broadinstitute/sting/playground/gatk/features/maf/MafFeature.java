/*
 * Copyright (c) 2010 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.features.maf;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 24, 2011
 * Time: 12:01:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class MafFeature implements Feature {
    private String contig;                      // our contig location
    private int start;                          // our starting location, zero based
    private int stop;                           // our stopping location

    private String refAllele = ".";             // the reference allele
    private String[] observedTumAlleles = null;          // The sequences of the observed alleles in tumor
    private String[] observedNormAlleles = null;          // The sequences of the observed alleles in normal
    private String tumorSampleId = null;
    private String normalSampleId = null;
    private String hugoSymbol = null;
    private Classification classification = null;

    public enum Type {
        UNKNOWN,SNP,MNP,INS,DEL
    };

    public enum Classification {
        Unclassified, Intergenic,Intron,Noncoding_transcript,UTR3,UTR5,Flank5,Silent,Missense, Nonsense, Splice_site, miRNA,
        Frameshift, Inframe, Stop_deletion, Promoter,De_novo_start
    }

    private Type type = Type.UNKNOWN;

    /**
     * create the dbSNP feature, given the following information:
     *
     * @param contig the contig rsID
     * @param start  the start position, one based
     * @param stop   the stop position, one based
     */
    MafFeature(String contig,
                 int start,
                 int stop) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
    }

    public void setVariantType(String t) {
        type=Type.valueOf(t);
    }

    public void setObservedTumor(String[] obs) {
        observedTumAlleles = obs;
    }

    public void setObservedTumor(String allele1, String allele2) {
        observedTumAlleles = new String[2];
        observedTumAlleles[0] = allele1;
        observedTumAlleles[1] = allele2;
    }

    public void setRefAllele(String ref) {
        this.refAllele = ref;
    }

    public void setTumorSample(String sampleId) {
        this.tumorSampleId = sampleId;
    }

    public void setNormalSample(String sampleId) {
        this.normalSampleId = sampleId;
    }

    public String getRefBases() {
        return refAllele;
    }

    public String getHugoGeneSymbol() {
        return hugoSymbol;
    }

    public String setHugoGeneSymbol(String genename) {
        return hugoSymbol = genename;
    }

    /**
     * Returns list of alleles (represented as strings) observed in Tumor. Returned alleles
     * could be redundant (e.g. if we have homozygous non-ref at ploidy 2+).
     * @return
     */
    public List<String> getObservedTumorAlleleList() {
        return Arrays.asList(observedTumAlleles);
    }

    /**
     * Returns list of alleles (represented as strings) observed in Tumor. Returned alleles
     * could be redundant (e.g. if we have homozygous non-ref at ploidy 2+).
     * @return
     */
    public List<String> getObservedNormalAlleleList() {
        if ( observedNormAlleles == null ) {
            // if we got no ref allele observations recorded in the maf, we assume its ref/ref (somatic event)
            List<String> l = new ArrayList<String>(2);
            l.add(refAllele);
            l.add(refAllele);
            return l;
        }
        else return Arrays.asList(observedTumAlleles);
    }

    /** Returns a (non-redundant) list of all distinct alleles
     * observed at the site, plus a reference allele (whether it
     * was actually observed or not). The reference allele is always returned as first
     * element of the list.
     * @return
     */
    public List<String> getAllAlleleList() {
        List<String> l = new ArrayList<String>();
        l.add(refAllele);
        for ( String a : observedTumAlleles ) {
            if ( l.contains(a) ) continue;
            l.add(a);
        }
        if ( observedNormAlleles != null ) {
            for ( String a : observedNormAlleles ) {
                if ( l.contains(a) ) continue;      // already have this allele
                l.add(a);
            }
        }
        return l;
    }

    /** Returns a (non-redundant) list of all distinct non-reference alleles
     * observed at the site
     * @return
     */
    public List<String> getAllNonRefAlleleList() {
        List<String> l = new ArrayList<String>();
        for ( String a : observedTumAlleles ) {
            if ( l.contains(a) ) continue;      // already have this allele
            if ( a.equals(refAllele)) continue; // allele is ref, we do not need it
            l.add(a); 
        }
        if ( observedNormAlleles != null ) {
            for ( String a : observedNormAlleles ) {
                if ( l.contains(a) ) continue;      // already have this allele
                if ( a.equals(refAllele)) continue; // allele is ref, we do not need it
                l.add(a);
            }
        }
        return l;
    }

    public String getTumorSampleId() { return tumorSampleId; }
    public String getNormalSampleId() { return normalSampleId; }

    public boolean isRefAllele(String a) { return refAllele.equals(a); }

    public Type getType() { return type; }

    public int lengthOnRef() {
        switch ( type ) {
            case SNP:
            case MNP:
            case DEL:
                return refAllele.length();
            case INS:
                return 0;
            default:
                throw new StingException("Unrecognized event type in Maf record: "+type);
        }
    }

    public boolean isSomatic() {
        if ( observedTumAlleles[0].equals(refAllele) && observedTumAlleles[1].equals(refAllele) ) return false; // tumor is ref
        // we get here only if tumor is non-ref
        if ( observedNormAlleles == null ) return true; // norm alleles are omitted from maf only if they are all ref
        if ( observedNormAlleles[0].equals(refAllele) && observedNormAlleles[1].equals(refAllele) ) return true;
        return false;
    }

    public void setVariantClassification(String s) {
        if ( s.equals("IGR") ) { classification = Classification.Intergenic ; return; }
        if ( s.equals("Intron") ) { classification = Classification.Intron ; return; }
        if ( s.equals("3'UTR") ) { classification = Classification.UTR3 ; return; }
        if ( s.equals("5'UTR") ) { classification = Classification.UTR5 ; return; }
        if ( s.equals("5'-Flank") ) { classification = Classification.Flank5 ; return; }
        if ( s.equals("Silent") ) { classification = Classification.Silent ; return; }
        if ( s.equals("Non-coding_Transcript")) { classification = Classification.Noncoding_transcript; return; }
        if ( s.equals("Missense") || s.equals("Missense_Mutation") ) { classification = Classification.Missense ; return; }
        if ( s.equals("Nonsense_Mutation") ) { classification = Classification.Nonsense ; return; }
        if ( s.equals("Splice_Site") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("miRNA") ) { classification = Classification.miRNA ; return; }
        if ( s.equals("Frame_Shift_Ins") ) { classification = Classification.Frameshift ; return; }
        if ( s.equals("Frame_Shift_Del") ) { classification = Classification.Frameshift ; return; }
        if ( s.equals("In_Frame_Ins") ) { classification = Classification.Inframe ; return; }
        if ( s.equals("In_Frame_Del") ) { classification = Classification.Inframe ; return; }
        if ( s.equals("Stop_Codon_Del") ) { classification = Classification.Stop_deletion ; return; }
        if ( s.equals("Splice_Site_Del") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("Splice_Site_Ins") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("Splice_Site_SNP") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("Promoter") ) { classification = Classification.Promoter ; return; }
        if ( s.equals("De_novo_Start") ) { classification = Classification.De_novo_start ; return; }
        if ( s.equals("TX-REF-MISMATCH") ) { classification = Classification.Unclassified ; return; }
        throw new UserException.MalformedFile("Unknown variant classification: " + s);
    }

    public Classification getVariantClassification() {
        return classification;
    }

   /*
     * the required getting and setter methods
     */

    public String getChr() {
        return contig;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return stop;
    }
    
}
