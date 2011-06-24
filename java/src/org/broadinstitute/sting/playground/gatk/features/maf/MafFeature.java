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
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
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
        Frameshift, Inframe, Stop_deletion, Promoter,De_novo_start, De_novo_start_out_of_frame
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

    public void setHugoGeneSymbol(String genename) {
        int pos = genename.indexOf('|');
        if ( pos < 0 ) {
            hugoSymbol = genename;
        } else {
            hugoSymbol = genename.substring(0,pos);
        }
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
        if ( s.equals("3'UTR") || s.equals("3'-UTR")) { classification = Classification.UTR3 ; return; }
        if ( s.equals("5'UTR") || s.equals("5'-UTR")) { classification = Classification.UTR5 ; return; }
        if ( s.equals("5'-Flank") ) { classification = Classification.Flank5 ; return; }
        if ( s.equals("Silent") || s.equals("Synonymous")) { classification = Classification.Silent ; return; }
        if ( s.equals("Non-coding_Transcript")) { classification = Classification.Noncoding_transcript; return; }
        if ( s.equals("Missense") || s.equals("Missense_Mutation") ) { classification = Classification.Missense ; return; }
        if ( s.equals("Nonsense_Mutation") || s.equals("Nonsense") ) { classification = Classification.Nonsense ; return; }
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
        if ( s.equals("De_novo_Start_OutOfFrame") ) { classification = Classification.De_novo_start_out_of_frame ; return; }
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

class MafAdaptor implements VariantContextAdaptors.VCAdaptor {
    /**
     * Converts Maf features to VariantContext.
     * @return MafFeature.
     */
    @Override
    public Class<? extends Feature> getAdaptableFeatureType() { return MafFeature.class; }

    /**
     * convert to a Variant Context, given:
     * @param name the name of the ROD
     * @param input the Rod object, in this case a MafFeature
     * @return a VariantContext object
     */
//        VariantContext convert(String name, Object input) {
//            return convert(name, input, null);
//        }

    /**
     * convert to a Variant Context, given:
     * @param name  the name of the ROD
     * @param input the Rod object, in this case a MafFeature
     * @param ref   the reference context
     * @return a VariantContext object
     */
    @Override
    public VariantContext convert(String name, Object input, ReferenceContext ref) {

        if ( ref == null )
            throw new UnsupportedOperationException("Conversion from MAF to VariantContext requires a reference context, null received");

        MafFeature maf = (MafFeature)input;
        if ( ! Allele.acceptableAlleleBases(maf.getRefBases()) )
            return null;

        List<Allele> alleles = new ArrayList<Allele>();

        Allele refAllele = Allele.create(maf.getRefBases(), true);
        // add the reference allele:
        alleles.add(refAllele);

        // add all of the alt alleles
        for ( String alt : maf.getAllNonRefAlleleList() ) {
            if ( ! Allele.acceptableAlleleBases(alt) ) {
                //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                return null;
            }
            alleles.add(Allele.create(alt, false));
        }

        // make a mapping from sample to genotype

        String normalSample = maf.getNormalSampleId();
        String tumorSample = maf.getTumorSampleId();

//                String[] genotypeStrings = hapmap.getGenotypes();

        Map<String, Genotype> genotypes = new HashMap<String, Genotype>(2);

        addGenotype(genotypes, normalSample, maf.getObservedNormalAlleleList(),maf.getRefBases());
        addGenotype(genotypes,tumorSample,maf.getObservedTumorAlleleList(),maf.getRefBases());


        HashMap<String, Object> attrs = new HashMap<String, Object>(10);
        // fill attributes:
        if ( maf.getHugoGeneSymbol() != null && ! maf.getHugoGeneSymbol().equals("Unknown"))
            attrs.put("Gene",maf.getHugoGeneSymbol());

        if ( maf.isSomatic() ) {
            attrs.put(VCFConstants.SOMATIC_KEY,true);
            attrs.put("SS","Somatic");
        } else {
            attrs.put("SS","Germline");
        }

        if ( maf.getVariantClassification() != null ) {
            switch(maf.getVariantClassification()) {
                case Intergenic: attrs.put("VC","Genomic"); break;
                case Intron: attrs.put("VC","Intron"); break;
                case Noncoding_transcript: attrs.put("VC","Noncoding_transcript"); break;
                case UTR3: attrs.put("VC","3'UTR"); break;
                case UTR5: attrs.put("VC","5'UTR"); break;
                case Flank5: attrs.put("VC","5'flank"); break;
                case Promoter: attrs.put("VC","5'flank"); break;
                case De_novo_start: attrs.put("VC","De_novo_start"); break;
                case De_novo_start_out_of_frame: attrs.put("VC","De_novo_start_out_of_frame"); break;
                case Silent: attrs.put("VC","Silent"); break;
                case Missense: attrs.put("VC","Missense"); break;
                case Nonsense: attrs.put("VC","Nonsense"); break;
                case Splice_site: attrs.put("VC","Splice_site"); break;
                case miRNA: attrs.put("VC","miRNA"); break;
                case Frameshift: attrs.put("VC","Frameshift"); break;
                case Inframe: attrs.put("VC","Inframe"); break;
                case Stop_deletion: attrs.put("VC","Stop_codon_deletion");
                case Unclassified: attrs.put("VC","Unclassified");
                default:
            }
        }

        attrs.put("VT",maf.getType());

//                attrs.put(VariantContext.ID_KEY, hapmap.getName());
        int end = maf.getEnd();
        VariantContext vc = new VariantContext(name, maf.getChr(), maf.getStart(), end, alleles,
                genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, attrs);
        return vc;
    }

    private void addGenotype(Map<String,Genotype> dest, String sampleId, List<String> alleles, String refAllele) {
        List<Allele> myAlleles = new ArrayList<Allele>(2);

        boolean success = true;

        for ( String a : alleles ) {
            if ( a.isEmpty() || a.contains("N") || a.contains(".")) return; // bad allele found
            myAlleles.add(Allele.create(a,refAllele.equals(a)));
        }
        dest.put(sampleId, new Genotype(sampleId,myAlleles));
    }

}