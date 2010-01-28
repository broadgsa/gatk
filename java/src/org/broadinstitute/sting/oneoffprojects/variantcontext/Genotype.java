package org.broadinstitute.sting.oneoffprojects.variantcontext;

import java.util.*;

/**
 * @author ebanks
 *         <p/>
 *         Class Genotype
 *         <p/>
 *         This class emcompasses all the basic information about a genotype
 */
public class Genotype extends AttributedObject {
    private List<Allele> alleles = new ArrayList<Allele>();
    private String sampleName = null;

    // todo -- do genotypes need to have locations?  Or is it sufficient to have an
    // todo -- associated VC with a location?  One nasty implication is that people will have to
    // todo -- pass around both a Variant Context and genotypes.  Although users can always just package up
    // the associated genotypes into the VC itself.

    public Genotype(VariantContext vc, List<String> alleles, String sampleName, double negLog10PError) {
        this(resolveAlleles(vc, alleles), sampleName, negLog10PError);
    }

    public Genotype(List<Allele> alleles, String sampleName, double negLog10PError) {
        setAlleles(alleles);
        setSampleName(sampleName);
        setNegLog10PError(negLog10PError);
    }

    /**
     * @return the alleles for this genotype
     */
    public List<Allele> getAlleles() {
        return alleles;
    }

    public List<Allele> getAlleles(Allele allele) {
        List<Allele> al = new ArrayList<Allele>();
        for ( Allele a : alleles )
            if ( a.equals(allele) )
                al.add(a);

        return al;
    }


    /**
     * @return the ploidy of this genotype
     */
    public int getPloidy() { return alleles.size(); }

    public enum Type {
        NO_CALL,
        HOM_REF,
        HET,
        HOM_VAR
    }

    public Type getType() {
        Allele firstAllele = alleles.get(0);

        if ( firstAllele.isNoCall() ) {
            return Type.NO_CALL;
        }

        for (Allele a : alleles) {
            if ( ! firstAllele.equals(a) )
                return Type.HET;
        }
        return firstAllele.isReference() ? Type.HOM_REF : Type.HOM_VAR;
    }

    /**
     * @return true if all observed alleles are the same (regardless of whether they are ref or alt)
     */
    public boolean isHom() { return getType() == Type.HOM_REF || getType() == Type.HOM_VAR; }

    /**
     * @return true if we're het (observed alleles differ)
     */
    public boolean isHet() { return getType() == Type.HET; }

    /**
     * @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF)
     */
    public boolean isNoCall() { return getType() == Type.NO_CALL; }

//    /**
//     * @return true if this is a variant genotype, false if it's reference
//     */
//    public boolean isVariant() {
//        for ( Allele allele : alleles ) {
//            if ( allele.isNonReference() )
//                return true;
//        }
//        return false;
//    }

    /**
     * @return the sample name
     */
    public String getSampleName() {
        return sampleName;
    }

    /**
     * Sets the sample name
     *
     * @param sampleName    the sample name
     */
    public void setSampleName(String sampleName) {
        if ( sampleName == null ) throw new IllegalArgumentException("Sample name cannot be null " + this);
        this.sampleName = sampleName;
    }

    /**
     *
     * @param alleles
     */
    public void setAlleles(List<Allele> alleles) {
        this.alleles.clear();

        // todo -- add validation checking here

        if ( alleles == null ) throw new IllegalArgumentException("BUG: alleles cannot be null in setAlleles");
        if ( alleles.size() == 0) throw new IllegalArgumentException("BUG: alleles cannot be of size 0 in setAlleles");

        int nNoCalls = 0;
        for ( Allele allele : alleles ) { nNoCalls += allele.isNoCall() ? 1 : 0; }
        if ( nNoCalls > 0 && nNoCalls != alleles.size() )
            throw new IllegalArgumentException("BUG: alleles include some No Calls and some Calls, an illegal state " + this);

        for ( Allele allele : alleles ) {
            addAllele(allele);
        }
    }

    public void addAllele(Allele allele) {
        if ( allele == null ) throw new IllegalArgumentException("BUG: Cannot add a null allele to a genotype");
        this.alleles.add(allele);
    }


    private static List<Allele> resolveAlleles(VariantContext vc, List<String> alleleStrings) {
        List<Allele> alleles = new ArrayList<Allele>();
        for ( String alleleString : alleleStrings ) {
            Allele allele = vc.getAllele(alleleString);

            if ( allele == null ) {
                if ( Allele.wouldBeNoCallAllele(alleleString.getBytes()) ) {
                    allele = new Allele(alleleString);
                } else {
                    throw new IllegalArgumentException("Allele " + alleleString + " not present in the list of alleles in VariantContext " + vc);
                }
            }

            alleles.add(allele);
        }

        return alleles;
    }

    public String toString() {
        return String.format("[GT: %s %s %s Q%.2f %s]", getSampleName(), getAlleles(), getType(), getNegLog10PError(), getAttributes());
    }
}