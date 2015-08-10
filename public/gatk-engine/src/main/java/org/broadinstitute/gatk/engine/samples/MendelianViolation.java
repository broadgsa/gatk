/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.samples;

import org.broadinstitute.gatk.engine.samples.Sample;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.MathUtils;

import java.util.*;

/**
 * User: carneiro / lfran
 * Date: 3/9/11
 * Time: 12:38 PM
 *
 * Class for the identification and tracking of mendelian violation. It can be used in 2 distinct ways:
 * - Either using an instance of the MendelianViolation class to track mendelian violations for each of the families while
 * walking over the variants
 * - Or using the static methods to directly get information about mendelian violation in a family at a given locus
 *
 */
public class MendelianViolation {
    //List of families with violations
    private List<String> violationFamilies;

    //Call information
    private int nocall = 0;
    private int familyCalled = 0;
    private int varFamilyCalled = 0;
    private int lowQual = 0;

    private boolean allCalledOnly = true;

    //Stores occurrences of inheritance
    private EnumMap<GenotypeType, EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> inheritance;

    private int violations_total=0;

    private double minGenotypeQuality;

    private boolean abortOnSampleNotFound;

    //Number of families with genotype information for all members
    public int getFamilyCalledCount(){
        return familyCalled;
    }

    //Number of families with genotype information for all members
    public int getVarFamilyCalledCount(){
        return varFamilyCalled;
    }

    //Number of families missing genotypes for one or more of their members
    public int getFamilyNoCallCount(){
        return nocall;
    }

    //Number of families with genotypes below the set quality threshold
    public int getFamilyLowQualsCount(){
        return lowQual;
    }

    public int getViolationsCount(){
        return violations_total;
    }

    //Count of alt alleles inherited from het parents (no violation)
    public int getParentHetInheritedVar(){
        return getParentsHetHetInheritedVar() + getParentsRefHetInheritedVar() + getParentsVarHetInheritedVar();
    }

    //Count of ref alleles inherited from het parents (no violation)
    public int getParentHetInheritedRef(){
        return getParentsHetHetInheritedRef() + getParentsRefHetInheritedRef() + getParentsVarHetInheritedRef();
    }

    //Count of HomRef/HomRef/HomRef trios
    public int getRefRefRef(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF);
    }

    //Count of HomVar/HomVar/HomVar trios
    public int getVarVarVar(){
        return inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR);
    }

    //Count of HomRef/HomVar/Het trios
    public int getRefVarHet(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR).get(GenotypeType.HET) +
                inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
    }

    //Count of Het/Het/Het trios
    public int getHetHetHet(){
        return inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HET);
    }

    //Count of Het/Het/HomRef trios
    public int getHetHetHomRef(){
        return inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_REF);
    }

    //Count of Het/Het/HomVar trios
    public int getHetHetHomVar(){
        return inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_VAR);
    }

    //Count of ref alleles inherited from Het/Het parents (no violation)
    public int getParentsHetHetInheritedRef(){
        return inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HET)
               + 2*inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_REF);
        //return parentsHetHet_childRef;
    }

    //Count of var alleles inherited from Het/Het parents (no violation)
    public int getParentsHetHetInheritedVar(){
        return inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HET)
               + 2*inheritance.get(GenotypeType.HET).get(GenotypeType.HET).get(GenotypeType.HOM_VAR);
        //return parentsHetHet_childVar;
    }

    //Count of ref alleles inherited from HomRef/Het parents (no violation)
    public int getParentsRefHetInheritedRef(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HET).get(GenotypeType.HOM_REF)
               + inheritance.get(GenotypeType.HET).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF);
        //return parentsHomRefHet_childRef;
    }

    //Count of var alleles inherited from HomRef/Het parents (no violation)
    public int getParentsRefHetInheritedVar(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HET).get(GenotypeType.HET)
               + inheritance.get(GenotypeType.HET).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
        //return parentsHomRefHet_childVar;
    }

    //Count of ref alleles inherited from HomVar/Het parents (no violation)
    public int getParentsVarHetInheritedRef(){
        return inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HET).get(GenotypeType.HET)
               + inheritance.get(GenotypeType.HET).get(GenotypeType.HOM_VAR).get(GenotypeType.HET);
        //return parentsHomVarHet_childRef;
    }

    //Count of var alleles inherited from HomVar/Het parents (no violation)
    public int getParentsVarHetInheritedVar(){
        return inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HET).get(GenotypeType.HOM_VAR)
               + inheritance.get(GenotypeType.HET).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR);
        //return parentsHomVarHet_childVar;
    }

    //Count of violations of the type HOM_REF/HOM_REF -> HOM_VAR
    public int getParentsRefRefChildVar(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR);
    }

    //Count of violations of the type HOM_REF/HOM_REF -> HET
    public int getParentsRefRefChildHet(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF).get(GenotypeType.HET);
    }

    //Count of violations of the type HOM_REF/HET -> HOM_VAR
    public int getParentsRefHetChildVar(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HET).get(GenotypeType.HOM_VAR)
                + inheritance.get(GenotypeType.HET).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR);
    }

    //Count of violations of the type HOM_REF/HOM_VAR -> HOM_VAR
    public int getParentsRefVarChildVar(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR)
                + inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR);
    }

    //Count of violations of the type HOM_REF/HOM_VAR -> HOM_REF
    public int getParentsRefVarChildRef(){
        return inheritance.get(GenotypeType.HOM_REF).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF)
                + inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF).get(GenotypeType.HOM_REF);
    }

    //Count of violations of the type HOM_VAR/HET -> HOM_REF
    public int getParentsVarHetChildRef(){
        return inheritance.get(GenotypeType.HET).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF)
                + inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HET).get(GenotypeType.HOM_REF);
    }

    //Count of violations of the type HOM_VAR/HOM_VAR -> HOM_REF
    public int getParentsVarVarChildRef(){
        return inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_REF);
    }

    //Count of violations of the type HOM_VAR/HOM_VAR -> HET
    public int getParentsVarVarChildHet(){
        return inheritance.get(GenotypeType.HOM_VAR).get(GenotypeType.HOM_VAR).get(GenotypeType.HET);
    }


    //Count of violations of the type HOM_VAR/? -> HOM_REF
    public int getParentVarChildRef(){
        return getParentsRefVarChildRef() + getParentsVarHetChildRef() +getParentsVarVarChildRef();
    }

    //Count of violations of the type HOM_REF/? -> HOM_VAR
    public int getParentRefChildVar(){
        return getParentsRefVarChildVar() + getParentsRefHetChildVar() +getParentsRefRefChildVar();
    }

    //Returns a String containing all trios where a Mendelian violation was observed.
    //The String is formatted "mom1+dad1=child1,mom2+dad2=child2,..."
    public String getViolationFamiliesString(){
        if(violationFamilies.isEmpty())
            return "";

        Iterator<String> it = violationFamilies.iterator();
        String violationFams = it.next();
        while(it.hasNext()){
            violationFams += ","+it.next();
        }
        return violationFams;
    }

    public List<String> getViolationFamilies(){
        return violationFamilies;
    }

    static final int[] mvOffsets = new int[] { 1,2,5,6,8,11,15,18,20,21,24,25 };
    static final int[] nonMVOffsets = new int[]{ 0,3,4,7,9,10,12,13,14,16,17,19,22,23,26 };

    public double getMinGenotypeQuality() {
        return minGenotypeQuality;
    }

   /**
     * Constructor
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     *
     */
    public MendelianViolation(double minGenotypeQualityP) {
        this(minGenotypeQualityP,true);
    }

    /**
     * Constructor
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     * @param abortOnSampleNotFound - Whether to stop execution if a family is passed but no relevant genotypes are found. If false, then the family is ignored.
     */
    public MendelianViolation(double minGenotypeQualityP, boolean abortOnSampleNotFound) {
        minGenotypeQuality = minGenotypeQualityP;
        this.abortOnSampleNotFound = abortOnSampleNotFound;
        violationFamilies = new ArrayList<String>();
        createInheritanceMap();
    }

    /**
     * Constructor
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     * @param abortOnSampleNotFound - Whether to stop execution if a family is passed but no relevant genotypes are found. If false, then the family is ignored.
     * @param completeTriosOnly - whether only complete trios are considered or parent/child pairs are too.
     */
    public MendelianViolation(double minGenotypeQualityP, boolean abortOnSampleNotFound, boolean completeTriosOnly) {
        minGenotypeQuality = minGenotypeQualityP;
        this.abortOnSampleNotFound = abortOnSampleNotFound;
        violationFamilies = new ArrayList<String>();
        createInheritanceMap();
        allCalledOnly = completeTriosOnly;
    }

    /**
     * @param families the families to be checked for Mendelian violations
     * @param vc the variant context to extract the genotypes and alleles for mom, dad and child.
     * @return whether or not there is a mendelian violation at the site.
     */
    public int countViolations(Map<String, Set<Sample>> families, VariantContext vc){

        //Reset counts
        nocall = 0;
        lowQual = 0;
        familyCalled = 0;
        varFamilyCalled = 0;
        violations_total=0;
        violationFamilies.clear();
        clearInheritanceMap();

        for(Set<Sample> family : families.values()){
            Iterator<Sample> sampleIterator = family.iterator();
            Sample sample;
            while(sampleIterator.hasNext()){
                sample = sampleIterator.next();
                if(sample.getParents().size() > 0)
                    updateViolations(sample.getFamilyID(),sample.getMaternalID(), sample.getPaternalID(), sample.getID() ,vc);
            }
        }
        return violations_total;
    }

    public boolean isViolation(Sample mother, Sample father, Sample child, VariantContext vc){

        //Reset counts
        nocall = 0;
        lowQual = 0;
        familyCalled = 0;
        varFamilyCalled = 0;
        violations_total=0;
        violationFamilies.clear();
        clearInheritanceMap();
        updateViolations(mother.getFamilyID(),mother.getID(),father.getID(),child.getID(),vc);
        return violations_total>0;
    }


    private void updateViolations(String familyId, String motherId, String fatherId, String childId, VariantContext vc){

            int count;
            Genotype gMom = vc.getGenotype(motherId);
            Genotype gDad = vc.getGenotype(fatherId);
            Genotype gChild = vc.getGenotype(childId);

            if (gMom == null || gDad == null || gChild == null){
                if(abortOnSampleNotFound)
                    throw new IllegalArgumentException(String.format("Variant %s:%d: Missing genotypes for family %s: mom=%s dad=%s family=%s", vc.getChr(), vc.getStart(), familyId, motherId, fatherId, childId));
                else
                    return;
            }
            //Count No calls
            if(allCalledOnly && (!gMom.isCalled() || !gDad.isCalled() || !gChild.isCalled())){
                nocall++;
            }
            else if (!gMom.isCalled() && !gDad.isCalled() || !gChild.isCalled()){
                nocall++;
            }
            //Count lowQual. Note that if min quality is set to 0, even values with no quality associated are returned
            else if (minGenotypeQuality>0 && (gMom.getPhredScaledQual()   < minGenotypeQuality ||
                gDad.getPhredScaledQual()   < minGenotypeQuality ||
                gChild.getPhredScaledQual() < minGenotypeQuality )) {
                lowQual++;
            }
            else{
                //Count all families per loci called
                familyCalled++;
                //If the family is all homref, not too interesting
                if(!(gMom.isHomRef() && gDad.isHomRef() && gChild.isHomRef()))
                {
                    varFamilyCalled++;
                    if(isViolation(gMom, gDad, gChild)){
                        violationFamilies.add(familyId);
                        violations_total++;
                    }
                }
                count = inheritance.get(gMom.getType()).get(gDad.getType()).get(gChild.getType());
                inheritance.get(gMom.getType()).get(gDad.getType()).put(gChild.getType(),count+1);

            }
    }

    /**
     * Evaluate the genotypes of mom, dad, and child to detect Mendelian violations
     *
     * @param gMom
     * @param gDad
     * @param gChild
     * @return true if the three genotypes represent a Mendelian violation; false otherwise
     */
    public static boolean isViolation(final Genotype gMom, final Genotype gDad, final Genotype gChild) {
        //1 parent is no "call
        if(!gMom.isCalled()){
            return (gDad.isHomRef() && gChild.isHomVar()) || (gDad.isHomVar() && gChild.isHomRef());
        }
        else if(!gDad.isCalled()){
            return (gMom.isHomRef() && gChild.isHomVar()) || (gMom.isHomVar() && gChild.isHomRef());
        }
        //Both parents have genotype information
        return !(gMom.getAlleles().contains(gChild.getAlleles().get(0)) && gDad.getAlleles().contains(gChild.getAlleles().get(1)) ||
            gMom.getAlleles().contains(gChild.getAlleles().get(1)) && gDad.getAlleles().contains(gChild.getAlleles().get(0)));
    }

    private void createInheritanceMap(){

        inheritance = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>>(GenotypeType.class);
        for(GenotypeType mType : GenotypeType.values()){
            inheritance.put(mType, new EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>(GenotypeType.class));
            for(GenotypeType dType : GenotypeType.values()){
                inheritance.get(mType).put(dType, new EnumMap<GenotypeType,Integer>(GenotypeType.class));
                for(GenotypeType cType : GenotypeType.values()){
                    inheritance.get(mType).get(dType).put(cType, 0);
                }
            }
        }

    }

    private void clearInheritanceMap(){
        for(GenotypeType mType : GenotypeType.values()){
            for(GenotypeType dType : GenotypeType.values()){
                for(GenotypeType cType : GenotypeType.values()){
                    inheritance.get(mType).get(dType).put(cType, 0);
                }
            }
        }
    }

    /**
     * @return the likelihood ratio for a mendelian violation
     */
    public double violationLikelihoodRatio(VariantContext vc, String motherId, String fatherId, String childId) {
        double[] logLikAssignments = new double[27];
        // the matrix to set up is
        // MOM   DAD    CHILD
        //                    |-  AA
        //   AA     AA    |    AB
        //                    |-   BB
        //                    |- AA
        //  AA     AB     |   AB
        //                    |- BB
        // etc. The leaves are counted as 0-11 for MVs and 0-14 for non-MVs
        double[] momGL = vc.getGenotype(motherId).getLikelihoods().getAsVector();
        double[] dadGL = vc.getGenotype(fatherId).getLikelihoods().getAsVector();
        double[] childGL = vc.getGenotype(childId).getLikelihoods().getAsVector();
        int offset = 0;
        for ( int oMom = 0; oMom < 3; oMom++ ) {
            for ( int oDad = 0; oDad < 3; oDad++ ) {
                for ( int oChild = 0; oChild < 3; oChild ++ ) {
                    logLikAssignments[offset++] = momGL[oMom] + dadGL[oDad] + childGL[oChild];
                }
            }
        }
        double[] mvLiks = new double[12];
        double[] nonMVLiks = new double[15];
        for ( int i = 0; i < 12; i ++ ) {
            mvLiks[i] = logLikAssignments[mvOffsets[i]];
        }

        for ( int i = 0; i < 15; i++) {
            nonMVLiks[i] = logLikAssignments[nonMVOffsets[i]];
        }

        return MathUtils.log10sumLog10(mvLiks) - MathUtils.log10sumLog10(nonMVLiks);
    }

}
