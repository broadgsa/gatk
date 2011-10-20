package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Phases a trio VCF (child phased by transmission, implied phase carried over to parents).  Given genotypes for a trio,
 * this walker modifies the genotypes (if necessary) to reflect the most likely configuration given the genotype
 * likelihoods and inheritance constraints, phases child by transmission and carries over implied phase to the parents
 * (their alleles in their genotypes are ordered as transmitted|untransmitted).  Computes probability that the
 * determined phase is correct given that the genotype configuration is correct (useful if you want to use this to
 * compare phasing accuracy, but want to break that comparison down by phasing confidence in the truth set).  Optionally
 * filters out sites where the phasing is indeterminate (site has no-calls), ambiguous (everyone is heterozygous), or
 * the genotypes exhibit a Mendelian violation.  This walker assumes there are only three samples in the VCF file to
 * begin.
 */
public class PhaseByTransmission extends RodWalker<HashMap<Byte,Integer>, HashMap<Byte,Integer>> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName = "mvf",required = false,fullName = "MendelianViolationsFile", doc="File to output the mendelian violation details.")
    private PrintStream mvFile = null;

    @Output
    protected VCFWriter vcfWriter = null;

    private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
    private final String SOURCE_NAME = "PhaseByTransmission";

    private final Double MENDELIAN_VIOLATION_PRIOR = 1e-8;

    private ArrayList<Sample> trios = new ArrayList<Sample>();

    //Matrix of priors for all genotype combinations
    private EnumMap<Genotype.Type,EnumMap<Genotype.Type,EnumMap<Genotype.Type,Integer>>> mvCountMatrix;

    //Matrix of allele transmission
    private EnumMap<Genotype.Type,EnumMap<Genotype.Type,EnumMap<Genotype.Type,TrioPhase>>> transmissionMatrix;

    //Metrics counters hashkeys
    private final Byte NUM_TRIO_GENOTYPES_CALLED = 0;
    private final Byte NUM_TRIO_GENOTYPES_NOCALL = 1;
    private final Byte NUM_TRIO_GENOTYPES_PHASED = 2;
    private final Byte NUM_HET = 3;
    private final Byte NUM_HET_HET_HET = 4;
    private final Byte NUM_VIOLATIONS = 5;

   private enum AlleleType {
        NO_CALL,
        REF,
        VAR,
        UNPHASED_REF,
        UNPHASED_VAR
    }

    //Stores a trio-genotype
    private class TrioPhase {

        private ArrayList<AlleleType> trioAlleles = new ArrayList<AlleleType>(6);

        private ArrayList<AlleleType> getAlleles(Genotype.Type genotype){
        ArrayList<AlleleType> alleles = new ArrayList<AlleleType>(2);
        if(genotype == Genotype.Type.HOM_REF){
            alleles.add(AlleleType.REF);
            alleles.add(AlleleType.REF);
        }
        else if(genotype == Genotype.Type.HET){
            alleles.add(AlleleType.REF);
            alleles.add(AlleleType.VAR);
        }
        else if(genotype == Genotype.Type.HOM_VAR){
            alleles.add(AlleleType.VAR);
            alleles.add(AlleleType.VAR);
        }
        else{
            alleles.add(AlleleType.NO_CALL);
            alleles.add(AlleleType.NO_CALL);
        }
        return alleles;
    }

    private ArrayList<AlleleType> phaseSingleIndividualAlleles(Genotype.Type genotype){
        if(genotype == Genotype.Type.HET){
            ArrayList<AlleleType> phasedAlleles = new ArrayList<AlleleType>(2);
            phasedAlleles.add(AlleleType.UNPHASED_REF);
            phasedAlleles.add(AlleleType.UNPHASED_VAR);
            return phasedAlleles;
        }
        else
            return  getAlleles(genotype);
    }

    private ArrayList<AlleleType> phaseMonoParentFamilyAlleles(Genotype.Type parent, Genotype.Type child){
        ArrayList<AlleleType> phasedAlleles = new ArrayList<AlleleType>(4);
        //Special case for Het/Het as it is ambiguous
        if(parent == Genotype.Type.HET && child == Genotype.Type.HET){
            phasedAlleles.add(AlleleType.UNPHASED_REF);
            phasedAlleles.add(AlleleType.UNPHASED_VAR);
            phasedAlleles.add(AlleleType.UNPHASED_REF);
            phasedAlleles.add(AlleleType.UNPHASED_VAR);
        }

        ArrayList<AlleleType> parentAlleles = getAlleles(parent);
        ArrayList<AlleleType> childAlleles = getAlleles(child);

        int childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(0));
        if(childTransmittedAlleleIndex > -1){
           phasedAlleles.add(parentAlleles.get(0));
           phasedAlleles.add(parentAlleles.get(1));
           phasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
           phasedAlleles.add(childAlleles.get(0));
        }
        else if((childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(1))) > -1){
           phasedAlleles.add(parentAlleles.get(1));
           phasedAlleles.add(parentAlleles.get(0));
           phasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
           phasedAlleles.add(childAlleles.get(0));
        }
        else{
            parentAlleles.addAll(childAlleles);
            for(AlleleType allele : parentAlleles){
                if(allele == AlleleType.REF){
                    phasedAlleles.add(AlleleType.UNPHASED_REF);
                }
                else if(allele == AlleleType.VAR){
                    phasedAlleles.add(AlleleType.UNPHASED_VAR);
                }
                else{
                    phasedAlleles.add(AlleleType.NO_CALL);
                }
            }
        }

        return phasedAlleles;
    }

    private ArrayList<AlleleType> phaseFamilyAlleles(Genotype.Type mother, Genotype.Type father, Genotype.Type child){
        ArrayList<AlleleType> phasedAlleles = new ArrayList<AlleleType>(6);

        Set<ArrayList<AlleleType>> possiblePhasedChildGenotypes = new HashSet<ArrayList<AlleleType>>();
        ArrayList<AlleleType> motherAlleles = getAlleles(mother);
        ArrayList<AlleleType> fatherAlleles = getAlleles(father);
        ArrayList<AlleleType> childAlleles = getAlleles(child);

        //Build all possible child genotypes for the given parent's genotypes
        for (AlleleType momAllele : motherAlleles) {
            for (AlleleType fatherAllele : fatherAlleles) {
                ArrayList<AlleleType> possiblePhasedChildAlleles = new ArrayList<AlleleType>(2);
                possiblePhasedChildAlleles.add(momAllele);
                possiblePhasedChildAlleles.add(fatherAllele);
                possiblePhasedChildGenotypes.add(possiblePhasedChildAlleles);
            }
        }

        for (ArrayList<AlleleType> phasedChildGenotype : possiblePhasedChildGenotypes) {
            int firstAlleleIndex = phasedChildGenotype.indexOf(childAlleles.get(0));
            int secondAlleleIndex = phasedChildGenotype.lastIndexOf(childAlleles.get(1));
            if (firstAlleleIndex != secondAlleleIndex && firstAlleleIndex > -1 && secondAlleleIndex > -1) {
                //Add mother's alleles
                phasedAlleles.add(phasedChildGenotype.get(0));
                if(motherAlleles.get(0) != phasedAlleles.get(0))
                    phasedAlleles.add(motherAlleles.get(0));
                else
                    phasedAlleles.add(motherAlleles.get(1));

                //Add father's alleles
                phasedAlleles.add(phasedChildGenotype.get(1));
                if(fatherAlleles.get(0) != phasedAlleles.get(2))
                    phasedAlleles.add(fatherAlleles.get(0));
                else
                    phasedAlleles.add(fatherAlleles.get(1));

                //Add child's alleles
                phasedAlleles.addAll(phasedChildGenotype);
                return phasedAlleles;
            }
        }

        //If this is reached then no phasing could be found
        motherAlleles.addAll(fatherAlleles);
        motherAlleles.addAll(childAlleles);
        for(AlleleType allele : motherAlleles){
            if(allele == AlleleType.REF){
                    phasedAlleles.add(AlleleType.UNPHASED_REF);
                }
                else if(allele == AlleleType.VAR){
                    phasedAlleles.add(AlleleType.UNPHASED_VAR);
                }
                else{
                    phasedAlleles.add(AlleleType.NO_CALL);
                }
        }
        return phasedAlleles;
    }


    public TrioPhase(Genotype.Type mother, Genotype.Type father, Genotype.Type child){

        //Take care of cases where one or more family members are no call
        if(child == Genotype.Type.NO_CALL){
            trioAlleles.addAll(phaseSingleIndividualAlleles(mother));
            trioAlleles.addAll(phaseSingleIndividualAlleles(father));
            trioAlleles.add(AlleleType.NO_CALL);
            trioAlleles.add(AlleleType.NO_CALL);
        }
        else if(mother == Genotype.Type.NO_CALL){
            trioAlleles.add(AlleleType.NO_CALL);
            trioAlleles.add(AlleleType.NO_CALL);
            if(father == Genotype.Type.NO_CALL){
                trioAlleles.add(AlleleType.NO_CALL);
                trioAlleles.add(AlleleType.NO_CALL);
                trioAlleles.addAll(phaseSingleIndividualAlleles(child));
            }
            else
                trioAlleles.addAll(phaseMonoParentFamilyAlleles(father, child));
        }
        else if(father == Genotype.Type.NO_CALL){
            trioAlleles.addAll(phaseMonoParentFamilyAlleles(mother, child));
            trioAlleles.add(2, AlleleType.NO_CALL);
            trioAlleles.add(3, AlleleType.NO_CALL);
        }
        //Special case for Het/Het/Het as it is ambiguous
        else if(mother == Genotype.Type.HET && father  == Genotype.Type.HET && child == Genotype.Type.HET){
            trioAlleles.add(AlleleType.UNPHASED_REF);
            trioAlleles.add(AlleleType.UNPHASED_VAR);
            trioAlleles.add(AlleleType.UNPHASED_REF);
            trioAlleles.add(AlleleType.UNPHASED_VAR);
            trioAlleles.add(AlleleType.UNPHASED_REF);
            trioAlleles.add(AlleleType.UNPHASED_VAR);
        }
        //All family members have genotypes and at least one of them is not Het
        else{
            trioAlleles = phaseFamilyAlleles(mother, father, child);
        }
    }

    public ArrayList<Genotype> getPhasedGenotypes(Allele ref, Allele alt, Genotype motherGenotype, Genotype fatherGenotype, Genotype childGenotype, int transmissionProb,ArrayList<Genotype> phasedGenotypes){
        phasedGenotypes.add(getPhasedGenotype(ref,alt,motherGenotype,transmissionProb,trioAlleles.subList(0,2)));
        phasedGenotypes.add(getPhasedGenotype(ref,alt,fatherGenotype,transmissionProb,trioAlleles.subList(2,4)));
        phasedGenotypes.add(getPhasedGenotype(ref,alt,childGenotype,transmissionProb,trioAlleles.subList(4,6)));
        return phasedGenotypes;
    }

   private Genotype getPhasedGenotype(Allele refAllele, Allele altAllele, Genotype genotype, int transmissionProb, List<AlleleType> phasedAlleles){

       //Add the transmission probability
       Map<String, Object> genotypeAttributes = new HashMap<String, Object>();
       genotypeAttributes.putAll(genotype.getAttributes());
       genotypeAttributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, transmissionProb);
       genotype = Genotype.modifyAttributes(genotype, genotypeAttributes);

       boolean isPhased = true;

       List<Allele> alleles = new ArrayList<Allele>(2);

        //If unphased, return original genotype
        for(AlleleType allele : phasedAlleles){
            if(allele == AlleleType.NO_CALL){
                return genotype;
            }
            //Otherwise add the appropriate allele
            else if(allele == AlleleType.UNPHASED_REF){
                isPhased = false;
                alleles.add(refAllele);
            }
            else if(allele == AlleleType.UNPHASED_VAR){
                isPhased = false;
                alleles.add(altAllele);
            }
            else if(allele == AlleleType.REF){
                alleles.add(refAllele);
            }
            else if(allele == AlleleType.VAR){
                alleles.add(altAllele);
            }
        }

       return new Genotype(genotype.getSampleName(), alleles, genotype.getLikelihoods().getLog10GQ(genotype.getType()), null, genotype.getAttributes(), isPhased);
   }


    }

    /**
     * Parse the familial relationship specification, and initialize VCF writer
     */
    public void initialize() {
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //Get the trios from the families passed as ped
        setTrios();
        if(trios.size()<1)
            throw new UserException.BadInput("No PED file passed or no trios found in PED file. Aborted.");


        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred score of the phase given that the genotypes are correct"));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        buildMatrices();

        if(mvFile != null)
            mvFile.println("#CHROM\tPOS\tFILTER\tAC\tFAMILY\tTP\tMOTHER_GT\tMOTHER_DP\tMOTHER_RAD\tMOTHER_AAD\tMOTHER_HRPL\tMOTHER_HETPL\tMOTHER_HAPL\tFATHER_GT\tFATHER_DP\tFATHER_RAD\tFATHER_AAD\tFATHER_HRPL\tFATHER_HETPL\tFATHER_HAPL\tCHILD_GT\tCHILD_DP\tCHILD_RAD\tCHILD_AAD\tCHILD_HRPL\tCHILD_HETPL\tCHILD_HAPL");

    }

    /**
     * Select Trios only
     */
    private void setTrios(){

        Map<String,Set<Sample>> families = this.getSampleDB().getFamilies();
        Set<Sample> family;
        ArrayList<Sample> parents;
        for(String familyID : families.keySet()){
            family = families.get(familyID);
            if(family.size()!=3){
                logger.info(String.format("Caution: Family %s has %d members; At the moment Phase By Transmission only supports trios. Family skipped.",familyID,family.size()));
            }
            else{
                for(Sample familyMember : family){
                     parents = familyMember.getParents();
                     if(parents.size()>0){
                        if(family.containsAll(parents))
                            this.trios.add(familyMember);
                        else
                            logger.info(String.format("Caution: Family %s is not a trio; At the moment Phase By Transmission only supports trios. Family skipped.",familyID));
                        break;
                     }
                }
            }

        }



    }

    private void buildMatrices(){
        mvCountMatrix = new EnumMap<Genotype.Type,EnumMap<Genotype.Type,EnumMap<Genotype.Type,Integer>>>(Genotype.Type.class);
        transmissionMatrix = new EnumMap<Genotype.Type,EnumMap<Genotype.Type,EnumMap<Genotype.Type,TrioPhase>>>(Genotype.Type.class);
        for(Genotype.Type mother : Genotype.Type.values()){
           mvCountMatrix.put(mother,new EnumMap<Genotype.Type,EnumMap<Genotype.Type,Integer>>(Genotype.Type.class));
           transmissionMatrix.put(mother,new EnumMap<Genotype.Type,EnumMap<Genotype.Type,TrioPhase>>(Genotype.Type.class));
           for(Genotype.Type father : Genotype.Type.values()){
               mvCountMatrix.get(mother).put(father,new EnumMap<Genotype.Type, Integer>(Genotype.Type.class));
               transmissionMatrix.get(mother).put(father,new EnumMap<Genotype.Type,TrioPhase>(Genotype.Type.class));
                for(Genotype.Type child : Genotype.Type.values()){
                    mvCountMatrix.get(mother).get(father).put(child, getCombinationMVCount(mother, father, child));
                    transmissionMatrix.get(mother).get(father).put(child,new TrioPhase(mother,father,child));
                }
            }
        }
    }

    private int getCombinationMVCount(Genotype.Type mother, Genotype.Type father, Genotype.Type child){

        //Child is no call => No MV
        if(child == Genotype.Type.NO_CALL || child == Genotype.Type.UNAVAILABLE)
            return 0;
        //Add parents with genotypes for the evaluation
        ArrayList<Genotype.Type> parents = new ArrayList<Genotype.Type>();
        if (!(mother == Genotype.Type.NO_CALL || mother == Genotype.Type.UNAVAILABLE))
            parents.add(mother);
        if (!(father == Genotype.Type.NO_CALL || father == Genotype.Type.UNAVAILABLE))
            parents.add(father);

        //Both parents no calls => No MV
        if (parents.isEmpty())
            return 0;

        //If at least one parent had a genotype, then count the number of ref and alt alleles that can be passed
        int parentsNumRefAlleles = 0;
        int parentsNumAltAlleles = 0;

        for(Genotype.Type parent : parents){
            if(parent == Genotype.Type.HOM_REF){
                parentsNumRefAlleles++;
            }
            else if(parent == Genotype.Type.HET){
                parentsNumRefAlleles++;
                parentsNumAltAlleles++;
            }
            else if(parent == Genotype.Type.HOM_VAR){
                parentsNumAltAlleles++;
            }
        }

        //Case Child is HomRef
        if(child == Genotype.Type.HOM_REF){
            if(parentsNumRefAlleles == parents.size())
                return 0;
            else return (parents.size()-parentsNumRefAlleles);
        }

        //Case child is HomVar
        if(child == Genotype.Type.HOM_VAR){
            if(parentsNumAltAlleles == parents.size())
                return 0;
            else return parents.size()-parentsNumAltAlleles;
        }

        //Case child is Het
        if(child == Genotype.Type.HET && ((parentsNumRefAlleles > 0 && parentsNumAltAlleles > 0) || parents.size()<2))
            return 0;

        //MV
        return 1;
    }

    private Double getCombinationPrior(Genotype.Type mother, Genotype.Type father, Genotype.Type child){

        double nonMVPrior = 1.0 - 12*MENDELIAN_VIOLATION_PRIOR;

        //Child is no call => No MV
        if(child == Genotype.Type.NO_CALL || child == Genotype.Type.UNAVAILABLE)
            return nonMVPrior;
        //Add parents with genotypes for the evaluation
        ArrayList<Genotype.Type> parents = new ArrayList<Genotype.Type>();
        if (!(mother == Genotype.Type.NO_CALL || mother == Genotype.Type.UNAVAILABLE))
            parents.add(mother);
        if (!(father == Genotype.Type.NO_CALL || father == Genotype.Type.UNAVAILABLE))
            parents.add(father);

        //Both parents no calls => No MV
        if (parents.isEmpty())
            return nonMVPrior;

        //If at least one parent had a genotype, then count the number of ref and alt alleles that can be passed
        int parentsNumRefAlleles = 0;
        int parentsNumAltAlleles = 0;

        for(Genotype.Type parent : parents){
            if(parent == Genotype.Type.HOM_REF){
                parentsNumRefAlleles++;
            }
            else if(parent == Genotype.Type.HET){
                parentsNumRefAlleles++;
                parentsNumAltAlleles++;
            }
            else if(parent == Genotype.Type.HOM_VAR){
                parentsNumAltAlleles++;
            }
        }

        //Case Child is HomRef
        if(child == Genotype.Type.HOM_REF){
            if(parentsNumRefAlleles == parents.size())
                return nonMVPrior;
            else return Math.pow(MENDELIAN_VIOLATION_PRIOR, parents.size()-parentsNumRefAlleles);
        }

        //Case child is HomVar
        if(child == Genotype.Type.HOM_VAR){
            if(parentsNumAltAlleles == parents.size())
                return nonMVPrior;
            else return Math.pow(MENDELIAN_VIOLATION_PRIOR, parents.size()-parentsNumAltAlleles);
        }

        //Case child is Het
        if(child == Genotype.Type.HET && ((parentsNumRefAlleles > 0 && parentsNumAltAlleles > 0) || parents.size()<2))
            return nonMVPrior;

        //MV
        return MENDELIAN_VIOLATION_PRIOR;
    }

    private int countFamilyGenotypeDiff(Genotype.Type motherOriginal,Genotype.Type fatherOriginal,Genotype.Type childOriginal,Genotype.Type motherNew,Genotype.Type fatherNew,Genotype.Type childNew){
        int count = 0;
        if(motherOriginal!=motherNew)
            count++;
        if(fatherOriginal!=fatherNew)
            count++;
        if(childOriginal!=childNew)
            count++;
        return count;
    }


    private boolean phaseTrioGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child,ArrayList<Genotype> finalGenotypes) {

        //For now, only consider trios with complete information
        //TODO: Phasing of trios with missing information
        if(!(mother.isCalled() && father.isCalled() && child.isCalled())) {
            finalGenotypes.add(mother);
            finalGenotypes.add(father);
            finalGenotypes.add(child);
            return false;
        }

        //Get the PL
        Map<Genotype.Type,Double> motherLikelihoods = mother.getLikelihoods().getAsMap(true);
        Map<Genotype.Type,Double> fatherLikelihoods = father.getLikelihoods().getAsMap(true);
        Map<Genotype.Type,Double> childLikelihoods = child.getLikelihoods().getAsMap(true);

        //Prior vars
        double bestConfigurationLikelihood = 0.0;
        double norm = 0.0;
        boolean isMV = false;
        int bestConfigurationGenotypeDiffs=4;
        Genotype.Type bestMotherGenotype = mother.getType();
        Genotype.Type bestFatherGenotype = father.getType();
        Genotype.Type bestChildGenotype = child.getType();

        //Get the most likely combination
        int mvCount;
        double configurationLikelihood;
        int configurationGenotypeDiffs;
        for(Map.Entry<Genotype.Type,Double> motherGenotype : motherLikelihoods.entrySet()){
            for(Map.Entry<Genotype.Type,Double> fatherGenotype : fatherLikelihoods.entrySet()){
                for(Map.Entry<Genotype.Type,Double> childGenotype : childLikelihoods.entrySet()){
                    mvCount = mvCountMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey());
                    configurationLikelihood =  mvCount>0 ? Math.pow(MENDELIAN_VIOLATION_PRIOR,mvCount)*motherGenotype.getValue()*fatherGenotype.getValue()*childGenotype.getValue() : (1.0-12*MENDELIAN_VIOLATION_PRIOR)*motherGenotype.getValue()*fatherGenotype.getValue()*childGenotype.getValue();
                    norm += configurationLikelihood;
                    configurationGenotypeDiffs = countFamilyGenotypeDiff(mother.getType(),father.getType(),child.getType(),motherGenotype.getKey(),fatherGenotype.getKey(),childGenotype.getKey());
                    //Keep this combination if
                    //It has a better likelihood
                    //Or it has the same likelihood but requires less changes from original genotypes
                    if ((configurationLikelihood > bestConfigurationLikelihood) ||
                            (configurationLikelihood == bestConfigurationLikelihood && configurationGenotypeDiffs < bestConfigurationGenotypeDiffs)) {
                        bestConfigurationLikelihood = configurationLikelihood;
                        bestMotherGenotype = motherGenotype.getKey();
                        bestFatherGenotype = fatherGenotype.getKey();
                        bestChildGenotype = childGenotype.getKey();
                        isMV = mvCount>0;
                        bestConfigurationGenotypeDiffs=configurationGenotypeDiffs;
                    }
                }
            }
        }

        //Get the phased alleles for the genotype configuration
        TrioPhase phasedTrioGenotypes = transmissionMatrix.get(bestMotherGenotype).get(bestFatherGenotype).get(bestChildGenotype);

        //Return the phased genotypes
        phasedTrioGenotypes.getPhasedGenotypes(ref,alt,mother,father,child,MathUtils.probabilityToPhredScale(1-(bestConfigurationLikelihood / norm)),finalGenotypes);
        return isMV;

    }

    /**
     * For each variant in the file, determine the phasing for the child and replace the child's genotype with the trio's genotype
     *
     * @param tracker  the reference meta-data tracker
     * @param ref      the reference context
     * @param context  the alignment context
     * @return null
     */
    @Override
    public HashMap<Byte,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        //Local cars to avoid lookups on increment
        int numTrioGenotypesCalled = 0;
        int numTrioGenotypesNoCall = 0;
        int numTrioGenotypesPhased = 0;
        int numHet = 0 ;
        int numHetHetHet = 0;
        int numMVs = 0;

        if (tracker != null) {
            VariantContext vc = tracker.getFirstValue(variantCollection.variants, context.getLocation());

            Map<String, Genotype> genotypeMap = vc.getGenotypes();

            boolean isMV;

            for (Sample sample : trios) {
                Genotype mother = vc.getGenotype(sample.getMaternalID());
                Genotype father = vc.getGenotype(sample.getPaternalID());
                Genotype child = vc.getGenotype(sample.getID());

                //Skip trios where any of the genotype is missing in the variant context
                if(mother == null || father == null | child == null)
                    continue;

                ArrayList<Genotype> trioGenotypes = new ArrayList<Genotype>(3);
                isMV = phaseTrioGenotypes(vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), mother, father, child,trioGenotypes);

                Genotype phasedMother = trioGenotypes.get(0);
                Genotype phasedFather = trioGenotypes.get(1);
                Genotype phasedChild = trioGenotypes.get(2);

                genotypeMap.put(phasedMother.getSampleName(), phasedMother);
                genotypeMap.put(phasedFather.getSampleName(), phasedFather);
                genotypeMap.put(phasedChild.getSampleName(), phasedChild);

                //Increment metrics counters
                if(phasedMother.isCalled() && phasedFather.isCalled() && phasedChild.isCalled()){
                   numTrioGenotypesCalled++;
                   if(phasedMother.isPhased())
                       numTrioGenotypesPhased++;

                    if(phasedMother.isHet() || phasedFather.isHet() || phasedChild.isHet()){
                        numHet++;
                        if(phasedMother.isHet() && phasedFather.isHet() && phasedChild.isHet()){
                            numHetHetHet++;
                        }else if(!phasedMother.isPhased()){
                            int x =9;
                        }
                    }

                    if(isMV){
                        numMVs++;
                        if(mvFile != null)
                            mvFile.println(String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s:%s:%s:%s\t%s:%s:%s:%s\t%s:%s:%s:%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),sample.toString(),phasedMother.getAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),phasedMother.getGenotypeString(),phasedMother.getAttribute(VCFConstants.DEPTH_KEY),phasedMother.getAttribute("AD"),phasedMother.getLikelihoods().toString(),phasedFather.getGenotypeString(),phasedFather.getAttribute(VCFConstants.DEPTH_KEY),phasedFather.getAttribute("AD"),phasedFather.getLikelihoods().toString(),phasedChild.getGenotypeString(),phasedChild.getAttribute(VCFConstants.DEPTH_KEY),phasedChild.getAttribute("AD"),phasedChild.getLikelihoods().toString()));
                    }
                }else{
                    numTrioGenotypesNoCall++;
                }

            }


            VariantContext newvc = VariantContext.modifyGenotypes(vc, genotypeMap);

            vcfWriter.add(newvc);
        }

        HashMap<Byte,Integer> metricsCounters = new HashMap<Byte, Integer>(5);
        metricsCounters.put(NUM_TRIO_GENOTYPES_CALLED,numTrioGenotypesCalled);
        metricsCounters.put(NUM_TRIO_GENOTYPES_NOCALL,numTrioGenotypesNoCall);
        metricsCounters.put(NUM_TRIO_GENOTYPES_PHASED,numTrioGenotypesPhased);
        metricsCounters.put(NUM_HET,numHet);
        metricsCounters.put(NUM_HET_HET_HET,numHetHetHet);
        metricsCounters.put(NUM_VIOLATIONS,numMVs);
        return metricsCounters;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public HashMap<Byte,Integer> reduceInit() {
        HashMap<Byte,Integer> metricsCounters = new HashMap<Byte, Integer>(5);
        metricsCounters.put(NUM_TRIO_GENOTYPES_CALLED,0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_NOCALL,0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_PHASED,0);
        metricsCounters.put(NUM_HET,0);
        metricsCounters.put(NUM_HET_HET_HET,0);
        metricsCounters.put(NUM_VIOLATIONS,0);
        return metricsCounters;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public HashMap<Byte,Integer> reduce(HashMap<Byte,Integer> value, HashMap<Byte,Integer> sum) {
        sum.put(NUM_TRIO_GENOTYPES_CALLED,value.get(NUM_TRIO_GENOTYPES_CALLED)+sum.get(NUM_TRIO_GENOTYPES_CALLED));
        sum.put(NUM_TRIO_GENOTYPES_NOCALL,value.get(NUM_TRIO_GENOTYPES_NOCALL)+sum.get(NUM_TRIO_GENOTYPES_NOCALL));
        sum.put(NUM_TRIO_GENOTYPES_PHASED,value.get(NUM_TRIO_GENOTYPES_PHASED)+sum.get(NUM_TRIO_GENOTYPES_PHASED));
        sum.put(NUM_HET,value.get(NUM_HET)+sum.get(NUM_HET));
        sum.put(NUM_HET_HET_HET,value.get(NUM_HET_HET_HET)+sum.get(NUM_HET_HET_HET));
        sum.put(NUM_VIOLATIONS,value.get(NUM_VIOLATIONS)+sum.get(NUM_VIOLATIONS));
        return sum;
    }


    @Override
    public void onTraversalDone(HashMap<Byte,Integer> result) {
        logger.info("Number of complete trio-genotypes: " + result.get(NUM_TRIO_GENOTYPES_CALLED));
        logger.info("Number of trio-genotypes containing no call(s): " + result.get(NUM_TRIO_GENOTYPES_NOCALL));
        logger.info("Number of trio-genotypes phased: " + result.get(NUM_TRIO_GENOTYPES_PHASED));
        logger.info("Number of resulting Hom/Hom/Hom trios: " + result.get(NUM_HET));
        logger.info("Number of resulting Het/Het/Het trios: " + result.get(NUM_HET_HET_HET));
        logger.info("Number of remaining mendelian violations: " + result.get(NUM_VIOLATIONS));
    }
}
