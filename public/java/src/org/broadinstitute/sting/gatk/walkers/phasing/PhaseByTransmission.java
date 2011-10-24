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
 * Computes the most likely genotype combination and phases trios and parent/child pairs
 *
 * <p>
 * PhaseByTransmission is a GATK tool that 1) computes the most likely genotype combination and phases trios and parent/child pairs given their genotype likelihoods and a mutation prior and 2) phases
 * all sites were parent/child transmission can be inferred unambiguously. It reports the genotype combination (and hence phasing) probability.
 * Ambiguous sites are:
 * <ul>
 *     <li>Sites where all individuals are heterozygous</li>
 *     <li>Sites where there is a Mendelian violation</li>
 * </ul>
 * Missing genotypes are handled as follows:
 * <ul>
 *     <li>In parent/child pairs: If an individual genotype is missing at one site, the other one is phased if it is homozygous. No phasing probability is emitted.</li>
 *     <li>In trios: If the child is missing, parents are treated as separate individuals and phased if homozygous. No phasing probability is emitted.</li>
 *     <li>In trios: If one of the parents is missing, it is handled like a parent/child pair. Phasing is done unless both the parent and child are heterozygous and a phasing probabilitt is emitted.</li>
 *     <li>In trios: If two individuals are missing, the remaining individual is phased if it is homozygous. No phasing probability is emitted.</li>
 * </ul>
 *
 * <h2>Input</h2>
 * <p>
 * <ul>
 *     <li>A VCF variant set containing trio(s) and/or parent/child pair(s).</li>
 *     <li>A PED pedigree file containing the description of the individuals relationships.</li>
 * </ul>
 * </p>
 *
 * <h2>Options</h2>
 * <p>
 *     <ul>
 *         <li>MendelianViolationsFile: An optional argument for reporting. If a file is specified, all sites that remain in mendelian violation after being assigned the most likely genotype
 *         combination will be reported there. Information reported: chromosome, position, filter, allele count in VCF, family, transmission probability,
 *         and each individual genotype, depth, allelic depth and likelihoods.</li>
 *         <li>DeNovoPrior: Mutation prio; default is 1e-8</li>
 *     </ul>
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An VCF with genotypes recalibrated as most likely under the familial constraint and phased by descent where non ambiguous..
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T PhaseByTransmission \
 *   -V input.vcf \
 *   -ped input.ped \
 *   -o output.vcf
 * </pre>
 *
 */
public class PhaseByTransmission extends RodWalker<HashMap<Byte,Integer>, HashMap<Byte,Integer>> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName = "mvf",required = false,fullName = "MendelianViolationsFile", doc="File to output the mendelian violation details.")
    private PrintStream mvFile = null;

    @Argument(shortName = "prior",required = false,fullName = "DeNovoPrior", doc="Prior for de novo mutations. Default: 1e-8")
    private double deNovoPrior=1e-8;

    @Output
    protected VCFWriter vcfWriter = null;

    private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
    private final String SOURCE_NAME = "PhaseByTransmission";

    public final double NO_TRANSMISSION_PROB = -1.0;

    private ArrayList<Sample> trios = new ArrayList<Sample>();

    //Matrix of priors for all genotype combinations
    private EnumMap<Genotype.Type,EnumMap<Genotype.Type,EnumMap<Genotype.Type,Integer>>> mvCountMatrix;

    //Matrix of allele transmission
    private EnumMap<Genotype.Type,EnumMap<Genotype.Type,EnumMap<Genotype.Type,TrioPhase>>> transmissionMatrix;

    //Metrics counters hash keys
    private final Byte NUM_TRIO_GENOTYPES_CALLED = 0;
    private final Byte NUM_TRIO_GENOTYPES_NOCALL = 1;
    private final Byte NUM_TRIO_GENOTYPES_PHASED = 2;
    private final Byte NUM_HET = 3;
    private final Byte NUM_HET_HET_HET = 4;
    private final Byte NUM_VIOLATIONS = 5;

    private enum FamilyMember {
        MOTHER,
        FATHER,
        CHILD
    }

    //Stores a conceptual trio or parent/child pair genotype combination along with its phasing.
    //This combination can then be "applied" to a given trio or pair using the getPhasedGenotypes method.
    private class TrioPhase {

        //Create 2 fake alleles
        //The actual bases will never be used but the Genotypes created using the alleles will be.
        private final Allele REF = Allele.create("A",true);
        private final Allele VAR = Allele.create("A",false);
        private final Allele NO_CALL = Allele.create(".",false);
        private final String DUMMY_NAME = "DummySample";

       private EnumMap<FamilyMember,Genotype> trioPhasedGenotypes = new EnumMap<FamilyMember, Genotype>(FamilyMember.class);

        private ArrayList<Allele> getAlleles(Genotype.Type genotype){
        ArrayList<Allele> alleles = new ArrayList<Allele>(2);
        if(genotype == Genotype.Type.HOM_REF){
            alleles.add(REF);
            alleles.add(REF);
        }
        else if(genotype == Genotype.Type.HET){
            alleles.add(REF);
            alleles.add(VAR);
        }
        else if(genotype == Genotype.Type.HOM_VAR){
            alleles.add(VAR);
            alleles.add(VAR);
        }
        else if(genotype == Genotype.Type.NO_CALL){
            alleles.add(NO_CALL);
            alleles.add(NO_CALL);
        }
        else{
            return null;
        }
        return alleles;
    }

    //Create a new Genotype based on information from a single individual
    //Homozygous genotypes will be set as phased, heterozygous won't be
    private void phaseSingleIndividualAlleles(Genotype.Type genotype, FamilyMember familyMember){
        if(genotype == Genotype.Type.HOM_REF || genotype == Genotype.Type.HOM_VAR){
            trioPhasedGenotypes.put(familyMember, new Genotype(DUMMY_NAME, getAlleles(genotype), Genotype.NO_NEG_LOG_10PERROR, null, null, true));
        }
        else
            trioPhasedGenotypes.put(familyMember, new Genotype(DUMMY_NAME,getAlleles(genotype),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
    }

    //Find the phase for a parent/child pair
    private void phasePairAlleles(Genotype.Type parentGenotype, Genotype.Type childGenotype, FamilyMember parent){

        //Special case for Het/Het as it is ambiguous
        if(parentGenotype == Genotype.Type.HET && childGenotype == Genotype.Type.HET){
            trioPhasedGenotypes.put(parent, new Genotype(DUMMY_NAME, getAlleles(parentGenotype), Genotype.NO_NEG_LOG_10PERROR, null, null, false));
            trioPhasedGenotypes.put(FamilyMember.CHILD, new Genotype(DUMMY_NAME,getAlleles(childGenotype),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
        }

        ArrayList<Allele> parentAlleles = getAlleles(parentGenotype);
        ArrayList<Allele> childAlleles = getAlleles(childGenotype);
        ArrayList<Allele> parentPhasedAlleles = new ArrayList<Allele>(2);
        ArrayList<Allele> childPhasedAlleles = new ArrayList<Allele>(2);

        //If there is a possible phasing between the mother and child => phase
        int childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(0));
        if(childTransmittedAlleleIndex > -1){
           trioPhasedGenotypes.put(parent, new Genotype(DUMMY_NAME, parentAlleles, Genotype.NO_NEG_LOG_10PERROR, null, null, true));
           childPhasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
           childPhasedAlleles.add(childAlleles.get(0));
           trioPhasedGenotypes.put(FamilyMember.CHILD, new Genotype(DUMMY_NAME, childPhasedAlleles, Genotype.NO_NEG_LOG_10PERROR, null, null, true));
        }
        else if((childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(1))) > -1){
           parentPhasedAlleles.add(parentAlleles.get(1));
           parentPhasedAlleles.add(parentAlleles.get(0));
           trioPhasedGenotypes.put(parent, new Genotype(DUMMY_NAME, parentPhasedAlleles, Genotype.NO_NEG_LOG_10PERROR, null, null, true));
           childPhasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
           childPhasedAlleles.add(childAlleles.get(0));
           trioPhasedGenotypes.put(FamilyMember.CHILD, new Genotype(DUMMY_NAME, childPhasedAlleles, Genotype.NO_NEG_LOG_10PERROR, null, null, true));
        }
        //This is a Mendelian Violation => Do not phase
        else{
            trioPhasedGenotypes.put(parent, new Genotype(DUMMY_NAME,getAlleles(parentGenotype),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
            trioPhasedGenotypes.put(FamilyMember.CHILD, new Genotype(DUMMY_NAME,getAlleles(childGenotype),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
        }
    }

    //Phases a family by transmission
    private void phaseFamilyAlleles(Genotype.Type mother, Genotype.Type father, Genotype.Type child){

        Set<ArrayList<Allele>> possiblePhasedChildGenotypes = new HashSet<ArrayList<Allele>>();
        ArrayList<Allele> motherAlleles = getAlleles(mother);
        ArrayList<Allele> fatherAlleles = getAlleles(father);
        ArrayList<Allele> childAlleles = getAlleles(child);

        //Build all possible child genotypes for the given parent's genotypes
        for (Allele momAllele : motherAlleles) {
            for (Allele fatherAllele : fatherAlleles) {
                ArrayList<Allele> possiblePhasedChildAlleles = new ArrayList<Allele>(2);
                possiblePhasedChildAlleles.add(momAllele);
                possiblePhasedChildAlleles.add(fatherAllele);
                possiblePhasedChildGenotypes.add(possiblePhasedChildAlleles);
            }
        }

        for (ArrayList<Allele> childPhasedAllelesAlleles : possiblePhasedChildGenotypes) {
            int firstAlleleIndex = childPhasedAllelesAlleles.indexOf(childAlleles.get(0));
            int secondAlleleIndex = childPhasedAllelesAlleles.lastIndexOf(childAlleles.get(1));
            //If a possible combination has been found, create the genotypes
            if (firstAlleleIndex != secondAlleleIndex && firstAlleleIndex > -1 && secondAlleleIndex > -1) {
                //Create mother's genotype
                ArrayList<Allele> motherPhasedAlleles = new ArrayList<Allele>(2);
                motherPhasedAlleles.add(childPhasedAllelesAlleles.get(0));
                if(motherAlleles.get(0) != motherPhasedAlleles.get(0))
                    motherPhasedAlleles.add(motherAlleles.get(0));
                else
                    motherPhasedAlleles.add(motherAlleles.get(1));
                trioPhasedGenotypes.put(FamilyMember.MOTHER, new Genotype(DUMMY_NAME,motherPhasedAlleles,Genotype.NO_NEG_LOG_10PERROR,null,null,true));

                //Create father's genotype
                ArrayList<Allele> fatherPhasedAlleles = new ArrayList<Allele>(2);
                fatherPhasedAlleles.add(childPhasedAllelesAlleles.get(1));
                if(fatherAlleles.get(0) != fatherPhasedAlleles.get(0))
                    fatherPhasedAlleles.add(fatherAlleles.get(0));
                else
                    fatherPhasedAlleles.add(fatherAlleles.get(1));
                trioPhasedGenotypes.put(FamilyMember.FATHER, new Genotype(DUMMY_NAME,fatherPhasedAlleles,Genotype.NO_NEG_LOG_10PERROR,null,null,true));

                //Create child's genotype
                trioPhasedGenotypes.put(FamilyMember.CHILD, new Genotype(DUMMY_NAME,childPhasedAllelesAlleles,Genotype.NO_NEG_LOG_10PERROR,null,null,true));

                //Once a phased combination is found; exit
                return;
            }
        }

        //If this is reached then no phasing could be found
        trioPhasedGenotypes.put(FamilyMember.MOTHER, new Genotype(DUMMY_NAME,getAlleles(mother),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
        trioPhasedGenotypes.put(FamilyMember.FATHER, new Genotype(DUMMY_NAME,getAlleles(father),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
        trioPhasedGenotypes.put(FamilyMember.CHILD, new Genotype(DUMMY_NAME,getAlleles(child),Genotype.NO_NEG_LOG_10PERROR,null,null,false));
    }

    /*  Constructor: Creates a conceptual trio genotype combination from the given genotypes.
        If one or more genotypes are set as NO_CALL or UNAVAILABLE, it will phase them like a pair
        or single individual.
    */
    public TrioPhase(Genotype.Type mother, Genotype.Type father, Genotype.Type child){

        //Take care of cases where one or more family members are no call
        if(child == Genotype.Type.NO_CALL || child == Genotype.Type.UNAVAILABLE){
            phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
            phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
            phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
        }
        else if(mother == Genotype.Type.NO_CALL || mother == Genotype.Type.UNAVAILABLE){
            phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
            if(father == Genotype.Type.NO_CALL || father == Genotype.Type.UNAVAILABLE){
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
            }
            else
                phasePairAlleles(father, child, FamilyMember.FATHER);
        }
        else if(father == Genotype.Type.NO_CALL || father == Genotype.Type.UNAVAILABLE){
            phasePairAlleles(mother, child, FamilyMember.MOTHER);
            phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
        }
        //Special case for Het/Het/Het as it is ambiguous
        else if(mother == Genotype.Type.HET && father  == Genotype.Type.HET && child == Genotype.Type.HET){
            phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
            phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
            phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
        }
        //All family members have genotypes and at least one of them is not Het
        else{
            phaseFamilyAlleles(mother, father, child);
        }
    }

        /**
         * Applies the trio genotype combination to the given trio.
         * @param ref: Reference allele
         * @param alt: Alternate allele
         * @param motherGenotype: Genotype of the mother to phase using this trio genotype combination
         * @param fatherGenotype: Genotype of the father to phase using this trio genotype combination
         * @param childGenotype: Genotype of the child to phase using this trio genotype combination
         * @param transmissionProb: Probability for this trio genotype combination to be correct (pass NO_TRANSMISSION_PROB if unavailable)
         * @param phasedGenotypes: An ArrayList<Genotype> to which the newly phased genotypes are added in the following order: Mother, Father, Child
         */
    public void getPhasedGenotypes(Allele ref, Allele alt, Genotype motherGenotype, Genotype fatherGenotype, Genotype childGenotype, double transmissionProb,ArrayList<Genotype> phasedGenotypes){
        phasedGenotypes.add(getPhasedGenotype(ref,alt,motherGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.MOTHER)));
        phasedGenotypes.add(getPhasedGenotype(ref,alt,fatherGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.FATHER)));
        phasedGenotypes.add(getPhasedGenotype(ref,alt,childGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.CHILD)));
    }

   private Genotype getPhasedGenotype(Allele refAllele, Allele altAllele, Genotype genotype, double transmissionProb, Genotype phasedGenotype){

       //Handle null, missing and unavailable genotypes
       //Note that only cases where a null/missing/unavailable genotype was passed in the first place can lead to a null/missing/unavailable
       //genotype so it is safe to return the original genotype in this case.
       if(genotype == null || !phasedGenotype.isAvailable() || phasedGenotype.isNoCall())
           return genotype;

       //Add the transmission probability
       Map<String, Object> genotypeAttributes = new HashMap<String, Object>();
       genotypeAttributes.putAll(genotype.getAttributes());
       if(transmissionProb>NO_TRANSMISSION_PROB)
            genotypeAttributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, MathUtils.probabilityToPhredScale(1-(transmissionProb)));

       ArrayList<Allele> phasedAlleles = new ArrayList<Allele>(2);
       for(Allele allele : phasedGenotype.getAlleles()){
           if(allele.isReference())
               phasedAlleles.add(refAllele);
           else if(allele.isNonReference())
               phasedAlleles.add(altAllele);
           //At this point there should not be any other alleles left
           else
               throw new UserException(String.format("BUG: Unexpected allele: %s. Please report.",allele.toString()));

       }

       //Compute the new Log10Error if the genotype is different from the original genotype
       double negLog10Error;
       if(genotype.getType() == phasedGenotype.getType())
           negLog10Error = genotype.getNegLog10PError();
       else
          negLog10Error =  genotype.getLikelihoods().getNegLog10GQ(phasedGenotype.getType());

       return new Genotype(genotype.getSampleName(), phasedAlleles, negLog10Error, null, genotypeAttributes, phasedGenotype.isPhased());
   }


    }

    /**
     * Parse the familial relationship specification, build the transmission matrices and initialize VCF writer
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
        headerLines.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred score of the genotype combination and phase given that the genotypes are correct"));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        buildMatrices();

        if(mvFile != null)
            mvFile.println("#CHROM\tPOS\tFILTER\tAC\tFAMILY\tTP\tMOTHER_GT\tMOTHER_DP\tMOTHER_RAD\tMOTHER_AAD\tMOTHER_HRPL\tMOTHER_HETPL\tMOTHER_HAPL\tFATHER_GT\tFATHER_DP\tFATHER_RAD\tFATHER_AAD\tFATHER_HRPL\tFATHER_HETPL\tFATHER_HAPL\tCHILD_GT\tCHILD_DP\tCHILD_RAD\tCHILD_AAD\tCHILD_HRPL\tCHILD_HETPL\tCHILD_HAPL");

    }

    /**
     * Select trios and parent/child pairs only
     */
    private void setTrios(){

        Map<String,Set<Sample>> families = this.getSampleDB().getFamilies();
        Set<Sample> family;
        ArrayList<Sample> parents;
        for(String familyID : families.keySet()){
            family = families.get(familyID);
            if(family.size()<2 || family.size()>3){
                logger.info(String.format("Caution: Family %s has %d members; At the moment Phase By Transmission only supports trios and parent/child pairs. Family skipped.",familyID,family.size()));
            }
            else{
                for(Sample familyMember : family){
                     parents = familyMember.getParents();
                     if(parents.size()>0){
                        if(family.containsAll(parents))
                            this.trios.add(familyMember);
                        else
                            logger.info(String.format("Caution: Family %s skipped as it is not a trio nor a parent/child pair; At the moment Phase By Transmission only supports trios and parent/child pairs. Family skipped.",familyID));
                        break;
                     }
                }
            }

        }



    }

    //Create the transmission matrices
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

    //Returns the number of Mendelian Violations for a given genotype combination.
    //If one of the parents genotype is missing, it will consider it as a parent/child pair
    //If the child genotype or both parents genotypes are missing, 0 is returned.
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

    //Given two trio genotypes combinations, returns the number of different genotypes between the two combinations.
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

    //Get a Map of genotype likelihoods. If the genotype is NO_CALL or UNAVAILABLE, the Map will contain a single
    //NO_CALL resp. UNAVAILABLE element with a likelihood of 1.0
    private EnumMap<Genotype.Type,Double> getLikelihoodsAsMapSafeNull(Genotype genotype){
        if(genotype == null || !genotype.isAvailable()){
            EnumMap<Genotype.Type,Double> likelihoods = new EnumMap<Genotype.Type, Double>(Genotype.Type.class);
            likelihoods.put(Genotype.Type.UNAVAILABLE,1.0);
            return likelihoods;
        }
        else if(genotype.isNoCall()){
            EnumMap<Genotype.Type,Double> likelihoods = new EnumMap<Genotype.Type, Double>(Genotype.Type.class);
            likelihoods.put(Genotype.Type.NO_CALL,1.0);
            return likelihoods;
        }
        return genotype.getLikelihoods().getAsMap(true);
    }

    //Returns the Genotype.Type; returns UNVAILABLE if given null
    private Genotype.Type getTypeSafeNull(Genotype genotype){
        if(genotype == null)
            return Genotype.Type.UNAVAILABLE;
        return genotype.getType();
    }


    /**
     * Phases the genotypes of the given trio. If one of the parents is null, it is considered a parent/child pair.
     * @param ref: Reference allele
     * @param alt: Alternative allele
     * @param mother: Mother's genotype
     * @param father: Father's genotype
     * @param child: Child's genotype
     * @param finalGenotypes: An ArrayList<Genotype> that will be added the genotypes phased by transmission in the following order: Mother, Father, Child
     * @return
     */
    private boolean phaseTrioGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child,ArrayList<Genotype> finalGenotypes) {

        //Get the PL
        Map<Genotype.Type,Double> motherLikelihoods = getLikelihoodsAsMapSafeNull(mother);
        Map<Genotype.Type,Double> fatherLikelihoods = getLikelihoodsAsMapSafeNull(father);
        Map<Genotype.Type,Double> childLikelihoods = getLikelihoodsAsMapSafeNull(child);

        //Prior vars
        double bestConfigurationLikelihood = 0.0;
        double norm = 0.0;
        boolean isMV = false;
        int bestConfigurationGenotypeDiffs=4;
        Genotype.Type bestMotherGenotype = getTypeSafeNull(mother);
        Genotype.Type bestFatherGenotype = getTypeSafeNull(father);
        Genotype.Type bestChildGenotype = getTypeSafeNull(child);

        //Get the most likely combination
        //Only check for most likely combination if at least a parent and the child have genotypes
        if(childLikelihoods.size()>2 && (motherLikelihoods.size() + fatherLikelihoods.size())>3){
            int mvCount;
            double configurationLikelihood;
            int configurationGenotypeDiffs;
            for(Map.Entry<Genotype.Type,Double> motherGenotype : motherLikelihoods.entrySet()){
                for(Map.Entry<Genotype.Type,Double> fatherGenotype : fatherLikelihoods.entrySet()){
                    for(Map.Entry<Genotype.Type,Double> childGenotype : childLikelihoods.entrySet()){
                        mvCount = mvCountMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey());
                        configurationLikelihood =  mvCount>0 ? Math.pow(deNovoPrior,mvCount)*motherGenotype.getValue()*fatherGenotype.getValue()*childGenotype.getValue() : (1.0-11*deNovoPrior)*motherGenotype.getValue()*fatherGenotype.getValue()*childGenotype.getValue();
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

            //normalize the best configuration probability
            bestConfigurationLikelihood = bestConfigurationLikelihood / norm;
        }
        else{
            bestConfigurationLikelihood = NO_TRANSMISSION_PROB;
        }

        //Get the phased alleles for the genotype configuration
        TrioPhase phasedTrioGenotypes = transmissionMatrix.get(bestMotherGenotype).get(bestFatherGenotype).get(bestChildGenotype);

        //Return the phased genotypes
        phasedTrioGenotypes.getPhasedGenotypes(ref,alt,mother,father,child,bestConfigurationLikelihood,finalGenotypes);
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

                //Keep only trios and parent/child pairs
                if(mother == null && father == null || child == null)
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
                        if(phasedMother.isHet() && phasedFather.isHet() && phasedChild.isHet())
                            numHetHetHet++;
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
     * Initializes the reporting counters.
     *
     * @return All counters initialized to 0
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
     * Adds the value of the site phased to the reporting counters.
     *
     * @param value Site values
     * @param sum   accumulator for the reporting counters
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


    /**
     * Reports statistics on the phasing by transmission process.
     * @param result Accumulator with all counters.
     */
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
