package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.HapMapAlleleFrequenciesROD;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.Formatter;
import static java.lang.Math.log10;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import net.sf.picard.PicardException;

public class PopPriorWalker extends LocusWalker<Integer, Integer> {

    @Argument(fullName = "quality_score_cutoff", doc="quality score cutoff", required = false) public Byte qualityScoreCutoff;

    public PopPriorWalker() {
    }

    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    protected class OldAndBustedGenotypeLikelihoods
    {

        public double[] likelihoods;
        public String[] genotypes;

        OldAndBustedGenotypeLikelihoods()
        {
            likelihoods = new double[10];
            genotypes   = new String[10];

            genotypes[0] = "AA";
            genotypes[1] = "AC";
            genotypes[2] = "AG";
            genotypes[3] = "AT";
            genotypes[4] = "CC";
            genotypes[5] = "CG";
            genotypes[6] = "CT";
            genotypes[7] = "GG";
            genotypes[8] = "GT";
            genotypes[9] = "TT";
        }

        void add(char ref, char read, byte qual)
        {
            double p_error = Math.pow(10.0, (double)qual / -10);
            for (int i = 0; i < genotypes.length; i++)
            {
                likelihoods[i] += AlleleLikelihood(ref, read, genotypes[i], p_error);
            }
        }

        double AlleleLikelihood(char ref, char read, String genotype, double p_error)
        {
            char h1 = genotype.charAt(0);
            char h2 = genotype.charAt(1);

            double p_base;

            if      ((h1 == h2) && (h1 == read))                 { p_base = Math.log10(1-p_error); }
            else if ((h1 != h2) && (h1 == read) || (h2 == read)) { p_base = Math.log10(0.5 - (p_error/2.0)); }
            else                                                 { p_base = Math.log10(p_error); }

            return p_base;
        }

        // TODO: horrible horrible idea... but it's ok for right now.  should really
        // implement this more cleanly so you can get sorted genotypes without losing the
        // ability to call add
        public void doSort() {
            Integer[] permutation       = Utils.SortPermutation(likelihoods);
            String[] sorted_genotypes   = Utils.PermuteArray(genotypes, permutation);
            double[] sorted_likelihoods = Utils.PermuteArray(likelihoods, permutation);

            this.genotypes = sorted_genotypes;
            this.likelihoods = sorted_likelihoods;
        }

    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        char upRef = Character.toUpperCase(ref.getBase());
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        // Look up hapmap and dbsnp priors
        String rodString = "";

        rodDbSNP dbsnpInfo = null;
        HapMapAlleleFrequenciesROD hapmap = null;

        for ( ReferenceOrderedDatum datum : tracker.getAllRods() )
        {
            if ( datum != null )
            {
                if ( datum instanceof rodDbSNP) {
                    dbsnpInfo = (rodDbSNP)datum;

                } else if (datum instanceof HapMapAlleleFrequenciesROD) {
                    hapmap = (HapMapAlleleFrequenciesROD) datum;
                }
            }
        }

        // Accumulate genotype likelihoods
        int aCount = 0;
        int cCount = 0;
        int gCount = 0;
        int tCount = 0;
        OldAndBustedGenotypeLikelihoods glAll = new OldAndBustedGenotypeLikelihoods();
        OldAndBustedGenotypeLikelihoods glForward = new OldAndBustedGenotypeLikelihoods();
        OldAndBustedGenotypeLikelihoods glReverse = new OldAndBustedGenotypeLikelihoods();
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            // TODO: could this be done better elsewhere?
            if (read.getNotPrimaryAlignmentFlag() || read.getDuplicateReadFlag() || read.getReadUnmappedFlag()) {
                continue;
            }

            int offset = offsets.get(i);
            char base = read.getReadString().charAt(offset);
            byte qual = read.getBaseQualities()[offset];

            if (qual == 0) { continue; }
            
            if (qualityScoreCutoff != null && qual < qualityScoreCutoff) { continue; }

            if (base == 'A') { aCount++; }
            if (base == 'C') { cCount++; }
            if (base == 'G') { gCount++; }
            if (base == 'T') { tCount++; }

            // System.out.println(read.getReadName() + " " + base + " " + qual );
            
            glAll.add(ref.getBase(), base, qual);

            if (read.getReadNegativeStrandFlag()) {
                glReverse.add(ref.getBase(), base, qual);
            } else {
                glForward.add(ref.getBase(), base, qual);
            }
        }

        // apply priors
        String priorType = "NONE";

        double[] priors;

        // hapmap overrides dbsnp which overrides "novel"
        if (hapmap != null) {
            List<String> knownAlleles = hapmap.getAllelesFWD();
            priorType = "HAPMAP";
            rodString = "[HAPMAP" + hapmap.toSimpleString() + "]";

            priors = getKnownSiteKnownFreqPriors(((byte)(upRef & 0xff)), knownAlleles, hapmap.getVarAlleleFreq());
        } else if (dbsnpInfo != null && dbsnpInfo.isSNP()) {
            List<String> knownAlleles = dbsnpInfo.getAllelesFWD();
            priorType = "DBSNP";
            rodString = "[DBSNP: " + dbsnpInfo.toMediumString() + "]";

            priors = getKnownSitePriors(((byte)(upRef & 0xff)), knownAlleles);
        } else {
            priors = getNovelSitePriors(((byte)(upRef & 0xff)));
        }



        StringBuilder priorString = new StringBuilder();
        priorString.append(priorType);
        for(int i=0; i<priors.length; i++) {

            // apply the prior!!!
            double log10prior = log10(priors[i]);
            glAll.likelihoods[i] += log10prior;
            glForward.likelihoods[i] += log10prior;
            glReverse.likelihoods[i] += log10prior;

            priorString.append(" ");
            priorString
                    .append(glAll.genotypes[i])
                    .append(":")
                    .append(new Formatter().format("%4.2f",log10prior));

        }


        // FIXME: horrible way to do this!!!!
        glAll.doSort();
        glForward.doSort();
        glReverse.doSort();


        String gt = glAll.genotypes[glAll.genotypes.length-1];
        double btnb = glAll.likelihoods[glAll.genotypes.length-1] - glAll.likelihoods[glAll.genotypes.length-2];
        double btr = 1.0d; // FIXME: do we need this?

        // just to mimic old behavior
        String type = "homozygous";
        if (gt.charAt(0) != gt.charAt(1)) { type = "heterozygous-SNP"; }
        else if (gt.charAt(0) != upRef) { type = "homozygous-SNP"; }


        out.print(String.format("%s %d %s %s %4.6f %4.6f %s A:%d C:%d G:%d T:%d %d ",
                      context.getContig(),
                      context.getPosition(),
                      upRef,
                      gt,
                      btnb,
                      btr,
                      type,
                      aCount,
                      cCount,
                      gCount,
                      tCount,
                      aCount + cCount + gCount + tCount));

        out.print(dumpTheories(glAll));

        out.print(String.format("%s %s ",
                      priorString,
                      rodString));

        out.print(dumpTheories(glForward));
        out.print(dumpTheories(glReverse));
        out.print("\n");

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() {
        return 0;
    }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
    }

    protected static final double NOVEL_HET_PRIOR = 2*1E-4;
    protected static final double NOVEL_OTHER_PRIOR = 1E-5;
    protected static final double NOVEL_REF_PRIOR = 1d - 4*NOVEL_HET_PRIOR - 6*NOVEL_OTHER_PRIOR;

    protected static final double KNOWN_HET_PRIOR = .05;
    protected static final double KNOWN_HNREF_PRIOR = .05*.05;
    protected static final double KNOWN_OTHER_PRIOR = 1E-5;
    protected static final double KNOWN_REF_PRIOR = 1d - KNOWN_HET_PRIOR - KNOWN_HNREF_PRIOR - 8* KNOWN_OTHER_PRIOR;

    protected static final double FREQ_OTHER_PRIOR = 1E-5;

    public double[] getNovelSitePriors(final byte referenceBase) {
        //TODO: convert everything to use DiploidGenotype
        DiploidGenotype[] theories = DiploidGenotype.values();
        double[] priors = new double[theories.length];

        DiploidGenotype refTheory = DiploidGenotype.fromBases(referenceBase, referenceBase);

        // use novel variant prior
        // ref (e.g. AA) -> 1 - 4*(2*10^-4) - 6*(10^-5)
        // het (e.g. AC, AG, AT) -> 2*10^-4 (AC, AG, AT,
        // other (e.g. CC, CG, CT, GG, GT, TT) -> 10^-5
        int i=0;
        for (DiploidGenotype theory : theories) {
            double prior = 0;
            if (theory == refTheory) {
                prior = NOVEL_REF_PRIOR;
            } else if (theory.getAllele1() == referenceBase || theory.getAllele2() == referenceBase) {
                prior = NOVEL_HET_PRIOR;
            } else { // not reference, and not het with reference base
                prior = NOVEL_OTHER_PRIOR;
            }
            priors[i++] = prior;
        }

        return priors;
    }

    public double[] getKnownSitePriors(final byte referenceBase,
                                       final List<String> alleles) {

        //TODO: convert everything to use DiploidGenotype
        DiploidGenotype[] theories = DiploidGenotype.values();
        double[] priors = new double[theories.length];

        DiploidGenotype refTheory = DiploidGenotype.fromBases(referenceBase, referenceBase);

        if (alleles != null && alleles.size() > 2 ) {
            // FIXME: handle multi-allele case!
            return getNovelSitePriors(referenceBase);
        }


        StringBuilder sb = new StringBuilder();
        if (alleles != null ) {
            for(String allele : alleles) { sb.append(allele); }
        }
        String allelesString = sb.toString();


        DiploidGenotype knownAlleles = DiploidGenotype.fromBases(allelesString.getBytes());
        byte altAllele = getAlternateAllele(knownAlleles, referenceBase);
        DiploidGenotype hnref = DiploidGenotype.fromBases(altAllele, altAllele);


        // do we not have a freq (in dbsnp, known alleles)
        // AA -
        // AG - .05
        // GG - .05^2
        // other - 10^-5
        int i=0;
        for (DiploidGenotype theory : theories) {
            double prior = 0;
            if (theory == refTheory) {
                prior = KNOWN_REF_PRIOR;
            } else if (theory == knownAlleles) {
                prior = KNOWN_HET_PRIOR;
            } else if (theory == hnref) {
                prior = KNOWN_HNREF_PRIOR;
            } else {
                prior = NOVEL_OTHER_PRIOR;
            }
            priors[i++] = prior;
        }

        return priors;

    }

    public double[] getKnownSiteKnownFreqPriors(final byte referenceBase,
                              final List<String> alleles,
                              final Double altAllelefrequency) {

        //TODO: convert everything to use DiploidGenotype
        DiploidGenotype[] theories = DiploidGenotype.values();
        double[] priors = new double[theories.length];

        DiploidGenotype refTheory = DiploidGenotype.fromBases(referenceBase, referenceBase);

        // if we're dealing with more than one allele here... make sure they all have no frequency data
        if (alleles != null && alleles.size() > 2 ) {
            throw new PicardException("Can't handle multiple known alleles with a single known frequency!");
        }


        StringBuilder sb = new StringBuilder();
        if (alleles != null ) {
            for(String allele : alleles) { sb.append(allele); }
        }
        String allelesString = sb.toString();


//            System.out.println("FREQ");
        // assert that we have just one set of frequencies...
//FIXME:            if (aaf.size() > 1) { throw new PicardException("Attempting to call with multiple allele frequencies!"); }
        DiploidGenotype knownAlleles = DiploidGenotype.fromBases(allelesString.getBytes());
        byte altAllele = getAlternateAllele(knownAlleles, referenceBase);
        DiploidGenotype hnref = DiploidGenotype.fromBases(altAllele, altAllele);

        double q = altAllelefrequency;
        double p = 1d - q;


        // do we have a freq (in hapmap, known alleles)
        // AA -
        // AG -
        // GG -
        // other - 10^-5
        int i=0;
        for (DiploidGenotype theory : theories) {
            double prior = 0;
            if (theory == refTheory) {
                prior = p*p;                            //TODO: should this really sum to 1?
            } else if (theory == knownAlleles) {
                prior = 2*p*q;
            } else if (theory == hnref) {
                prior = q*q;
            } else {
                prior = FREQ_OTHER_PRIOR;
            }
            priors[i++] = prior;
        }


        return priors;
    }

    protected byte getAlternateAllele(final DiploidGenotype gt, final byte refAllele) {
        return (gt.getAllele1() == refAllele)?gt.getAllele2():gt.getAllele1();
    }

    protected String dumpTheories(OldAndBustedGenotypeLikelihoods gl) {
        StringBuilder sb = new StringBuilder();
        for (int i = gl.genotypes.length-1; i >= 0; i--)
        {
            sb.append(String.format("%s:%4.2f ", gl.genotypes[i], gl.likelihoods[i]));
        }
        return sb.toString();
    }


}
