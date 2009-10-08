package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ListUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.playground.utils.PoolUtils;
import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Oct 8, 2009
 * Time: 9:44:35 AM
 * To change this template use File | Settings | File Templates.
 */
public class PowerBelowFrequencyWalker extends LocusWalker<Integer,Integer> {
    @Argument(fullName="lodThreshold", shortName="lod", doc="Threshold for log likelihood ratio to be called a SNP. Defaults to 3.0", required = false)
    public double lodThresh = 3.0;

    @Argument(fullName="minimumQScore", shortName="qm", doc="Use bases whose phred (Q) score meets or exceeds this number. Defaults to 0", required = false)
    public byte minQ = 0;

    @Argument(fullName="poolSize", shortName="ps", doc="Number of individuals in the pool", required = true)
    public int numIndividuals = 0;

    @Argument(fullName="alleleFrequency", shortName="af", doc="Calculate power for all allele frequencies below this. Defaults to 4", required=false)
    public int alleleFreq = 4;

    @Argument(fullName="useMeanProb", doc="Use the mean probability as the \"average quality\" rather than median Q-score")
    boolean useMean = false;

    public void initialize() {
        if ( alleleFreq < 1 ) {
            String err = "Allele frequency (-af) must be greater than or equal to one.";
            throw new StingException(err);
        }
    }

    public Integer reduceInit() {
        out.print(makeHeader());
        return 0;
    }

    public Integer reduce(Integer mapint, Integer prevint) {
        return 0;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String output = String.format("%s", context.getLocation().toString());

        // threshold reads if necessary

        if ( minQ > 0 ) {
            Pair<List<SAMRecord>, List<Integer>> thresh = PoolUtils.thresholdReadsByQuality(context.getReads(),context.getOffsets(),minQ);
            context = new AlignmentContext(context.getLocation(), thresh.getFirst(), thresh.getSecond());
        }

        // calculate powers and put into output string

        for ( int i = 1; i <= alleleFreq; i ++ ) {
            output = String.format("%s\t%f",output,calculatePowerAtFrequency(context,i));
        }

        // print the output string

        out.printf("%s%n",output);

        return 0;
    }

    public double calculatePowerAtFrequency( AlignmentContext context, int alleles ) {
        return theoreticalPower( context.numReads(), getMeanQ(context), alleles );
    }

    public byte getMeanQ( AlignmentContext context ) {
        byte meanQ;
        if ( useMean ) {
            meanQ = QualityUtils.probToQual(expectedMatchRate(context));
        } else {
            meanQ = ListUtils.getQScoreMedian(context.getReads(),context.getOffsets());
        }

        return meanQ;
    }

    public double expectedMatchRate(AlignmentContext context) {
        int nReads = context.numReads();
        double matches = 0.0;
        for ( int r = 0; r < nReads; r ++ ) {
            matches += QualityUtils.qualToProb(context.getReads().get(r).getBaseQualities()[context.getOffsets().get(r)]);
        }

        return matches/nReads;
    }

    public String makeHeader() {
        // create the header
        String header = "chrm:pos";
        for ( int i = 1; i <= alleleFreq; i ++ ) {
            header = header + "\tPower_at_"+Integer.toString(i);
        }

        return String.format("%s%n", header);
    }

    public double theoreticalPower( int depth, byte q, int alleles ) {
        double power;
        if( depth <= 0 ) {
            power = 0.0;
        } else {
            double p_error = QualityUtils.qualToErrorProb(q);
            double snpProp = getSNPProportion(alleles);
            double kterm_num = Math.log10( snpProp * (1 - p_error) + (1 - snpProp) * (p_error/3) );
            double kterm_denom = Math.log10( p_error/3 );
            double dkterm_num = Math.log10( snpProp * (p_error/3) + (1 - snpProp) * (1 - p_error) );
            double dkterm_denom = Math.log10( 1 - p_error);

            int kaccept = (int) Math.ceil( ( lodThresh - ( (double) depth ) * ( dkterm_num - dkterm_denom ) ) /
                    ( kterm_num - kterm_denom- dkterm_num + dkterm_denom ) );

            // we will reject the null hypothesis if we see kaccept or more SNPs, the power is the probability that this occurs
            // we can optimize this by checking to see which sum is smaller

            if ( depth - kaccept < kaccept ) {// kaccept > depth/2 - calculate power as P[hits between kaccept and depth]
                power = MathUtils.cumBinomialProbLog(kaccept, depth, depth, snpProp);
            } else { // kaccept < depth/2 - calculate power as 1-P[hits between 0 and kaccept]
                power = 1-MathUtils.cumBinomialProbLog(0,kaccept,depth,snpProp);
            }
        }

        return power;
    }

    private double getSNPProportion(int alleles) {
        return ((double)alleles)/numIndividuals;
    }


}
