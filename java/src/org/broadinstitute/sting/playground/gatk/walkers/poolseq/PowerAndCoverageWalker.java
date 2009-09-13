package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.playground.utils.SQuad;
import org.broadinstitute.sting.playground.utils.ReadOffsetQuad;
import org.broadinstitute.sting.playground.utils.PoolUtils;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import net.sf.samtools.SAMRecord;

import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Sep 12, 2009
 * Time: 12:20:42 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
public class PowerAndCoverageWalker extends LocusWalker<SQuad<Integer>, SQuad<Long>> {

    @Argument(fullName="suppressLocusPrinting", doc="Walker does not output locus data to any file. Only total coverage information is given at the end", required=false)
    public boolean suppressPrinting = false;

    @Argument(fullName="outputFile", shortName="of", doc="Locus information is printed to this file instead of standard out", required = false)
    public String outputFile = null;

    @Argument(fullName="bootStrapIterations", shortName="bs", doc="Use a bootstrap method with this many iterations", required=false)
    public int bootstrapIterations = 0;

    @Argument(fullName="lodThreshold", shortName="lod", doc="Threshold for log likelihood ratio to be called a SNP. Defaults to 3.0", required = false)
    public double lodThresh = 3.0;

    @Argument(fullName="minimumQScore", shortName="qm", doc="Use bases whose phred (Q) score meets or exceeds this number. Defaults to 0", required = false)
    public byte minQ = 0;

    @Argument(fullName="aboveQScoreOutputOnly", doc="Output only the information for reads that exceed the q-score threshold", required=false)
    public boolean aboveQScoreOutputOnly = false;

    @Argument(fullName="outputRawStatistics", shortName="rs", doc="Walker outputs the median quality score and power for the un-thresholded reads as well as those that are thresholded", required=false)
    public boolean outputRawStatistics = false;

    @Argument(fullName="poolSize", shortName="ps", doc="Number of individuals in the pool", required = true)
    public int numIndividuals = 0;

    protected PrintStream outputWriter = null;

    public void initialize() {
        if(numIndividuals <= 0) {
            throw new StingException("Pool size must be greater than 1. You input "+numIndividuals);
        }
        if( outputFile != null ) {
            try {
                outputWriter = new PrintStream(outputFile);
            } catch (FileNotFoundException e) {
                String errMsg = "The filepath input as an argument to -of, "+outputFile+" was a directory or could not be found or created.";
                throw new StingException(errMsg, e);
            }
                outputWriter.printf("%s%n", generateHeaderString());
        } else if( ! suppressPrinting ) {
            out.printf("%s%n", generateHeaderString());
        }
    }

    public SQuad<Long> reduceInit() {
        return new SQuad<Long>(0l,0l,0l,0l);
    }

    /* REDUCE
     * Add up filtered and unfiltered coverage in the forward
     * and reverse directions.
     */
    public SQuad<Long> reduce( SQuad<Integer> newCvg, SQuad<Long> oldCvg) {
        return new SQuad<Long>(newCvg.getFirst() + oldCvg.getFirst(), newCvg.getSecond() + oldCvg.getSecond(),
                               newCvg.getThird() + oldCvg.getThird(), newCvg.getFourth() + oldCvg.getFourth());
    }
    /* MAP
     * Get the coverage in each direction above and below threshold and calculate the estimated power
     * if printing is not suppressed print this information to the file
     * pass coverage information on to Reduce
     */
    public SQuad<Integer> map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ReadOffsetQuad readsByDirection = splitReadsByDirection(context.getReads(),context.getOffsets());
        ReadOffsetQuad readsByDirectionAboveThresh;

        if(minQ > 0) {
            readsByDirectionAboveThresh = thresholdReadsByQuality(readsByDirection);
        } else {
            readsByDirectionAboveThresh = readsByDirection;
        }


        SQuad<Byte> medianQScoresAboveThresh = getMedianQScores(readsByDirectionAboveThresh);
        SQuad<Double> powerAboveThresh = calculatePower(readsByDirectionAboveThresh, medianQScoresAboveThresh);
        SQuad<Byte> medianQScoresUnthresh;
        SQuad<Double> powerUnthresh;

        if(! outputRawStatistics) {
            medianQScoresUnthresh = null;
            powerUnthresh = null;
        } else {
            medianQScoresUnthresh = getMedianQScores(readsByDirection);
            powerUnthresh= calculatePower(readsByDirection, medianQScoresUnthresh);
        }

        if ( ! suppressPrinting && outputWriter == null) {
            out.printf("%s%n",outputString(readsByDirection,readsByDirectionAboveThresh,powerAboveThresh, powerUnthresh,
                              medianQScoresAboveThresh, medianQScoresUnthresh, context.getLocation()));
        } else if (! suppressPrinting && outputWriter != null) {
            outputWriter.printf("%s%n",outputString(readsByDirection,readsByDirectionAboveThresh,powerAboveThresh, powerUnthresh,
                                       medianQScoresAboveThresh, medianQScoresUnthresh, context.getLocation()));
        }

        return new SQuad<Integer>(readsByDirection.numReadsFirst(), readsByDirection.numReadsSecond(),
                                  readsByDirectionAboveThresh.numReadsFirst(), readsByDirectionAboveThresh.numReadsSecond());
        
    }

    // private helper methods

    /*
     * Make different headers for the output file depending on printing-related command line arguments
     */
    private String generateHeaderString() {
        String header;
        if (minQ <= 0) {
            header = "Chrom:Pos  CoverageFwd  CoverageRev  CoverageCom  MedianQFwd  MedianQRev  MedianQCom  PowerFwd  PowerRev  PowerCom";
        } else if (! aboveQScoreOutputOnly && ! outputRawStatistics){
            String m = Byte.toString(minQ);
            header = "Chrom:Pos  CoverageFwd  CoverageRev  CoverageCom  CoverageAboveQ"+m+"Fwd  CoverageAboveQ"+m+"Rev  CoverageAboveQ"+m+"Com  MedianQFwd  MedianQRev  MedianQCom  PowerFwd  PowerRev  PowerCom";
        } else if (! outputRawStatistics){
            header = "Chrom:Pos  CoverageFwd  CoverageRev  CoverageCom  MedianQFwd  MedianQRev  MedianQCom  PowerFwd  PowerRev  PowerCom";
        } else {
            String m = Byte.toString(minQ);
            header = "Chrom:Pos  CoverageFwd  CoverageRev  CoverageCom  CoverageAboveQ"+m+"Fwd  CoverageAboveQ"+m+"Rev  "+
                     "CoverageAboveQ"+m+"Com  RawMedianQFwd  RawMedianQRev  RawMedianQCom  MedianQAboveQ"+m+"Fwd  "+
                     "MedianQAboveQ"+m+"Rev  MedianQAboveQ"+m+"Com  RawPowerFwd  RawPowerRev  RawPowerCom  "+
                     "PowerAboveQ"+m+"Fwd  PowerAboveQ"+m+"Rev  PowerAboveQ"+m+"Com";
        }

        return header;
    }
    /*
     * Make different output strings depending on printing-related command-line options
     */
    private String outputString(ReadOffsetQuad reads, ReadOffsetQuad threshReads, SQuad<Double> pow, SQuad<Double> rawPow, SQuad<Byte> medQ, SQuad<Byte> rawQ, GenomeLoc loc) {
        String outputString;
        if(minQ <= 0) {
            outputString = String.format("%s %d %d %d %d %d %d %f %f %f", loc.toString(), reads.numReadsFirst(), reads.numReadsSecond(), reads.numReads(),
                                          medQ.getFirst(), medQ.getSecond(), medQ.getThird(), pow.getFirst(), pow.getSecond(), pow.getThird());
        } else if (! aboveQScoreOutputOnly && ! outputRawStatistics ){
            outputString = String.format("%s %d %d %d %d %d %d %d %d %d %f %f %f", loc.toString(), reads.numReadsFirst(), reads.numReadsSecond(), reads.numReads(),
                                         threshReads.numReadsFirst(), threshReads.numReadsSecond(), threshReads.numReads(), medQ.getFirst(),
                                         medQ.getSecond(), medQ.getThird(), pow.getFirst(), pow.getSecond(), pow.getThird());
        } else if (! outputRawStatistics ){
            outputString = String.format("%s %d %d %d %d %d %d %f %f %f", loc.toString(), threshReads.numReadsFirst(), threshReads.numReadsSecond(),
                                          threshReads.numReads(), medQ.getFirst(), medQ.getSecond(), medQ.getThird(), pow.getFirst(),
                                          pow.getSecond(), pow.getThird());
        } else {
            outputString = String.format("%s %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f", loc.toString(), reads.numReadsFirst(),
                                          reads.numReadsSecond(), reads.numReads(), threshReads.numReadsFirst(), threshReads.numReadsSecond(),
                                          threshReads.numReads(), rawQ.getFirst(), rawQ.getSecond(), rawQ.getThird(), medQ.getFirst(),
                                          medQ.getSecond(), medQ.getThird(), rawPow.getFirst(), rawPow.getSecond(), rawPow.getThird(),
                                          pow.getFirst(), pow.getSecond(), pow.getThird());
        }
        return outputString;
    }

    private ReadOffsetQuad splitReadsByDirection(List<SAMRecord> reads, List<Integer> offsets) {
        return PoolUtils.splitReadsByReadDirection(reads,offsets);
    }

    private ReadOffsetQuad thresholdReadsByQuality(ReadOffsetQuad readQuad) {
        return new ReadOffsetQuad(PoolUtils.thresholdReadsByQuality(readQuad.getFirstReadOffsetPair(), minQ),
                                  PoolUtils.thresholdReadsByQuality(readQuad.getSecondReadOffsetPair(), minQ));
    }

    private SQuad<Byte> getMedianQScores(ReadOffsetQuad reads) {
        return new SQuad<Byte>(ListUtils.getQScoreMedian(reads.getFirstReads(), reads.getFirstOffsets()),
                               ListUtils.getQScoreMedian(reads.getSecondReads(), reads.getSecondOffsets()),
                               ListUtils.getQScoreMedian(reads.getReadsCombined(), reads.getOffsetsCombined()),
                               (byte) 0);
    }

    private SQuad<Double> calculatePower(ReadOffsetQuad readQuad, SQuad<Byte> medianQs) {
        SQuad<Double> power;
        if( bootstrapIterations > 0 ) {
            power = bootstrapPowerByDirection(readQuad);
        } else {
            power = theoreticalPowerByDirection(readQuad,medianQs);
        }
        return power;
    }

    private SQuad<Double> theoreticalPowerByDirection(ReadOffsetQuad readQuad, SQuad<Byte> medianQs) {
        return new SQuad<Double>( theoreticalPower(readQuad.numReadsFirst(), medianQs.getFirst()),
                                  theoreticalPower(readQuad.numReadsSecond(), medianQs.getSecond()),
                                  theoreticalPower(readQuad.numReads(), medianQs.getThird()),
                                  0.0 );
    }

    private SQuad<Double> bootstrapPowerByDirection(ReadOffsetQuad readQuad) {
        return new SQuad<Double> ( bootstrapPower(readQuad.getFirstReads(), readQuad.getFirstOffsets()),
                                   bootstrapPower(readQuad.getSecondReads(), readQuad.getSecondOffsets()),
                                   bootstrapPower(readQuad.getReadsCombined(), readQuad.getOffsetsCombined()),
                                   0.0 );
    }

    private double theoreticalPower(int depth, byte q) {
        double power;
        if( depth <= 0 ) {
            power = 0.0;
        } else {
            double p_error = QualityUtils.qualToErrorProb(q);
            double snpProp = getSNPProportion();
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

    private double bootstrapPower(List<SAMRecord> reads, List<Integer> offsets) {
        double power = 0.0;
        double snpProp = getSNPProportion();
        for ( int boot = 0; boot < bootstrapIterations; boot ++) {
            ReadOffsetQuad partition = PoolUtils.coinTossPartition(reads, offsets, getSNPProportion());
            if ( PoolUtils.calculateLogLikelihoodOfSample(partition, snpProp) >= lodThresh ) {
                power += 1.0/bootstrapIterations;
            }
        }
        return power;
    }

    private double getSNPProportion() {
        return 1/(2.0*numIndividuals);
    }

}
