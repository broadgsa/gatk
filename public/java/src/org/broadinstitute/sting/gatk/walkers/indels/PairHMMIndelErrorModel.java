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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

/*import org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates.Covariate;
import org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates.RecalDataManager;
import org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates.RecalDatum;
import org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates.RecalibrationArgumentCollection;
*/


public class PairHMMIndelErrorModel {


    public static final int BASE_QUAL_THRESHOLD = 20;


    private static final int MATCH_OFFSET = 0;
    private static final int X_OFFSET = 1;
    private static final int Y_OFFSET = 2;

    private static final int DIAG = 0;
    private static final int UP = 1;
    private static final int LEFT = 2;

    private static final int DIAG_GOTO_M = 0;
    private static final int DIAG_GOTO_X = 1;
    private static final int DIAG_GOTO_Y = 2;

    private static final int UP_GOTO_M = 4;
    private static final int UP_GOTO_X = 5;
    private static final int UP_GOTO_Y = 6;

    private static final int LEFT_GOTO_M = 8;
    private static final int LEFT_GOTO_X = 9;
    private static final int LEFT_GOTO_Y = 10;

    private static final int[] ACTIONS_M = {DIAG_GOTO_M, DIAG_GOTO_X, DIAG_GOTO_Y};
    private static final int[] ACTIONS_X = {UP_GOTO_M, UP_GOTO_X, UP_GOTO_Y};
    private static final int[] ACTIONS_Y = {LEFT_GOTO_M, LEFT_GOTO_X, LEFT_GOTO_Y};


    private final double logGapOpenProbability;
    private final double logGapContinuationProbability;

    private boolean DEBUG = false;

    private static final int MAX_CACHED_QUAL = 127;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    private final static double LOG_ONE_HALF;
    private final static double END_GAP_COST;

    private static final int START_HRUN_GAP_IDX = 4;
    private static final int MAX_HRUN_GAP_IDX = 20;

    private static final double MIN_GAP_OPEN_PENALTY = 30.0;
    private static final double MIN_GAP_CONT_PENALTY = 10.0;
    private static final double GAP_PENALTY_HRUN_STEP = 1.0; // each increase in hrun decreases gap penalty by this.


    private boolean doViterbi = false;

    private final boolean useAffineGapModel = true;
    private boolean doContextDependentPenalties = false;

    private final double[] GAP_OPEN_PROB_TABLE;
    private final double[] GAP_CONT_PROB_TABLE;

    private boolean getGapPenaltiesFromFile = false;

    private int SMOOTHING = 1;
    private int MAX_QUALITY_SCORE = 50;
    private int PRESERVE_QSCORES_LESS_THAN = 5;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
//copy+
/*    private RecalDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>(); // List of covariates to be used in this calculation
    private static final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private static final Pattern OLD_RECALIBRATOR_HEADER = Pattern.compile("^rg,.*");
    private static final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");
    protected static final String EOF_MARKER = "EOF";
    private long numReadsWithMalformedColorSpace = 0;
    private RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
    private NestedHashMap qualityScoreByFullCovariateKey = new NestedHashMap(); // Caches the result of performSequentialQualityCalculation(..) for all sets of covariate values.
  */
//copy-
    static {
        LOG_ONE_HALF= -Math.log10(2.0);
        END_GAP_COST = LOG_ONE_HALF;

        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10, -k/10.);


            baseMatchArray[k] =  Math.log10(1-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public  PairHMMIndelErrorModel(double indelGOP, double indelGCP, boolean deb, boolean doCDP, boolean dovit,boolean gpf, File RECAL_FILE) {

        this(indelGOP, indelGCP, deb, doCDP, dovit);
        this.getGapPenaltiesFromFile = gpf;

        // read data from recal file
        // gdebug - start copy from TableRecalibrationWalker
/*        if (gpf) {
            boolean sawEOF = false;
            boolean REQUIRE_EOF = false;

            int lineNumber = 0;
            boolean foundAllCovariates = false;
            // Get a list of all available covariates
            final List<Class<? extends Covariate>> classes = new PluginManager<Covariate>(Covariate.class).getPlugins();

            try {
                for ( String line : new XReadLines(RECAL_FILE) ) {
                    lineNumber++;
                    if ( EOF_MARKER.equals(line) ) {
                        sawEOF = true;
                    } else if( COMMENT_PATTERN.matcher(line).matches() || OLD_RECALIBRATOR_HEADER.matcher(line).matches() )  {
                        ; // Skip over the comment lines, (which start with '#')
                    }
                    // Read in the covariates that were used from the input file
                    else if( COVARIATE_PATTERN.matcher(line).matches() ) { // The line string is either specifying a covariate or is giving csv data
                        if( foundAllCovariates ) {
                            throw new UserException.MalformedFile( RECAL_FILE, "Malformed input recalibration file. Found covariate names intermingled with data in file: " + RECAL_FILE );
                        } else { // Found the covariate list in input file, loop through all of them and instantiate them
                            String[] vals = line.split(",");
                            for( int iii = 0; iii < vals.length - 3; iii++ ) { // There are n-3 covariates. The last three items are nObservations, nMismatch, and Qempirical
                                boolean foundClass = false;
                                for( Class<?> covClass : classes ) {
                                    if( (vals[iii] + "Covariate").equalsIgnoreCase( covClass.getSimpleName() ) ) {
                                        foundClass = true;
                                        try {
                                            Covariate covariate = (Covariate)covClass.newInstance();
                                            requestedCovariates.add( covariate );
                                        } catch (Exception e) {
                                            throw new DynamicClassResolutionException(covClass, e);
                                        }

                                    }
                                }

                                if( !foundClass ) {
                                    throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration file. The requested covariate type (" + (vals[iii] + "Covariate") + ") isn't a valid covariate option." );
                                }
                            }
                        }

                    } else { // Found a line of data
                        if( !foundAllCovariates ) {
                            foundAllCovariates = true;

                            // At this point all the covariates should have been found and initialized
                            if( requestedCovariates.size() < 2 ) {
                                throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Covariate names can't be found in file: " + RECAL_FILE );
                            }

                            final boolean createCollapsedTables = true;

                            // Initialize any covariate member variables using the shared argument collection
                            for( Covariate cov : requestedCovariates ) {
                                cov.initialize( RAC );
                            }
                            // Initialize the data hashMaps
                            dataManager = new RecalDataManager( createCollapsedTables, requestedCovariates.size() );

                        }
                        addCSVData(RECAL_FILE, line); // Parse the line and add the data to the HashMap
                    }
                }

            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotReadInputFile(RECAL_FILE, "Can not find input file", e);
            } catch ( NumberFormatException e ) {
                throw new UserException.MalformedFile(RECAL_FILE, "Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
            }

            if ( !sawEOF ) {
                final String errorMessage = "No EOF marker was present in the recal covariates table; this could mean that the file is corrupted or was generated with an old version of the CountCovariates tool.";
                if ( REQUIRE_EOF )
                    throw new UserException.MalformedFile(RECAL_FILE, errorMessage);
            }

            if( dataManager == null ) {
                throw new UserException.MalformedFile(RECAL_FILE, "Can't initialize the data manager. Perhaps the recal csv file contains no data?");
            }

            // Create the tables of empirical quality scores that will be used in the sequential calculation
            dataManager.generateEmpiricalQualities( SMOOTHING, MAX_QUALITY_SCORE );
        }
        // debug end copy
  */
    }
    /**
     * For each covariate read in a value and parse it. Associate those values with the data itself (num observation and num mismatches)
     */
 /*
    private void addCSVData(final File file, final String line) {
        final String[] vals = line.split(",");

        // Check if the data line is malformed, for example if the read group string contains a comma then it won't be parsed correctly
        if( vals.length != requestedCovariates.size() + 3 ) { // +3 because of nObservations, nMismatch, and Qempirical
            throw new UserException.MalformedFile(file, "Malformed input recalibration file. Found data line with too many fields: " + line +
                    " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
        }

        final Object[] key = new Object[requestedCovariates.size()];
        Covariate cov;
        int iii;
        for( iii = 0; iii < requestedCovariates.size(); iii++ ) {
            cov = requestedCovariates.get( iii );
            key[iii] = cov.getValue( vals[iii] );
        }

        // Create a new datum using the number of observations, number of mismatches, and reported quality score
        final RecalDatum datum = new RecalDatum( Long.parseLong( vals[iii] ), Long.parseLong( vals[iii + 1] ), Double.parseDouble( vals[1] ), 0.0 );
        // Add that datum to all the collapsed tables which will be used in the sequential calculation
        dataManager.addToAllTables( key, datum, PRESERVE_QSCORES_LESS_THAN );
    }

*/
    public  PairHMMIndelErrorModel(double indelGOP, double indelGCP, boolean deb, boolean doCDP, boolean dovit) {
        this(indelGOP, indelGCP, deb, doCDP);
        this.doViterbi = dovit;
    }

    public PairHMMIndelErrorModel(double indelGOP, double indelGCP, boolean deb, boolean doCDP) {


        this.logGapOpenProbability = -indelGOP/10.0; // QUAL to log prob
        this.logGapContinuationProbability = -indelGCP/10.0; // QUAL to log prob
        this.doContextDependentPenalties = doCDP;
        this.DEBUG = deb;


        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new double[MAX_HRUN_GAP_IDX];

        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = logGapOpenProbability;
            GAP_CONT_PROB_TABLE[i] = logGapContinuationProbability;
        }

        double gop = logGapOpenProbability;
        double gcp = logGapContinuationProbability;
        double step = GAP_PENALTY_HRUN_STEP/10.0;

        double maxGOP = -MIN_GAP_OPEN_PENALTY/10.0;  // phred to log prob
        double maxGCP = -MIN_GAP_CONT_PENALTY/10.0;  // phred to log prob

        for (int i=START_HRUN_GAP_IDX; i < MAX_HRUN_GAP_IDX; i++) {
            gop += step;
            if (gop > maxGOP)
                gop = maxGOP;

            gcp += step;
            if(gcp > maxGCP)
                gcp = maxGCP;
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

    }

    private double computeReadLikelihoodGivenHaplotype(byte[] haplotypeBases, byte[] readBases, byte[] readQuals) {
        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        // initialize path metric and traceback memories for likelihood computation
        double[][] pathMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestMetricArray = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        pathMetricArray[0][0]= 0;//Double.NEGATIVE_INFINITY;

        for (int i=1; i < X_METRIC_LENGTH; i++) {
            pathMetricArray[i][0] = 0;
            bestMetricArray[i][0] = UP;
        }

        for (int j=1; j < Y_METRIC_LENGTH; j++) {
            pathMetricArray[0][j] = 0;//logGapOpenProbability + (j-1) * logGapContinuationProbability;
            bestMetricArray[0][j] = LEFT;
        }

        for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
            for (int indJ=1; indJ < Y_METRIC_LENGTH; indJ++) {

                byte x = readBases[indI-1];
                byte y = haplotypeBases[indJ-1];
                byte qual = readQuals[indI-1];

                double bestMetric = 0.0;
                int bestMetricIdx = 0;

                // compute metric for match/mismatch
                // workaround for reads whose bases quality = 0,
                if (qual < 1)
                    qual = 1;

                if (qual > MAX_CACHED_QUAL)
                    qual = MAX_CACHED_QUAL;

                double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];
                double[] metrics = new double[3];

                metrics[DIAG] = pathMetricArray[indI-1][indJ-1] + pBaseRead;
                metrics[UP] = pathMetricArray[indI-1][indJ] + logGapOpenProbability;//(end?0.0:logGapOpenProbability);
                metrics[LEFT] = pathMetricArray[indI][indJ-1] + logGapOpenProbability;//(end?0.0:logGapOpenProbability);

                if (doViterbi) {
                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(metrics);

                pathMetricArray[indI][indJ] = bestMetric;
                bestMetricArray[indI][indJ] = bestMetricIdx;

            }
        }


        double bestMetric=0.0;
        int bestMetricIdx=0,bestI=X_METRIC_LENGTH - 1, bestJ=Y_METRIC_LENGTH - 1;

        for (int i=0; i < X_METRIC_LENGTH; i ++ ) {
            int j= Y_METRIC_LENGTH-1;

            if (pathMetricArray[i][j] > bestMetric) {
                bestMetric = pathMetricArray[i][j];
                bestI = i;
                bestJ = j;
            }
        }
        for (int j=0; j < Y_METRIC_LENGTH; j++ ) {
            int i= X_METRIC_LENGTH-1;
            if (pathMetricArray[i][j] >= bestMetric) {
                bestMetric = pathMetricArray[i][j];
                bestI = i;
                bestJ = j;
            }
        }

        if (DEBUG && doViterbi) {

            String haplotypeString = new String (haplotypeBases);
            String readString = new String(readBases);


            int i = bestI;
            int j = bestJ;


            System.out.println("Simple NW");

            while (i >0 || j >0) {
                bestMetricIdx = bestMetricArray[i][j];
                System.out.print(bestMetricIdx);
                if (bestMetricIdx == UP) {
                    // insert gap in Y
                    haplotypeString = haplotypeString.substring(0,j)+"-"+haplotypeString.substring(j);
                    i--;
                } else if (bestMetricIdx == LEFT) {
                    readString = readString.substring(0,i)+"-"+readString.substring(i);
                    j--;
                }
                else {
                    i--; j--;
                }
            }




            System.out.println("\nAlignment: ");
            System.out.println("R:"+readString);
            System.out.println("H:"+haplotypeString);
            System.out.println();


        }
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);

        return bestMetric;


    }

    static private void getContextHomopolymerLength(final byte[] refBytes, int[] hrunArray) {
        // compute forward hrun length, example:
        // AGGTGACCCCCCTGAGAG
        // 001000012345000000
        int runCount = 0;
        hrunArray[0] = 0;
        int[] hforward = new int[hrunArray.length];
        int[] hreverse = new int[hrunArray.length];

        for (int i = 1; i < refBytes.length; i++) {
            if (refBytes[i] == refBytes[i-1])
                hforward[i] = hforward[i-1]+1;
            else
                hforward[i] = 0;
        }

        // do similar thing for reverse length, example:
        // AGGTGACCCCCCTGAGAG
        // 021000543210000000
        // and then accumulate with forward values.
        // Total:
        // AGGTGACCCCCCTGAGAG
        // 022000555555000000
        for (int i=refBytes.length-1; i > 0; i--) {
            if (refBytes[i-1] == refBytes[i])
                hreverse[i-1] += hreverse[i]+1;
        }

        for (int i = 1; i < refBytes.length; i++)
            hrunArray[i] = hforward[i]+hreverse[i];
    }


    private double computeReadLikelihoodGivenHaplotypeAffineGaps(byte[] haplotypeBases, byte[] readBases, byte[] readQuals,
                                                                 double[] currentGOP, double[] currentGCP) {

        final int X_METRIC_LENGTH = readBases.length+1;
        final int Y_METRIC_LENGTH = haplotypeBases.length+1;

        // initialize path metric and traceback memories for likelihood computation
        double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayM = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayX = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        int[][] bestActionArrayY = new int[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        double c,d;
        matchMetricArray[0][0]= END_GAP_COST;//Double.NEGATIVE_INFINITY;

        for (int i=1; i < X_METRIC_LENGTH; i++) {
            //initialize first column
            matchMetricArray[i][0]  = Double.NEGATIVE_INFINITY;
            YMetricArray[i][0]      = Double.NEGATIVE_INFINITY;
            XMetricArray[i][0]      = END_GAP_COST*(i);//logGapOpenProbability + (i-1)*logGapContinuationProbability;

            bestActionArrayX[i][0] = bestActionArrayY[i][0] = bestActionArrayM[i][0] = UP_GOTO_X;
        }

        for (int j=1; j < Y_METRIC_LENGTH; j++) {
            // initialize first row
            matchMetricArray[0][j]  = Double.NEGATIVE_INFINITY;
            XMetricArray[0][j]      = Double.NEGATIVE_INFINITY;
            YMetricArray[0][j]      = END_GAP_COST*(j);//logGapOpenProbability + (j-1) * logGapContinuationProbability;

            bestActionArrayY[0][j] = bestActionArrayM[0][j] = bestActionArrayX[0][j] = LEFT_GOTO_Y;
        }

        for (int indI=1; indI < X_METRIC_LENGTH; indI++) {
            int im1 = indI-1;
            for (int indJ=1; indJ < Y_METRIC_LENGTH; indJ++) {
                int jm1 = indJ-1;
                byte x = readBases[im1];
                byte y = haplotypeBases[jm1];
                byte qual = readQuals[im1];

                double bestMetric = 0.0;
                int bestMetricIdx = 0;

                // compute metric for match/mismatch
                // workaround for reads whose bases quality = 0,
                if (qual < 1)
                    qual = 1;

                if (qual > MAX_CACHED_QUAL)
                    qual = MAX_CACHED_QUAL;

                double pBaseRead =  (x == y)? baseMatchArray[(int)qual]:baseMismatchArray[(int)qual];


                double[] metrics = new double[3];


                if (doViterbi) {
                    // update match array
                    metrics[MATCH_OFFSET] = matchMetricArray[im1][jm1] + pBaseRead;
                    metrics[X_OFFSET] = XMetricArray[im1][jm1] + pBaseRead;
                    metrics[Y_OFFSET] = YMetricArray[im1][jm1] + pBaseRead;

                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(matchMetricArray[im1][jm1] + pBaseRead, XMetricArray[im1][jm1] + pBaseRead,
                            YMetricArray[im1][jm1] + pBaseRead);

                matchMetricArray[indI][indJ] = bestMetric;
                bestActionArrayM[indI][indJ] = ACTIONS_M[bestMetricIdx];

                // update X array
                // State X(i,j): X(1:i) aligned to a gap in Y(1:j).
                // When in last column of X, ie X(1:i) aligned to full Y, we don't want to penalize gaps

                //c = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: currentGOP[jm1]);
                //d = (indJ==Y_METRIC_LENGTH-1? END_GAP_COST: currentGCP[jm1]);
                if (getGapPenaltiesFromFile) {
                    c = currentGOP[im1];
                    d = logGapContinuationProbability;

                } else {
                    c = currentGOP[jm1];
                    d = currentGCP[jm1];
                }
                if (indJ == Y_METRIC_LENGTH-1)
                    c = d = END_GAP_COST;

                if (doViterbi) {
                    metrics[MATCH_OFFSET] = matchMetricArray[im1][indJ] + c;
                    metrics[X_OFFSET] = XMetricArray[im1][indJ] + d;
                    metrics[Y_OFFSET] = Double.NEGATIVE_INFINITY; //YMetricArray[indI-1][indJ] + logGapOpenProbability;

                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(matchMetricArray[im1][indJ] + c, XMetricArray[im1][indJ] + d);

                XMetricArray[indI][indJ] = bestMetric;
                bestActionArrayX[indI][indJ] = ACTIONS_X[bestMetricIdx];

                // update Y array
                //c = (indI==X_METRIC_LENGTH-1? END_GAP_COST: currentGOP[jm1]);
                //d = (indI==X_METRIC_LENGTH-1? END_GAP_COST: currentGCP[jm1]);
                if (getGapPenaltiesFromFile) {
                    c = currentGOP[im1];
                    d = logGapContinuationProbability;
                }
                else {
                    c = currentGOP[jm1];
                    d = currentGCP[jm1];                        
                }
                if (indI == X_METRIC_LENGTH-1)
                    c = d = END_GAP_COST;



                if (doViterbi) {
                    metrics[MATCH_OFFSET] = matchMetricArray[indI][jm1] + c;
                    metrics[X_OFFSET] = Double.NEGATIVE_INFINITY; //XMetricArray[indI][indJ-1] + logGapOpenProbability;
                    metrics[Y_OFFSET] = YMetricArray[indI][jm1] + d;

                    bestMetricIdx = MathUtils.maxElementIndex(metrics);
                    bestMetric = metrics[bestMetricIdx];
                }
                else
                    bestMetric = MathUtils.softMax(matchMetricArray[indI][jm1] + c, YMetricArray[indI][jm1] + d);

                YMetricArray[indI][indJ] = bestMetric;
                bestActionArrayY[indI][indJ] = ACTIONS_Y[bestMetricIdx];



            }
        }

        double bestMetric;
        double metrics[] = new double[3];
        int bestTable=0, bestI=X_METRIC_LENGTH - 1, bestJ=Y_METRIC_LENGTH - 1;
        metrics[MATCH_OFFSET] = matchMetricArray[bestI][bestJ];
        metrics[X_OFFSET] = XMetricArray[bestI][bestJ];
        metrics[Y_OFFSET] = YMetricArray[bestI][bestJ];
        if (doViterbi) {
            bestTable = MathUtils.maxElementIndex(metrics);
            bestMetric = metrics[bestTable];
        }
        else
            bestMetric = MathUtils.softMax(metrics);

        // Do traceback (needed only for debugging!)
        if (DEBUG && doViterbi) {

            int bestAction;
            int i = bestI;
            int j = bestJ;


            System.out.println("Affine gap NW");


            String haplotypeString = new String (haplotypeBases);
            String readString = new String(readBases);


            while (i >0 || j >0) {
                if (bestTable == X_OFFSET) {
                    // insert gap in Y
                    haplotypeString = haplotypeString.substring(0,j)+"-"+haplotypeString.substring(j);
                    bestAction = bestActionArrayX[i][j];
                }
                else if (bestTable == Y_OFFSET) {
                    readString = readString.substring(0,i)+"-"+readString.substring(i);
                    bestAction = bestActionArrayY[i][j];

                }
                else {
                    bestAction = bestActionArrayM[i][j];
                }
                System.out.print(bestAction);


                // bestAction contains action to take at next step
                // encoding of bestAction: upper 2 bits = direction, lower 2 bits = next table

                // bestTable and nextDirection for next step
                bestTable = bestAction & 0x3;
                int nextDirection = bestAction >> 2;
                if (nextDirection == UP) {
                    i--;
                } else if (nextDirection == LEFT) {
                    j--;
                } else { //  if (nextDirection == DIAG)
                    i--; j--;
                }

            }




            System.out.println("\nAlignment: ");
            System.out.println("R:"+readString);
            System.out.println("H:"+haplotypeString);
            System.out.println();


        }
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);

        return bestMetric;

    }

    private void fillGapProbabilities(int[] hrunProfile,
                                      double[] contextLogGapOpenProbabilities, double[] contextLogGapContinuationProbabilities) {
        // fill based on lookup table
        for (int i = 0; i < hrunProfile.length; i++) {
            if (hrunProfile[i] >= MAX_HRUN_GAP_IDX) {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
            }
            else {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[hrunProfile[i]];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[hrunProfile[i]];
            }
        }
    }
    public synchronized double[] computeReadHaplotypeLikelihoods(ReadBackedPileup pileup, LinkedHashMap<Allele,Haplotype> haplotypeMap,
                                                                   ReferenceContext ref, int eventLength,
                                                                   HashMap<PileupElement, LinkedHashMap<Allele,Double>> indelLikelihoodMap){

        int numHaplotypes = haplotypeMap.size();
        double[][] haplotypeLikehoodMatrix = new double[numHaplotypes][numHaplotypes];
        double readLikelihoods[][] = new double[pileup.getReads().size()][numHaplotypes];
        int readIdx=0;

        LinkedHashMap<Allele,double[]> gapOpenProbabilityMap = new LinkedHashMap<Allele,double[]>();
        LinkedHashMap<Allele,double[]> gapContProbabilityMap = new LinkedHashMap<Allele,double[]>();

        if (DEBUG) {
            System.out.println("Reference bases:");
            System.out.println(new String(ref.getBases()));
        }

        if (doContextDependentPenalties && !getGapPenaltiesFromFile)   {
            // will context dependent probabilities based on homopolymer run. Probabilities are filled based on total complete haplotypes.


            for (Allele a: haplotypeMap.keySet()) {
                Haplotype haplotype = haplotypeMap.get(a);
                byte[] haplotypeBases = haplotype.getBasesAsBytes();
                double[] contextLogGapOpenProbabilities = new double[haplotypeBases.length];
                double[] contextLogGapContinuationProbabilities = new double[haplotypeBases.length];

                // get homopolymer length profile for current haplotype
                int[] hrunProfile = new int[haplotypeBases.length];
                getContextHomopolymerLength(haplotypeBases,hrunProfile);
                if (DEBUG) {
                    System.out.println("Haplotype bases:");
                    System.out.println(new String(haplotypeBases));
                    for (int i=0; i < hrunProfile.length; i++)
                        System.out.format("%d",hrunProfile[i]);
                    System.out.println();
                }
                fillGapProbabilities(hrunProfile, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities);

                gapOpenProbabilityMap.put(a,contextLogGapOpenProbabilities);
                gapContProbabilityMap.put(a,contextLogGapContinuationProbabilities);

            }
        }
        for (PileupElement p: pileup) {

            // check if we've already computed likelihoods for this pileup element (i.e. for this read at this location)
            if (indelLikelihoodMap.containsKey(p)) {
                HashMap<Allele,Double> el = indelLikelihoodMap.get(p);
                int j=0;
                for (Allele a: haplotypeMap.keySet()) {
                    readLikelihoods[readIdx][j++] = el.get(a);
                }
            }
            else {
                //System.out.format("%d %s\n",p.getRead().getAlignmentStart(), p.getRead().getClass().getName());
                GATKSAMRecord read = ReadUtils.hardClipAdaptorSequence(p.getRead());
                if (read == null)
                    continue;

                if(ReadUtils.is454Read(read) && !getGapPenaltiesFromFile) {
                    continue;
                }

                double[] recalQuals = null;

 /*
                if (getGapPenaltiesFromFile) {
                    RecalDataManager.parseSAMRecord( read, RAC );


                    recalQuals = new double[read.getReadLength()];

                    //compute all covariate values for this read
                    final Comparable[][] covariateValues_offset_x_covar =
                            RecalDataManager.computeCovariates((GATKSAMRecord) read, requestedCovariates);
                    // For each base in the read
                    for( int offset = 0; offset < read.getReadLength(); offset++ ) {

                        final Object[] fullCovariateKey = covariateValues_offset_x_covar[offset];

                        Byte qualityScore = (Byte) qualityScoreByFullCovariateKey.get(fullCovariateKey);
                        if(qualityScore == null)
                        {
                            qualityScore = performSequentialQualityCalculation( fullCovariateKey );
                            qualityScoreByFullCovariateKey.put(qualityScore, fullCovariateKey);
                        }

                        recalQuals[offset] = -((double)qualityScore)/10.0;
                    }

                    // for each read/haplotype combination, compute likelihoods, ie -10*log10(Pr(R | Hi))
                    // = sum_j(-10*log10(Pr(R_j | Hi) since reads are assumed to be independent
                    if (DEBUG)  {
                        System.out.format("\n\nStarting read:%s S:%d US:%d E:%d UE:%d C:%s\n",read.getReadName(),
                                read.getAlignmentStart(),
                                read.getUnclippedStart(), read.getAlignmentEnd(), read.getUnclippedEnd(),
                                read.getCigarString());

                        byte[] bases = read.getReadBases();
                        for (int k = 0; k < recalQuals.length; k++) {
                            System.out.format("%c",bases[k]);
                        }
                        System.out.println();

                        for (int k = 0; k < recalQuals.length; k++) {
                            System.out.format("%.0f ",recalQuals[k]);
                        }
                        System.out.println();
                    }
                }        */
                // get bases of candidate haplotypes that overlap with reads
                final int trailingBases = 3;

                long readStart = read.getUnclippedStart();
                long readEnd = read.getUnclippedEnd();

                int numStartSoftClippedBases, numEndSoftClippedBases;

                // see if we want to use soft clipped bases. Aligners may soft clip all bases at insertions because they don't match,
                // but they're actually consistent with the insertion!
                // Rule: if a read starts in interval [eventStart-eventLength,eventStart+1] and we are at an insertion, we'll use all soft clipped bases at the beginning.
                // Conversely, if a read ends at [eventStart,eventStart+eventLength] we'll use all soft clipped bases in the end of the read.
                long eventStartPos = ref.getLocus().getStart();

                // compute total number of clipped bases (soft or hard clipped)
                numStartSoftClippedBases = read.getAlignmentStart()- read.getUnclippedStart();
                numEndSoftClippedBases = read.getUnclippedEnd()- read.getAlignmentEnd();

                // check for hard clips (never consider these bases):
                Cigar c = read.getCigar();
                CigarElement first = c.getCigarElement(0);
                CigarElement last = c.getCigarElement(c.numCigarElements()-1);
                int numStartHardClippedBases = 0, numEndHardClippedBases = 0;

                if (first.getOperator() == CigarOperator.H) {
                    numStartHardClippedBases = first.getLength();
                }

                if (last.getOperator() == CigarOperator.H) {
                    numEndHardClippedBases = last.getLength();
                }

                // correct for hard clips
                numStartSoftClippedBases -= numStartHardClippedBases;
                numEndSoftClippedBases -= numEndHardClippedBases;
                readStart += numStartHardClippedBases;
                readEnd -= numEndHardClippedBases;

                // remove soft clips if necessary
                if ((read.getAlignmentStart()>=eventStartPos-eventLength && read.getAlignmentStart() <= eventStartPos+1) ||
                        (read.getAlignmentEnd() >= eventStartPos && read.getAlignmentEnd() <= eventStartPos + eventLength)) {
                    numStartSoftClippedBases = 0;
                    numEndSoftClippedBases = 0;
                }



                byte[] unclippedReadBases, unclippedReadQuals;

                int numStartClippedBases = numStartSoftClippedBases;
                int numEndClippedBases = numEndSoftClippedBases;
                unclippedReadBases = read.getReadBases();
                unclippedReadQuals = read.getBaseQualities();

                // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
                // and may leave a string of Q2 bases still hanging off the reads.
                for (int i=numStartSoftClippedBases; i < unclippedReadBases.length; i++) {
                    if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                        numStartClippedBases++;
                    else
                        break;

                }
                for (int i=unclippedReadBases.length-numEndSoftClippedBases-1; i >= 0; i-- ){
                    if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                        numEndClippedBases++;
                    else
                        break;
                }

                int extraOffset = Math.abs(eventLength);

                long start = Math.max(readStart + numStartClippedBases - trailingBases - ReadUtils.getFirstInsertionOffset(read)-extraOffset, 0);
                long stop =  readEnd -numEndClippedBases  + trailingBases + ReadUtils.getLastInsertionOffset(read)+extraOffset;

                // Variables start and stop are coordinates (inclusive) where we want to get the haplotype from.
                int readLength = read.getReadLength()-numStartSoftClippedBases-numEndSoftClippedBases;
                // check if start of read will be before start of reference context
                if (start < ref.getWindow().getStart())// read starts before haplotype: read will have to be cut
                    start = ref.getWindow().getStart();

                // check also if end of read will go beyond reference context
                if (stop > ref.getWindow().getStop())
                    stop = ref.getWindow().getStop();

                // if there's an insertion in the read, the read stop position will be less than start + read legnth,
                // but we want to compute likelihoods in the whole region that a read might overlap
                if (stop <= start + readLength) {
                    stop = start + readLength-1;
                }

                // ok, we now figured out total number of clipped bases on both ends.
                // Figure out where we want to place the haplotype to score read against
                if (DEBUG)
                    System.out.format("numStartClippedBases: %d numEndClippedBases: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d\n",
                            numStartClippedBases, numEndClippedBases, ref.getWindow().getStart(), ref.getWindow().getStop(), start, stop, read.getReadLength());



                LinkedHashMap<Allele,Double> readEl = new LinkedHashMap<Allele,Double>();

                if (numStartClippedBases + numEndClippedBases >= unclippedReadBases.length) {
                    if (DEBUG)
                        System.out.println("BAD READ!!");

                    int j=0;
                    for (Allele a: haplotypeMap.keySet()) {
                        readEl.put(a,0.0);
                        readLikelihoods[readIdx][j++] = 0.0;
                    }

                }
                else {
                    byte[] readBases = Arrays.copyOfRange(unclippedReadBases,numStartClippedBases,
                            unclippedReadBases.length-numEndClippedBases);

                    byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,numStartClippedBases,
                            unclippedReadBases.length-numEndClippedBases);

                    double[] recalCDP = null;
                    if (getGapPenaltiesFromFile) {
                        recalCDP = Arrays.copyOfRange(recalQuals,numStartClippedBases,
                                unclippedReadBases.length-numEndClippedBases);

                    }

                    if (DEBUG) {
                        System.out.println("Read bases:");
                        System.out.println(new String(readBases));
                    }

                    int j=0;
                    for (Allele a: haplotypeMap.keySet()) {


                        Haplotype haplotype = haplotypeMap.get(a);
                        if (stop > haplotype.getStopPosition())
                            stop = haplotype.getStopPosition();

                        if (start < haplotype.getStartPosition())
                            start = haplotype.getStartPosition();

                        // cut haplotype bases
                        long indStart = start - haplotype.getStartPosition();
                        long indStop =  stop - haplotype.getStartPosition();

                        byte[] haplotypeBases = Arrays.copyOfRange(haplotype.getBasesAsBytes(),
                                (int)indStart, (int)indStop);

                        if (DEBUG) {
                            System.out.println("Haplotype to test:");
                            System.out.println(new String(haplotypeBases));
                        }

                        Double readLikelihood = 0.0;
                        if (useAffineGapModel) {

                            double[] currentContextGOP = null;
                            double[] currentContextGCP = null;

                            if (doContextDependentPenalties) {

                               if (getGapPenaltiesFromFile) {
                                   readLikelihood = computeReadLikelihoodGivenHaplotypeAffineGaps(haplotypeBases, readBases, readQuals, recalCDP, null);

                               }  else {
                                   currentContextGOP = Arrays.copyOfRange(gapOpenProbabilityMap.get(a), (int)indStart, (int)indStop);
                                   currentContextGCP = Arrays.copyOfRange(gapContProbabilityMap.get(a), (int)indStart, (int)indStop);
                                   readLikelihood = computeReadLikelihoodGivenHaplotypeAffineGaps(haplotypeBases, readBases, readQuals, currentContextGOP, currentContextGCP);
                               }
                            }

                        }
                        else
                            readLikelihood = computeReadLikelihoodGivenHaplotype(haplotypeBases, readBases, readQuals);

                        readEl.put(a,readLikelihood);
                        readLikelihoods[readIdx][j++] = readLikelihood;
                    }
                }
                indelLikelihoodMap.put(p,readEl);
            }
            readIdx++;
        }

        if (DEBUG) {
            System.out.println("\nLikelihood summary");
            for (readIdx=0; readIdx < pileup.getReads().size(); readIdx++) {
                System.out.format("Read Index: %d ",readIdx);
                for (int i=0; i < readLikelihoods[readIdx].length; i++)
                    System.out.format("L%d: %f ",i,readLikelihoods[readIdx][i]);
                System.out.println();
            }

        }
        for (int i=0; i < numHaplotypes; i++) {
            for (int j=i; j < numHaplotypes; j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                 for (readIdx=0; readIdx < pileup.getReads().size(); readIdx++) {

                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    if (Double.isInfinite(readLikelihoods[readIdx][i]) && Double.isInfinite(readLikelihoods[readIdx][j]))
                        continue;
                    haplotypeLikehoodMatrix[i][j] += ( MathUtils.softMax(readLikelihoods[readIdx][i],
                            readLikelihoods[readIdx][j]) + LOG_ONE_HALF);

                }


            }
        }

        return getHaplotypeLikelihoods(haplotypeLikehoodMatrix);

    }

    public static double[] getHaplotypeLikelihoods(double[][] haplotypeLikehoodMatrix) {
        int hSize = haplotypeLikehoodMatrix.length;
        double[] genotypeLikelihoods = new double[hSize*(hSize+1)/2];

        int k=0;
        double maxElement = Double.NEGATIVE_INFINITY;
        for (int j=0; j < hSize; j++) {
            for (int i=0; i <= j; i++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
                if (haplotypeLikehoodMatrix[i][j] > maxElement)
                    maxElement = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize
        for (int i=0; i < genotypeLikelihoods.length; i++)
            genotypeLikelihoods[i] -= maxElement;

        return genotypeLikelihoods;
    }

    /**
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     *   - calculate the global quality score shift across all data [DeltaQ]
     *   - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     *      -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     *   - The final shift equation is:
     *
     *      Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     * @param key The list of Comparables that were calculated from the covariates
     * @return A recalibrated quality score as a byte
     */
 /*
    private byte performSequentialQualityCalculation( final Object... key ) {

        final byte qualFromRead = (byte)Integer.parseInt(key[1].toString());
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] qualityScoreCollapsedKey = new Object[2];
        final Object[] covariateCollapsedKey = new Object[3];

        // The global quality shift (over the read group only)
        readGroupCollapsedKey[0] = key[0];
        final RecalDatum globalRecalDatum = ((RecalDatum)dataManager.getCollapsedTable(0).get( readGroupCollapsedKey ));
        double globalDeltaQ = 0.0;
        if( globalRecalDatum != null ) {
            final double globalDeltaQEmpirical = globalRecalDatum.getEmpiricalQuality();
            final double aggregrateQReported = globalRecalDatum.getEstimatedQReported();
            globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
        }

        // The shift in quality between reported and empirical
        qualityScoreCollapsedKey[0] = key[0];
        qualityScoreCollapsedKey[1] = key[1];
        final RecalDatum qReportedRecalDatum = ((RecalDatum)dataManager.getCollapsedTable(1).get( qualityScoreCollapsedKey ));
        double deltaQReported = 0.0;
        if( qReportedRecalDatum != null ) {
            final double deltaQReportedEmpirical = qReportedRecalDatum.getEmpiricalQuality();
            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        // The shift in quality due to each covariate by itself in turn
        double deltaQCovariates = 0.0;
        double deltaQCovariateEmpirical;
        covariateCollapsedKey[0] = key[0];
        covariateCollapsedKey[1] = key[1];
        for( int iii = 2; iii < key.length; iii++ ) {
            covariateCollapsedKey[2] =  key[iii]; // The given covariate
            final RecalDatum covariateRecalDatum = ((RecalDatum)dataManager.getCollapsedTable(iii).get( covariateCollapsedKey ));
            if( covariateRecalDatum != null ) {
                deltaQCovariateEmpirical = covariateRecalDatum.getEmpiricalQuality();
                deltaQCovariates += ( deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported) );
            }
        }

        final double newQuality = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
        return QualityUtils.boundQual( (int)Math.round(newQuality), (byte)MAX_QUALITY_SCORE );

        // Verbose printouts used to validate with old recalibrator
        //if(key.contains(null)) {
        //    System.out.println( key  + String.format(" => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte));
        //}
        //else {
        //    System.out.println( String.format("%s %s %s %s => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 key.get(0).toString(), key.get(3).toString(), key.get(2).toString(), key.get(1).toString(), qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte) );
        //}

        //return newQualityByte;

    }
*/
}
