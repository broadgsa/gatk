package org.broadinstitute.sting.oneoffprojects.firehosesummary;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 12, 2010
 */
public class DepthStatisticsCalculator extends SummaryStatisticsCalculator {
    private int numBasesAbove0x; // % of bases with ANY coverage
    private int numBasesAbove3x; // % of bases with 4x+
    private int numBasesAbove9x; // % of bases with 10x+
    private int numBasesAbove24x; // % of bases with 25x+
    private int numBasesAbove49x; // % of bases with 50x+
    private int numBasesAbove99x; // % of bases with 100x+

    private int targetsAbove0x;
    private int targetsAbove3x;
    private int targetsAbove9x;
    private int targetsAbove24x;
    private int targetsAbove49x;
    private int targetsAbove99x;
    private int numTargets;

    public static double[] DEPTH_CUTOFFS = {1,4,10,25,50,100};

    public DepthStatisticsCalculator(String name) {
        super(name);
        numBasesAbove0x = 0;
        numBasesAbove3x = 0;
        numBasesAbove9x = 0;
        numBasesAbove24x = 0;
        numBasesAbove49x = 0;
        numBasesAbove99x = 0;
        targetsAbove99x = 0;
        targetsAbove49x = 0;
        targetsAbove24x = 0;
        targetsAbove9x = 0;
        targetsAbove3x = 0;
        targetsAbove0x = 0;
        numTargets = 0;
    }

    public void updateLocus(int depth) {
        super.update(depth);

        if ( depth > 99 ) {
            numBasesAbove99x++;
        }
        if ( depth > 49 ) {
            numBasesAbove49x++;
        }
        if ( depth > 24 ) {
            numBasesAbove24x++;
        }
        if ( depth > 9 ) {
            numBasesAbove9x++;
        }
        if ( depth > 3 ) {
            numBasesAbove3x++;
        }
        if ( depth > 0 ) {
            numBasesAbove0x++;
        }
    }

    public void updateTargets(int targetLength, int totalCoverage) {
        double avgCvg = ( (double) totalCoverage ) / ( (double) targetLength);

        if ( avgCvg >= 100 ) {
            targetsAbove99x ++;
        }
        if ( avgCvg >= 50 ) {
            targetsAbove49x++;
        }
        if ( avgCvg >= 25 ) {
            targetsAbove24x++;
        }
        if ( avgCvg >= 10 ) {
            targetsAbove9x++;
        }
        if ( avgCvg >= 4 ) {
            targetsAbove3x++;
        }
        if ( avgCvg >= 1 ) {
            targetsAbove0x++;
        }

        numTargets++;
    }

    public double[] getLocusProportions() {
        double[] proportions = new double[6];

        proportions[0] = ( (double) numBasesAbove0x )/( (double) getSampleSize() );
        proportions[2] = ( (double) numBasesAbove9x )/( (double) getSampleSize() );
        proportions[3] = ( (double) numBasesAbove24x )/( (double) getSampleSize() );
        proportions[4] = ( (double) numBasesAbove49x )/( (double) getSampleSize() );
        proportions[5] = ( (double) numBasesAbove99x )/( (double) getSampleSize() );
        proportions[1] = ( (double) numBasesAbove3x )/( (double) getSampleSize() );

        return proportions;
    }

    public double[] getTargetProportions() {
        double[] proportions = new double[6];

        proportions[0] = ( (double) targetsAbove0x )/( (double) numTargets );
        proportions[2] = ( (double) targetsAbove9x )/( (double) numTargets );
        proportions[3] = ( (double) targetsAbove24x )/( (double) numTargets );
        proportions[4] = ( (double) targetsAbove49x )/( (double) numTargets );
        proportions[5] = ( (double) targetsAbove99x )/( (double) numTargets );
        proportions[1] = ( (double) targetsAbove3x )/( (double) numTargets );

        return proportions;
    }

    public double getPercentWellCoveredLoci() {
        return 10*( (double) numBasesAbove9x )/( (double) getSampleSize() );
    }

    public double getPercentWellCoveredTargets() {
        return 10*( (double) targetsAbove9x )/( (double) numTargets );
    }
}
