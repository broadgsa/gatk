package org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 12:39 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NumericFeatureAggregator extends ReadFeatureAggregator<Integer> {

    private int min;
    private int max;

    public NumericFeatureAggregator(RFAArgumentCollection col) {
        super(col);
        min = -1;
        max = -1;
    }

    protected void aggregate(Integer datum) {
        if ( min == -1 ) {
            min = datum;
        } else if ( max == -1 ) {
            if ( datum > min ) {
                max = datum;
            } else {
                max = min;
                min = datum;
            }
        } else if ( datum > max ) {
            update(max);
            max = datum;
        } else if ( datum < min ) {
            update(min);
            min = datum;
        } else {
            update(datum);
        }
    }

    protected void update(Integer datum) {
        double oldMean = mean;
        mean += (datum - mean)/(1+nReads);
        var = ((nReads*var) + (datum - oldMean)*(datum-mean))/++nReads;
    }
}
