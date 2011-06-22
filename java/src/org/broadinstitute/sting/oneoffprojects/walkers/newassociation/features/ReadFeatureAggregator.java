package org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 12:33 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReadFeatureAggregator<X> {

    protected double mean;
    protected double var;
    protected int nReads;

    public ReadFeatureAggregator(RFAArgumentCollection collection) {
        init(collection);
        mean = 0.0;
        var = 0.0;
        nReads = 0;
    }

    public void aggregate(SAMRecord record) {
        if ( featureDefined(record) ) {
            aggregate(extractFeature(record));
        }
    }

    protected abstract void aggregate(X feature);

    protected abstract boolean featureDefined(SAMRecord record);

    protected abstract X extractFeature(SAMRecord record);

    public double getMean() { return mean; }
    public double getVar() { return var; }
    public double getUnbiasedVar() { return var*( (double) nReads)/(nReads-1); }
    public int getnReads() { return nReads; }

    public void init(RFAArgumentCollection collection) { }

    public X parse(SAMRecord read) {
        if ( featureDefined(read) ) {
            return extractFeature(read);
        } else {
            return null;
        }
    }

    public String parseStr(SAMRecord read) {
        if ( featureDefined(read) ) {
            return extractFeature(read).toString();
        } else {
            return "undefined";
        }
    }

}
