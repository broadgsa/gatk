package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;

/**
 * Test routine for new VariantContext object
 */
public class TestVariantContextWalker extends RodWalker<Integer, Integer> {

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc cur = context.getLocation();

        if ( ref == null )
            return 0;
        else {
            RODRecordList<ReferenceOrderedDatum> dbsnpList = tracker.getTrackData("dbsnp", null);

            if (dbsnpList == null)
                return 0;
            else {
                int n = 0;
                for (ReferenceOrderedDatum d : dbsnpList) {
                    rodDbSNP dbsnpRecord = (rodDbSNP)d;
                    VariantContext vc = VariantContextAdaptors.dbsnp2VariantContext(dbsnpRecord);
                    if ( vc != null ) {
                        n++;
                        System.out.printf("%s%n", vc);
                    }
                }

                return n;
            }
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer point, Integer sum) {
        return point + sum;
    }
}