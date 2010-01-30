package org.broadinstitute.sting.oneoffprojects.walkers.vcftools;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import java.lang.Long;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Jan 29, 2010
 */
public abstract class SimpleVCFIntersectWalker extends RodWalker<VCFRecordPair,Long>{

    public void initialize() {

    }

    public Long reduceInit() {
        return 0l;
    }

    public VCFRecordPair map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return null;
    }

    public Long reduce(VCFRecordPair records, long prevReduce) {
        return 0l;
    }

    public void onTraversalDone(long finalReduce) {
        return;
    }
}

class VCFRecordPair {
    public VCFRecord rec1;
    public VCFRecord rec2;
}
