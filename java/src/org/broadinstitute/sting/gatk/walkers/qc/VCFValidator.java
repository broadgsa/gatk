package org.broadinstitute.sting.gatk.walkers.qc;


import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.StingException;

/**
 * A light-weight validator for a VCF file.
 */
@Requires(value={},referenceMetaData=@RMD(name="vcf",type= RodVCF.class))
public class VCFValidator extends RodWalker<Integer, Integer> {

    /**
     * It's about as simple as things come right now.  We let the rod system process all of the
     * entries in the file, and if no errors pop up in processing, then it validates!
     */

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker != null ) {
            RODRecordList<ReferenceOrderedDatum> rodlist = tracker.getTrackData("vcf", null);
            if ( rodlist != null ) {
                RodVCF rod = (RodVCF)rodlist.getRecords().get(0);
                if ( (rod.isSNP() || rod.isReference()) &&  rod.getReference().charAt(0) != ref.getBase() )
                    throw new StingException("The reference base (" + ref.getBase() + ") does not match the base from the VCF record (" + rod.getReference() + ")");
            }
        }
        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        out.println("The input file is a valid VCF");
    }
}
