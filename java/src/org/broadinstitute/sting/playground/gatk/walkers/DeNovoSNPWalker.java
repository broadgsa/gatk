package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenotypeUtils;

import java.util.List;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Sep 1, 2009
 * Time: 11:04:55 AM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="child",type= AllelicVariant.class), @RMD(name="parent1",type= AllelicVariant.class), @RMD(name="parent2",type= AllelicVariant.class)})
//@Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="eval",type=AllelicVariant.class), @RMD(name="dbsnp",type=AllelicVariant.class),@RMD(name="hapmap-chip",type=AllelicVariant.class), @RMD(name="interval",type=IntervalRod.class)})
//@Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="dbsnp",type=AllelicVariant.class)})
public class DeNovoSNPWalker extends RefWalker<String, Integer>{

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        AllelicVariant child = (AllelicVariant)tracker.lookup("child", null);
        AllelicVariant parent1 = (AllelicVariant)tracker.lookup("parent1", null);
        AllelicVariant parent2 = (AllelicVariant)tracker.lookup("parent2", null);
        AllelicVariant dbsnp = (AllelicVariant)tracker.lookup("dbSNP", null);
        if (child != null && parent1 != null && parent2 != null) {
            if (child.isSNP() &&
                child.getVariationConfidence() > 5 && // BTR > 5
                parent1.isReference() &&
                parent1.getConsensusConfidence() > 5 &&
                parent2.isReference() && 
                parent2.getConsensusConfidence() > 5
                ) {

                double sumConfidences = 0.5 * (0.5 * child.getVariationConfidence() + Math.min(parent1.getConsensusConfidence(), parent2.getConsensusConfidence()));

                out.format("%s\t", child.getLocation().getContig());
                out.format("%s\t", child.getLocation().getStart());
                out.format("%.4f\t", sumConfidences);
                out.format("%.4f\t", child.getVariationConfidence());
                out.format("%.4f\t", parent1.getConsensusConfidence());
                out.format("%.4f\t", parent2.getConsensusConfidence());
                out.format("%s\t", dbsnp != null);

                out.format ("%s\t", child.toString());
                out.format ("%s\t", parent1.toString());
                out.format ("%s", parent2.toSimpleString());
                if (dbsnp != null)
                    out.format ("\tDBSNP\t:%s", dbsnp.toString());
                out.println();
            }
        }

        return "";
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce(String line, Integer a) {
        return 1;
    }


}
