package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariationRod;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Sep 1, 2009
 * Time: 11:04:55 AM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="child",type= VariationRod.class), @RMD(name="parent1",type= VariationRod.class), @RMD(name="parent2",type= VariationRod.class)})
//@Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="eval",type=AllelicVariant.class), @RMD(name="dbsnp",type=AllelicVariant.class),@RMD(name="hapmap-chip",type=AllelicVariant.class), @RMD(name="interval",type=IntervalRod.class)})
//@Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="dbsnp",type=AllelicVariant.class)})
public class DeNovoSNPWalker extends RefWalker<String, Integer>{

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Variation child = (Variation)tracker.lookup("child", null);
        Variation parent1 = (Variation)tracker.lookup("parent1", null);
        Variation parent2 = (Variation)tracker.lookup("parent2", null);
        Variation dbsnp = (Variation)tracker.lookup("dbSNP", null);
        if (child != null && parent1 != null && parent2 != null) {
            if (!(parent1 instanceof VariantBackedByGenotype) || !(parent2 instanceof VariantBackedByGenotype))
                        throw new StingException("Both parents ROD tracks must be backed by genotype data. Ensure that your parent rod(s) contain genotyping information");
            if (child.isSNP() &&
                child.getNegLog10PError() > 5 && // BTR > 5
                parent1.isReference() &&
                ((VariantBackedByGenotype)parent1).getCalledGenotype().getNegLog10PError() > 5 &&
                parent2.isReference() && 
                ((VariantBackedByGenotype)parent2).getCalledGenotype().getNegLog10PError() > 5
                ) {

                double sumConfidences = 0.5 * (0.5 * child.getNegLog10PError() +
                        Math.min(((VariantBackedByGenotype)parent1).getCalledGenotype().getNegLog10PError(),
                                 ((VariantBackedByGenotype)parent2).getCalledGenotype().getNegLog10PError()));

                out.format("%s\t", child.getLocation().getContig());
                out.format("%s\t", child.getLocation().getStart());
                out.format("%.4f\t", sumConfidences);
                out.format("%.4f\t", child.getNegLog10PError());
                out.format("%.4f\t", ((VariantBackedByGenotype)parent1).getCalledGenotype().getNegLog10PError());
                out.format("%.4f\t", ((VariantBackedByGenotype)parent2).getCalledGenotype().getNegLog10PError());
                out.format("%s\t", dbsnp != null);

                out.format ("%s\t", child.toString());
                out.format ("%s\t", parent1.toString());
                out.format ("%s", parent2.toString());
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
