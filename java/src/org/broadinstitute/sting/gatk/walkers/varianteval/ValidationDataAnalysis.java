package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Sep 30, 2009
 * Time: 5:28:57 PM
 * To change this template use File | Settings | File Templates.
 */

//  @Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="eval",type=VariationRod.class), @RMD(name="dbsnp",type= rodDbSNP.class),@RMD(name="hapmap-chip",type=RodGenotypeChipAsGFF.class), @RMD(name="interval",type=IntervalRod.class), @RMD(name="validation",type=RodGenotypeChipAsGFF.class)})

public class ValidationDataAnalysis extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    private int calls_at_validated_sites = 0;
    private int calls_at_sites_validated_true = 0;
    private int validated_sites = 0;

    public ValidationDataAnalysis() {
        super("validation_data_analysis");
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {

        validated_sites++;
        Variation val_data = (Variation) tracker.lookup("validation", null);
        Variation dbsnp = (Variation) tracker.lookup("dbsnp",null);

        if (eval != null) {
            calls_at_sites_validated_true++;
            //out.format("Has validaiton data: %s\n", val_data.getLocation());
            if (val_data != null) {
                //out.format("Validated true: %s\n", val_data.getLocation());
                calls_at_validated_sites++;
            }
        }
        //out.println(context.getLocation());

        return null;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        if (calls_at_validated_sites > 0) { // only output info if there were any validation sites encountered
            s.add(String.format("validated sites               %d", validated_sites));
            s.add(String.format("calls at validated sites      %d", calls_at_validated_sites));
            s.add(String.format("calls at sites validated true %d", calls_at_sites_validated_true));
            s.add(String.format("%% validated true             %f", (float) calls_at_validated_sites / calls_at_sites_validated_true));
        }else{
            s.add("No validation data encountered");
        }
        return s;
    }
}
