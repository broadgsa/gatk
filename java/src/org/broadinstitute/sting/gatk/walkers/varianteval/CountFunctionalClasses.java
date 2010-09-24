package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

@Analysis(name = "Count Functional Classes", description = "Counts instances of different functional variant classes (provided the variants are annotated with that information)")
public class CountFunctionalClasses extends VariantEvaluator {
    // the following fields are in output order:
    @DataPoint(description = "miRNA")
    long nMiRNA= 0;

    @DataPoint(description = "3'-UTR")
    long nUTR3 = 0;

    @DataPoint(description = "Intron")
    long nIntron = 0;

    @DataPoint(description = "Splice-site")
    long nSpliceSite= 0;

    @DataPoint(description = "Read-through")
    long nReadThrough = 0;

    @DataPoint(description = "Nonsense")
    long nNonsense = 0;

    @DataPoint(description = "Missense")
    long nMissense = 0;

    @DataPoint(description = "Synonymous")
    long nSynonymous = 0;

    @DataPoint(description = "5'-UTR")
    long nUTR5= 0;

    @DataPoint(description = "Promoter")
    long nPromoter = 0;

    public CountFunctionalClasses(VariantEvalWalker parent) {
        super(parent);
    }

    public String getName() {
        return "functionalclasses";
    }

    public boolean enabled() {
        return false;
    }

    public int getComparisonOrder() {
        return 1;
    }

    public String update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String type = vc1.getAttributeAsString("type");

        if (type != null) {
            if      (type.equalsIgnoreCase("miRNA")) { nMiRNA++; }
            else if (type.equalsIgnoreCase("3'-UTR")) { nUTR3++; }
            else if (type.equalsIgnoreCase("Intron")) { nIntron++; }
            else if (type.equalsIgnoreCase("Splice_site")) { nSpliceSite++; }
            else if (type.equalsIgnoreCase("Read-through")) { nReadThrough++; }
            else if (type.equalsIgnoreCase("Nonsense")) { nNonsense++; }
            else if (type.equalsIgnoreCase("Missense")) { nMissense++; }
            else if (type.equalsIgnoreCase("Synonymous")) { nSynonymous++; }
            else if (type.equalsIgnoreCase("5'-UTR")) { nUTR5++; }
            else if (type.equalsIgnoreCase("Promoter")) { nPromoter++; }
        }

        return null;
    }
}
