package org.broadinstitute.sting.oneoffprojects.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.utils.report.tags.Analysis;

import java.util.ArrayList;
import java.util.List;

/**
 * Extends the Per-sample variant evaluator class and returns, for each sample, the number of variants, the Ti/Tv, and
 * the comp overlap. It does this only on sites where the sample is identified as hom var, or het.
 */
@org.broadinstitute.sting.utils.report.tags.Analysis(name="FunctionalClassBySample",description="Count of SNPs by functional class by sample")
public class FunctionalClassBySample extends VariantEvaluatorBySample {

    public FunctionalClassBySample(VariantEvalWalker parent) { super(parent); }

    public List<SampleDataPoint> getDataPoints() {
        List<SampleDataPoint> points = new ArrayList<SampleDataPoint>(10);
        points.add(new FCPoint("miRNA","miRNA"));
        points.add(new FCPoint("3'-UTR","3'-UTR"));
        points.add(new FCPoint("Intron","Intron"));
        points.add(new FCPoint("Splice_site","Splice_site"));
        points.add(new FCPoint("Read-through","Read-through"));
        points.add(new FCPoint("Nonsense","Nonsense"));
        points.add(new FCPoint("Missense","Missense"));
        points.add(new FCPoint("Synonymous","Synonymous"));
        points.add(new FCPoint("5'-UTR","5'-UTR"));
        points.add(new FCPoint("Promoter","Promoter"));

        return points;
    }

    public String getTableName() { return "Functional Class Counts by Sample"; }

    public String getName() { return "Functional Class Counts by Sample"; }

    public int getComparisonOrder() { return 1; }

    public boolean enabled() { return true; }

    public boolean includeGenotype(Genotype g) { return ( ! g.isFiltered() ) && ( g.isHet() || g.isHomVar() ); }

}

class FCPoint extends SampleDataPoint {
    private String matchStr;
    private int count;

    public FCPoint(String fcName, String fcMatch) {
        super(fcName);
        matchStr = fcMatch;
        count = 0;
    }

    public void update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc == null ) {
            return;
        } else {
            String type = vc.getAttributeAsString("type","none");
            if ( type.equalsIgnoreCase(matchStr) ) {
                count++;
            }
        }
    }

    public String toString() {
        return String.format("%d",count);
    }
}
