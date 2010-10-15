package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.utils.TableType;

import java.util.List;

/**
 *
 */

@Analysis(name = "Allele Frequency Comparison", description = "Compare allele frequency and counts between eval and comp")
public class AlleleFrequencyComparison extends VariantEvaluator {
    private static int MAX_AC_COUNT = 100; // todo -- command line argument?

    @DataPoint(description="Counts of eval frequency versus comp frequency")
    AFTable afTable = new AFTable();

    @DataPoint(description="Counts of eval AC versus comp AC")
    ACTable acTable = new ACTable(MAX_AC_COUNT);

    public boolean enabled() { return true; }

    public int getComparisonOrder() { return 2; }

    public String getName() { return "Allele Frequency Comparison"; }

    public AlleleFrequencyComparison(VariantEvalWalker parent) {
        super(parent);
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        if ( ! (isValidVC(eval) && isValidVC(comp))  ) {
            return null;
        } else {
            afTable.update(getAF(eval),getAF(comp));
            acTable.update(getAC(eval),getAC(comp));
        }

        return null; // there is nothing interesting
    }

    private static boolean isValidVC(final VariantContext vc) {
        return (vc != null && !vc.isFiltered() && vc.getAlternateAlleles().size() == 1 && vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) && vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
    }

    private static double getAF(VariantContext vc) {
        return ((List<Double>) vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)).get(0);
    }

    private static int getAC(VariantContext vc) {
        return ((List<Integer>) vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY)).get(0);
    }
}

class AFTable implements TableType {

    protected int[][] afCounts = new int[101][101];

    public Object[] getRowKeys() {
        String[] afKeys = new String[101];
        for ( int f = 0; f < 101; f ++ ) {
            afKeys[f] = String.format("%.2f",(f+0.0)/100.0);
        }

        return afKeys;
    }

    public Object[] getColumnKeys() {
        return getRowKeys(); // nice thing about symmetric tables
    }

    public Object getCell(int i, int j) {
        return afCounts[i][j];
    }

    public String getName() {
        return "Allele Frequency Concordance";
    }

    public void update(double eval, double comp) {
        afCounts[af2index(eval)][af2index(comp)]++;
    }

    private int af2index(double d) {
        return (int) Math.round(100*d);
    }
}

class ACTable implements TableType {
    protected int[][] acCounts;
    protected int maxAC;

    public ACTable(int acMaximum) {
        maxAC = acMaximum;
        acCounts = new int[acMaximum+1][acMaximum+1];
    }

    public Object[] getRowKeys() {
        String[] acKeys = new String[maxAC+1];
        for ( int i = 0 ; i <= maxAC ; i ++ ) {
            acKeys[i] = String.format("%d",i);
        }

        return acKeys;
    }

    public Object[] getColumnKeys() {
        return getRowKeys();
    }

    public Object getCell(int i, int j) {
        return acCounts[i][j];
    }

    public String getName() {
        return "Allele Counts Concordance";
    }

    public void update(int eval, int comp) {
        eval = eval > maxAC ? maxAC : eval;
        comp = comp > maxAC ? maxAC : comp;

        acCounts[eval][comp]++;
    }

}
