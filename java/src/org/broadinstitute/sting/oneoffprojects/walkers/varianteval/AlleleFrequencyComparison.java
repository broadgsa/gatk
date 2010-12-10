package org.broadinstitute.sting.oneoffprojects.walkers.varianteval;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvaluator;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.utils.TableType;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
            // todo -- this is a godawful hack. The "right way" isn't working, so do it the unsafe way for now. Note that
            // todo -- this precludes getting the AC/AF values from the info field because some may not be there...
            /*if ( missingField(eval) ) {
                recalculateCounts(eval);
            }
            if ( missingField(comp) ) {
                recalculateCounts(comp);
            }*/
            HashMap<String,Object> evalCounts = new HashMap<String,Object>(2);
            HashMap<String,Object> compCounts = new HashMap<String,Object>(2);

            VariantContextUtils.calculateChromosomeCounts(eval,evalCounts,false);
            VariantContextUtils.calculateChromosomeCounts(comp,compCounts,false);
            afTable.update(((List<Double>)evalCounts.get("AF")).get(0),((List<Double>)compCounts.get("AF")).get(0));
            acTable.update(((List<Integer>)evalCounts.get("AC")).get(0),((List<Integer>)compCounts.get("AC")).get(0));
        }

        return null; // there is nothing interesting
    }

    private static boolean missingField(final VariantContext vc) {
        return ! ( vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) && vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) );
    }

    private void recalculateCounts(VariantContext vc) {
        Map<String,Object> attributes = new HashMap<String,Object>();
        VariantContextUtils.calculateChromosomeCounts(vc,attributes,false);
        vc = VariantContext.modifyAttributes(vc,attributes);
        getLogger().debug(String.format("%s %s | %s %s",attributes.get("AC"),attributes.get("AF"),vc.getAttribute("AC"),vc.getAttribute("AF")));
        if ( attributes.size() == 2 && missingField(vc) ) {
            throw new org.broadinstitute.sting.utils.exceptions.StingException("VariantContext should have had attributes modified but did not");
        }
    }

    private static boolean isValidVC(final VariantContext vc) {
        return (vc != null && !vc.isFiltered() && vc.getAlternateAlleles().size() == 1);
    }

    private static double getAF(VariantContext vc) {
        Object af = vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY);
        if ( af == null ) {
	    //throw new UserException("Variant context "+vc.getName()+" does not have allele frequency entry which is required for this walker");
            // still none after being re-computed; this is 0.00
            return 0.00;
        } else if ( List.class.isAssignableFrom(af.getClass())) {
            return ( (List<Double>) af ).get(0);
        } else if ( String.class.isAssignableFrom(af.getClass())) {
            // two possibilities
            String s = (String) af;
            try {
                if ( s.startsWith("[") ) {
                    return Double.parseDouble(s.replace("\\[","").replace("\\]",""));
                } else {
                    return Double.parseDouble(s);
                }
            } catch (NumberFormatException e) {
                throw new UserException("Allele frequency field may be improperly formatted, found AF="+s,e);
            }
        } else if ( Double.class.isAssignableFrom(vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY).getClass())) {
            return (Double) af;
        } else {
            throw new UserException(String.format("Class of Allele Frequency does not appear to be formated, had AF=%s, of class %s",af.toString(),af.getClass()));
        }
    }

    private static int getAC(VariantContext vc) {
        Object ac = vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY);
        if ( ac == null ) {
            // still none after being re computed; this is 0
            return 0;
        } else if ( List.class.isAssignableFrom(ac.getClass())) {
            return ( (List<Integer>) ac ).get(0);
        } else if ( String.class.isAssignableFrom(ac.getClass())) {
            // two possibilities
            String s = (String) ac;
            try {
                if ( s.startsWith("[") ) {
                    return Integer.parseInt(s.replace("\\[","").replace("\\]",""));
                } else {
                    return Integer.parseInt(s);
                }
            } catch (NumberFormatException e) {
                throw new UserException(String.format("Allele count field may be improperly formatted, found AC=%s for record %s:%d",ac,vc.getChr(),vc.getStart()),e);
            }
        } else if ( Integer.class.isAssignableFrom(ac.getClass())) {
            return (Integer) ac;
        } else {
            throw new UserException(String.format("Class of Allele Frequency does not appear to be formated, had AF=%s, of class %s",ac.toString(),ac.getClass()));
        }
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
