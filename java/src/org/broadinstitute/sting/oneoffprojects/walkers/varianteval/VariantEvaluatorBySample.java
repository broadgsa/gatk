package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvaluator;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;

import java.util.*;

/**
 * An abstract way to break variant analyses down by sample. SampleDataPoint objects (e.g. its inheritors) are propagated
 * into a per-sample table, which is updated only when a specific sample's genotype is such that the module-defined
 * includeGenotype(G) returns true.
 * @Author chartl
 */
public abstract class VariantEvaluatorBySample extends VariantEvaluator {
    @DataPoint(name="VariantEvaluatorBySample",description="Evaluation broken down by sample")
    EvalBySample evalBySample;

    public VariantEvaluatorBySample(VariantEvalWalker parent) {
        super(parent);
        evalBySample = initializeTable();
    }

    public abstract String getTableName();

    public abstract List<SampleDataPoint> getDataPoints();

    public abstract boolean includeGenotype(Genotype g);

    public EvalBySample initializeTable() {
        if ( enabled() ) {
            EvalBySample ebs = new EvalBySample(getTableName(),getDataPoints());
            return ebs;
        } else {
            return null;
        }
    }

    // note -- this only updates at all sites after the first site where a sample has been identified containing a variant genotype
    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for ( Map.Entry<String,List<SampleDataPoint>> entry : evalBySample.sampleAndEvalResults.entrySet() ) {
            for ( SampleDataPoint dp : entry.getValue() ) {
                dp.update0(tracker,ref,context);
            }
        }
    }

    public String update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for ( String sample : vc1.getSampleNames() ) {
            if ( includeGenotype(vc1.getGenotype(sample)) ) {
                if ( ! evalBySample.containsKey(sample) ) {
                    evalBySample.put(sample,getDataPoints());
                }

                for ( SampleDataPoint dp : evalBySample.sampleAndEvalResults.get(sample) ) {
                    dp.update1(vc1,tracker,ref,context);
                }
            }
        }

        return null; // don't return interesting sites
    }

    public String update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc1 == null ) {
            return null; // cannot update by sample if there are no samples
        }
        for ( String sample : vc1.getSampleNames() ) {
            if ( includeGenotype(vc1.getGenotype(sample)) ) {
                if ( ! evalBySample.containsKey(sample) ) {
                    evalBySample.put(sample,getDataPoints());
                }

                for ( SampleDataPoint dp : evalBySample.sampleAndEvalResults.get(sample) ) {
                    dp.update2(vc1,vc2,tracker,ref,context);
                }
            }
        }

        return null; // don't return interesting sites
    }

    @Override
    public void finalizeEvaluation() {
        evalBySample.finalizeTable();
    }

}

abstract class SampleDataPoint {
    public String name;
    public String sampleName;

    public SampleDataPoint(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setSampleName(String sName) {
        sampleName = sName;
    }

    public abstract String toString();

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {}

    public void update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {}

    public void update2(VariantContext eval, VariantContext comp,RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {}

    public void finalizeCalculation() {}
}

class EvalBySample implements TableType {
    public String[] evalNames;
    public TreeMap<String, List<SampleDataPoint>> sampleAndEvalResults;
    public String name;
    private HashMap<String,Integer> nameToDataPointOffset;

    private Object[][] finalizedResults;

    public EvalBySample(String name, Collection<SampleDataPoint> evals) {
        int i = 0;
        this.evalNames = new String[evals.size()];
        this.nameToDataPointOffset = new HashMap<String,Integer>(evals.size());
        for ( SampleDataPoint s : evals ) {
            this.evalNames[i] = s.getName();
            this.nameToDataPointOffset.put(s.getName(),i);
            i++;
        }

        this.name = name;
        this.sampleAndEvalResults = new TreeMap<String,List<SampleDataPoint>>();
    }

    public Object[] getColumnKeys() {
        //System.out.printf("%s%n","Call to column keys");
        return evalNames;
    }

    public String getCell(int x, int y) {
        return finalizedResults[x][y].toString();
    }

    public String getName() {
        return name;
    }

    public Object[] getRowKeys() {
        String[] rowNames = new String[sampleAndEvalResults.size()];
        int i = 0;
        for ( Map.Entry<String,List<SampleDataPoint>> e : sampleAndEvalResults.entrySet() ) {
            rowNames[i] = e.getKey();
            i++;
        }

        //System.out.printf("%s%n","Call to row keys");

        return rowNames;
    }

    public void finalizeTable() {
        if ( sampleAndEvalResults == null || sampleAndEvalResults.size() == 0 ) {
            finalizedResults = new Object[0][0];
            return; // todo -- early return is hacky
        }
        finalizedResults = new Object[sampleAndEvalResults.size()][sampleAndEvalResults.firstEntry().getValue().size()];
        int i = 0;
        for ( Map.Entry<String,List<SampleDataPoint>> evalBySample : sampleAndEvalResults.entrySet() ) {
            int j = 0;
            for ( SampleDataPoint o : evalBySample.getValue() ) {
                o.finalizeCalculation();
                finalizedResults[i][j] = o;
                j++;
            }
            i++;
        }
    }

    public boolean containsKey(String sample) {
        return sampleAndEvalResults.containsKey(sample);
    }

    public void put(String sample, List<SampleDataPoint> dataPoints) {
        for ( SampleDataPoint dp : dataPoints ) {
            dp.setSampleName(sample);
        }
        sampleAndEvalResults.put(sample,dataPoints);
    }

}
