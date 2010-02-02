package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.List;
import java.util.ArrayList;
import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 20, 2009
 * Time: 4:03:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class PooledFrequencyAnalysis extends BasicPoolVariantAnalysis implements PoolAnalysis, GenotypeAnalysis {
    /*
        Maintains an array of basic variant analyses broken down by estimated frequency
         */
    private VariantDBCoverage[] coverageAnalysisByFrequency;
    private VariantCounter[] variantCounterByFrequency;
    private TransitionTranversionAnalysis[] transitionTransversionByFrequency;

    public PooledFrequencyAnalysis(int poolSize, String knownDBSNPName ) {
        super("Pooled_Frequency_Analysis",poolSize);
        if ( poolSize > 0 ) {
            coverageAnalysisByFrequency = new VariantDBCoverage[getNumberOfAllelesInPool()+1];
            variantCounterByFrequency = new VariantCounter[getNumberOfAllelesInPool()+1];
            transitionTransversionByFrequency = new TransitionTranversionAnalysis[getNumberOfAllelesInPool()+1];
            for ( int j = 0; j < getNumberOfAllelesInPool()+1; j ++ ) {
                coverageAnalysisByFrequency[j] = new VariantDBCoverage(knownDBSNPName);
                variantCounterByFrequency[j] = new VariantCounter();
                transitionTransversionByFrequency[j] = new TransitionTranversionAnalysis();
            }
        }
    }      

    public void initialize(VariantEvalWalker master, PrintStream out1, PrintStream out2, String name) {
        super.initialize(master,out1,out2,name);

    }

    public String update( Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context ) {
        int f = calculateFrequencyFromVariation(eval);
        String interest1 = coverageAnalysisByFrequency[f].update(eval,tracker,ref,context);
        String interest2 = variantCounterByFrequency[f].update(eval,tracker,ref,context);
        String interest3 = transitionTransversionByFrequency[f].update(eval,tracker,ref,context);

        if ( interest1 != null ) {
            return "coverageAnalysis_Frequency_"+f+"\t"+interest1;
        }

        if ( interest2 != null ) {
            return "VariantCounter_Frequency_"+f+"\t"+interest2;
        }

        if ( interest3 != null ) {
            return "transitionTransversion_Frequency_"+f+"\t"+interest3;
        }

        return null;
    }

    public int calculateFrequencyFromVariation(Variation eval) {
        double freq = eval.getNonRefAlleleFrequency();
        return (int) Math.round(getNumberOfAllelesInPool()*freq);
    }

    @Override
    public List<String> done() {
        List<String> s = new ArrayList<String>();
        for ( int j = 0; j < getNumberOfAllelesInPool()+1; j ++ ) {
            s.add("Pool_Frequency_Analysis_Frequency:\t"+j);
            s.addAll(coverageAnalysisByFrequency[j].done());
            s.addAll(variantCounterByFrequency[j].done());
            s.addAll(transitionTransversionByFrequency[j].done());
        }

        return s;
    }

}
