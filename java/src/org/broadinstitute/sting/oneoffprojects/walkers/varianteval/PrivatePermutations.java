package org.broadinstitute.sting.oneoffprojects.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvaluator;
import org.broadinstitute.sting.utils.report.tags.Analysis;
import org.broadinstitute.sting.utils.report.tags.DataPoint;
import org.broadinstitute.sting.utils.report.utils.TableType;

import java.util.Arrays;
import java.util.Collection;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 22, 2010
 * Time: 12:22:08 PM
 * To change this template use File | Settings | File Templates.
 */
@Analysis(name = "PrivatePermutations", description = "Number of additional mutations from each new sample; random permutations")
public class PrivatePermutations extends VariantEvaluator {
    private final int NUM_PERMUTATIONS = 50;
    private final double LOW_GQ_PCT = 0.95;
    private final double LOW_GQ_THRSH = 30.0;
    private boolean initialized = false;
    private long skipped = 0l;

    @DataPoint(name="Marginal Number of Mutations",description="Number of additional mutations from each new sample; random permutations")
    AdditionalBySample permuteCounts = null;

    String[][] permutations;

    public boolean enabled() {
        return true;
    }

    public int getComparisonOrder() {
        return 2;
    }

    public String getName() {
        return "PrivatePermutations";
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval != null && ! initialized ) {
            this.veWalker.getLogger().warn("Initializing...");
            initialize(eval);
            initialized = true;
        }

        if ( isGood(eval) ) {
            if ( comp != null && ! comp.isFiltered() ) {
                return null;
            }

            int order_offset = 0;
            for ( String[] ordering : permutations ) {
                int sample_offset = 0;
                for ( String sample : ordering ) {
                    if ( eval.getGenotype(sample).isHet() || eval.getGenotype(sample).isHomVar() ) {
                        break;
                    }
                    sample_offset ++;
                }

                permuteCounts.additionalValue[order_offset][sample_offset]++;
                order_offset++;
            }
        } else {
            skipped++;    
        }

        return null;
    }

    private boolean isGood(VariantContext vc) {
        if ( vc == null || vc.isFiltered() || (vc.getHetCount() + vc.getHomVarCount() == 0) ) { // todo -- should be is variant, but need to ensure no alt alleles at ref sites
            return false;
        } else {
            Collection<Genotype> gtypes = vc.getGenotypes().values();
            int ngood = 0;
            for ( Genotype g : gtypes) {
                if ( g.getPhredScaledQual() >= LOW_GQ_THRSH ) {
                    ngood ++;
                }
            }

            return ( (0.0+ngood)/(0.0+gtypes.size()) >= LOW_GQ_PCT );
        }
    }

    public PrivatePermutations(VariantEvalWalker parent) {
        super(parent);
    }

    public void initialize(VariantContext vc) {
        Set<String> permuteSamples = vc.getSampleNames();
        permutations = new String[NUM_PERMUTATIONS][permuteSamples.size()];
        veWalker.getLogger().warn(String.format("Num samples: %d",permuteSamples.size()));
        int offset = 0;
        for ( String s : permuteSamples ) {
            permutations[0][offset] = s;
            offset ++;
        }
        
        for ( int p = 1; p < NUM_PERMUTATIONS ; p++ ) {
            permutations[p] = permutations[0].clone();
            for ( int o = 0; o < permutations[p].length; o ++ ) {
                int r = (int) Math.floor(Math.random()*(o+1));
                String swap = permutations[p][r];
                permutations[p][r] = permutations[p][o];
                permutations[p][o] = swap;
            }
        }

        permuteCounts = new AdditionalBySample();
        permuteCounts.additionalValue = new int[NUM_PERMUTATIONS][permuteSamples.size()];
    }

    class AdditionalBySample implements TableType {
        int[][] additionalValue;
        //String[][] permutationNames;
        String[] rowKeys = null;
        String[] colKeys = null;

        public Object[] getRowKeys() {
            if ( rowKeys == null ) {
                rowKeys = new String[additionalValue.length];
                for ( int i = 0; i < additionalValue.length; i ++ ) {
                    rowKeys[i] = String.format("%s%d","P",i);
                }
            }


            return rowKeys;
        }

        public String getCell(int x, int y) {
            return String.format("%d",additionalValue[x][y]);
        }

        public String getName() { return "Marginal Number of Mutations"; }

        public Object[] getColumnKeys() {
            if ( colKeys == null ) {
                colKeys = new String[additionalValue[0].length];
                for ( int i = 0; i < additionalValue[0].length; i ++ ) {
                    colKeys[i] = String.format("%s%d","S",i);
                }
            }

            return colKeys;
        }
    }

    public void finalizeEvaluation() {
        veWalker.getLogger().info(String.format("Skipped: %d",skipped));    
    }
}
