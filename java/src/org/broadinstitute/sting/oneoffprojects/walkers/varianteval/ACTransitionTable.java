package org.broadinstitute.sting.oneoffprojects.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvaluator;
import org.broadinstitute.sting.utils.collections.Pair;
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
@Analysis(name = "ACTransitionMatrix", description = "Number of additional genotypes from each new sample; random permutations")
public class ACTransitionTable extends VariantEvaluator {
    private final int NUM_PERMUTATIONS = 50;
    private final double LOW_GQ_PCT = 0.95;
    private final double LOW_GQ_THRSH = 30.0;
    private boolean initialized = false;
    private long skipped = 0l;

    @DataPoint(name="Het transitions",description="AC[s] = AC[s-1]+1 and AC[s] = AC[s-1]+2 transitions")
    TransitionTable transitions = null;
    @DataPoint(name="Private permutations",description="Marginal increase in number of sites per sample")
    PermutationCounts privatePermutations;
    @DataPoint(name="AC2 Permutations",description="Marginal increase in number of AC=2 sites, per sample")
    PermutationCounts doubletonPermutations;
    @DataPoint(name="AC3 Permutations",description="Marginal increase in number of tripleton sites, per sample")
    PermutationCounts tripletonPermutations;

    String[][] permutations;

    public boolean enabled() {
        return true;
    }

    public int getComparisonOrder() {
        return 2;
    }

    public String getName() {
        return "ACTransitionTable";
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
                int variant_ac = 0;
                for ( String sample : ordering ) {
                    if ( eval.getGenotype(sample).isHet() ) {
                        variant_ac++;
                        transitions.hetTransitionCounts[order_offset][variant_ac-1][sample_offset]++;
                    } else if ( eval.getGenotype(sample).isHomVar() ) {
                        variant_ac += 2;
                        transitions.homTransitionCounts[order_offset][variant_ac-1][sample_offset]++;
                    } else {
                        // todo -- note, unclear how to treat no calls. Is the hom in het,ref,ref,nocall,hom sample 4 or 5?
                        // todo -- do we want to tabulate P[sample i is not variant | some variant]? This is just combinatorics so i left it out
                        if ( variant_ac > 0 ) {
                            transitions.stationaryCounts[order_offset][variant_ac-1][sample_offset]++;
                        }
                    }
                    sample_offset ++;
                }
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
                if ( g.isCalled() && g.getPhredScaledQual() >= LOW_GQ_THRSH ) {
                    ngood ++;
                }
            }

            return ( (0.0+ngood)/(0.0+gtypes.size()) >= LOW_GQ_PCT );
        }
    }

    public ACTransitionTable(VariantEvalWalker parent) {
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

        transitions = new TransitionTable();
        transitions.hetTransitionCounts = new int[NUM_PERMUTATIONS][permuteSamples.size()*2][permuteSamples.size()];
        transitions.homTransitionCounts = new int[NUM_PERMUTATIONS][permuteSamples.size()*2][permuteSamples.size()];
        transitions.stationaryCounts = new int[NUM_PERMUTATIONS][permuteSamples.size()*2][permuteSamples.size()];
        privatePermutations = new PermutationCounts(1,transitions);
        doubletonPermutations = new PermutationCounts(2,transitions);
        tripletonPermutations = new PermutationCounts(3,transitions);
    }

    public void finalizeEvaluation() { // note: data points are null when this is called (wtf?)
        veWalker.getLogger().info(String.format("Skipped: %d",skipped));
    }

    class TransitionTable implements TableType {
        int[][][] hetTransitionCounts;
        int[][][] homTransitionCounts;
        int[][][] stationaryCounts;
        String[][] countAverages;
        String[] rowKeys = null;
        String[] colKeys = null;

        public Object[] getRowKeys() {
            if ( rowKeys == null ) {
                rowKeys = new String[3*hetTransitionCounts[0].length];
                for ( int i = 0; i < hetTransitionCounts[0].length; i ++ ) {
                    rowKeys[i] = String.format("%s%d%s","AC_",i,"_(het)");
                }
                for ( int i = 0; i < hetTransitionCounts[0].length; i ++ ) {
                    rowKeys[hetTransitionCounts[0].length+i] = String.format("%s%d%s","AC_",i,"_(hom)");
                }
                for ( int i = 0; i < hetTransitionCounts[0].length; i ++ ) {
                    rowKeys[2*hetTransitionCounts[0].length+i] = String.format("%s%d%s","AC_",i,"_(ref)");
                }
            }


            return rowKeys;
        }

        public String getCell(int x, int y) {
            if ( countAverages == null ) {
                countAverages = new String[hetTransitionCounts[0].length*3][hetTransitionCounts[0][0].length];
                for ( int sam = 0; sam < hetTransitionCounts[0][0].length; sam ++) {
                    for ( int idx = 0 ; idx < hetTransitionCounts[0].length; idx ++ ) {
                        int totalTimesAtACSample = 0;
                        int totalStationary = 0;
                        int totalAC1Shift = 0;
                        int totalAC2Shift = 0;
                        for ( int p = 0; p < hetTransitionCounts.length; p++ ) {
                            totalStationary += stationaryCounts[p][idx][sam];
                            totalAC2Shift += (idx+2 >= hetTransitionCounts[0][0].length) ? 0 : homTransitionCounts[p][idx+2][sam];
                            totalAC1Shift += (idx+1 >= hetTransitionCounts[0][0].length) ? 0 : hetTransitionCounts[p][idx+1][sam];
                        }
                        totalTimesAtACSample = totalStationary+totalAC1Shift+totalAC2Shift;
                        countAverages[idx][sam] = formatProp(totalAC1Shift,totalTimesAtACSample);
                        countAverages[hetTransitionCounts[0].length+idx][sam] = formatProp(totalAC2Shift,totalTimesAtACSample);
                        countAverages[hetTransitionCounts[0].length*2+idx][sam] = formatProp(totalStationary,totalTimesAtACSample);
                    }
                }
            }

            return countAverages[x][y] == null ? "0.00" : countAverages[x][y];
        }

        private String formatProp(int num, int denom) {
            return (denom != 0) ? String.format("%.4f", ((double) num)/denom) : "0.0";
        }

        public String getName() { return "AC Transition Tables"; }

        public Object[] getColumnKeys() {
            if ( colKeys == null ) {
                colKeys = new String[hetTransitionCounts[0][0].length];
                for ( int ac = 0; ac < hetTransitionCounts[0][0].length; ac ++ ) {
                    colKeys[ac] = String.format("Sample_%d",ac);
                }
            }

            return colKeys;
        }
    }


    class PermutationCounts implements TableType {
        int acToExtract;
        TransitionTable table;
        String[] rowNames;
        String[] colNames;

        public PermutationCounts(int ac, TransitionTable tTable) {
            acToExtract = ac;
            table = tTable;
        }

        public String[] getRowKeys() {
            //System.out.printf("%s%n",table);
            if ( rowNames == null ) {
                rowNames = new String[table.stationaryCounts.length];
                for ( int p = 0 ; p < rowNames.length; p ++ ) {
                    rowNames[p] = String.format("Perm%d",p+1);
                }
            }

            return rowNames;
        }

        public String[] getColumnKeys() {
            if ( colNames == null ) {
                colNames = new String[table.stationaryCounts[0][0].length];
                for ( int s = 0 ; s < colNames.length; s ++ ) {
                    colNames[s] = String.format("Sample%d",s+1);
                }
            }

            return colNames;
        }

        public Integer getCell(int x, int y) {
            return table.hetTransitionCounts[x][acToExtract-1][y] +
                    ( (acToExtract > table.homTransitionCounts[0][0].length) ? 0 : table.homTransitionCounts[x][acToExtract-1][y]);
        }

        public String getName() {
            return String.format("PermutationCountsAC%d",acToExtract);
        }

        public void init() {
            getRowKeys();
            getColumnKeys();
            getCell(1,1);
        }
    }


}

