package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.io.*;

@By(DataSource.REFERENCE)
@Requires(DataSource.REFERENCE)
@Allows(DataSource.REFERENCE)
public class VariantEvalWalker extends RefWalker<Integer, Integer> {
    @Argument(shortName="minDiscoveryQ", doc="Phred-scaled minimum LOD to consider an evaluation SNP a variant", required=false)
    public int minDiscoveryQ = -1;

    @Argument(shortName="printVariants", doc="If true, prints the variants in all of the variant tracks that are examined", required=false)
    public boolean printVariants = false;

    @Argument(shortName="badHWEThreshold", doc="XXX", required=false)
    public double badHWEThreshold = 0.001;

    String analysisFilenameBase = null;

    String COMMENT_STRING = "";

    HashMap<String, ArrayList<VariantAnalysis>> analysisSets;

    final String ALL_SNPS = "all";
    final String SINGLETON_SNPS = "singletons";
    final String TWOHIT_SNPS = "2plus_hit";
    final String[] ALL_ANALYSIS_NAMES = { ALL_SNPS, SINGLETON_SNPS, TWOHIT_SNPS };

    public void initialize() {
        // setup the path to the analysis
        if ( this.getToolkit().getArguments().outFileName != null ) {
            analysisFilenameBase = this.getToolkit().getArguments().outFileName + ".analysis.";
        }

        analysisSets = new HashMap<String, ArrayList<VariantAnalysis>>();
        for ( String setName : ALL_ANALYSIS_NAMES ) {
            analysisSets.put(setName, initializeAnalysisSet(setName));
        }
    }

    private ArrayList<VariantAnalysis> getAnalysisSet(final String name) {
        return analysisSets.get(name);
    }

    private ArrayList<VariantAnalysis> initializeAnalysisSet(final String setName) {
        ArrayList<VariantAnalysis> analyses = new ArrayList<VariantAnalysis>();

        //
        // Add new analyzes here!
        //
        analyses.add(new VariantCounter());
        analyses.add(new VariantDBCoverage("dbSNP"));
        analyses.add(new TransitionTranversionAnalysis());
        analyses.add(new PairwiseDistanceAnalysis());
        analyses.add(new HardyWeinbergEquilibrium(badHWEThreshold));


        if ( printVariants ) analyses.add(new VariantMatcher("dbSNP"));

        for ( VariantAnalysis analysis : analyses ) {
            analysis.initialize(this, openAnalysisOutputStream(setName, analysis));
        }

        return analyses;
    }

    /**
     * Returns the filename of the analysis output file where output for an analysis with
     * @param name
     * @param params
     * @return
     */
    public String getAnalysisFilename(final String name, final List<String> params) {
        if ( analysisFilenameBase == null )
            return null;
        else
            return analysisFilenameBase + Utils.join(".", Utils.cons(name, params));
    }

    public PrintStream openAnalysisOutputStream(final String setName, VariantAnalysis analysis) {
        final String filename = getAnalysisFilename(setName + "." + analysis.getName(), analysis.getParams());
        if ( filename == null )
            return out;
        else {
            File file = new File(filename);
            try {
                return new PrintStream(new FileOutputStream(file));
            } catch (FileNotFoundException e) {
                throw new StingException("Couldn't open analysis output file " + filename, e);
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        // Iterate over each analysis, and update it
        AllelicVariant eval = (AllelicVariant)tracker.lookup("eval", null);

        //if ( eval!=null ) System.out.printf("Eval: %f %d %b%n", eval.getVariationConfidence(), minDiscoveryQ, eval.getVariationConfidence() < minDiscoveryQ);
        if ( eval != null && eval.getVariationConfidence() < minDiscoveryQ )
            eval = null;

        updateAnalysisSet(ALL_SNPS, eval, tracker, ref, context);

        if ( eval instanceof SNPCallFromGenotypes ) {
            SNPCallFromGenotypes call = (SNPCallFromGenotypes)eval;
            int nVarGenotypes = call.nHetGenotypes() + call.nHomVarGenotypes();
            //System.out.printf("%d variant genotypes at %s%n", nVarGenotypes, calls);
            final String s = nVarGenotypes == 1 ? SINGLETON_SNPS : TWOHIT_SNPS;
            updateAnalysisSet(s, eval, tracker, ref, context);
        }

        return 1;
    }


    public void updateAnalysisSet(final String analysisSetName, AllelicVariant eval,
                                  RefMetaDataTracker tracker, char ref, LocusContext context) {
        // Iterate over each analysis, and update it
        for ( VariantAnalysis analysis : getAnalysisSet(analysisSetName) ) {
            String s = analysis.update(eval, tracker, ref, context);
            if ( s != null ) analysis.getPrintStream().println(s);
        }
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return treeReduce(sum,value);
    }
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer result) {
        for ( String analysisSetName : ALL_ANALYSIS_NAMES ) {
            printAnalysisSet(analysisSetName);
        }
    }

    private void printAnalysisSet( final String analysisSetName ) {
        Date now = new Date();
        for ( VariantAnalysis analysis : getAnalysisSet(analysisSetName) ) {
            PrintStream stream = analysis.getPrintStream(); // getAnalysisOutputStream(analysisSetName + "." + analysis.getName(), analysis.getParams());
            stream.printf("%s%s%n", COMMENT_STRING, Utils.dupString('-', 78));
            stream.printf("%sAnalysis set       %s%n", COMMENT_STRING, analysisSetName);
            stream.printf("%sAnalysis name      %s%n", COMMENT_STRING, analysis.getName());
            stream.printf("%sAnalysis params    %s%n", COMMENT_STRING, Utils.join(" ", analysis.getParams()));
            stream.printf("%sAnalysis class     %s%n", COMMENT_STRING, analysis );
            stream.printf("%sAnalysis time      %s%n", COMMENT_STRING, now );
            for ( String line : analysis.done()) {
                stream.printf("%s  %s%n", COMMENT_STRING, line);
            }
        }
    }
}
