package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;
import org.broadinstitute.sting.utils.genotype.glf.SinglePointCall;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.io.*;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
@By(DataSource.REFERENCE)
@Requires(DataSource.REFERENCE)
@Allows(DataSource.REFERENCE)
public class VariantEvalWalker extends RefWalker<Integer, Integer> {
    @Argument(shortName="minDiscoveryQ", doc="Phred-scaled minimum LOD to consider an evaluation SNP a variant", required=false)
    public int minDiscoveryQ = -1;

    @Argument(shortName="printVariants", doc="If true, prints the variants in all of the variant tracks that are examined", required=false)
    public boolean printVariants = false;

    @Argument(shortName="badHWEThreshold", doc="XXX", required=false)
    public double badHWEThreshold = 1e-3;

    @Argument(shortName="evalContainsGenotypes", doc="If true, the input list of variants will be treated as a genotyping file, containing assertions of actual genotype values for a particular person.  Analyses that only make sense on at the population level will be disabled, while those operating on genotypes will be enabled", required=false)
    public boolean evalContainsGenotypes = false;

    String analysisFilenameBase = null;

    String COMMENT_STRING = "";

    final String knownSNPDBName = "dbSNP";
    final String genotypeChipName = "hapmap-chip";

    HashMap<String, ArrayList<VariantAnalysis>> analysisSets;

    long nSites = 0;

    final String ALL_SNPS = "all";
    final String SINGLETON_SNPS = "singletons";
    final String TWOHIT_SNPS = "2plus_hit";
    final String KNOWN_SNPS = "known";
    final String NOVEL_SNPS = "novel";
    final String[] POPULATION_ANALYSIS_NAMES = { ALL_SNPS, SINGLETON_SNPS, TWOHIT_SNPS, KNOWN_SNPS, NOVEL_SNPS };
    final String[] GENOTYPE_ANALYSIS_NAMES = { ALL_SNPS, KNOWN_SNPS, NOVEL_SNPS };
    String[] ALL_ANALYSIS_NAMES = null;

    public void initialize() {
        ALL_ANALYSIS_NAMES = evalContainsGenotypes ? GENOTYPE_ANALYSIS_NAMES : POPULATION_ANALYSIS_NAMES;

        // setup the path to the analysis
        if ( this.getToolkit().getArguments().outFileName != null ) {
            analysisFilenameBase = this.getToolkit().getArguments().outFileName + "."; // + ".analysis.";
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
        analyses.add(new VariantDBCoverage(knownSNPDBName));
        analyses.add(new GenotypeConcordance(genotypeChipName));
        analyses.add(new TransitionTranversionAnalysis());
        analyses.add(new NeighborDistanceAnalysis());
        analyses.add(new HardyWeinbergEquilibrium(badHWEThreshold));
        analyses.add(new ClusterCounterAnalysis());
        analyses.add(new CallableBasesAnalysis());

        //
        // Filter out analyzes inappropriate for our evaluation type Population or Genotype
        //
        Iterator<VariantAnalysis> iter = analyses.iterator();
        while ( iter.hasNext() ) {
            VariantAnalysis analysis = iter.next();
            boolean disableForGenotyping = evalContainsGenotypes && ! (analysis instanceof GenotypeAnalysis);
            boolean disableForPopulation = ! evalContainsGenotypes && ! (analysis instanceof PopulationAnalysis);
            boolean disable = disableForGenotyping | disableForPopulation;
            String causeName = disableForGenotyping ? "population" : (disableForPopulation ? "genotype" : null);
            if ( disable ) {
                logger.info(String.format("Disabling %s-only analysis %s in set %s", causeName, analysis, setName));
                iter.remove();
            }
        }


        if ( printVariants ) analyses.add(new VariantMatcher(knownSNPDBName));

        for ( VariantAnalysis analysis : analyses ) {
            initializeAnalysisOutputStream(setName, analysis);
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

    public void initializeAnalysisOutputStream(final String setName, VariantAnalysis analysis) {
        final String filename = getAnalysisFilename(setName + "." + analysis.getName(), analysis.getParams());
        if ( filename == null )
            analysis.initialize(this, out, filename);
        else {
            File file = new File(filename);
            try {
                analysis.initialize(this, new PrintStream(new FileOutputStream(file)), filename);
            } catch (FileNotFoundException e) {
                throw new StingException("Couldn't open analysis output file " + filename, e);
            }
        }
    }

    // OMG this is painful
    private rodVariants glf2geli(char ref, RodGLF glf) {
        SinglePointCall rec = (SinglePointCall)glf.mRecord;
        // contig pos refBase depth maxMappingQ bestGenotype lodBtr lodBtnb genotypes
        Integer[] sorted = Utils.SortPermutation(rec.getLikelihoods());
        int bestIndex = sorted[0];
        char[] refs = {ref, ref};
        String homRef = new String(refs);
        int refIndex = rodVariants.Genotype.valueOf(homRef).ordinal();
        double refLikelihood = rec.getLikelihoods()[refIndex];
        double bestLikelihood = rec.getLikelihoods()[bestIndex];
        double secondBestLikelihood = rec.getLikelihoods()[sorted[1]];

        rodVariants var = new rodVariants("eval");
        var.loc = glf.getLocation();
        var.refBase = ref;
        var.depth = rec.getReadDepth();
        var.maxMappingQuality = rec.getRmsMapQ();
        var.bestGenotype = rodVariants.Genotype.values()[bestIndex].toString();
        var.lodBtr = Math.abs((bestLikelihood - refLikelihood) / GLFRecord.LIKELIHOOD_SCALE_FACTOR);
        var.lodBtnb = Math.abs((bestLikelihood - secondBestLikelihood) / GLFRecord.LIKELIHOOD_SCALE_FACTOR);
        var.genotypeLikelihoods = rec.getLikelihoods();
        for ( int i = 0; i < var.genotypeLikelihoods.length; i++ )
            var.genotypeLikelihoods[i] /= GLFRecord.LIKELIHOOD_SCALE_FACTOR;

        if ( false ) {
            System.out.printf("Converting : %s%n", glf);
            System.out.printf("    homRef: %s%n", homRef);
            System.out.printf("    refindex : %d%n", refIndex);
            System.out.printf("    bestIndex : %d%n", sorted[0]);
            System.out.printf("    2ndindex  : %d%n", sorted[1]);
            System.out.printf("    => %s%n", var);
        }
        
        return var;
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        nSites++;

        if ( tracker.lookup("eval", null) instanceof RodGLF) {
            tracker.bind("eval", glf2geli(ref, (RodGLF)tracker.lookup("eval", null)));
        }

        // Iterate over each analysis, and update it
        AllelicVariant eval = (AllelicVariant)tracker.lookup("eval", null);

        //if ( eval!=null ) System.out.printf("Eval: %f %d %b%n", eval.getVariationConfidence(), minDiscoveryQ, eval.getVariationConfidence() < minDiscoveryQ);
        if ( eval != null && eval.getVariationConfidence() < minDiscoveryQ )
            eval = null;

        // update stats about all of the SNPs
        updateAnalysisSet(ALL_SNPS, eval, tracker, ref, context);

        // update the known / novel set by checking whether the knownSNPDBName track has an entry here
        if ( eval != null ) {
            AllelicVariant dbsnp = (AllelicVariant)tracker.lookup(knownSNPDBName, null);
            String noveltySet = dbsnp == null ? NOVEL_SNPS : KNOWN_SNPS;
            updateAnalysisSet(noveltySet, eval, tracker, ref, context);
        }

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
            if ( s != null ) analysis.getCallPrintStream().println(s);
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
        out.printf("Writing analysis set %s", analysisSetName);
        Date now = new Date();
        for ( VariantAnalysis analysis : getAnalysisSet(analysisSetName) ) {
            analysis.finalize(nSites);
            PrintStream stream = analysis.getSummaryPrintStream();
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
