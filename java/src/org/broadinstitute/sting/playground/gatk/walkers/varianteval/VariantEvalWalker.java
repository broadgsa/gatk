package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

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
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="eval",type=VariationRod.class)}) // right now we have no base variant class for rods, this should change
@Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="eval",type=VariationRod.class), @RMD(name="dbsnp",type=rodDbSNP.class),@RMD(name="hapmap-chip",type=RodGenotypeChipAsGFF.class), @RMD(name="interval",type=IntervalRod.class)})
public class VariantEvalWalker extends RefWalker<Integer, Integer> {
    @Argument(shortName="minConfidenceScore", doc="Minimum confidence score to consider an evaluation SNP a variant", required=false)
    public int minConfidenceScore = -1;

    @Argument(shortName="printVariants", doc="If true, prints the variants in all of the variant tracks that are examined", required=false)
    public boolean printVariants = false;

    @Argument(shortName="badHWEThreshold", doc="XXX", required=false)
    public double badHWEThreshold = 1e-3;

    @Argument(fullName="evalContainsGenotypes", shortName = "G", doc="If true, the input list of variants will be treated as a genotyping file, containing assertions of actual genotype values for a particular person.  Analyses that only make sense on at the population level will be disabled, while those operating on genotypes will be enabled", required=false)
    public boolean evalContainsGenotypes = false;

    @Argument(fullName="explode", shortName = "E", doc="Old style formatting, with each analysis split into separate files.", required=false)
    public boolean explode = false;

    @Argument(fullName="includeViolations", shortName = "V", doc="If provided, violations will be written out along with summary information", required=false)
    public boolean includeViolations = false;

    @Argument(fullName="extensiveSubsets", shortName = "A", doc="If provided, output will be calculated over a lot of subsets, by default we only operate over all variants", required=false)
    public boolean extensiveSubsets = false;

    @Argument(fullName="supressDateInformation", doc="This flag indicates that we want to suppress the date information from the output, so that if can be diff'ed against previous evals.", required=false)
    public boolean supressDateInformation = false;

    @Argument(fullName = "numPeopleInPool", shortName="PS", doc="If using a variant file from a pooled caller, this field provides the number of individuals in each pool", required=false)
    public int numPeopleInPool = 1;

    @Argument(fullName = "pathToHapmapPoolFile", shortName="HPF", doc="If using a variant file from a pooled caller on pools of hapmap individuals, this field provides a filepath to the pool construction file listing which hapmap individuals are in which pool", required=false)
    public String pathToHapmapPoolFile = null;

    String analysisFilenameBase = null;

    final String knownSNPDBName = "dbSNP";
    final String genotypeChipName = "hapmap-chip";

    HashMap<String, ArrayList<VariantAnalysis>> analysisSets;

    PrintStream perLocusStream = null;

    long nSites = 0;

    final String ALL_SNPS = "all";
    final String SINGLETON_SNPS = "singletons";
    final String TWOHIT_SNPS = "2plus_hit";
    final String KNOWN_SNPS = "known";
    final String NOVEL_SNPS = "novel";
    final String[] POPULATION_ANALYSIS_NAMES = { ALL_SNPS, SINGLETON_SNPS, TWOHIT_SNPS, KNOWN_SNPS, NOVEL_SNPS };
    final String[] GENOTYPE_ANALYSIS_NAMES = { ALL_SNPS, KNOWN_SNPS, NOVEL_SNPS };
    final String[] SIMPLE_ANALYSIS_NAMES = { ALL_SNPS };
    String[] ALL_ANALYSIS_NAMES = null;


    public void initialize() {
        ALL_ANALYSIS_NAMES = SIMPLE_ANALYSIS_NAMES;
        if ( extensiveSubsets )
            ALL_ANALYSIS_NAMES = evalContainsGenotypes ? GENOTYPE_ANALYSIS_NAMES : POPULATION_ANALYSIS_NAMES;

        // setup the path to the analysis
        if ( this.getToolkit().getArguments().outFileName != null ) {
            analysisFilenameBase = this.getToolkit().getArguments().outFileName + "."; // + ".analysis.";
        }

        analysisSets = new HashMap<String, ArrayList<VariantAnalysis>>();
        for ( String setName : ALL_ANALYSIS_NAMES ) {
            analysisSets.put(setName, initializeAnalysisSet(setName));
        }
        // THIS IS A HACK required in order to reproduce the behavior of old (and imperfect) RODIterator and
        // hence to pass the integration test. The new iterator this code is now using does see ALL the SNPs,
        // whether masked by overlapping indels/other events or not.
        //TODO process correctly all the returned dbSNP rods at each location
        BrokenRODSimulator.attach("dbSNP");
    }

    private ArrayList<VariantAnalysis> getAnalysisSet(final String name) {
        return analysisSets.containsKey(name) ? analysisSets.get(name) : null;
    }

    private ArrayList<VariantAnalysis> initializeAnalysisSet(final String setName) {
        ArrayList<VariantAnalysis> analyses = new ArrayList<VariantAnalysis>();

        //
        // Add new analyses here!
        //
        analyses.add(new PooledGenotypeConcordance(pathToHapmapPoolFile));
        analyses.add(new VariantCounter());
        analyses.add(new VariantDBCoverage(knownSNPDBName));
        analyses.add(new GenotypeConcordance(genotypeChipName));
        analyses.add(new TransitionTranversionAnalysis());
        analyses.add(new NeighborDistanceAnalysis());
        analyses.add(new HardyWeinbergEquilibrium(badHWEThreshold));
        analyses.add(new ClusterCounterAnalysis());
        analyses.add(new CallableBasesAnalysis());
        analyses.add(new IndelMetricsAnalysis());

        //
        // Filter out analyses inappropriate for our evaluation type Population or Genotype
        //
        Iterator<VariantAnalysis> iter = analyses.iterator();
        while ( iter.hasNext() ) {
            VariantAnalysis analysis = iter.next();
            boolean disableForGenotyping = evalContainsGenotypes && ! (analysis instanceof GenotypeAnalysis);
            boolean disableForPopulation = ! evalContainsGenotypes && ! (analysis instanceof PopulationAnalysis);
            boolean disableForPools = pathToHapmapPoolFile == null;
            boolean disable = disableForGenotyping | disableForPopulation | disableForPools;
            String causeName = disableForGenotyping ? "population" : (disableForPopulation ? "genotype" : ( disableForPools ? "pool" : null ));
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

        try {
            if ( perLocusStream == null )
                perLocusStream = filename == null ? out : new PrintStream(new File(analysisFilenameBase + "interesting_sites"));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }

        if ( filename == null || ! explode )
            analysis.initialize(this, out, perLocusStream, filename);
        else {
            File file = new File(filename);
            try {
                analysis.initialize(this, new PrintStream(new FileOutputStream(file)), perLocusStream, filename);
            } catch (FileNotFoundException e) {
                throw new StingException("Couldn't open analysis output file " + filename, e);
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        nSites++;
        // Iterate over each analysis, and update it
        Variation eval = (Variation)tracker.lookup("eval", null);

        if ( eval != null )
               if ( eval.getNegLog10PError() < minConfidenceScore ) eval = null;

        // update stats about all of the SNPs
        updateAnalysisSet(ALL_SNPS, eval, tracker, ref.getBase(), context);

        // update the known / novel set by checking whether the knownSNPDBName track has an entry here
        if ( eval != null ) {
//            if ( ref.getLocus().getStart() >= 10168704 && ref.getLocus().getStop() <= 10168728) System.out.println("###DbSNP from MAP: ");
            Variation dbsnp = (Variation)BrokenRODSimulator.simulate_lookup("dbSNP",ref.getLocus(),tracker);
//            if ( ref.getLocus().getStart() >= 10168704 && ref.getLocus().getStop() <= 10168728) System.out.println("###\n");

//            RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData("dbSNP",null);

            //
            //TODO process correctly all the returned dbSNP rods at each location
//            if ( last_interval.containsP(ref.getLocus()) ) dbsnp = last_rod; // old RODIterator kept returning the same ROD until we completely walk out of it
//            else {
//                if ( rods != null && rods.size() > 0 ) dbsnp = (Variation)rods.getRecords().get(0);
//                if ( dbsnp != null ) {
//                     last_rod = dbsnp;
//                     last_interval = dbsnp.getLocation(); // remember what we just read
//                }
//            }

//            Variation dbsnp = (Variation)tracker.lookup(knownSNPDBName, null);
            String noveltySet = dbsnp == null ? NOVEL_SNPS : KNOWN_SNPS;
//            if ( dbsnp != null ) out.println(ref.getLocus()+" DBSNP RECORD "+dbsnp.getLocation());
            updateAnalysisSet(noveltySet, eval, tracker, ref.getBase(), context);
        }

        // are we a population backed call? then update
        if ( eval instanceof SNPCallFromGenotypes) {
            SNPCallFromGenotypes call = (SNPCallFromGenotypes)eval;
            int nVarGenotypes = call.nHetGenotypes() + call.nHomVarGenotypes();
            //System.out.printf("%d variant genotypes at %s%n", nVarGenotypes, calls);
            final String s = nVarGenotypes == 1 ? SINGLETON_SNPS : TWOHIT_SNPS;
            updateAnalysisSet(s, eval, tracker, ref.getBase(), context);
        }
        return 1;
    }


    public void updateAnalysisSet(final String analysisSetName, Variation eval,
                                  RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        // Iterate over each analysis, and update it
        if ( getAnalysisSet(analysisSetName) != null ) {
            for ( VariantAnalysis analysis : getAnalysisSet(analysisSetName) ) {
                String s = analysis.update(eval, tracker, ref, context);
                if ( s != null && includeViolations ) {
                    analysis.getCallPrintStream().println(getLineHeader(analysisSetName, "flagged", analysis.getName()) + s);
                }
            }
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

    private String getLineHeader( final String analysisSetName, final String keyword, final String analysis) {
        String s = Utils.join(",", Arrays.asList(analysisSetName, keyword, analysis));
        return s + Utils.dupString(' ', 50 - s.length());
    }

    private void printAnalysisSet( final String analysisSetName ) {
        //out.printf("Writing analysis set %s", analysisSetName);
        Date now = new Date();
        for ( VariantAnalysis analysis : getAnalysisSet(analysisSetName) ) {
            String header = getLineHeader(analysisSetName, "summary", analysis.getName());
            analysis.finalize(nSites);
            PrintStream stream = analysis.getSummaryPrintStream();
            stream.printf("%s%s%n", header, Utils.dupString('-', 78));
            //stream.printf("%s Analysis set       %s%n", analysisSetName, , analysisSetName);
            stream.printf("%sAnalysis name      %s%n", header, analysis.getName());
            stream.printf("%sAnalysis params    %s%n", header, Utils.join(" ", analysis.getParams()));
            stream.printf("%sAnalysis class     %s%n", header, analysis.getClass().getName());
            if (!supressDateInformation) stream.printf("%sAnalysis time      %s%n", header, now);
            for ( String line : analysis.done()) {
                stream.printf("%s%s%n", header, line);
            }
        }
    }
}
