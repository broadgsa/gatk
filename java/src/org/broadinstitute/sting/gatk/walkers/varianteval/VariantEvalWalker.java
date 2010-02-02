package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeEncoding;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import java.util.regex.Pattern;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * A robust and general purpose tool for characterizing the quality of SNPs, Indels, and other variants that includes basic
 * counting, ti/tv, dbSNP% (if -D is provided), concordance to chip or validation data, and will show interesting sites (-V)
 * that are FNs, FP, etc.
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="eval",type=ReferenceOrderedDatum.class)}) // right now we have no base variant class for rods, this should change
//@Allows(value={DataSource.REFERENCE},referenceMetaData = {@RMD(name="eval",type=ReferenceOrderedDatum.class), @RMD(name="dbsnp",type=rodDbSNP.class),@RMD(name="hapmap-chip",type=ReferenceOrderedDatum.class), @RMD(name="interval",type=IntervalRod.class), @RMD(name="validation",type=RodGenotypeChipAsGFF.class)})
//public class VariantEvalWalker extends RefWalker<Integer, Integer> {
public class VariantEvalWalker extends RodWalker<Integer, Integer> {
    @Argument(shortName="minPhredConfidenceScore", doc="Minimum confidence score to consider an evaluation SNP a variant", required=false)
    public double minConfidenceScore = -1.0;

    @Argument(shortName="printVariants", doc="If true, prints the variants in all of the variant tracks that are examined", required=false)
    public boolean printVariants = false;

    @Argument(shortName="badHWEThreshold", doc="Only sites with deviations froim Hardy-Weinberg equilibrium with P-values < than this threshold are flagged", required=false)
    public double badHWEThreshold = 1e-3;

    @Argument(fullName="evalContainsGenotypes", shortName = "G", doc="If true, the input list of variants will be treated as a genotyping file, containing assertions of actual genotype values for a particular person.  Analyses that only make sense on at the population level will be disabled, while those operating on genotypes will be enabled", required=false)
    public boolean evalContainsGenotypes = false;

    @Argument(fullName="explode", shortName = "E", doc="Old style formatting, with each analysis split into separate files.", required=false)
    public boolean explode = false;

    @Argument(fullName="includeViolations", shortName = "V", doc="If provided, violations will be written out along with summary information", required=false)
    public boolean mIncludeViolations = false;

    @Argument(fullName="extensiveSubsets", shortName = "A", doc="If provided, output will be calculated over a lot of subsets, by default we only operate over all variants", required=false)
    public boolean extensiveSubsets = false;

    @Argument(fullName="supressDateInformation", doc="This flag indicates that we want to suppress the date information from the output, so that if can be diff'ed against previous evals.", required=false)
    public boolean supressDateInformation = false;

    @Argument(fullName="includeFilteredRecords", doc="If true, variation record with filter fields at are true will be included in the analysis", required=false)
    public boolean includeFilteredRecords = false;


    @Argument(fullName = "samplesFile", shortName="samples", doc="When running an analysis on one or more individuals with truth data, this field provides a filepath to the listing of which samples are used (and are used to name corresponding rods with -B)", required=false)
    public String samplesFile = null;

    @Argument(fullName = "sampleName", shortName="sampleName", doc="When running an analysis on one or more individuals with truth data, provide this parameter to only analyze the genotype of this sample", required=false)
    public String sampleName = null;
    @Argument(fullName = "vcfInfoSelector", shortName="vcfInfoSelector", doc="When running an analysis on one or more individuals with truth data, provide this parameter to only analyze the genotype of this sample", required=false)
    public String vcfInfoSelector = null;

    String analysisFilenameBase = null;

    final String knownSNPDBName = "dbSNP";
    final String One1KGSNPNames = "1kg";
    final String genotypeChipName = "hapmap-chip";

    HashMap<ANALYSIS_TYPE, ArrayList<VariantAnalysis>> analysisSets;

    PrintStream perLocusStream = null;

    long nMappedSites = 0;

    // the types of analysis we support, and the string tags we associate with the enumerated value
    enum ANALYSIS_TYPE {
	// todo -- differeniate into three classes -- snps at known snp sites, snps not at known snp site but covered by known indel, and novel 
        ALL_SNPS("all"),
        SINGLETON_SNPS("singletons"),
        TWOHIT_SNPS("2plus_hit"),
        KNOWN_SNPS("known"),
        SNPS_AT_NON_SNP_SITES("snp_at_known_non_snps"),
        NOVEL_SNPS("novel"),
        FILTERED_SNPS("filtered");

        private final String value;
        ANALYSIS_TYPE(String value) { this.value = value;}

        public String toString() { return value; }

    }

//    final ANALYSIS_TYPE[] POPULATION_ANALYSIS_NAMES = {ANALYSIS_TYPE.ALL_SNPS,
//            ANALYSIS_TYPE.SINGLETON_SNPS,
//            ANALYSIS_TYPE.TWOHIT_SNPS,
//            ANALYSIS_TYPE.KNOWN_SNPS,
//            ANALYSIS_TYPE.SNPS_AT_NON_SNP_SITES,
//            ANALYSIS_TYPE.NOVEL_SNPS};
//    final ANALYSIS_TYPE[] GENOTYPE_ANALYSIS_NAMES = {ANALYSIS_TYPE.ALL_SNPS,
//            ANALYSIS_TYPE.KNOWN_SNPS,
//            ANALYSIS_TYPE.SNPS_AT_NON_SNP_SITES,
//            ANALYSIS_TYPE.NOVEL_SNPS};

    final ANALYSIS_TYPE[] POPULATION_ANALYSIS_NAMES = {ANALYSIS_TYPE.ALL_SNPS,
            ANALYSIS_TYPE.SINGLETON_SNPS,
            ANALYSIS_TYPE.TWOHIT_SNPS,
            ANALYSIS_TYPE.KNOWN_SNPS,
            ANALYSIS_TYPE.SNPS_AT_NON_SNP_SITES,
            ANALYSIS_TYPE.NOVEL_SNPS,
            ANALYSIS_TYPE.FILTERED_SNPS};
    final ANALYSIS_TYPE[] GENOTYPE_ANALYSIS_NAMES = {ANALYSIS_TYPE.ALL_SNPS,
            ANALYSIS_TYPE.KNOWN_SNPS,
            ANALYSIS_TYPE.SNPS_AT_NON_SNP_SITES,
            ANALYSIS_TYPE.NOVEL_SNPS,
            ANALYSIS_TYPE.FILTERED_SNPS};

    final ANALYSIS_TYPE[] SIMPLE_ANALYSIS_NAMES = {ANALYSIS_TYPE.ALL_SNPS};
    ANALYSIS_TYPE[] ALL_ANALYSIS_NAMES = null;

    public void initialize() {
        ALL_ANALYSIS_NAMES = SIMPLE_ANALYSIS_NAMES;
        if (extensiveSubsets)
            ALL_ANALYSIS_NAMES = evalContainsGenotypes ? GENOTYPE_ANALYSIS_NAMES : POPULATION_ANALYSIS_NAMES;

        // setup the path to the analysis
        if (this.getToolkit().getArguments().outFileName != null) {
            analysisFilenameBase = this.getToolkit().getArguments().outFileName + "."; // + ".analysis.";
        }

        analysisSets = new HashMap<ANALYSIS_TYPE, ArrayList<VariantAnalysis>>();
        for (ANALYSIS_TYPE setName : ALL_ANALYSIS_NAMES) {
            analysisSets.put(setName, initializeAnalysisSet(setName));
        }
        // THIS IS A HACK required in order to reproduce the behavior of old (and imperfect) RODIterator and
        // hence to pass the integration test. The new iterator this code is now using does see ALL the SNPs,
        // whether masked by overlapping indels/other events or not.
        //TODO process correctly all the returned dbSNP rods at each location
        //BrokenRODSimulator.attach("dbSNP");
    }


    public long getNMappedSites() {
        return nMappedSites;
    }

    private ArrayList<VariantAnalysis> getAnalysisSet(final ANALYSIS_TYPE name) {
        return analysisSets.containsKey(name) ? analysisSets.get(name) : null;
    }

    private ArrayList<VariantAnalysis> initializeAnalysisSet(final ANALYSIS_TYPE setName) {
        ArrayList<VariantAnalysis> analyses = new ArrayList<VariantAnalysis>();

        //
        // Add new analyses here!
        //
        analyses.add(new ValidationDataAnalysis());

        analyses.add(new VariantCounter());
        analyses.add(new VariantDBCoverage(knownSNPDBName));
        analyses.add(new VariantDBCoverage(One1KGSNPNames));

        if ( samplesFile != null ) {
            //if ( numPeopleInPool < 1 )
            //    analyses.add(new GenotypeConcordance(samplesFile, true));
            //else
            analyses.add(new PooledConcordance(samplesFile, true));
        } else {
            analyses.add(new GenotypeConcordance(genotypeChipName, false));
        }

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
        while (iter.hasNext()) {
            VariantAnalysis analysis = iter.next();
            VariantAnalysis removed = null;
            if ( evalContainsGenotypes ) {
                if ( ! (analysis instanceof GenotypeAnalysis ) ) {
                    removed = analysis;
                    iter.remove();
                }
            } else if ( ! (analysis instanceof PopulationAnalysis || analysis instanceof PoolAnalysis) ) {
                removed = analysis;
                iter.remove();
//            } else if ( numPeopleInPool > 1 && ! ( analysis instanceof PoolAnalysis ) ) {
            } else if ( analysis instanceof PoolAnalysis && samplesFile == null ) {
                removed = analysis;
                iter.remove();
            }

            if ( removed != null ) {
                logger.info(String.format("Disabling analysis %s in set %s", removed, setName));
            }
        }


        if (printVariants) analyses.add(new VariantMatcher(knownSNPDBName));

        for (VariantAnalysis analysis : analyses) {
            initializeAnalysisOutputStream(setName, analysis);
        }

        return analyses;
    }

    /**
     * Returns the filename of the analysis output file where output for an analysis with
     *
     * @param name
     * @param params
     *
     * @return
     */
    public String getAnalysisFilename(final String name, final List<String> params) {
        if (analysisFilenameBase == null)
            return null;
        else
            return analysisFilenameBase + Utils.join(".", Utils.cons(name, params));
    }

    public void initializeAnalysisOutputStream(final ANALYSIS_TYPE setName, VariantAnalysis analysis) {
        final String filename = getAnalysisFilename(setName + "." + analysis.getName(), analysis.getParams());

        try {
            if (perLocusStream == null)
                perLocusStream = filename == null ? out : new PrintStream(new File(analysisFilenameBase + "interesting_sites"));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }

        if (filename == null || !explode)
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
        //System.out.printf("Tracker at %s is %s, ref is %s, skip is %d mapped is %d%n", context.getLocation(), tracker, ref, context.getSkippedBases(), nMappedSites);
        nMappedSites += context.getSkippedBases();

        //System.out.printf("Tracker at %s is %s, ref is %s%n", context.getLocation(), tracker, ref);
        //if ( ref == null )
        //    out.printf("Last position was %s: skipping %d bases%n",
        //        context.getLocation(), context.getSkippedBases() );
        if ( ref == null ) { // we are seeing the last site
            return 0;
        }

        nMappedSites++;

        int nBoundGoodRods = tracker.getNBoundRodTracks("interval");
        if (nBoundGoodRods > 0) {
            //System.out.printf("%s: n = %d%n", context.getLocation(), nBoundGoodRods );

            // Iterate over each analysis, and update it
            Variation eval = (Variation) tracker.lookup("eval", null);
            Variation evalForFilter = null;

            // ensure that the variation we're looking at is bi-allelic
            if ( eval != null && ! eval.isBiallelic() )
                eval = null;  

            if (eval != null)
                if (eval.getNegLog10PError() * 10.0 < minConfidenceScore) eval = null;

            if ( eval != null && (eval instanceof RodVCF) ) {
                if ( sampleName != null ) {
                    // code to grab a particular sample from a VCF
                    Variation evalOld = eval;
//                    System.out.printf("original is %s%n", evalOld);
                    eval = fakeVCFForSample((RodVCF)eval, ref, sampleName);
                    if ( eval != null && eval.isSNP() ) {
//                        System.out.printf("sample   is %s%n", eval);
//                        System.out.printf("Replacing %s with %s%n", evalOld, eval);
//                        System.out.printf("%n");
                    } else {
                        eval = null;
                    }
                }

                if ( vcfInfoSelector != null && eval != null ) {
                    String[] keyValue = vcfInfoSelector.split("=");
                    Map<String, String> map = ((RodVCF)eval).getRecord().getInfoValues();
                    if ( map.containsKey(keyValue[0]) && ! Pattern.matches(keyValue[1], map.get(keyValue[0]) ) )
                        eval = null;
                }

                if ( eval != null && ((RodVCF)eval).mCurrentRecord.isFiltered() ) {
                    evalForFilter = eval;
                    if ( ! includeFilteredRecords ) {
                        eval = null;    // we are not including filtered records, so set eval to null
                    }
                    //System.out.printf("Rejecting filtered record %s%n", eval);
                }
            }

            // update stats about all of the SNPs
            updateAnalysisSet(ANALYSIS_TYPE.ALL_SNPS, eval, tracker, ref.getBase(), context);

            // update the known / novel set by checking whether the knownSNPDBName track has an entry here
            if (eval != null) {
                ANALYSIS_TYPE noveltySet = getNovelAnalysisType(tracker);
                updateAnalysisSet(noveltySet, eval, tracker, ref.getBase(), context);
            }

            if (evalForFilter != null) {
                updateAnalysisSet(ANALYSIS_TYPE.FILTERED_SNPS, evalForFilter, tracker, ref.getBase(), context);
            }


            // are we a population backed call? then update
            if (eval instanceof SNPCallFromGenotypes) {
                SNPCallFromGenotypes call = (SNPCallFromGenotypes) eval;
                int nVarGenotypes = call.nHetGenotypes() + call.nHomVarGenotypes();
                //System.out.printf("%d variant genotypes at %s%n", nVarGenotypes, calls);
                final ANALYSIS_TYPE s = nVarGenotypes == 1 ? ANALYSIS_TYPE.SINGLETON_SNPS : ANALYSIS_TYPE.TWOHIT_SNPS;
                updateAnalysisSet(s, eval, tracker, ref.getBase(), context);
            }
        }

        return 1;
    }

    private RodVCF fakeVCFForSample(RodVCF eval, ReferenceContext ref, final String sampleName) {
        VCFGenotypeRecord genotype = (VCFGenotypeRecord)eval.getGenotype(sampleName);
        if ( genotype.getNegLog10PError() > 0 ) {
            VCFRecord record = new VCFRecord(ref.getBase(), ref.getLocus(), "GT");
            record.setAlternateBases(eval.getRecord().getAlternateAlleles());
            record.addGenotypeRecord(genotype);
            record.setQual(10*genotype.getNegLog10PError());
            record.setFilterString(eval.getFilterString());
            record.addInfoFields(eval.getInfoValues());
            return new RodVCF("fakeVCFForSample", record, eval.getReader());
        } else {
            return null;
        }
    }

    private ANALYSIS_TYPE getNovelAnalysisType(RefMetaDataTracker tracker) {
        RODRecordList<ReferenceOrderedDatum> dbsnpList = tracker.getTrackData("dbsnp", null);

        if (dbsnpList == null)
            return ANALYSIS_TYPE.NOVEL_SNPS;

        for (ReferenceOrderedDatum d : dbsnpList) {
            if (((rodDbSNP) d).isSNP()) {
                return ANALYSIS_TYPE.KNOWN_SNPS;
            }
        }

        return ANALYSIS_TYPE.SNPS_AT_NON_SNP_SITES;

        // old and busted way of doing this
//        Variation dbsnp = (Variation) BrokenRODSimulator.simulate_lookup("dbSNP", ref.getLocus(), tracker);
//
//        ANALYSIS_TYPE noveltySet = null;
//        if ( dbsnp == null ) noveltySet = ANALYSIS_TYPE.NOVEL_SNPS;
//        else if ( dbsnp.isSNP() ) noveltySet = ANALYSIS_TYPE.KNOWN_SNPS;
//        else {  // if ( dbsnp.isIndel() )
//            noveltySet = ANALYSIS_TYPE.SNPS_AT_NON_SNP_SITES;
    }

    public boolean includeViolations() { return mIncludeViolations; }

    public void updateAnalysisSet(final ANALYSIS_TYPE analysisSetName, Variation eval,
                                  RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        // Iterate over each analysis, and update it
        ArrayList<VariantAnalysis> set = getAnalysisSet(analysisSetName);
        if ( set != null ) {
            for ( VariantAnalysis analysis : set ) {
                String s = analysis.update(eval, tracker, ref, context);
                if ( s != null && includeViolations() ) {
                    analysis.getCallPrintStream().println(getLineHeader(analysisSetName, "flagged", analysis.getName()) + s);
                }
            }
        }
    }

    // Given result of map function
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return treeReduce(sum, value);
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer result) {
        for (ANALYSIS_TYPE analysisSetName : ALL_ANALYSIS_NAMES) {
            printAnalysisSet(analysisSetName);
        }
    }

    private String getLineHeader(final ANALYSIS_TYPE analysisSetName, final String keyword, final String analysis) {
        String s = Utils.join(",", Arrays.asList(analysisSetName, keyword, analysis));
        return s + Utils.dupString(' ', Math.max(50 - s.length(), 1));
    }

    private void printAnalysisSet(final ANALYSIS_TYPE analysisSetName) {
        //out.printf("Writing analysis set %s", analysisSetName);
        Date now = new Date();
        for (VariantAnalysis analysis : getAnalysisSet(analysisSetName)) {
            String header = getLineHeader(analysisSetName, "summary", analysis.getName());
            analysis.finalize(getNMappedSites());
            PrintStream stream = analysis.getSummaryPrintStream();
            stream.printf("%s%s%n", header, Utils.dupString('-', 78));
            //stream.printf("%s Analysis set       %s%n", analysisSetName, , analysisSetName);
            stream.printf("%sAnalysis name      %s%n", header, analysis.getName());
            stream.printf("%sAnalysis params    %s%n", header, Utils.join(" ", analysis.getParams()));
            stream.printf("%sAnalysis class     %s%n", header, analysis.getClass().getName());
            if (!supressDateInformation) stream.printf("%sAnalysis time      %s%n", header, now);
            for (String line : analysis.done()) {
                stream.printf("%s%s%n", header, line);
            }
        }
    }
}
