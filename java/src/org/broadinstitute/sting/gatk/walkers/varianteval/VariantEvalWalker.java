package org.broadinstitute.sting.gatk.walkers.varianteval;

import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.gatk.walkers.varianteval.tags.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.*;
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.Tranche;
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import net.sf.picard.reference.ReferenceSequence;
import java.io.FileNotFoundException;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.*;

/**
 * General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)
 */
@Reference(window=@Window(start=-50, stop=50))
public class VariantEvalWalker extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    // Output arguments
    @Output
    protected PrintStream out;

    // Help arguments
    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit")
    protected Boolean LIST = false;

    // Partitioning the data arguments
    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", required=false)
    protected ArrayList<String> SELECT_EXPS = new ArrayList<String>();

    @Argument(shortName="selectName", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", required=false)
    protected ArrayList<String> SELECT_NAMES = new ArrayList<String>();

    @Argument(fullName="sample", shortName="sn", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", required=false)
    protected Set<String> SAMPLE_EXPRESSIONS;

    @Argument(shortName="knownName", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected String[] KNOWN_NAMES = {DbSNPHelper.STANDARD_DBSNP_TRACK_NAME};

    // Stratification arguments
    @Argument(fullName="stratificationModule", shortName="ST", doc="One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless -noS is specified)", required=false)
    protected String[] STRATIFICATIONS_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardStratifications", shortName="noST", doc="Do not use the standard stratification modules by default (instead, only those that are specified with the -S option)")
    protected Boolean NO_STANDARD_STRATIFICATIONS = false;

    @Argument(fullName="onlyVariantsOfType", shortName="VT", doc="If provided, only variants of these types will be considered during the evaluation, in ", required=false)
    protected Set<VariantContext.Type> typesToUse = null;

    // Evaluator arguments
    @Argument(fullName="evalModule", shortName="EV", doc="One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless -noE is specified)", required=false)
    protected String[] MODULES_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardModules", shortName="noEV", doc="Do not use the standard modules by default (instead, only those that are specified with the -E option)")
    protected Boolean NO_STANDARD_MODULES = false;

    // Other arguments
    @Argument(fullName="numSamples", shortName="ns", doc="Number of samples (used if no samples are available in the VCF file", required=false)
    protected Integer NUM_SAMPLES = 0;

    @Argument(fullName="minPhaseQuality", shortName="mpq", doc="Minimum phasing quality", required=false)
    protected double MIN_PHASE_QUALITY = 10.0;

    @Argument(shortName="family", doc="If provided, genotypes in will be examined for mendelian violations: this argument is a string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
    protected String FAMILY_STRUCTURE;

    @Argument(shortName="mvq", fullName="mendelianViolationQualThreshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 50;

    @Argument(fullName="tranchesFile", shortName="tf", doc="The input tranches file describing where to cut the data", required=false)
    private String TRANCHE_FILENAME = null;

    @Argument(fullName="ancestralAlignments", shortName="aa", doc="Fasta file with ancestral alleles", required=false)
    private File ancestralAlignmentsFile = null;

    // Variables
    private Set<SortableJexlVCMatchExp> jexlExpressions = new TreeSet<SortableJexlVCMatchExp>();
    private Set<String> compNames = new TreeSet<String>();
    private Set<String> knownNames = new TreeSet<String>();
    private Set<String> evalNames = new TreeSet<String>();

    private Set<String> sampleNamesForEvaluation = new TreeSet<String>();
    private Set<String> sampleNamesForStratification = new TreeSet<String>();
    private int numSamples = 0;

    // The list of stratifiers and evaluators to use
    private TreeSet<VariantStratifier> stratificationObjects = null;

    // The set of all possible evaluation contexts
    private HashMap<StateKey, NewEvaluationContext> evaluationContexts = null;

    // Output report
    private GATKReport report = null;

    // Public constants
    private static String ALL_SAMPLE_NAME = "all";

    // Utility class
    private final VariantEvalUtils variantEvalUtils = new VariantEvalUtils(this);

    // Ancestral alignments
    private IndexedFastaSequenceFile ancestralAlignments = null;

    /**
     * Initialize the stratifications, evaluations, evaluation contexts, and reporting object
     */
    public void initialize() {
        // Just list the modules, and exit quickly.
        if (LIST) { variantEvalUtils.listModulesAndExit(); }

        // Categorize each rod as an eval or a comp rod.
        for ( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if ( d.getName().startsWith("eval") ) {
                evalNames.add(d.getName());
            } else if ( d.getName().startsWith("comp") || d.getName().startsWith(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                compNames.add(d.getName());
            } else {
                logger.info(String.format("Not evaluating ROD binding '%s' because the name did not start with %s, comp, or eval", d.getName(), Utils.join(", ", KNOWN_NAMES)));
            }
        }

        // Barf if we don't have any eval tracks.
        if (evalNames.size() == 0) {
            throw new UserException("No evaluation tracks were specified.  Please bind one or more callsets to evaluate using the -B argument with a trackname that starts with the word 'eval'.");
        }

        // Add a dummy comp track if none exists
        if (compNames.size() == 0) {
            compNames.add("none");
        }

        // Set up set of known names
        knownNames.addAll(Arrays.asList(KNOWN_NAMES));

        // Now that we have all the rods categorized, determine the sample list from the eval rods.
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), evalNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        // Load the sample list
        sampleNamesForEvaluation.addAll(SampleUtils.getSamplesFromCommandLineInput(vcfSamples, SAMPLE_EXPRESSIONS));
        numSamples = NUM_SAMPLES > 0 ? NUM_SAMPLES : sampleNamesForEvaluation.size();

        if (Arrays.asList(STRATIFICATIONS_TO_USE).contains("Sample")) {
            sampleNamesForStratification.addAll(sampleNamesForEvaluation);
        }
        sampleNamesForStratification.add(ALL_SAMPLE_NAME);

        // Initialize select expressions
        for (VariantContextUtils.JexlVCMatchExp jexl : VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS)) {
            SortableJexlVCMatchExp sjexl = new SortableJexlVCMatchExp(jexl.name, jexl.exp);
            jexlExpressions.add(sjexl);
        }

        // Add select expressions for anything in the tranches file
        if ( TRANCHE_FILENAME != null ) {
            // we are going to build a few select names automatically from the tranches file
            for ( Tranche t : Tranche.readTranches(new File(TRANCHE_FILENAME)) ) {
                logger.info("Adding select for all variant above the pCut of : " + t);
                SELECT_EXPS.add(String.format(VariantRecalibrator.VQS_LOD_KEY + " >= %.2f", t.minVQSLod));
                SELECT_NAMES.add(String.format("TS-%.2f", t.ts));
            }
        }

        // Initialize the set of stratifications and evaluations to use
        stratificationObjects = variantEvalUtils.initializeStratificationObjects(this, NO_STANDARD_STRATIFICATIONS, STRATIFICATIONS_TO_USE);
        Set<Class<? extends VariantEvaluator>> evaluationObjects = variantEvalUtils.initializeEvaluationObjects(NO_STANDARD_MODULES, MODULES_TO_USE);

        // Initialize the evaluation contexts
        evaluationContexts = variantEvalUtils.initializeEvaluationContexts(stratificationObjects, evaluationObjects, null, null);

        // Initialize report table
        report = variantEvalUtils.initializeGATKReport(stratificationObjects, evaluationObjects);

        // Load ancestral alignments
        if (ancestralAlignmentsFile != null) {
            try {
                ancestralAlignments = new IndexedFastaSequenceFile(ancestralAlignmentsFile);
            } catch (FileNotFoundException e) {
                throw new ReviewedStingException(String.format("The ancestral alignments file, '%s', could not be found", ancestralAlignmentsFile.getAbsolutePath()));
            }
        }
    }

    /**
     * Collect relevant information from each variant in the supplied VCFs
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for ( NewEvaluationContext nec : evaluationContexts.values() ) {
            synchronized (nec) {
                nec.update0(tracker, ref, context);
            }
        }

        if (tracker != null) {
            String aastr = (ancestralAlignments == null) ? null : new String(ancestralAlignments.getSubsequenceAt(ref.getLocus().getContig(), ref.getLocus().getStart(), ref.getLocus().getStop()).getBases());

            //      track           sample  vc
            HashMap<String, HashMap<String, VariantContext>> vcs = variantEvalUtils.getVariantContexts(tracker, ref, compNames, evalNames, typesToUse != null);

            for ( String compName : compNames ) {
                VariantContext comp = vcs.containsKey(compName) && vcs.get(compName) != null && vcs.get(compName).containsKey(ALL_SAMPLE_NAME) ? vcs.get(compName).get(ALL_SAMPLE_NAME) : null;

                for ( String evalName : evalNames ) {
                    for ( String sampleName : sampleNamesForStratification ) {
                        VariantContext eval = vcs.containsKey(evalName) && vcs.get(evalName) != null ? vcs.get(evalName).get(sampleName) : null;

                        if ( typesToUse != null ) {
                            if ( eval != null && ! typesToUse.contains(eval.getType()) ) eval = null;
                            if ( comp != null && ! typesToUse.contains(comp.getType()) ) comp = null;
//                            if ( eval != null ) logger.info("Keeping " + eval);
                        }

                        if (eval != null && aastr != null) {
                            HashMap<String, Object> newAts = new HashMap<String, Object>(eval.getAttributes());
                            newAts.put("ANCESTRALALLELE", aastr);

                            eval = VariantContext.modifyAttributes(eval, newAts);
                        }

                        HashMap<VariantStratifier, ArrayList<String>> stateMap = new HashMap<VariantStratifier, ArrayList<String>>();
                        for ( VariantStratifier vs : stratificationObjects ) {
                            ArrayList<String> states = vs.getRelevantStates(ref, tracker, comp, compName, eval, evalName, sampleName);
                            stateMap.put(vs, states);
                        }

                        ArrayList<StateKey> stateKeys = new ArrayList<StateKey>();
                        variantEvalUtils.initializeStateKeys(stateMap, null, null, stateKeys);

                        HashSet<StateKey> stateKeysHash = new HashSet<StateKey>(stateKeys);

                        for ( StateKey stateKey : stateKeysHash ) {
                            NewEvaluationContext nec = evaluationContexts.get(stateKey);

                            synchronized (nec) {
                                nec.apply(tracker, ref, context, comp, eval);
                            }
                        }
                    }
                }
            }
        }

        return null;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) { return null; }

    @Override
    public Integer reduceInit() { return null; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return null; }

    /**
     * Output the finalized report
     *
     * @param result  an integer that doesn't get used for anything
     */
    public void onTraversalDone(Integer result) {
        logger.info("Finalizing variant report");

        for ( StateKey stateKey : evaluationContexts.keySet() ) {
            NewEvaluationContext nec = evaluationContexts.get(stateKey);

            for ( VariantEvaluator ve : nec.getEvaluationClassList().values() ) {
                ve.finalizeEvaluation();

                AnalysisModuleScanner scanner = new AnalysisModuleScanner(ve);
                Map<Field, DataPoint> datamap = scanner.getData();

                for (Field field : datamap.keySet()) {
                    try {
                        field.setAccessible(true);

                        if (field.get(ve) instanceof TableType) {
                            TableType t = (TableType) field.get(ve);

                            String subTableName = ve.getClass().getSimpleName() + "." + field.getName();
                            String subTableDesc = datamap.get(field).description();

                            GATKReportTable table;
                            if (!report.hasTable(subTableName)) {
                                report.addTable(subTableName, subTableDesc);
                                table = report.getTable(subTableName);

                                table.addPrimaryKey("entry", false);
                                table.addColumn(subTableName, subTableName);

                                for ( VariantStratifier vs : stratificationObjects ) {
                                    String columnName = vs.getClass().getSimpleName();

                                    table.addColumn(columnName, "unknown");
                                }

                                table.addColumn("row", "unknown");

                                for ( Object o : t.getColumnKeys() ) {
                                    String c;

                                    if (o instanceof String) {
                                        c = (String) o;
                                    } else {
                                        c = o.toString();
                                    }

                                    table.addColumn(c, 0.0);
                                }
                            } else {
                                table = report.getTable(subTableName);
                            }

                            for (int row = 0; row < t.getRowKeys().length; row++) {
                                String r = (String) t.getRowKeys()[row];

                                for ( VariantStratifier vs : stratificationObjects ) {
                                    String columnName = vs.getClass().getSimpleName();

                                    table.set(stateKey.toString() + r, columnName, stateKey.get(vs.getClass().getSimpleName()));
                                }

                                for (int col = 0; col < t.getColumnKeys().length; col++) {
                                    String c;
                                    if (t.getColumnKeys()[col] instanceof String) {
                                        c = (String) t.getColumnKeys()[col];
                                    } else {
                                        c = t.getColumnKeys()[col].toString();
                                    }

                                    String newStateKey = stateKey.toString() + r;
                                    table.set(newStateKey, c, t.getCell(row, col));

                                    table.set(newStateKey, "row", r);
                                }
                            }
                        } else {
                            GATKReportTable table = report.getTable(ve.getClass().getSimpleName());

                            for ( VariantStratifier vs : stratificationObjects ) {
                                String columnName = vs.getClass().getSimpleName();

                                table.set(stateKey.toString(), columnName, stateKey.get(vs.getClass().getSimpleName()));
                            }

                            table.set(stateKey.toString(), field.getName(), field.get(ve));
                        }
                    } catch (IllegalAccessException e) {
                        throw new StingException("IllegalAccessException: " + e);
                    }
                }
            }
        }

        report.print(out);
    }

    // Accessors
    public Logger getLogger() { return logger; }

    public int getNumSamples() { return numSamples; }

    public double getMinPhaseQuality() { return MIN_PHASE_QUALITY; }

    public String getFamilyStructure() { return FAMILY_STRUCTURE; }

    public double getMendelianViolationQualThreshold() { return MENDELIAN_VIOLATION_QUAL_THRESHOLD; }

    public TreeSet<VariantStratifier> getStratificationObjects() { return stratificationObjects; }

    public static String getAllSampleName() { return ALL_SAMPLE_NAME; }

    public Set<String> getKnownNames() { return knownNames; }

    public Set<String> getEvalNames() { return evalNames; }

    public Set<String> getSampleNamesForEvaluation() { return sampleNamesForEvaluation; }

    public Set<String> getSampleNamesForStratification() { return sampleNamesForStratification; }

    public Set<String> getCompNames() { return compNames; }

    public Set<SortableJexlVCMatchExp> getJexlExpressions() { return jexlExpressions; }

    public Set<String> getContigNames() {
        final TreeSet<String> contigs = new TreeSet<String>();
        for( final SAMSequenceRecord r :  getToolkit().getReferenceDataSource().getReference().getSequenceDictionary().getSequences()) {
            contigs.add(r.getSequenceName());
        }
        return contigs;
    }

}