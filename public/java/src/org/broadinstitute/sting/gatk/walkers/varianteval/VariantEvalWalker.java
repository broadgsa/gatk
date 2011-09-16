package org.broadinstitute.sting.gatk.walkers.varianteval;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.JexlExpression;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.*;
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.Tranche;
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.*;

/**
 * General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)
 *
 * <p>
 * Given a variant callset, it is common to calculate various quality control metrics. These metrics include the number of
 * raw or filtered SNP counts; ratio of transition mutations to transversions; concordance of a particular sample's calls
 * to a genotyping chip; number of singletons per sample; etc. Furthermore, it is often useful to stratify these metrics
 * by various criteria like functional class (missense, nonsense, silent), whether the site is CpG site, the amino acid
 * degeneracy of the site, etc. VariantEval facilitates these calculations in two ways: by providing several built-in
 * evaluation and stratification modules, and by providing a framework that permits the easy development of new evaluation
 * and stratification modules.
 *
 * <h2>Input</h2>
 * <p>
 * One or more variant sets to evaluate plus any number of comparison sets.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Evaluation tables detailing the results of the eval modules which were applied.
 * For example:
 * <pre>
 * output.eval.gatkreport:
 * ##:GATKReport.v0.1 CountVariants : Counts different classes of variants in the sample
 * CountVariants  CompRod   CpG      EvalRod  JexlExpression  Novelty  nProcessedLoci  nCalledLoci  nRefLoci  nVariantLoci  variantRate ...
 * CountVariants  dbsnp     CpG      eval     none            all      65900028        135770       0         135770        0.00206024  ...
 * CountVariants  dbsnp     CpG      eval     none            known    65900028        47068        0         47068         0.00071423  ...
 * CountVariants  dbsnp     CpG      eval     none            novel    65900028        88702        0         88702         0.00134601  ...
 * CountVariants  dbsnp     all      eval     none            all      65900028        330818       0         330818        0.00502000  ...
 * CountVariants  dbsnp     all      eval     none            known    65900028        120685       0         120685        0.00183133  ...
 * CountVariants  dbsnp     all      eval     none            novel    65900028        210133       0         210133        0.00318866  ...
 * CountVariants  dbsnp     non_CpG  eval     none            all      65900028        195048       0         195048        0.00295976  ...
 * CountVariants  dbsnp     non_CpG  eval     none            known    65900028        73617        0         73617         0.00111710  ...
 * CountVariants  dbsnp     non_CpG  eval     none            novel    65900028        121431       0         121431        0.00184265  ...
 * ...
 * </pre>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantEval \
 *   -o output.eval.gatkreport \
 *   --eval:set1 set1.vcf \
 *   --eval:set2 set2.vcf \
 *   [--comp comp.vcf]
 * </pre>
 *
 */
@Reference(window=@Window(start=-50, stop=50))
public class VariantEvalWalker extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

    @Output
    protected PrintStream out;

    /**
     * The variant file(s) to evaluate.
     */
    @Input(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", required=true)
    public List<RodBinding<VariantContext>> evals;

    /**
     * The variant file(s) to compare against.
     */
    @Input(fullName="comp", shortName = "comp", doc="Input comparison file(s)", required=false)
    public List<RodBinding<VariantContext>> compsProvided = Collections.emptyList();
    private List<RodBinding<VariantContext>> comps = new ArrayList<RodBinding<VariantContext>>();

    /**
     * dbSNP comparison VCF.  By default, the dbSNP file is used to specify the set of "known" variants.
     * Other sets can be specified with the -knownName (--known_names) argument.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // Help arguments
    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit", required=false)
    protected Boolean LIST = false;

    // Partitioning the data arguments
    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", required=false)
    protected ArrayList<String> SELECT_EXPS = new ArrayList<String>();

    @Argument(shortName="selectName", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", required=false)
    protected ArrayList<String> SELECT_NAMES = new ArrayList<String>();

    @Argument(fullName="sample", shortName="sn", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", required=false)
    protected Set<String> SAMPLE_EXPRESSIONS;

    /**
     * List of rod tracks to be used for specifying "known" variants other than dbSNP.
     */
    @Argument(shortName="knownName", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected HashSet<String> KNOWN_NAMES = new HashSet<String>();
    List<RodBinding<VariantContext>> knowns = new ArrayList<RodBinding<VariantContext>>();

    // Stratification arguments
    @Argument(fullName="stratificationModule", shortName="ST", doc="One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless -noS is specified)", required=false)
    protected String[] STRATIFICATIONS_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardStratifications", shortName="noST", doc="Do not use the standard stratification modules by default (instead, only those that are specified with the -S option)", required=false)
    protected Boolean NO_STANDARD_STRATIFICATIONS = false;

    /**
     * See the -list argument to view available modules.
     */
    @Argument(fullName="evalModule", shortName="EV", doc="One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless -noE is specified)", required=false)
    protected String[] MODULES_TO_USE = {};

    @Argument(fullName="doNotUseAllStandardModules", shortName="noEV", doc="Do not use the standard modules by default (instead, only those that are specified with the -E option)", required=false)
    protected Boolean NO_STANDARD_MODULES = false;

    // Other arguments
    @Argument(fullName="numSamples", shortName="ns", doc="Number of samples (used if no samples are available in the VCF file", required=false)
    protected Integer NUM_SAMPLES = 0;

    @Argument(fullName="minPhaseQuality", shortName="mpq", doc="Minimum phasing quality", required=false)
    protected double MIN_PHASE_QUALITY = 10.0;

    /**
     * This argument is a string formatted as dad+mom=child where these parameters determine which sample names are examined.
     */
    @Argument(shortName="family", doc="If provided, genotypes in will be examined for mendelian violations", required=false)
    protected String FAMILY_STRUCTURE;

    @Argument(shortName="mvq", fullName="mendelianViolationQualThreshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 50;

    @Argument(fullName="ancestralAlignments", shortName="aa", doc="Fasta file with ancestral alleles", required=false)
    private File ancestralAlignmentsFile = null;

    @Argument(fullName="requireStrictAlleleMatch", shortName="strict", doc="If provided only comp and eval tracks with exactly matching reference and alternate alleles will be counted as overlapping", required=false)
    private boolean requireStrictAlleleMatch = false;

    // Variables
    private Set<SortableJexlVCMatchExp> jexlExpressions = new TreeSet<SortableJexlVCMatchExp>();

    private Set<String> sampleNamesForEvaluation = new TreeSet<String>();
    private Set<String> sampleNamesForStratification = new TreeSet<String>();
    private int numSamples = 0;

    // The list of stratifiers and evaluators to use
    private TreeSet<VariantStratifier> stratificationObjects = null;

    // The set of all possible evaluation contexts
    private HashMap<StateKey, NewEvaluationContext> evaluationContexts = null;

    // important stratifications
    private boolean byFilterIsEnabled = false;
    private boolean perSampleIsEnabled = false;

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

        // maintain the full list of comps
        comps.addAll(compsProvided);
        if ( dbsnp.dbsnp.isBound() ) {
            comps.add(dbsnp.dbsnp);
            knowns.add(dbsnp.dbsnp);
        }

        // Add a dummy comp track if none exists
        if ( comps.size() == 0 )
            comps.add(new RodBinding<VariantContext>(VariantContext.class, "none", "UNBOUND", "", new Tags()));

        // Set up set of additional knowns
        for ( RodBinding<VariantContext> compRod : comps ) {
            if ( KNOWN_NAMES.contains(compRod.getName()) )
                knowns.add(compRod);
        }

        // Collect the eval rod names
        Set<String> evalNames = new TreeSet<String>();
        for ( RodBinding<VariantContext> evalRod : evals )
            evalNames.add(evalRod.getName());

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

        // Initialize the set of stratifications and evaluations to use
        stratificationObjects = variantEvalUtils.initializeStratificationObjects(this, NO_STANDARD_STRATIFICATIONS, STRATIFICATIONS_TO_USE);
        Set<Class<? extends VariantEvaluator>> evaluationObjects = variantEvalUtils.initializeEvaluationObjects(NO_STANDARD_MODULES, MODULES_TO_USE);
        for ( VariantStratifier vs : getStratificationObjects() ) {
            if ( vs.getClass().getSimpleName().equals("Filter") )
                byFilterIsEnabled = true;
            else if ( vs.getClass().getSimpleName().equals("Sample") )
                perSampleIsEnabled = true;
        }

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

            //      --------- track ---------           sample  - VariantContexts -
            HashMap<RodBinding<VariantContext>, HashMap<String, Set<VariantContext>>> evalVCs = variantEvalUtils.bindVariantContexts(tracker, ref, evals, byFilterIsEnabled, true, perSampleIsEnabled);
            HashMap<RodBinding<VariantContext>, HashMap<String, Set<VariantContext>>> compVCs = variantEvalUtils.bindVariantContexts(tracker, ref, comps, byFilterIsEnabled, false, false);

            // for each eval track
            for ( final RodBinding<VariantContext> evalRod : evals ) {
                final HashMap<String, Set<VariantContext>> evalSet = evalVCs.containsKey(evalRod) ? evalVCs.get(evalRod) : new HashMap<String, Set<VariantContext>>(0);

                // for each sample stratifier
                for ( final String sampleName : sampleNamesForStratification ) {
                    Set<VariantContext> evalSetBySample = evalSet.get(sampleName);
                    if ( evalSetBySample == null ) {
                        evalSetBySample = new HashSet<VariantContext>(1);
                        evalSetBySample.add(null);
                    }

                    // for each eval in the track
                    for ( VariantContext eval : evalSetBySample ) {
                        // deal with ancestral alleles if requested
                        if ( eval != null && aastr != null ) {
                            HashMap<String, Object> newAts = new HashMap<String, Object>(eval.getAttributes());
                            newAts.put("ANCESTRALALLELE", aastr);
                            eval = VariantContext.modifyAttributes(eval, newAts);
                        }

                        // for each comp track
                        for ( final RodBinding<VariantContext> compRod : comps ) {
                            // no sample stratification for comps
                            final HashMap<String, Set<VariantContext>> compSetHash = compVCs.get(compRod);
                            final Set<VariantContext> compSet = (compSetHash == null || compSetHash.size() == 0) ? new HashSet<VariantContext>(0) : compVCs.get(compRod).values().iterator().next();

                            // find the comp
                            final VariantContext comp = findMatchingComp(eval, compSet);

                            HashMap<VariantStratifier, List<String>> stateMap = new HashMap<VariantStratifier, List<String>>();
                            for ( VariantStratifier vs : stratificationObjects ) {
                                List<String> states = vs.getRelevantStates(ref, tracker, comp, compRod.getName(), eval, evalRod.getName(), sampleName);
                                stateMap.put(vs, states);
                            }

                            ArrayList<StateKey> stateKeys = new ArrayList<StateKey>();
                            variantEvalUtils.initializeStateKeys(stateMap, null, null, stateKeys);

                            HashSet<StateKey> stateKeysHash = new HashSet<StateKey>(stateKeys);

                            for ( StateKey stateKey : stateKeysHash ) {
                                NewEvaluationContext nec = evaluationContexts.get(stateKey);

                                // eval against the comp
                                synchronized (nec) {
                                    nec.apply(tracker, ref, context, comp, eval);
                                }

                                // eval=null against all comps of different type
                                for ( VariantContext otherComp : compSet ) {
                                    if ( otherComp != comp ) {
                                        synchronized (nec) {
                                            nec.apply(tracker, ref, context, otherComp, null);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return null;
    }

    private VariantContext findMatchingComp(final VariantContext eval, final Set<VariantContext> comps) {
        // if no comps, return null
        if ( comps == null || comps.isEmpty() )
            return null;

        // if no eval, return any comp
        if ( eval == null )
            return comps.iterator().next();

        // find all of the matching comps
        List<VariantContext> matchingComps = new ArrayList<VariantContext>(comps.size());
        for ( VariantContext comp : comps ) {
            if ( comp.getType() == eval.getType() )
                matchingComps.add(comp);
        }

        // if no matching comp, return null
        if ( matchingComps.size() == 0 )
            return null;

        // find the comp which matches both the reference allele and alternate allele from eval
        Allele altEval = eval.getAlternateAlleles().size() == 0 ? null : eval.getAlternateAllele(0);
        for ( VariantContext comp : matchingComps ) {
            Allele altComp = comp.getAlternateAlleles().size() == 0 ? null : comp.getAlternateAllele(0);
            if ( (altEval == null && altComp == null) || (altEval != null && altEval.equals(altComp) && eval.getReference().equals(comp.getReference())) )
                return comp;
        }

        // if none match, just return the first one unless we require a strict match
        return (requireStrictAlleleMatch ? null : matchingComps.get(0));
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

    public List<RodBinding<VariantContext>> getKnowns() { return knowns; }

    public List<RodBinding<VariantContext>> getEvals() { return evals; }

    public Set<String> getSampleNamesForEvaluation() { return sampleNamesForEvaluation; }

    public Set<String> getSampleNamesForStratification() { return sampleNamesForStratification; }

    public List<RodBinding<VariantContext>> getComps() { return comps; }

    public Set<SortableJexlVCMatchExp> getJexlExpressions() { return jexlExpressions; }

    public Set<String> getContigNames() {
        final TreeSet<String> contigs = new TreeSet<String>();
        for( final SAMSequenceRecord r :  getToolkit().getReferenceDataSource().getReference().getSequenceDictionary().getSequences()) {
            contigs.add(r.getSequenceName());
        }
        return contigs;
    }

}