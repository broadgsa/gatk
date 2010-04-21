/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableVariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.playground.utils.report.ReportMarshaller;
import org.broadinstitute.sting.playground.utils.report.VE2ReportFactory;
import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.classloader.PackageUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;

// todo -- evalations should support comment lines
// todo -- add Mendelian variable explanations (nDeNovo and nMissingTransmissions)

// todo -- site frequency spectrum eval (freq. of variants in eval as a function of their AC and AN numbers)
// todo -- clustered SNP counter
// todo -- HWEs
// todo -- indel metrics [count of sizes in/del should be in CountVariants]
// todo -- synonymous / non-synonmous ratio, or really just comparison of observed vs. expected biological annotation values

// todo -- Performance:
// todo -- deal with performance issues with variant contexts

// todo -- port over SNP density walker:
// todo -- see walker for WG calc but will need to make it work with intervals correctly

// todo -- counts of snps per target [target name, gene, etc]

// todo -- add subgroup of known variants as to those at hapmap sites [it's in the dbSNP record]

// Todo -- should really include argument parsing @annotations from subclass in this walker.  Very
// todo -- useful general capability.  Right now you need to add arguments to VariantEval2 to handle new
// todo -- evaluation arguments (which is better than passing a string!)

// todo -- these really should be implemented as default select expression
// todo Extend VariantEval, our general-purpose tool for SNP evaluation, to differentiate Ti/Tv at CpG islands and also
// todo classify (and count) variants into coding, non-coding, synonomous/non-symonomous, 2/4 fold degenerate sites, etc.
// todo Assume that the incoming VCF has the annotations (you don't need to do this) but VE2 should split up results by
// todo these catogies automatically (using the default selects)

// todo -- this is really more a documentation issue.  Really would be nice to have a pre-defined argument packet that
// todo -- can be provided to the system
// todo -- We agreed to report two standard values for variant evaluation from here out. One, we will continue to report
// todo -- the dbSNP 129 rate. Additionally, we will start to report the % of variants found that have already been seen in
// todo -- 1000 Genomes. This should be implemented as another standard comp_1kg binding, pointing to only variants
// todo -- discovered and released by 1KG.  Might need to make this data set ourselves and keep it in GATK/data like
// todo -- dbsnp rod
//
// todo -- aux. plotting routines for VE2
//
// todo -- implement as select statment, but it's hard for multi-sample calls.
// todo -- Provide separate dbsnp rates for het only calls and any call where there is at least one hom-var genotype,
// todo -- since hets are much more likely to be errors
//
// todo -- Add Heng's hom run metrics -- single sample haplotype block lengths


/**
 * Test routine for new VariantContext object
 */
public class VariantEvalWalker extends RodWalker<Integer, Integer> {
    // --------------------------------------------------------------------------------------------------------------
    //
    // walker arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", required=false)
    protected String[] SELECT_EXPS = {};

    @Argument(shortName="selectName", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", required=false)
    protected String[] SELECT_NAMES = {};

    @Argument(shortName="known", doc="Name of ROD bindings containing variant sites that should be treated as known when splitting eval rods into known and novel subsets", required=false)
    protected String[] KNOWN_NAMES = {rodDbSNP.STANDARD_DBSNP_TRACK_NAME};

    @Argument(shortName="sample", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", required=false)
    protected String[] SAMPLES = {};
    private List<String> SAMPLES_LIST = null;

    //
    // Arguments for choosing which modules to run
    //
    @Argument(fullName="evalModule", shortName="E", doc="One or more specific eval modules to apply to the eval track(s)", required=false)
    protected String[] modulesToUse = {};

    @Argument(fullName="useNoModules", shortName="none", doc="Use no eval modules", required=false)
    protected Boolean USE_NO_MODULES = false;

    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit")
    protected Boolean LIST = false;

    //
    // Arguments for Mendelian Violation calculations
    //
    @Argument(shortName="family", doc="If provided, genotypes in will be examined for mendelian violations: this argument is a string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
    protected String FAMILY_STRUCTURE;

    @Argument(shortName="MVQ", fullName="MendelianViolationQualThreshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 50;

    @Argument(shortName="outputVCF", fullName="InterestingSitesVCF", doc="If provided, interesting sites emitted to this vcf and the INFO field annotated as to why they are interesting", required=false)
    protected String outputVCF = null;

    private static double NO_MIN_QUAL_SCORE = -1.0;
    @Argument(shortName = "Q", fullName="minPhredConfidenceScore", doc="Minimum confidence score to consider an evaluation SNP a variant", required=false)
    public double minQualScore = NO_MIN_QUAL_SCORE;
    @Argument(shortName = "Qcomp", fullName="minPhredConfidenceScoreForComp", doc="Minimum confidence score to consider a comp SNP a variant", required=false)
    public double minCompQualScore = NO_MIN_QUAL_SCORE;

    // Right now we will only be looking at SNPS
    EnumSet<VariantContext.Type> ALLOW_VARIANT_CONTEXT_TYPES = EnumSet.of(VariantContext.Type.SNP, VariantContext.Type.NO_VARIATION);

    @Argument(shortName="rsID", fullName="rsID", doc="If provided, list of rsID and build number for capping known snps by their build date", required=false)
    protected String rsIDFile = null;

    @Argument(shortName="maxRsIDBuild", fullName="maxRsIDBuild", doc="If provided, only variants with rsIDs <= maxRsIDBuild will be included in the set of known snps", required=false)
    protected int maxRsIDBuild = Integer.MAX_VALUE;

    @Argument(shortName="reportType", fullName="reportType", doc="If provided, set the template type", required=false)
    protected VE2ReportFactory.VE2TemplateType reportType = VE2ReportFactory.defaultReportFormat;


    Set<String> rsIDsToExclude = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // private walker data
    //
    // --------------------------------------------------------------------------------------------------------------

    /** private class holding all of the information about a single evaluation group (e.g., for eval ROD) */
    private class EvaluationContext implements Comparable<EvaluationContext> {
        // useful for typing
        public String evalTrackName, compTrackName, novelty, filtered;
        public boolean enableInterestingSiteCaptures = false;
        VariantContextUtils.JexlVCMatchExp selectExp;
        Set<VariantEvaluator> evaluations;

        public boolean isIgnoringFilters()      { return filtered.equals(RAW_SET_NAME); }
        public boolean requiresFiltered()       { return filtered.equals(FILTERED_SET_NAME); }
        public boolean requiresNotFiltered()    { return filtered.equals(RETAINED_SET_NAME); }
        public boolean isIgnoringNovelty()      { return novelty.equals(ALL_SET_NAME); }
        public boolean requiresNovel()          { return novelty.equals(NOVEL_SET_NAME); }
        public boolean requiresKnown()          { return novelty.equals(KNOWN_SET_NAME); }

        public boolean isSelected() { return selectExp == null; }

        public String getDisplayName() {
            return Utils.join(".", Arrays.asList(evalTrackName, compTrackName, selectExp == null ? "all" : selectExp.name, filtered, novelty));
        }

        public int compareTo(EvaluationContext other) {
            return this.getDisplayName().compareTo(other.getDisplayName());
        }

        public EvaluationContext( String evalName, String compName, String novelty, String filtered, VariantContextUtils.JexlVCMatchExp selectExp ) {
            this.evalTrackName = evalName;
            this.compTrackName = compName;
            this.novelty = novelty;
            this.filtered = filtered;
            this.selectExp = selectExp;
            this.enableInterestingSiteCaptures = selectExp == null;
            this.evaluations = instantiateEvalationsSet();
        }
    }

    private List<EvaluationContext> contexts = null;

    // lists of all comp and eval ROD track names
    private Set<String> compNames = new HashSet<String>();
    private Set<String> evalNames = new HashSet<String>();

    private List<String> variantEvaluationNames = new ArrayList<String>();

    private static String RAW_SET_NAME      = "raw";
    private static String RETAINED_SET_NAME = "called";
    private static String FILTERED_SET_NAME = "filtered";
    private static String ALL_SET_NAME      = "all";
    private static String KNOWN_SET_NAME    = "known";
    private static String NOVEL_SET_NAME    = "novel";

    private static String NO_COMP_NAME = "N/A";

    private final static String CONTEXT_HEADER = "eval.comp.select.filter.novelty";
    private final static int N_CONTEXT_NAME_PARTS = CONTEXT_HEADER.split("\\.").length;
    private static int[] nameSizes = new int[N_CONTEXT_NAME_PARTS];
    static {
        int i = 0;
        for ( String elt : CONTEXT_HEADER.split("\\.") )
            nameSizes[i++] = elt.length();
    }

    // Dynamically determined variantEvaluation classes
    private List<Class<? extends VariantEvaluator>> evaluationClasses = null;

    /** output writer for interesting sites */
    private VCFWriter writer = null;
    private boolean wroteHeader = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    // --------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if ( LIST )
            listModulesAndExit();

        SAMPLES_LIST = Arrays.asList(SAMPLES);

        determineEvalations();

        List<VariantContextUtils.JexlVCMatchExp> selectExps = VariantContextUtils.initializeMatchExps(SELECT_NAMES, SELECT_EXPS);

        for ( ReferenceOrderedDataSource d : this.getToolkit().getRodDataSources() ) {
            if ( d.getName().startsWith("eval") ) {
                evalNames.add(d.getName());
            } else if ( d.getName().startsWith(rodDbSNP.STANDARD_DBSNP_TRACK_NAME) ) {
                // it feels like overkill (i.e. too much output) to include the dbsnp track as
                //  truth given that it's also used in the known/novel stratification.  If this
                //  becomes useful for some reason, uncomment this line...
                //compNames.add(d.getName());
            } else if ( d.getName().startsWith("hapmap") || d.getName().startsWith("comp") ) {
                compNames.add(d.getName());
            } else {
                logger.info("Not evaluating ROD binding " + d.getName());
            }
        }

        // if no comp rod was provided, we still want to be able to do evaluations, so use a default comp name
        if ( compNames.size() == 0 )
            compNames.add(NO_COMP_NAME);

        contexts = initializeEvaluationContexts(evalNames, compNames, selectExps);
        determineContextNamePartSizes();

        if ( outputVCF != null )
            writer = new VCFWriter(new File(outputVCF));

        if ( rsIDFile != null ) {
            if ( maxRsIDBuild == Integer.MAX_VALUE )
                throw new IllegalArgumentException("rsIDFile " + rsIDFile + " was given but associated max RSID build parameter wasn't available");
            rsIDsToExclude = getrsIDsToExclude(new File(rsIDFile), maxRsIDBuild);
        }
    }

    private void listModulesAndExit() {
        List<Class<? extends VariantEvaluator>> veClasses = PackageUtils.getClassesImplementingInterface(VariantEvaluator.class);
        out.println("\nAvailable eval modules:");
        for (int i = 0; i < veClasses.size(); i++)
            out.println("\t" + veClasses.get(i).getSimpleName());
        out.println();
        System.exit(0);
    }

    private static Set<String> getrsIDsToExclude(File rsIDFile, int maxRsIDBuild) {
        List<String> toExclude = new LinkedList<String>();

        int n = 1;
        try {
            for ( String line : new XReadLines(rsIDFile) ) {
                String parts[] = line.split(" ");
                if ( parts.length != 2 )
                    throw new StingException("Invalid rsID / build pair at " + n + " line = " + line );
                //System.out.printf("line %s %s %s%n", line, parts[0], parts[1]);
                if ( Integer.valueOf(parts[1]) > maxRsIDBuild ) {
                    //System.out.printf("Excluding %s%n", line);
                    toExclude.add("rs"+parts[0]);
                }
                n++;

                if ( n % 1000000 == 0 )
                    logger.info(String.format("Read %d rsIDs from rsID -> build file", n));
            }
        } catch (FileNotFoundException e) {
            throw new StingException(e.getMessage());
        }

        logger.info(String.format("Excluding %d of %d (%.2f%%) rsIDs found from builds > %d",
                toExclude.size(), n, ((100.0 * toExclude.size())/n), maxRsIDBuild));

        return new HashSet<String>(toExclude);
    }

    private final static String ID = "ID";
    private boolean excludeComp(VariantContext vc) {
        String id = vc != null && vc.hasAttribute(ID) ? vc.getAttributeAsString(ID) : null;
        boolean ex = rsIDsToExclude != null && id != null && rsIDsToExclude.contains(id);
        //System.out.printf("Testing id %s ex=%b against %s%n", id, ex, vc);
        return ex;
    }

    private void determineEvalations() {
        // create a map for all eval modules for easy lookup
        HashMap<String, Class<? extends VariantEvaluator>> classMap = new HashMap<String, Class<? extends VariantEvaluator>>();
        for ( Class<? extends VariantEvaluator> c : PackageUtils.getClassesImplementingInterface(VariantEvaluator.class) )
            classMap.put(c.getSimpleName(), c);

        if ( USE_NO_MODULES ) {
            evaluationClasses = new ArrayList<Class<? extends VariantEvaluator>>(0);            
        } else if ( modulesToUse.length == 0 ) {
            evaluationClasses = new ArrayList<Class<? extends VariantEvaluator>>(classMap.values());
        } else {
            // get the specific classes provided
            evaluationClasses = new ArrayList<Class<? extends VariantEvaluator>>(modulesToUse.length);
            for ( String module : modulesToUse ) {
                Class<? extends VariantEvaluator> moduleClass = classMap.get(module);
                if ( moduleClass == null )
                    throw new StingException("Class " + module + " is not found; please check that you have specified the class name correctly");
                evaluationClasses.add(moduleClass);
            }
        }

        for ( VariantEvaluator e : instantiateEvalationsSet() ) {
            // for collecting purposes
            variantEvaluationNames.add(e.getName());
            logger.debug("Including VariantEvaluator " + e.getName() + " of class " + e.getClass());
        }

        Collections.sort(variantEvaluationNames);
    }

    private <T> List<T> append(List<T> selectExps, T elt) {
        List<T> l = new ArrayList<T>(selectExps);
        l.add(elt);
        return l;
    }

    private List<EvaluationContext> initializeEvaluationContexts(Set<String> evalNames, Set<String> compNames, List<VariantContextUtils.JexlVCMatchExp> selectExps) {
        List<EvaluationContext> contexts = new ArrayList<EvaluationContext>();

        selectExps = append(selectExps, null);
        for ( String evalName : evalNames ) {
            for ( String compName : compNames ) {
                for ( VariantContextUtils.JexlVCMatchExp e : selectExps ) {
                    for ( String filteredName : Arrays.asList(RAW_SET_NAME, RETAINED_SET_NAME, FILTERED_SET_NAME) ) {
                        for ( String novelty : Arrays.asList(ALL_SET_NAME, KNOWN_SET_NAME, NOVEL_SET_NAME) ) {
                            EvaluationContext context = new EvaluationContext(evalName, compName, novelty, filteredName, e);
                            contexts.add(context);
                        }
                    }
                }
            }
        }

        Collections.sort(contexts);
        return contexts;
    }

    private Set<VariantEvaluator> instantiateEvalationsSet() {
        Set<VariantEvaluator> evals = new HashSet<VariantEvaluator>();
        Object[] args = new Object[]{this};
        Class[] argTypes = new Class[]{this.getClass()};

        for ( Class c : evaluationClasses ) {
            try {
                Constructor constructor = c.getConstructor(argTypes);
                VariantEvaluator eval = (VariantEvaluator)constructor.newInstance(args);
                evals.add(eval);
            } catch (InstantiationException e) {
                throw new StingException(String.format("Cannot instantiate class '%s': must be concrete class", c.getSimpleName()));
            } catch (NoSuchMethodException e) {
                throw new StingException(String.format("Cannot find expected constructor for class '%s': must have constructor accepting a single VariantEval2Walker object", c.getSimpleName()));
            } catch (IllegalAccessException e) {
                throw new StingException(String.format("Cannot instantiate class '%s':", c.getSimpleName()));
            } catch (InvocationTargetException e) {
                throw new StingException(String.format("Cannot instantiate class '%s':", c.getSimpleName()));
            }
        }

        return evals;
    }



    private boolean captureInterestingSitesOfEvalSet(EvaluationContext group) {
        //System.out.printf("checking %s%n", name);
        return group.requiresNotFiltered() && group.isIgnoringNovelty();
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    // --------------------------------------------------------------------------------------------------------------

    // todo -- call a single function to build a map from track name -> variant context / null for all
    //      -- eval + comp names.  Use this data structure to get data throughout rest of the loops here
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //System.out.printf("map at %s with %d skipped%n", context.getLocation(), context.getSkippedBases());

        Map<String, VariantContext> vcs = getVariantContexts(tracker, context);
        //Collection<VariantContext> comps = getCompVariantContexts(tracker, context);

        // to enable walking over pairs where eval or comps have no elements
        for ( EvaluationContext group : contexts ) {
            VariantContext vc = vcs.get(group.evalTrackName);

            //logger.debug(String.format("Updating %s with variant", vc));
            Set<VariantEvaluator> evaluations = group.evaluations;
            boolean evalWantsVC = applyVCtoEvaluation(vc, vcs, group);
            List<String> interestingReasons = new ArrayList<String>();

            for ( VariantEvaluator evaluation : evaluations ) {
                if ( evaluation.enabled() ) {
                    // we always call update0 in case the evaluation tracks things like number of bases covered
                    evaluation.update0(tracker, ref, context);

                    // the other updateN methods don't see a null context
                    if ( tracker == null )
                        continue;

                    // now call the single or paired update function
                    switch ( evaluation.getComparisonOrder() ) {
                        case 1:
                            if ( evalWantsVC && vc != null ) {
                                String interesting = evaluation.update1(vc, tracker, ref, context);
                                if ( interesting != null ) interestingReasons.add(interesting);
                            }
                            break;
                        case 2:
                            VariantContext comp = vcs.get(group.compTrackName);
                            if ( comp != null && minCompQualScore != NO_MIN_QUAL_SCORE && comp.getNegLog10PError() < (minCompQualScore / 10.0))
                                comp = null;                            
                            String interesting = evaluation.update2( evalWantsVC ? vc : null, comp, tracker, ref, context );
                            if ( interesting != null ) interestingReasons.add(interesting);
                            break;
                        default:
                            throw new StingException("BUG: Unexpected evaluation order " + evaluation);
                    }
                }
            }

            if ( tracker != null && group.enableInterestingSiteCaptures && captureInterestingSitesOfEvalSet(group) )
                writeInterestingSite(interestingReasons, vc, ref.getBase());
        }

        return 0;
    }

    private void writeInterestingSite(List<String> interestingReasons, VariantContext vc, char ref) {
        if ( vc != null && writer != null && interestingReasons.size() > 0 ) {
            // todo -- the vc == null check is because you can be interesting because you are a FN, and so VC == null
            MutableVariantContext mvc = new MutableVariantContext(vc);

            for ( String why : interestingReasons ) {
                String key, value;
                String[] parts = why.split("=");

                switch ( parts.length ) {
                    case 1:
                        key = parts[0];
                        value = "1";
                        break;
                    case 2:
                        key = parts[0];
                        value = parts[1];
                        break;
                    default:
                        throw new IllegalStateException("BUG: saw a interesting site reason sting with multiple = signs " + why);
                }

                mvc.putAttribute(key, value);
            }

            if ( ! wroteHeader ) {
                writer.writeHeader(VariantContextAdaptors.createVCFHeader(null, vc));
                wroteHeader = true;
            }

            writer.addRecord(VariantContextAdaptors.toVCF(mvc, ref));
            //interestingReasons.clear();
        }
    }

    private boolean applyVCtoEvaluation(VariantContext vc, Map<String, VariantContext> vcs, EvaluationContext group) {
        if ( vc == null )
            return true;

        if ( minQualScore != NO_MIN_QUAL_SCORE && vc.getNegLog10PError() < (minQualScore / 10.0)) {
            //System.out.printf("exclude %s%n", vc);
            return false;
        }

        if ( group.requiresFiltered() && vc.isNotFiltered() )
            return false;

        if ( group.requiresNotFiltered() && vc.isFiltered() )
            return false;

        boolean vcKnown = vcIsKnown(vc, vcs, KNOWN_NAMES);
        if ( group.requiresKnown() && ! vcKnown )
            return false;
        else if ( group.requiresNovel() && vcKnown )
            return false;

        if ( group.selectExp != null && ! VariantContextUtils.match(vc, group.selectExp) )
            return false;

        // nothing invalidated our membership in this set
        return true;
    }

    private boolean vcIsKnown(VariantContext vc, Map<String, VariantContext> vcs, String[] knownNames ) {
        for ( String knownName : knownNames ) {
            VariantContext known = vcs.get(knownName);
            if ( known != null && known.isNotFiltered() && known.getType() == vc.getType() ) {
                return true;
            }
        }

        return false;
    }

// can't handle this situation
// todo -- warning, this leads to some missing SNPs at complex loci, such as:
// todo -- 591     1       841619  841620  rs4970464       0       -       A       A       -/C/T   genomic mixed   unknown 0       0       near-gene-3     exact   1
// todo -- 591     1       841619  841620  rs62677860      0       +       A       A       C/T     genomic single  unknown 0       0       near-gene-3     exact   1
//
//logger.info(String.format("Ignore second+ events at locus %s in rod %s => rec is %s", context.getLocation(), rodList.getName(), rec));

    private Map<String, VariantContext> getVariantContexts(RefMetaDataTracker tracker, AlignmentContext context) {
        // todo -- we need to deal with dbSNP where there can be multiple records at the same start site.  A potential solution is to
        // todo -- allow the variant evaluation to specify the type of variants it wants to see and only take the first such record at a site
        Map<String, VariantContext> bindings = new HashMap<String, VariantContext>();
        if ( tracker != null ) {
            bindVariantContexts(bindings, evalNames, tracker, context, false);
            bindVariantContexts(bindings, compNames, tracker, context, true);
        }
        return bindings;
    }

    private void bindVariantContexts(Map<String, VariantContext> map, Collection<String> names,
                                     RefMetaDataTracker tracker, AlignmentContext context, boolean allowExcludes ) {
        for ( String name : names ) {
            Collection<VariantContext> contexts = tracker.getVariantContexts(name, ALLOW_VARIANT_CONTEXT_TYPES, context.getLocation(), true, true);
            if ( contexts.size() > 1 )
                throw new StingException("Found multiple variant contexts at " + context.getLocation());

            VariantContext vc = contexts.size() == 1 ? contexts.iterator().next() : null;

            if ( vc != null && vc.hasGenotypes(SAMPLES_LIST) && SAMPLES_LIST.size() > 0 ) {
                //if ( ! name.equals("eval") ) logger.info(String.format("subsetting VC %s", vc));
                vc = vc.subContextFromGenotypes(vc.getGenotypes(SAMPLES_LIST).values());
                //if ( ! name.equals("eval") ) logger.info(String.format("  => VC %s", vc));
            }

            map.put(name, allowExcludes && excludeComp(vc) ? null : vc);
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    // --------------------------------------------------------------------------------------------------------------
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer point, Integer sum) {
        return point + sum;
    }

    public VariantEvaluator getEvalByName(String name, Set<VariantEvaluator> s) {
        for ( VariantEvaluator e : s )
            if ( e.getName().equals(name) )
                return e;
        return null;
    }

    private void determineContextNamePartSizes() {
        for ( EvaluationContext group : contexts ) {
            String[] parts = group.getDisplayName().split("\\.");
            if ( parts.length != N_CONTEXT_NAME_PARTS ) {
                throw new StingException("Unexpected number of eval name parts " + group.getDisplayName() + " length = " + parts.length + ", expected " + N_CONTEXT_NAME_PARTS);
            } else {
                for ( int i = 0; i < parts.length; i++ )
                    nameSizes[i] = Math.max(nameSizes[i], parts[i].length());
            }
        }
    }

    private String formatKeyword(String keyWord) {
        //System.out.printf("keyword %s%n", keyWord);

        StringBuilder s = new StringBuilder();
        int i = 0;
        for ( String part : keyWord.split("\\.") ) {
            //System.out.printf("part %s %d%n", part, nameSizes[i]);
            s.append(String.format("%" + nameSizes[i] + "s ", part));
            i++;
        }

        return s.toString();
    }

    public void onTraversalDone(Integer result) {
        ReportMarshaller marshaller = VE2ReportFactory.createMarhsaller(out,reportType,createExtraOutputTags());
        for ( String evalName : variantEvaluationNames ) {
            for ( EvaluationContext group : contexts ) {
                VariantEvaluator eval = getEvalByName(evalName, group.evaluations);
                // finalize the evaluation
                eval.finalizeEvaluation();

                if ( eval.enabled() )
                    marshaller.write(createPrependNodeList(group),eval);
            }
        }
        marshaller.close();
    }

    /**
     * create some additional output lines about the analysis
     * @return a list of nodes to attach to the report as tags
     */
    private List<Node> createExtraOutputTags() {
        List<Node> list = new ArrayList();
        list.add(new Node("reference file",getToolkit().getArguments().referenceFile.getName(),"The reference sequence file"));
        for (String binding : getToolkit().getArguments().RODBindings)
            list.add(new Node("ROD binding",binding,"The reference sequence file"));
        return list;
    }


    /**
     * given the evaluation name, and the context, create the list of pre-pended nodes for the output system.
     * Currently it expects the the following list: jexl_expression, evaluation_name, comparison_name, filter_name,
     * novelty_name
     * @param group the evaluation context
     * @return a list of Nodes to prepend the analysis module output with
     */
    private List<Node> createPrependNodeList(EvaluationContext group) {
        // add the branching nodes: jexl expression, comparison track, eval track etc
        Node jexlNode = new Node("jexl_expression",(group.selectExp != null) ? group.selectExp.name : "none","The jexl filtering expression");
        Node compNode = new Node("comparison_name",group.compTrackName,"The comparison track name");
        Node evalNode = new Node("evaluation_name",group.evalTrackName,"The evaluation name");
        Node filterNode = new Node("filter_name",group.filtered,"The filter name");
        Node noveltyNode = new Node("novelty_name",group.novelty,"The novelty name");
        // the ordering is important below, this is the order the columns will appear in any output format
        return Arrays.asList(evalNode,compNode,jexlNode,filterNode,noveltyNode);
    }

    protected Logger getLogger() { return logger; }
}
