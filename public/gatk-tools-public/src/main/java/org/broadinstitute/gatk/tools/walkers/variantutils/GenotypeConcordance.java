/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeader;

import java.io.PrintStream;
import java.util.*;

/**
 * Genotype concordance between two callsets
 *
 * <p>
 *  This tool takes in two callsets (vcfs) and tabulates the number of sites which overlap and share alleles,
 *  and for each sample, the genotype-by-genotype counts (e.g. the number of sites at which a sample was
 *  called homozygous-reference in the EVAL callset, but homozygous-variant in the COMP callset). It outputs these
 *  counts as well as convenient proportions (such as the proportion of het calls in the EVAL which were called REF in
 *  the COMP) and metrics (such as NRD and NRS).
 *  </p>
 *
 *  <h3>Input</h3>
 *  <p>
 *  Genotype concordance requires two callsets (as it does a comparison): an EVAL and a COMP callset, specified via
 *  the -eval and -comp arguments. Typically, the EVAL callset is an experimental set you want to evaluate, while the
 *  COMP callset is a previously existing set used as a standard for comparison (taken to represent "truth").
 *  </p>
 *  <p>
 *  (Optional) Jexl expressions for genotype-level filtering of EVAL or COMP genotypes, specified via the -gfe and
 *  -cfe arguments, respectively.
 *  </p>
 *
 *  <h3>Output</h3>
 *  <p>
 *  Genotype Concordance writes a GATK report to the specified file (via -o), consisting of multiple tables of counts
 *  and proportions. These tables are constructed on a per-sample basis, and include counts of EVAL vs COMP genotype
 *  states.
 *  </p>
 *  <h4>Tables</h4>
 *  <p>
 *  Headers for the (non-moltenized -- see below) GenotypeConcordance counts and proportions tables give the genotype of
 *  the EVAL callset followed by the genotype of the COMP callset. For example the value corresponding to HOM_REF_HET
 *  reflects variants called HOM_REF in the EVAL callset and HET in the COMP callset. Variants for which the alternate
 *  alleles between the EVAL and COMP sample did not match are excluded from genotype comparisons and given in the
 *  "Mismatching_Alleles" field.
 *  </p>
 *  <p>
 *  It may be informative to reshape rows of the GenotypeConcordance counts and proportions tables into separate row-major tables
 *  where the columns indicate the COMP genotype and the rows indicate the EVAL genotype for easy comparison between the
 *  two callsets. This can be done with the gsa.reshape.concordance.table function in the gsalib R library.
 *  In Excel this can be accomplished using the OFFSET function.
 *  </p>
 *  <ul>
 *      <li><i>GenotypeConcordance_CompProportions</i>: Gives the proportion of variants in each category normalized to the total number of called genotypes in the COMP callset</li>
 *      <li><i>GenotypeConcordance_Counts</i>: Gives the counts for number of genotypes in each category</li>
 *      <li><i>GenotypeConcordance_EvalProportions</i>: Gives the proportion of genotypes in each category normalized to the total number of called genotypes in the EVAL callset</li>
 *      <li><i>GenotypeConcordance_Summary</i>: Summary statistics for the sum of all samples and each sample individually. See below for definitions.</li>
 *      <li><i>SiteConcordance_Summary</i>: Gives comparison counts of called genotypes and their alleles between the two callsets. See below for definitions.</li>
 *  </ul>
 *  </p>
 *
 *  <h4>Term and metrics definitions</h4>
 * <p>
 * <ul>
 *          <li><i>GenotypeConcordance_CompProportions</i>, <i>GenotypeConcordance_Counts</i>, and <i>GenotypeConcordance_EvalProportions</i></li>
 *          <ul>
 *               <li>NO_CALL: reported genotype is ./., indicating not enough data to call</li>
 *               <li>HET: heterozygous</li>
 *               <li>HOM_REF: homozygous reference</li>
 *               <li>HOM_VAR: homozygous variant</li>
 *               <li>UNAVAILABLE: variant is not called in this callset</li>
 *               <li>MIXED: something like ./1</li>
 *          </ul>
 *          <li><i>GenotypeConcordance_Summary</i></li>
 *          <ul>
 *               <li>Non-Reference_Sensitivity (NRS): sensitivity of the EVAL calls to polymorphic calls in the COMP set, calculated by (# true positive)/(# true polymorphic)</li>
 *               <li>Non-Reference_Discrepancy (NRD): genotype discordance excluding concordant reference sites, calculated by (# discordant sites)/(total excluding # HOM_REF_HOM_REF) = 1.0-(# HOM_VAR_HOM_VAR + # HET_HET)/(total excluding # HOM_REF_HOM_REF)</li>
 *               <li>Overall_Genotype_Concordance: overall concordance calculated by (# concordant genotypes)/(# genotypes)</li>
 *          </ul>
 *          <li><i>SiteConcordance_Summary</i></li>
 *          <ul>
 *               <li>ALLELES_MATCH: counts of calls at the same site where the alleles match</li>
 *               <li>ALLELES_DO_NOT_MATCH: counts of calls at the same location with different alleles, such as the EVAL set calling a 'G' alternate allele, and the comp set calling a 'T' alternate allele</li>
 *               <li>EVAL_SUBSET_TRUTH: (multi-alleleic sites only) ALT alleles for EVAL are a subset of ALT alleles for COMP. See also below.</li>
 *               <li>EVAL_SUPERSET_TRUTH: (multi-allelic sites only) ALT alleles for COMP are a subset of ALT alleles for EVAL. See also below.</li>
 *               <li>EVAL_ONLY: counts of sites present only in the EVAL set, not in the COMP set</li>
 *               <li>TRUTH_ONLY: counts of sites present only in the COMP set, not in the EVAL set</li>
 *          </ul>
 * </ul>
 * </p>
 *
 * <h4>Site-level allelic concordance</h4>
 *
 * <p>
 *  For strictly bi-allelic VCFs, only the ALLELES_MATCH, EVAL_ONLY, TRUTH_ONLY fields will be populated,
 *  but where multi-allelic sites are involved counts for EVAL_SUBSET_TRUTH and EVAL_SUPERSET_TRUTH will be generated.
 * </p>
 * <p>
 *  For example, in the following situation
 *  <pre>
 *    eval:  ref - A   alt - C
 *    comp:  ref - A   alt - C,T
 *  </pre>
 *  then the site is tabulated as EVAL_SUBSET_TRUTH. Were the situation reversed, it would be EVAL_SUPERSET_TRUTH.
 *  However, in the case where EVAL has both C and T alternate alleles, both must be observed in the genotypes
 *  (that is, there must be at least one of (0/1,1/1) and at least one of (0/2,1/2,2/2) in the genotype field). If
 *  one of the alleles has no observations in the genotype fields of the EVAL, the site-level concordance is
 *  tabulated as though that allele were not present in the record.
 * </p>
 *
 *  <h4>Monomorphic Records</h4>
 *  <p>
 *  A site which has an alternate allele, but which is monomorphic in samples, is treated as not having been
 *  discovered, and will be recorded in the TRUTH_ONLY column (if a record exists in the COMP set), or not at all
 *  (if no record exists in the COMP set).
 *  </p>
 *  <p>
 *  That is, in the situation
 *  <pre>
 *   eval:  ref - A   alt - C   genotypes - 0/0  0/0  0/0 ... 0/0
 *   comp:  ref - A   alt - C   ...         0/0  0/0  ...
 *  </pre>
 *  is equivalent to
 *  <pre>
 *   eval:  ref - A   alt - .   genotypes - 0/0  0/0  0/0 ... 0/0
 *   comp:  ref - A   alt - C   ...         0/0  0/0  ...
 *  </pre>
 *  </p>
 *  <p>
 *  When a record is present in the COMP set the *genotypes* for the monomorphic site will still be used to evaluate
 *  per-sample genotype concordance counts.
 * </p>
 *
 *  <h4>Filtered Records</h4>
 *  Filtered records are treated as though they were not present in the VCF, unless -ignoreSiteFilters is provided,
 *  in which case all records are used. There is currently no way to assess concordance metrics on filtered sites
 *  exclusively. SelectVariants can be used to extract filtered sites, and VariantFiltration used to un-filter them.
 *
 * <h4>Moltenized tables</h4>
 *
 * <p>These tables may be optionally moltenized via the -moltenize argument. That is, the standard table
 *
 *  <pre>
 *  Sample   NO_CALL_HOM_REF  NO_CALL_HET  NO_CALL_HOM_VAR   (...)
 *  NA12878       0.003        0.001            0.000        (...)
 *  NA12891       0.005        0.000            0.000        (...)
 *  </pre>
 *
 *  would instead be displayed
 *
 * <pre>
 *  NA12878  NO_CALL_HOM_REF   0.003
 *  NA12878  NO_CALL_HET       0.001
 *  NA12878  NO_CALL_HOM_VAR   0.000
 *  NA12891  NO_CALL_HOM_REF   0.005
 *  NA12891  NO_CALL_HET       0.000
 *  NA12891  NO_CALL_HOM_VAR   0.000
 *  (...)
 *  </pre>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T GenotypeConcordance \
 *   -R reference.fasta \
 *   -eval test_set.vcf \
 *   -comp truth_set.vcf \
 *   -o output.grp
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class GenotypeConcordance extends RodWalker<List<Pair<VariantContext,VariantContext>>,ConcordanceMetrics> {

    /**
     * The callset you want to evaluate, typically this is where you'd put 'unassessed' callsets.
     */
    @Input(fullName="eval",shortName="eval",doc="The variants and genotypes to evaluate",required=true)
    RodBinding<VariantContext> evalBinding;

    /**
     * The callset you want to treat as 'truth'. Can also be of unknown quality for the sake of callset comparisons.
     */
    @Input(fullName="comp",shortName="comp",doc="The variants and genotypes to compare against",required=true)
    RodBinding<VariantContext> compBinding;

    /**
     * The FILTER field of the eval and comp VCFs will be ignored. If this flag is not included, all FILTER sites will
     * be treated as not being present in the VCF. (That is, the genotypes will be assigned UNAVAILABLE, as distinct
     * from NO_CALL).
     */
    @Argument(fullName="ignoreFilters",doc="Filters will be ignored",required=false)
    boolean ignoreFilters = false;

    /**
     * A genotype level JEXL expression to apply to eval genotypes. Genotypes filtered in this way will be replaced by NO_CALL.
     * For instance: -gfe 'GQ<20' will set to no-call any genotype with genotype quality less than 20.
     */
    @Argument(shortName="gfe", fullName="genotypeFilterExpressionEval", doc="One or more criteria to use to set EVAL genotypes to no-call. "+
            "These genotype-level filters are only applied to the EVAL rod.", required=false)
    public ArrayList<String> genotypeFilterExpressionsEval = new ArrayList<String>();

    /**
     * Identical to -gfe except the filter is applied to genotypes in the comp rod.
     */
    @Argument(shortName="gfc", fullName="genotypeFilterExpressionComp", doc="One or more criteria to use to set COMP genotypes to no-call. "+
            "These genotype-level filters are only applied to the COMP rod.", required=false)
    public ArrayList<String> genotypeFilterExpressionsComp = new ArrayList<String>();

    /**
     * Moltenize the count and proportion tables. Rather than moltenizing per-sample data into a 2x2 table, it is fully
     * moltenized into elements. That is, WITHOUT this argument, each row of the table begins with the sample name and
     * proceeds directly with counts/proportions of eval/comp counts (for instance HOM_REF/HOM_REF, HOM_REF/NO_CALL).
     *
     * If the Moltenize argument is given, the output will begin with a sample name, followed by the contrastive genotype
     * type (such as HOM_REF/HOM_REF), followed by the count or proportion. This will significantly increase the number of
     * rows.
     */
    @Argument(shortName="moltenize",fullName="moltenize",doc="Molten rather than tabular output")
    public boolean moltenize = false;

    /**
     * Print sites where genotypes are mismatched between callsets along with annotations giving the genotype of each callset
     * to the given filename
     *
     */

    @Argument(shortName = "sites",required = false,fullName = "printInterestingSites", doc="File to output the discordant sites and genotypes.")
    private PrintStream sitesFile = null;

    @Output
    PrintStream out;

    private List<String> evalSamples;
    private List<String> compSamples;
    private List<VariantContextUtils.JexlVCMatchExp> evalJexls = null;
    private List<VariantContextUtils.JexlVCMatchExp> compJexls = null;

    // todo -- table with "proportion of overlapping sites" (not just eval/comp margins) [e.g. drop no-calls]
    //  (this will break all the integration tests of course, due to new formatting)

    public void initialize() {
        evalJexls = initializeJexl(genotypeFilterExpressionsEval);
        compJexls = initializeJexl(genotypeFilterExpressionsComp);
    }

    private List<VariantContextUtils.JexlVCMatchExp> initializeJexl(ArrayList<String> genotypeFilterExpressions) {
        ArrayList<String> dummyNames = new ArrayList<String>(genotypeFilterExpressions.size());
        int expCount = 1;
        for ( String exp : genotypeFilterExpressions ) {
            dummyNames.add(String.format("gfe%d",expCount++));
        }
        return VariantContextUtils.initializeMatchExps(dummyNames, genotypeFilterExpressions);
    }

    public ConcordanceMetrics reduceInit() {
        Map<String,VCFHeader> headerMap = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(evalBinding,compBinding));
        VCFHeader evalHeader = headerMap.get(evalBinding.getName());
        evalSamples = evalHeader.getGenotypeSamples();
        VCFHeader compHeader = headerMap.get(compBinding.getName());
        compSamples = compHeader.getGenotypeSamples();
        return new ConcordanceMetrics(evalHeader,compHeader, sitesFile);
    }


    public List<Pair<VariantContext,VariantContext>> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<Pair<VariantContext,VariantContext>> evalCompPair = new ArrayList<Pair<VariantContext,VariantContext>>(3);
        if ( tracker != null && (
                tracker.getValues(evalBinding,ref.getLocus()).size() > 0 ||
                        tracker.getValues(compBinding,ref.getLocus()).size() > 0 ) ) {

            List<VariantContext> eval = tracker.getValues(evalBinding,ref.getLocus());
            List<VariantContext> comp = tracker.getValues(compBinding,ref.getLocus());
            if ( eval.size() > 1 || comp.size() > 1 ) {
                if ( noDuplicateTypes(eval) && noDuplicateTypes(comp) ) {
                    logger.info("Eval or Comp Rod at position " + ref.getLocus().toString() + " has multiple records. Resolving.");
                    evalCompPair = resolveMultipleRecords(eval,comp);
                } else {
                    logger.warn("Eval or Comp Rod at position "+ref.getLocus().toString()+" has multiple records of the same type. This locus will be skipped.");
                }
            } else {
                // if a rod is missing, explicitly create a variant context with 'missing' genotypes. Slow, but correct.
                // note that if there is no eval rod there must be a comp rod, and also the reverse
                VariantContext evalContext = eval.size() == 1 ? eval.get(0) : createEmptyContext(comp.get(0),evalSamples);
                VariantContext compContext = comp.size() == 1 ? comp.get(0) : createEmptyContext(eval.get(0),compSamples);
                evalContext = filterGenotypes(evalContext,ignoreFilters,evalJexls);
                compContext = filterGenotypes(compContext,ignoreFilters,compJexls);
                evalCompPair.add(new Pair<VariantContext, VariantContext>(evalContext,compContext));
            }
        }

        return evalCompPair;
    }

    private boolean noDuplicateTypes(List<VariantContext> vcList) {
        HashSet<VariantContext.Type> types = new HashSet<VariantContext.Type>(vcList.size());
        for ( VariantContext vc : vcList ) {
            VariantContext.Type type = vc.getType();
            if ( types.contains(type) )
                return false;
            types.add(type);
        }

        return true;
    }

    /**
     * The point of this method is to match up pairs of evals and comps by their type (or alternate alleles for mixed).
     * Basically multiple records could exist for a site such as:
     * Eval: 20   4000     A     C
     * Eval: 20   4000     A    AC
     * Comp: 20   4000     A     C
     * So for each eval, loop through the comps. If the types match, or for mixed types if eval alleles (non-emptily)
     * intersect the comp alleles, pair them up and remove that comp records.
     * Continue until we're out of evals or comps. This is n^2, but should rarely actually happen.
     *
     * The remaining unpaired records get paird with an empty contexts. So in the example above we'd get a list of:
     *  1 - (20,4000,A/C  |  20,4000,A/C)
     *  2 - (20,4000,A/AC |    Empty    )
     * @param evalList - list of eval variant contexts
     * @param compList - list of comp variant contexts
     * @return resolved pairs of the input lists
     */
    private List<Pair<VariantContext,VariantContext>> resolveMultipleRecords(List<VariantContext> evalList, List<VariantContext> compList) {
        List<Pair<VariantContext,VariantContext>> resolvedPairs = new ArrayList<Pair<VariantContext,VariantContext>>(evalList.size()+compList.size()); // oversized but w/e
        List<VariantContext> pairedEval = new ArrayList<VariantContext>(evalList.size());
        for ( VariantContext eval : evalList ) {
            VariantContext.Type evalType = eval.getType();
            Set<Allele> evalAlleles = new HashSet<Allele>(eval.getAlternateAlleles());
            VariantContext pairedComp = null;
            for ( VariantContext comp : compList ) {
                if ( evalType.equals(comp.getType()) ) {
                    pairedComp = comp;
                    break;
                } else if ( eval.isMixed() || comp.isMixed() ) {
                    for ( Allele compAllele : comp.getAlternateAlleles() ) {
                        if ( evalAlleles.contains(compAllele) ) {
                            pairedComp = comp;
                            break;
                        }
                    }
                }
            }
            if ( pairedComp != null ) {
                compList.remove(pairedComp);
                resolvedPairs.add(new Pair<VariantContext, VariantContext>(filterGenotypes(eval,ignoreFilters,evalJexls),filterGenotypes(pairedComp,ignoreFilters,compJexls)));
                pairedEval.add(eval);
                if ( compList.size() < 1 )
                    break;
            }
        }
        evalList.removeAll(pairedEval);
        for ( VariantContext unpairedEval : evalList ) {
            resolvedPairs.add(new Pair<VariantContext, VariantContext>(filterGenotypes(unpairedEval,ignoreFilters,evalJexls),createEmptyContext(unpairedEval,compSamples)));
        }

        for ( VariantContext unpairedComp : compList ) {
            resolvedPairs.add(new Pair<VariantContext, VariantContext>(createEmptyContext(unpairedComp,evalSamples),filterGenotypes(unpairedComp,ignoreFilters,compJexls)));
        }

        return resolvedPairs;
    }

    public ConcordanceMetrics reduce(List<Pair<VariantContext,VariantContext>> evalCompList, ConcordanceMetrics metrics) {
        for ( Pair<VariantContext,VariantContext> evalComp : evalCompList){
            metrics.update(evalComp.getFirst(),evalComp.getSecond());

        }
        return metrics;
    }

    private static double repairNaN(double d) {
     if ( Double.isNaN(d) ) {
      return 0.0;
     }
     return d;
    }

    public void onTraversalDone(ConcordanceMetrics metrics) {
        // todo -- this is over 200 lines of code just to format the output and could use some serious cleanup
        GATKReport report = new GATKReport();
        GATKReportTable concordanceCounts = new GATKReportTable("GenotypeConcordance_Counts","Per-sample concordance tables: comparison counts",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceEvalProportions = new GATKReportTable("GenotypeConcordance_EvalProportions", "Per-sample concordance tables: proportions of genotypes called in eval",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceCompProportions = new GATKReportTable("GenotypeConcordance_CompProportions", "Per-sample concordance tables: proportions of genotypes called in comp",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceSummary = new GATKReportTable("GenotypeConcordance_Summary","Per-sample summary statistics: NRS, NRD, and OGC",2);
        GATKReportTable siteConcordance = new GATKReportTable("SiteConcordance_Summary","Site-level summary statistics",ConcordanceMetrics.SiteConcordanceType.values().length);
        if ( moltenize ) {
            concordanceCompProportions.addColumn("Sample","%s");
            concordanceCounts.addColumn("Sample","%s");
            concordanceEvalProportions.addColumn("Sample","%s");
            concordanceSummary.addColumn("Sample","%s");

            concordanceCompProportions.addColumn("Eval_Genotype","%s");
            concordanceCounts.addColumn("Eval_Genotype","%s");
            concordanceEvalProportions.addColumn("Eval_Genotype","%s");
            concordanceSummary.addColumn("Non-Reference_Discrepancy","%.3f");

            concordanceCompProportions.addColumn("Comp_Genotype","%s");
            concordanceCounts.addColumn("Comp_Genotype","%s");
            concordanceEvalProportions.addColumn("Comp_Genotype","%s");
            concordanceSummary.addColumn("Non-Reference_Sensitivity","%.3f");

            concordanceCompProportions.addColumn("Proportion","%.3f");
            concordanceCounts.addColumn("Count","%d");
            concordanceEvalProportions.addColumn("Proportion","%.3f");
            concordanceSummary.addColumn("Overall_Genotype_Concordance","%.3f");

            for ( Map.Entry<String,ConcordanceMetrics.GenotypeConcordanceTable> entry : metrics.getPerSampleGenotypeConcordance().entrySet() ) {
                ConcordanceMetrics.GenotypeConcordanceTable table = entry.getValue();
                for ( GenotypeType evalType : GenotypeType.values() ) {
                    for ( GenotypeType compType : GenotypeType.values() ) {
                        String rowKey = String.format("%s_%s_%s",entry.getKey(),evalType.toString(),compType.toString());
                        concordanceCounts.set(rowKey,"Sample",entry.getKey());
                        concordanceCounts.set(rowKey,"Eval_Genotype",evalType.toString());
                        concordanceCounts.set(rowKey,"Comp_Genotype",compType.toString());
                        int count = table.get(evalType, compType);
                        concordanceCounts.set(rowKey,"Count",count);
                        if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR) {
                            concordanceEvalProportions.set(rowKey,"Sample",entry.getKey());
                            concordanceEvalProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                            concordanceEvalProportions.set(rowKey,"Comp_Genotype",compType.toString());
                            concordanceEvalProportions.set(rowKey,"Proportion",repairNaN(( (double) count)/table.getnEvalGenotypes(evalType)));
                        }
                        if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF ) {
                            concordanceCompProportions.set(rowKey,"Sample",entry.getKey());
                            concordanceCompProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                            concordanceCompProportions.set(rowKey,"Comp_Genotype",compType.toString());
                            concordanceCompProportions.set(rowKey,"Proportion",repairNaN(( (double) count)/table.getnCompGenotypes(compType)));
                        }
                    }
                }
                String mismatchKey = String.format("%s_%s",entry.getKey(),"Mismatching");
                concordanceCounts.set(mismatchKey,"Sample",entry.getKey());
                concordanceCounts.set(mismatchKey,"Eval_Genotype","Mismatching_Alleles");
                concordanceCounts.set(mismatchKey,"Comp_Genotype","Mismatching_Alleles");
                concordanceEvalProportions.set(mismatchKey,"Sample",entry.getKey());
                concordanceEvalProportions.set(mismatchKey,"Eval_Genotype","Mismatching_Alleles");
                concordanceEvalProportions.set(mismatchKey,"Comp_Genotype","Mismatching_Alleles");
                concordanceCompProportions.set(mismatchKey,"Sample",entry.getKey());
                concordanceCompProportions.set(mismatchKey,"Eval_Genotype","Mismatching_Alleles");
                concordanceCompProportions.set(mismatchKey,"Comp_Genotype","Mismatching_Alleles");
                concordanceEvalProportions.set(mismatchKey,"Proportion", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes()));
                concordanceCompProportions.set(mismatchKey,"Proportion", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes()));
                concordanceCounts.set(mismatchKey,"Count",table.getnMismatchingAlt());
            }

            String sampleKey = "ALL";
            ConcordanceMetrics.GenotypeConcordanceTable table = metrics.getOverallGenotypeConcordance();
            for ( GenotypeType evalType : GenotypeType.values() ) {
                for ( GenotypeType compType : GenotypeType.values() ) {
                    String rowKey = String.format("%s_%s_%s",sampleKey,evalType.toString(),compType.toString());
                    concordanceCounts.set(rowKey,"Sample",sampleKey);
                    concordanceCounts.set(rowKey,"Eval_Genotype",evalType.toString());
                    concordanceCounts.set(rowKey,"Comp_Genotype",compType.toString());
                    int count = table.get(evalType, compType);
                    concordanceCounts.set(rowKey,"Count",count);
                    if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR) {
                        concordanceEvalProportions.set(rowKey,"Sample",sampleKey);
                        concordanceEvalProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                        concordanceEvalProportions.set(rowKey,"Comp_Genotype",compType.toString());
                        concordanceEvalProportions.set(rowKey,"Proportion",repairNaN(( (double) count)/table.getnEvalGenotypes(evalType)));
                    }
                    if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF ) {
                        concordanceCompProportions.set(rowKey,"Sample",sampleKey);
                        concordanceCompProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                        concordanceCompProportions.set(rowKey,"Comp_Genotype",compType.toString());
                        concordanceCompProportions.set(rowKey,"Proportion",repairNaN(( (double) count)/table.getnCompGenotypes(compType)));
                    }
                }
            }
            String rowKey = String.format("%s_%s",sampleKey,"Mismatching");
            concordanceCounts.set(rowKey,"Sample",sampleKey);
            concordanceCounts.set(rowKey,"Eval_Genotype","Mismatching_Alleles");
            concordanceCounts.set(rowKey,"Comp_Genotype","Mismatching_Alleles");
            concordanceEvalProportions.set(rowKey,"Sample",sampleKey);
            concordanceEvalProportions.set(rowKey,"Eval_Genotype","Mismatching_Alleles");
            concordanceEvalProportions.set(rowKey,"Comp_Genotype","Mismatching_Alleles");
            concordanceCompProportions.set(rowKey,"Sample",sampleKey);
            concordanceCompProportions.set(rowKey,"Eval_Genotype","Mismatching_Alleles");
            concordanceCompProportions.set(rowKey,"Comp_Genotype","Mismatching_Alleles");
            concordanceEvalProportions.set(rowKey,"Proportion", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes()));
            concordanceCompProportions.set(rowKey,"Proportion", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes()));
            concordanceCounts.set(rowKey,"Count",table.getnMismatchingAlt());

            for ( Map.Entry<String,Double> nrsEntry : metrics.getPerSampleNRS().entrySet() ) {
                concordanceSummary.set(nrsEntry.getKey(),"Sample",nrsEntry.getKey());
                concordanceSummary.set(nrsEntry.getKey(),"Non-Reference_Sensitivity",nrsEntry.getValue());
            }
            for ( Map.Entry<String,Double> nrdEntry : metrics.getPerSampleNRD().entrySet() ) {
                concordanceSummary.set(nrdEntry.getKey(),"Non-Reference_Discrepancy",nrdEntry.getValue());
            }
            for ( Map.Entry<String,Double> ogcEntry : metrics.getPerSampleOGC().entrySet() ) {
                concordanceSummary.set(ogcEntry.getKey(),"Overall_Genotype_Concordance",ogcEntry.getValue());
            }
            concordanceSummary.set("ALL_NRS_NRD","Sample","ALL");
            concordanceSummary.set("ALL_NRS_NRD","Non-Reference_Sensitivity",metrics.getOverallNRS());
            concordanceSummary.set("ALL_NRS_NRD","Non-Reference_Discrepancy",metrics.getOverallNRD());
            concordanceSummary.set("ALL_NRS_NRD","Overall_Genotype_Concordance",metrics.getOverallOGC());


            for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
                siteConcordance.addColumn(type.toString(),"%d");
            }

            for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
                siteConcordance.set("Comparison",type.toString(),metrics.getOverallSiteConcordance().get(type));
            }

        } else {
            concordanceCompProportions.addColumn("Sample","%s");
            concordanceCounts.addColumn("Sample","%s");
            concordanceEvalProportions.addColumn("Sample","%s");
            concordanceSummary.addColumn("Sample","%s");
            for ( GenotypeType evalType : GenotypeType.values() ) {
                for ( GenotypeType compType : GenotypeType.values() ) {
                    String colKey = String.format("%s_%s", evalType.toString(), compType.toString());
                    concordanceCounts.addColumn(colKey,"%d");
                    if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                        concordanceEvalProportions.addColumn(colKey,"%.3f");
                    if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                        concordanceCompProportions.addColumn(colKey,"%.3f");
                }
            }
            concordanceEvalProportions.addColumn("Mismatching_Alleles","%.3f");
            concordanceCompProportions.addColumn("Mismatching_Alleles","%.3f");
            concordanceCounts.addColumn("Mismatching_Alleles","%d");
            concordanceSummary.addColumn("Non-Reference Sensitivity","%.3f");
            concordanceSummary.addColumn("Non-Reference Discrepancy","%.3f");
            concordanceSummary.addColumn("Overall_Genotype_Concordance","%.3f");
            for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
                siteConcordance.addColumn(type.toString(),"%d");
            }

            for ( Map.Entry<String,ConcordanceMetrics.GenotypeConcordanceTable> entry : metrics.getPerSampleGenotypeConcordance().entrySet() ) {
                ConcordanceMetrics.GenotypeConcordanceTable table = entry.getValue();
                concordanceEvalProportions.set(entry.getKey(),"Sample",entry.getKey());
                concordanceCompProportions.set(entry.getKey(),"Sample",entry.getKey());
                concordanceCounts.set(entry.getKey(),"Sample",entry.getKey());
                for ( GenotypeType evalType : GenotypeType.values() ) {
                    for ( GenotypeType compType : GenotypeType.values() ) {
                        String colKey = String.format("%s_%s",evalType.toString(),compType.toString());
                        int count = table.get(evalType, compType);
                        concordanceCounts.set(entry.getKey(),colKey,count);
                        if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                            concordanceEvalProportions.set(entry.getKey(),colKey,repairNaN(( (double) count)/table.getnEvalGenotypes(evalType)));
                        if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                            concordanceCompProportions.set(entry.getKey(),colKey,repairNaN(( (double) count)/table.getnCompGenotypes(compType)));
                    }
                }
                concordanceEvalProportions.set(entry.getKey(),"Mismatching_Alleles", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes()));
                concordanceCompProportions.set(entry.getKey(),"Mismatching_Alleles", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes()));
                concordanceCounts.set(entry.getKey(),"Mismatching_Alleles",table.getnMismatchingAlt());
            }

            String rowKey = "ALL";
            concordanceCompProportions.set(rowKey,"Sample",rowKey);
            concordanceEvalProportions.set(rowKey,"Sample",rowKey);
            concordanceCounts.set(rowKey,"Sample",rowKey);
            ConcordanceMetrics.GenotypeConcordanceTable table = metrics.getOverallGenotypeConcordance();
            for ( GenotypeType evalType : GenotypeType.values() ) {
                for ( GenotypeType compType : GenotypeType.values() ) {
                    String colKey = String.format("%s_%s",evalType.toString(),compType.toString());
                    int count = table.get(evalType,compType);
                    concordanceCounts.set(rowKey,colKey,count);
                    if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                        concordanceEvalProportions.set(rowKey,colKey,repairNaN(( (double) count)/table.getnEvalGenotypes(evalType)));
                    if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                        concordanceCompProportions.set(rowKey,colKey,repairNaN(( (double) count)/table.getnCompGenotypes(compType)));
                }
            }
            concordanceEvalProportions.set(rowKey,"Mismatching_Alleles", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes()));
            concordanceCompProportions.set(rowKey,"Mismatching_Alleles", repairNaN(( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes()));
            concordanceCounts.set(rowKey,"Mismatching_Alleles",table.getnMismatchingAlt());

            for ( Map.Entry<String,Double> nrsEntry : metrics.getPerSampleNRS().entrySet() ) {
                concordanceSummary.set(nrsEntry.getKey(),"Sample",nrsEntry.getKey());
                concordanceSummary.set(nrsEntry.getKey(),"Non-Reference Sensitivity",nrsEntry.getValue());
            }
            for ( Map.Entry<String,Double> nrdEntry : metrics.getPerSampleNRD().entrySet() ) {
                concordanceSummary.set(nrdEntry.getKey(),"Non-Reference Discrepancy",nrdEntry.getValue());
            }
            for ( Map.Entry<String,Double> ogcEntry : metrics.getPerSampleOGC().entrySet() ) {
                concordanceSummary.set(ogcEntry.getKey(),"Overall_Genotype_Concordance",ogcEntry.getValue());
            }
            concordanceSummary.set("ALL","Sample","ALL");
            concordanceSummary.set("ALL","Non-Reference Sensitivity",metrics.getOverallNRS());
            concordanceSummary.set("ALL","Non-Reference Discrepancy",metrics.getOverallNRD());
            concordanceSummary.set("ALL","Overall_Genotype_Concordance",metrics.getOverallOGC());

            for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
                siteConcordance.set("Comparison",type.toString(),metrics.getOverallSiteConcordance().get(type));
            }
        }

        report.addTable(concordanceCompProportions);
        report.addTable(concordanceEvalProportions);
        report.addTable(concordanceCounts);
        report.addTable(concordanceSummary);
        report.addTable(siteConcordance);

        report.print(out);
    }

    public VariantContext createEmptyContext(VariantContext other, List<String> samples) {
        VariantContextBuilder builder = new VariantContextBuilder();
        // set the alleles to be the same
        builder.alleles(other.getAlleles());
        builder.loc(other.getChr(),other.getStart(),other.getEnd());
        // set all genotypes to empty
        List<Genotype> genotypes = new ArrayList<Genotype>(samples.size());
        for ( String sample : samples )
            genotypes.add(GenotypeBuilder.create(sample, new ArrayList<Allele>(0)));
        builder.genotypes(genotypes);
        return builder.make();
    }

    public VariantContext filterGenotypes(VariantContext context, boolean ignoreSiteFilter, List<VariantContextUtils.JexlVCMatchExp> exps) {
        if ( ! context.isFiltered() || ignoreSiteFilter ) {
            List<Genotype> filteredGenotypes = new ArrayList<Genotype>(context.getNSamples());
            for ( Genotype g : context.getGenotypes() ) {
                Map<VariantContextUtils.JexlVCMatchExp, Boolean> matchMap = VariantContextUtils.match(context, g, exps);
                boolean filtered = false;
                for ( Boolean b : matchMap.values() ) {
                    if ( b ) {
                        filtered = true;
                        break;
                    }
                }
                if ( filtered ) {
                    filteredGenotypes.add(GenotypeBuilder.create(g.getSampleName(),Arrays.asList(Allele.NO_CALL,Allele.NO_CALL),g.getExtendedAttributes()));
                } else {
                    filteredGenotypes.add(g);
                }
            }
            VariantContextBuilder builder = new VariantContextBuilder(context);
            builder.genotypes(filteredGenotypes);
            return builder.make();
        }

        VariantContextBuilder builder = new VariantContextBuilder();
        builder.alleles(Arrays.asList(context.getReference()));
        builder.loc(context.getChr(),context.getStart(),context.getEnd());
        List<Genotype> newGeno = new ArrayList<Genotype>(context.getNSamples());
        for ( Genotype g : context.getGenotypes().iterateInSampleNameOrder() ) {
            newGeno.add(GenotypeBuilder.create(g.getSampleName(),new ArrayList<Allele>()));
        }
        builder.genotypes(newGeno);
        return builder.make();
    }
}
