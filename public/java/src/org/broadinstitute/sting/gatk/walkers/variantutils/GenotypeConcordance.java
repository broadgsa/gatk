/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.PrintStream;
import java.util.*;

/**
 * A simple walker for performing genotype concordance calculations between two callsets. Outputs a GATK table with
 * per-sample and aggregate counts and frequencies, a summary table for NRD/NRS, and a table for site allele overlaps.
 *
 * <p>
 *  Genotype concordance takes in two callsets (vcfs) and tabulates the number of sites which overlap and share alleles,
 *  and for each sample, the genotype-by-genotype counts (for instance, the number of sites at which a sample was
 *  called homozygous reference in the EVAL callset, but homozygous variant in the COMP callset). It outputs these
 *  counts as well as convenient proportions (such as the proportion of het calls in the EVAL which were called REF in
 *  the COMP) and metrics (such as NRD and NRS).
 *
 *  <h2> INPUT </h2>
 *  <p>
 *  Genotype concordance requires two callsets (as it does a comparison): an EVAL and a COMP callset, specified via
 *  the -eval and -comp arguments
 *  <p>
 *  (Optional) Jexl expressions for genotype-level filtering of EVAL or COMP genotypes, specified via the -gfe and
 *  -cfe arguments, respectively.
 *
 *  <h2> OUTPUT </h2>
 *  Genotype Concordance writes a GATK report to the specified (via -o) file, consisting of multiple tables of counts
 *  and proportions. These tables may be optionally moltenized via the -moltenize argument.
 *
 */
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
        return new ConcordanceMetrics(evalHeader,compHeader);
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
        for ( Pair<VariantContext,VariantContext> evalComp : evalCompList)
            metrics.update(evalComp.getFirst(),evalComp.getSecond());
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
        GATKReportTable concordanceSummary = new GATKReportTable("GenotypeConcordance_Summary","Per-sample summary statistics: NRS and NRD",2);
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

            for ( Map.Entry<String,ConcordanceMetrics.GenotypeConcordanceTable> entry : metrics.getPerSampleGenotypeConcordance().entrySet() ) {
                ConcordanceMetrics.GenotypeConcordanceTable table = entry.getValue();
                for ( GenotypeType evalType : GenotypeType.values() ) {
                    for ( GenotypeType compType : GenotypeType.values() ) {
                        String rowKey = String.format("%s_%s_%s",entry.getKey(),evalType.toString(),compType.toString());
                        concordanceCounts.set(rowKey,"Sample",entry.getKey());
                        concordanceCounts.set(rowKey,"Eval_Genotype",evalType.toString());
                        concordanceCounts.set(rowKey,"Comp_Genotype",evalType.toString());
                        int count = table.get(evalType, compType);
                        concordanceCounts.set(rowKey,"Count",count);
                        if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR) {
                            concordanceEvalProportions.set(rowKey,"Sample",entry.getKey());
                            concordanceEvalProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                            concordanceEvalProportions.set(rowKey,"Comp_Genotype",evalType.toString());
                            concordanceEvalProportions.set(rowKey,"Proportion",repairNaN(( (double) count)/table.getnEvalGenotypes(evalType)));
                        }
                        if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF ) {
                            concordanceCompProportions.set(rowKey,"Sample",entry.getKey());
                            concordanceCompProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                            concordanceCompProportions.set(rowKey,"Comp_Genotype",evalType.toString());
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
                    concordanceCounts.set(rowKey,"Comp_Genotype",evalType.toString());
                    int count = table.get(evalType, compType);
                    concordanceCounts.set(rowKey,"Count",count);
                    if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR) {
                        concordanceEvalProportions.set(rowKey,"Sample",sampleKey);
                        concordanceEvalProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                        concordanceEvalProportions.set(rowKey,"Comp_Genotype",evalType.toString());
                        concordanceEvalProportions.set(rowKey,"Proportion",repairNaN(( (double) count)/table.getnEvalGenotypes(evalType)));
                    }
                    if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF ) {
                        concordanceCompProportions.set(rowKey,"Sample",sampleKey);
                        concordanceCompProportions.set(rowKey,"Eval_Genotype",evalType.toString());
                        concordanceCompProportions.set(rowKey,"Comp_Genotype",evalType.toString());
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
            concordanceSummary.set("ALL_NRS_NRD","Sample","ALL");
            concordanceSummary.set("ALL_NRS_NRD","Non-Reference_Sensitivity",metrics.getOverallNRS());
            concordanceSummary.set("ALL_NRS_NRD","Non-Reference_Discrepancy",metrics.getOverallNRD());


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
            concordanceSummary.set("ALL","Sample","ALL");
            concordanceSummary.set("ALL","Non-Reference Sensitivity",metrics.getOverallNRS());
            concordanceSummary.set("ALL","Non-Reference Discrepancy",metrics.getOverallNRD());

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
