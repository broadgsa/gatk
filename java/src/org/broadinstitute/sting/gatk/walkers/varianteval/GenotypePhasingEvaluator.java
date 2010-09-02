package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.playground.gatk.walkers.Biallele;
import org.broadinstitute.sting.playground.gatk.walkers.BialleleSNP;
import org.broadinstitute.sting.playground.gatk.walkers.ReadBackedPhasingWalker;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.*;

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

@Analysis(name = "Genotype Phasing Evaluation", description = "Evaluates the phasing of genotypes in different tracks")
public class GenotypePhasingEvaluator extends VariantEvaluator {
    protected final static Logger logger = Logger.getLogger(GenotypePhasingEvaluator.class);

    //
    //@Argument(fullName = "phaseQualityThresh", shortName = "phaseThresh", doc = "The minimum phasing quality score required to consider eval track phasing; [default:20.0]", required = false)
    //
    protected Double phaseQualityThresh = 20.0; // PQ = 20.0 <=> P(error) = 10^(-20/10) = 0.01, P(correct) = 0.99
    //
    //
    
    private VariantEvalWalker.EvaluationContext group = null;

    // a mapping from sample to stats
    @DataPoint(name = "samples", description = "the phasing statistics for each sample")
    SamplePhasingStatistics samplePhasingStatistics = null;

    SamplePreviousGenotypes samplePrevGenotypes = null;

    //
    //
    //
    HashMap<String, Genotype> prevCompGenotypes = new HashMap<String, Genotype>();
    HashMap<String, Genotype> prevEvalGenotypes = new HashMap<String, Genotype>();
    //
    //
    //

    public GenotypePhasingEvaluator(VariantEvalWalker parent) {
        super(parent);
        this.samplePhasingStatistics = new SamplePhasingStatistics(phaseQualityThresh);
        this.samplePrevGenotypes = new SamplePreviousGenotypes();
    }

    public String getName() {
        return "GenotypePhasingEvaluator";
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see pairs of (comp, eval)
    }

    public boolean enabled() {
        return true;
    }

    public String toString() {
        return getName() + ": <table>";
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        String interesting = null;
        if (ref == null)
            return interesting;
        GenomeLoc curLocus = ref.getLocus();
        this.group = group;

        logger.debug("update2() locus: " + curLocus);
        logger.debug("comp = " + comp + " eval = " + eval);

        if (comp == null || eval == null || comp.isFiltered()) // ignore existence of filtered variants in comp [NOTE that eval will be checked BOTH filtered and unfiltered (by VariantEvalWalker)]
            return interesting;

        if (!comp.isBiallelic() || !eval.isBiallelic() || !comp.getAlternateAllele(0).equals(eval.getAlternateAllele(0))) // these are not the same biallelic variants
            return interesting;

        Map<String, Genotype> compSampGenotypes = comp.getGenotypes();
        Map<String, Genotype> evalSampGenotypes = eval.getGenotypes();

        Set<String> allSamples = new HashSet<String>();
        allSamples.addAll(comp.getSampleNames());
        allSamples.addAll(eval.getSampleNames());

        for (String samp : allSamples) {
            logger.debug("sample = " + samp);

            Genotype compSampGt = compSampGenotypes.get(samp);
            Genotype evalSampGt = evalSampGenotypes.get(samp);
            if (compSampGt != null && evalSampGt != null && compSampGt.isHet() && evalSampGt.isHet()) {
                CompEvalGenotypes prevCompAndEval = samplePrevGenotypes.get(samp);
                if (prevCompAndEval != null && !prevCompAndEval.getLocus().onSameContig(curLocus)) // exclude curLocus if it is "phased" relative to a different chromosome
                    prevCompAndEval = null;

                // Replace the previous with current:
                samplePrevGenotypes.put(samp, curLocus, compSampGt, evalSampGt);
                if (prevCompAndEval != null) {
                    logger.debug("Potentially phaseable locus: " + curLocus);
                    Genotype prevCompSampGt = prevCompAndEval.getCompGenotpye();
                    Genotype prevEvalSampGt = prevCompAndEval.getEvalGenotype();

                    PhaseStats ps = samplePhasingStatistics.ensureSampleStats(samp);

                    boolean compSampIsPhased = genotypesArePhasedAboveThreshold(compSampGt);
                    boolean evalSampIsPhased = genotypesArePhasedAboveThreshold(evalSampGt);
                    if (compSampIsPhased || evalSampIsPhased) {
                        if (!evalSampIsPhased)
                            ps.onlyCompPhased++;
                        else if (!compSampIsPhased)
                            ps.onlyEvalPhased++;
                        else {
                            Biallele prevCompBiallele = new Biallele(prevCompSampGt);
                            Biallele compBiallele = new Biallele(compSampGt);

                            Biallele prevEvalBiallele = new Biallele(prevEvalSampGt);
                            Biallele evalBiallele = new Biallele(evalSampGt);

                            boolean topsMatch = (prevCompBiallele.getTopAllele().equals(prevEvalBiallele.getTopAllele()) && compBiallele.getTopAllele().equals(evalBiallele.getTopAllele()));
                            boolean topMatchesBottom = (prevCompBiallele.getTopAllele().equals(prevEvalBiallele.getBottomAllele()) && compBiallele.getTopAllele().equals(evalBiallele.getBottomAllele()));

                            if (topsMatch || topMatchesBottom) {
                                ps.phasesAgree++;
                            }
                            else {
                                ps.phasesDisagree++;
                            }
                        }
                    }
                    else {
                        ps.neitherPhased++;
                    }
                }
            }
            else { // This site (either non-existent or homozygous for at least one of them) breaks the phase of the next position:
                samplePrevGenotypes.put(samp, null);
            }
        }
        logger.debug("\n" + samplePhasingStatistics + "\n");

        return interesting;
    }

    public boolean genotypesArePhasedAboveThreshold(Genotype gt) {
        if (!gt.genotypesArePhased())
            return false;

        Object pq = gt.getAttributes().get("PQ");
        return (pq == null || (new Double(pq.toString()) >= phaseQualityThresh));
    }

    //
    //
    //
    public String update3(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        this.group = group;

        String interesting = null;
        if (comp == null && eval == null)
            return interesting;

        if (comp != null && comp.isFiltered()) // ignore existence of filtered variants in comp [NOTE that eval will be checked BOTH filtered and unfiltered (by VariantEvalWalker)]
            return interesting;

        Set<String> allSamples = new HashSet<String>();

        boolean compIsUsable = (comp != null && comp.isBiallelic());
        boolean evalIsUsable = (eval != null && eval.isBiallelic());
        if (compIsUsable && evalIsUsable && !comp.getAlternateAllele(0).equals(eval.getAlternateAllele(0))) { // these are not the same biallelic variants
            compIsUsable = false;
            evalIsUsable = false;
        }

        Map<String, Genotype> compSampGenotypes = null;
        if (compIsUsable) {
            compSampGenotypes = comp.getGenotypes();
            allSamples.addAll(comp.getSampleNames());
        }

        Map<String, Genotype> evalSampGenotypes = null;
        if (evalIsUsable) {
            evalSampGenotypes = eval.getGenotypes();
            allSamples.addAll(eval.getSampleNames());
        }

        // Note that missing or incompatible VariantContext (e.g., comp) might break the phase of the other VariantContext (e.g., eval):
        for (String samp : allSamples) {
            Genotype compSampGt = null;
            if (compSampGenotypes != null)
                compSampGt = compSampGenotypes.get(samp);
            boolean compIsHetForSample = (compIsUsable && compSampGt != null && compSampGt.isHet());

            Genotype evalSampGt = null;
            if (evalSampGenotypes != null)
                evalSampGt = evalSampGenotypes.get(samp);
            boolean evalIsHetForSample = (evalIsUsable && evalSampGt != null && evalSampGt.isHet());

            // Only compare phasing at het positions after het positions:
            if (compIsHetForSample && evalIsHetForSample) {
                

                // Put het genotypes into previous Genotypes:
                prevCompGenotypes.put(samp, compSampGt);
                prevEvalGenotypes.put(samp, evalSampGt);
            }
            else { // at least one is not a biallelic het site for samp:
                breakPhaseForUnphasedOrHomozygousSite(evalSampGt, samp, prevEvalGenotypes);
                breakPhaseForUnphasedOrHomozygousSite(compSampGt, samp, prevCompGenotypes);
            }
        }

        /*
        Biallele compBiallele = new Biallele();
        Biallele evalBiallele = null;
        */

        //interesting = determineStats(eval, validation);

        return interesting; // we don't capture any interesting sites
    }
    //
    //
    //

    private void breakPhaseForUnphasedOrHomozygousSite(Genotype gt, String samp, HashMap<String, Genotype> prevVcGenotypes) {
        if (gt == null)
            return;

        // Break the phase of gt [i.e., set the previous Genotype to null] since it is unphased (or homozygous):
        if (!genotypesArePhasedAboveThreshold(gt) || gt.isHom()) // otherwise, a phased het site lets the phase pass through
            prevVcGenotypes.put(samp, null);
    }

    public void finalizeEvaluation() {
    }
}



class CompEvalGenotypes {
    private GenomeLoc loc;
    private Genotype compGt;
    private Genotype evalGt;

    public CompEvalGenotypes(GenomeLoc loc, Genotype compGt, Genotype evalGt) {
        this.loc = loc;
        this.compGt = compGt;
        this.evalGt = evalGt;
    }

    public GenomeLoc getLocus() {
        return loc;
    }

    public Genotype getCompGenotpye() {
        return compGt;
    }
    public Genotype getEvalGenotype() {
        return evalGt;
    }

    public void setCompGenotype(Genotype compGt) {
        this.compGt = compGt;
    }

    public void setEvalGenotype(Genotype evalGt) {
        this.evalGt = evalGt;
    }
}

class SamplePreviousGenotypes {
    private HashMap<String, CompEvalGenotypes> sampleGenotypes = null;

    public SamplePreviousGenotypes() {
        this.sampleGenotypes = new HashMap<String, CompEvalGenotypes>();
    }

    public CompEvalGenotypes get(String sample) {
        return sampleGenotypes.get(sample);
    }

    public void put(String sample, CompEvalGenotypes compEvalGts) {
        sampleGenotypes.put(sample, compEvalGts);
    }

    public void put(String sample, GenomeLoc locus, Genotype compGt, Genotype evalGt) {
        sampleGenotypes.put(sample, new CompEvalGenotypes(locus, compGt, evalGt));
    }
}

class PhaseStats {
    public int neitherPhased;
    public int onlyCompPhased;
    public int onlyEvalPhased;
    public int phasesAgree;
    public int phasesDisagree;

    public PhaseStats() {
        this.neitherPhased = 0;
        this.onlyCompPhased = 0;
        this.onlyEvalPhased = 0;
        this.phasesAgree = 0;
        this.phasesDisagree = 0;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Neither phased: " + neitherPhased + "\tOnly Comp: " + onlyCompPhased + "\tOnly Eval: " + onlyEvalPhased + "\tSame phase: " + phasesAgree + "\tOpposite phase: " + phasesDisagree);
        return sb.toString();
    }

    public static String[] getFieldNamesArray() {
        return new String[]{"total", "neither", "only_comp", "only_eval", "both", "match", "switch", "switch_rate"};
    }

    public Object getField(int index) {
        switch (index) {
            case (0):
                return (neitherPhased + onlyCompPhased + onlyEvalPhased + phasesAgree + phasesDisagree);
            case (1):
                return neitherPhased;
            case (2):
                return onlyCompPhased;
            case (3):
                return onlyEvalPhased;
            case (4):
                return (phasesAgree + phasesDisagree);
            case (5):
                return phasesAgree;
            case (6):
                return phasesDisagree;
            case (7):
                return ((phasesDisagree == 0) ? 0 : ((double) phasesDisagree) / (phasesAgree + phasesDisagree));
            default:
                return -1;
        }
    }
}

/**
 * a table of sample names to genotype phasing statistics
 */
class SamplePhasingStatistics implements TableType {
    private HashMap<String, PhaseStats> sampleStats = null;
    private double phaseQualityThresh;

    public SamplePhasingStatistics(double phaseQualityThresh) {
        this.sampleStats = new HashMap<String, PhaseStats>();
        this.phaseQualityThresh = phaseQualityThresh;
    }

    public PhaseStats ensureSampleStats(String samp) {
        PhaseStats ps = sampleStats.get(samp);
        if (ps == null) {
            ps = new PhaseStats();
            sampleStats.put(samp, ps);
        }
        return ps;
    }

    /**
     * @return one row per sample
     */
    public String[] getRowKeys() {
        return sampleStats.keySet().toArray(new String[sampleStats.size()]);
    }

    /**
     * get the column keys
     *
     * @return a list of objects, in this case strings, that are the column names
     */
    public String[] getColumnKeys() {
        return PhaseStats.getFieldNamesArray();
    }

    public Object getCell(int x, int y) {
        String[] rowKeys = getRowKeys();
        PhaseStats ps = sampleStats.get(rowKeys[x]);
        return ps.getField(y);
    }

    public String getName() {
        return "Sample Phasing Statistics (for PQ >= " + phaseQualityThresh + ")";
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, PhaseStats> sampPhaseStatsEnt : sampleStats.entrySet()) {
            String sample = sampPhaseStatsEnt.getKey();
            PhaseStats ps = sampPhaseStatsEnt.getValue();

            sb.append(sample + "\t" + ps);
        }
        return sb.toString();
    }
}