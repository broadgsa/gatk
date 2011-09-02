package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.phasing.AllelePair;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.NewEvaluationContext;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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

    // a mapping from sample to stats
    @DataPoint(description = "the phasing statistics for each sample")
    SamplePhasingStatistics samplePhasingStatistics = null;

    SamplePreviousGenotypes samplePrevGenotypes = null;

    double minPhaseQuality = 10.0;

    public void initialize(VariantEvalWalker walker) {
        this.samplePhasingStatistics = new SamplePhasingStatistics(walker.getMinPhaseQuality());
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

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, NewEvaluationContext group) {
    //public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        Reasons interesting = new Reasons();
        if (ref == null)
            return interesting.toString();
        GenomeLoc curLocus = ref.getLocus();

        logger.debug("update2() locus: " + curLocus);
        logger.debug("comp = " + comp + " eval = " + eval);

        Set<String> allSamples = new HashSet<String>();

        Map<String, Genotype> compSampGenotypes = null;
        if (isRelevantToPhasing(comp)) {
            allSamples.addAll(comp.getSampleNames());
            compSampGenotypes = comp.getGenotypes();
        }

        Map<String, Genotype> evalSampGenotypes = null;
        if (isRelevantToPhasing(eval)) {
            allSamples.addAll(eval.getSampleNames());
            evalSampGenotypes = eval.getGenotypes();
        }

        for (String samp : allSamples) {
            logger.debug("sample = " + samp);

            Genotype compSampGt = null;
            if (compSampGenotypes != null)
                compSampGt = compSampGenotypes.get(samp);

            Genotype evalSampGt = null;
            if (evalSampGenotypes != null)
                evalSampGt = evalSampGenotypes.get(samp);

            if (compSampGt == null || evalSampGt == null) { // Since either comp or eval (or both) are missing the site, the best we can do is hope to preserve phase [if the non-missing one preserves phase]
                // Having an unphased site breaks the phasing for the sample [does NOT permit "transitive phasing"] - hence, must reset phasing knowledge for both comp and eval [put a null CompEvalGenotypes]:
                if (isNonNullButUnphased(compSampGt) || isNonNullButUnphased(evalSampGt))
                    samplePrevGenotypes.put(samp, null);
            }
            else { // Both comp and eval have a non-null Genotype at this site:
                AllelePair compAllelePair = new AllelePair(compSampGt);
                AllelePair evalAllelePair = new AllelePair(evalSampGt);

                boolean breakPhasing = false;
                if (compSampGt.isHet() != evalSampGt.isHet() || compSampGt.isHom() != evalSampGt.isHom())
                    breakPhasing = true; // since they are not both het or both hom
                else { // both are het, or both are hom:
                    boolean topMatchesTopAndBottomMatchesBottom = (topMatchesTop(compAllelePair, evalAllelePair) && bottomMatchesBottom(compAllelePair, evalAllelePair));
                    boolean topMatchesBottomAndBottomMatchesTop = (topMatchesBottom(compAllelePair, evalAllelePair) && bottomMatchesTop(compAllelePair, evalAllelePair));
                    if (!topMatchesTopAndBottomMatchesBottom && !topMatchesBottomAndBottomMatchesTop)
                        breakPhasing = true; // since the 2 VCFs have different diploid genotypes for this sample
                }

                if (breakPhasing) {
                    samplePrevGenotypes.put(samp, null); // nothing to do for this site, AND must remove any history for the future
                }
                else if (compSampGt.isHet() && evalSampGt.isHet()) {
                    /* comp and eval have the HET same Genotype at this site:
                       [Note that if both are hom, then nothing is done here, but the het history IS preserved].
                     */
                    CompEvalGenotypes prevCompAndEval = samplePrevGenotypes.get(samp);
                    if (prevCompAndEval != null && !prevCompAndEval.getLocus().onSameContig(curLocus)) // exclude curLocus if it is "phased" relative to a different chromosome
                        prevCompAndEval = null;

                    // Replace the previous hets with the current hets:
                    samplePrevGenotypes.put(samp, curLocus, compSampGt, evalSampGt);

                    if (prevCompAndEval != null) {
                        GenomeLoc prevLocus = prevCompAndEval.getLocus();
                        logger.debug("Potentially phaseable het locus: " + curLocus + " [relative to previous het locus: " + prevLocus + "]");
                        PhaseStats ps = samplePhasingStatistics.ensureSampleStats(samp);

                        boolean compSampIsPhased = genotypesArePhasedAboveThreshold(compSampGt);
                        boolean evalSampIsPhased = genotypesArePhasedAboveThreshold(evalSampGt);
                        if (compSampIsPhased || evalSampIsPhased) {
                            if (!evalSampIsPhased) {
                                ps.onlyCompPhased++;
                                //interesting.addReason("ONLY_COMP", samp, group, prevLocus, "");
                            }
                            else if (!compSampIsPhased) {
                                ps.onlyEvalPhased++;
                                //interesting.addReason("ONLY_EVAL", samp, group, prevLocus, "");
                            }
                            else { // both comp and eval are phased:                        
                                AllelePair prevCompAllelePair = new AllelePair(prevCompAndEval.getCompGenotpye());
                                AllelePair prevEvalAllelePair = new AllelePair(prevCompAndEval.getEvalGenotype());

                                // Sufficient to check only the top of comp, since we ensured that comp and eval have the same diploid genotypes for this sample:
                                boolean topsMatch = (topMatchesTop(prevCompAllelePair, prevEvalAllelePair) && topMatchesTop(compAllelePair, evalAllelePair));
                                boolean topMatchesBottom = (topMatchesBottom(prevCompAllelePair, prevEvalAllelePair) && topMatchesBottom(compAllelePair, evalAllelePair));

                                if (topsMatch || topMatchesBottom) {
                                    ps.phasesAgree++;

                                    Double compPQ = getPQ(compSampGt);
                                    Double evalPQ = getPQ(evalSampGt);
                                    if (compPQ != null && evalPQ != null && MathUtils.compareDoubles(compPQ, evalPQ) != 0) {
                                        //interesting.addReason("PQ_CHANGE", samp, group, prevLocus, compPQ + " -> " + evalPQ);
                                    }
                                }
                                else {
                                    ps.phasesDisagree++;
                                    logger.debug("SWITCHED locus: " + curLocus);
                                    //interesting.addReason("SWITCH", samp, group, prevLocus, toString(prevCompAllelePair, compAllelePair) + " -> " + toString(prevEvalAllelePair, evalAllelePair));
                                }
                            }
                        }
                        else {
                            ps.neitherPhased++;
                        }
                    }
                }
            }
        }
        logger.debug("\n" + samplePhasingStatistics + "\n");

        return interesting.toString();
    }

    public static boolean isRelevantToPhasing(VariantContext vc) {
        return (vc != null && !vc.isFiltered());
    }

    public boolean isNonNullButUnphased(Genotype gt) {
        return (gt != null && !genotypesArePhasedAboveThreshold(gt));
    }

    public boolean genotypesArePhasedAboveThreshold(Genotype gt) {
        if (gt.isHom()) // Can always consider a hom site to be phased to its predecessor, since its successor will only be phased to it if it's hom or "truly" phased
            return true;

        if (!gt.isPhased())
            return false;

        Double pq = getPQ(gt);
        return (pq == null || pq >= minPhaseQuality);
    }

    public static Double getPQ(Genotype gt) {
        Double d = gt.getAttributeAsDouble(ReadBackedPhasingWalker.PQ_KEY, -1);
        return d == -1 ? null : d;
    }

    public static boolean topMatchesTop(AllelePair b1, AllelePair b2) {
        return b1.getTopAllele().equals(b2.getTopAllele());
    }

    public static boolean topMatchesBottom(AllelePair b1, AllelePair b2) {
        return b1.getTopAllele().equals(b2.getBottomAllele());
    }

    public static boolean bottomMatchesTop(AllelePair b1, AllelePair b2) {
        return topMatchesBottom(b2, b1);
    }

    public static boolean bottomMatchesBottom(AllelePair b1, AllelePair b2) {
        return b1.getBottomAllele().equals(b2.getBottomAllele());
    }

    public String toString(AllelePair prev, AllelePair cur) {
        return prev.getTopAllele().getBaseString() + "+" + cur.getTopAllele().getBaseString() + "|" + prev.getBottomAllele().getBaseString() + "+" + cur.getBottomAllele().getBaseString();
    }

    public void finalizeEvaluation() {
    }

    private static class Reasons {
        private StringBuilder sb;

        public Reasons() {
            sb = new StringBuilder();
        }

//        public void addReason(String category, String sample, VariantEvalWalker.EvaluationContext evalGroup, GenomeLoc prevLoc, String reason) {
//             sb.append(category + "(" + sample + ", previous: " + prevLoc + " [" + evalGroup.compTrackName + ", " + evalGroup.evalTrackName + "]): " + reason + ";");
//        }

        public String toString() {
            if (sb.length() == 0)
                return null;

            return "reasons=" + sb.toString();
        }
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
    private double minPhaseQuality;

    public SamplePhasingStatistics(double minPhaseQuality) {
        this.sampleStats = new HashMap<String, PhaseStats>();
        this.minPhaseQuality = minPhaseQuality;
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
        return "Sample Phasing Statistics (for PQ >= " + minPhaseQuality + ")";
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
