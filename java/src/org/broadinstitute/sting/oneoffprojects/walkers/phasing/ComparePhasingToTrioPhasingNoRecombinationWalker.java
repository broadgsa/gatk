/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers.phasing;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.phasing.AllelePair;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
import org.broadinstitute.sting.gatk.walkers.phasing.WriteVCF;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.GenotypePhasingEvaluator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

import static org.broadinstitute.sting.utils.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = ComparePhasingToTrioPhasingNoRecombinationWalker.TRIO_ROD_NAME, type = ReferenceOrderedDatum.class), @RMD(name = ComparePhasingToTrioPhasingNoRecombinationWalker.PHASING_ROD_NAME, type = ReferenceOrderedDatum.class)})

@ReadFilters({ZeroMappingQualityReadFilter.class})
// Filter out all reads with zero mapping quality

public class ComparePhasingToTrioPhasingNoRecombinationWalker extends RodWalker<CompareResult, CompareToTrioPhasingStats> {
    public final static String TRIO_ROD_NAME = "trio";
    public final static String PHASING_ROD_NAME = "phasing";

    private final static int NUM_IN_TRIO = 3;

    private final static int DIPLOID = 2;

    @Output
    protected PrintStream out;

    @Argument(fullName = "trioAugmentedPhasing", shortName = "trioAugmentedPhasing", doc = "File to which trio-phased variants should be written", required = false)
    protected VCFWriter writer = null;
    
    @Argument(fullName = "diffTrioAndPhasingTracks", shortName = "diffTrioAndPhasingTracks", doc = "File to which comparisons of phasing information in 'trio' and 'phasing' tracks should be written", required = false)
    protected PrintStream diffTrioAndPhasingTracks = null;

    private CompareTrioAndPhasingTracks diffTrioAndPhasingCounts = null;

    private String phasingSample = null;

    private enum TrioStatus {
        PRESENT, MISSING, TRIPLE_HET
    }

    private GenomeLoc prevLoc = null;
    private VariantContext prevTrioVc = null;
    private TrioStatus prevTrioStatus = TrioStatus.MISSING;

    private Genotype prevPhasingGt = null;


    public void initialize() {
        initializeVcfWriter();

        // Will compare the phasing ALREADY present in the trio track [without regards to what this trio phasing mechanism (without recombination) would do]:
        if (diffTrioAndPhasingTracks != null)
            diffTrioAndPhasingCounts = new CompareTrioAndPhasingTracks();
    }

    private void initializeVcfWriter() {
        if (writer == null)
            return;

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        List<String> rodNames = new LinkedList<String>();
        rodNames.add(PHASING_ROD_NAME);
        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> samples = new TreeSet<String>(rodNameToHeader.get(PHASING_ROD_NAME).getGenotypeSamples());
        writer.writeHeader(new VCFHeader(hInfo, samples));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public CompareToTrioPhasingStats reduceInit() {
        return new CompareToTrioPhasingStats();
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public CompareResult map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLoc = ref.getLocus();
        VariantContext phasingVc = tracker.getVariantContext(ref, PHASING_ROD_NAME, curLoc);

        CompareToTrioPhasingStats stats = new CompareToTrioPhasingStats();
        CompareResult result = new CompareResult(phasingVc, stats);

        if (phasingVc == null || phasingVc.isFiltered())
            return result;

        Map<String, Genotype> phasingSampleToGt = phasingVc.getGenotypes();
        if (phasingSampleToGt.size() != 1)
            throw new UserException("Must provide EXACTLY one sample in " + PHASING_ROD_NAME + " track!");
        Map.Entry<String, Genotype> phasingSampGt = phasingSampleToGt.entrySet().iterator().next();
        String sample = phasingSampGt.getKey();
        if (phasingSample == null)
            phasingSample = sample;
        if (!sample.equals(phasingSample))
            throw new UserException("Must provide EXACTLY one sample!");
        Genotype curPhasingGt = phasingSampGt.getValue();
        if (!curPhasingGt.isHet())
            return result;

        VariantContext curTrioVc = tracker.getVariantContext(ref, TRIO_ROD_NAME, curLoc);
        boolean useTrioVc = (curTrioVc != null && !curTrioVc.isFiltered());

        Genotype sampleCurGtInTrio = null;
        if (useTrioVc) {
            sampleCurGtInTrio = curTrioVc.getGenotype(phasingSample);

            if (curTrioVc.getNSamples() > NUM_IN_TRIO || sampleCurGtInTrio == null)
                throw new UserException("Must provide trio data for sample: " + phasingSample);

            if (!new TreeSet<Allele>(curPhasingGt.getAlleles()).equals(new TreeSet<Allele>(sampleCurGtInTrio.getAlleles()))) {
                logger.warn("Locus " + curLoc + " breaks phase, since " + PHASING_ROD_NAME + " and " + TRIO_ROD_NAME + " tracks have different genotypes for " + phasingSample + "!");
                prevLoc = null;
                return result;
            }
        }

        // Now, we have a [trio-consistent] het genotype that may be phased or not [and we want to know if it could be phased based on trio information]:
        int processed = 1;

        TrioStatus currentTrioStatus = TrioStatus.MISSING;
        if (useTrioVc && curTrioVc.getNSamples() == NUM_IN_TRIO) {
            boolean allHet = true;
            for (int i = 0; i < NUM_IN_TRIO; i++) {
                if (!curTrioVc.getGenotype(i).isHet()) {
                    allHet = false;
                    break;
                }
            }

            if (allHet)
                currentTrioStatus = TrioStatus.TRIPLE_HET;
            else
                currentTrioStatus = TrioStatus.PRESENT;
        }

        if (prevLoc != null && curLoc.onSameContig(prevLoc)) {
            String trioPhaseStatus;
            stats.comparedSites++;
            String addToOutput = "";

            if (prevTrioStatus == TrioStatus.TRIPLE_HET || currentTrioStatus == TrioStatus.TRIPLE_HET) {
                trioPhaseStatus = "Het3";
            }
            else if (prevTrioStatus == TrioStatus.MISSING || currentTrioStatus == TrioStatus.MISSING) {
                trioPhaseStatus = "Missing";
            }
            else {
                if (prevTrioStatus != TrioStatus.PRESENT || currentTrioStatus != TrioStatus.PRESENT)
                    throw new ReviewedStingException("LOGICAL error: prevTrioStatus != TrioStatus.PRESENT || currentTrioStatus != TrioStatus.PRESENT");

                trioPhaseStatus = "trio_phased";
                stats.trioPhaseableSites++;

                if (writer != null) { // Phase the genotypes using the trio information:
                    String parent1 = null;
                    String parent2 = null;
                    for (Map.Entry<String, Genotype> trioEntry : curTrioVc.getGenotypes().entrySet()) {
                        String trioSample = trioEntry.getKey();
                        if (trioEntry.getValue().getPloidy() != DIPLOID)
                            throw new UserException("Each sample in trio must be diploid!");
                        if (trioSample.equals(phasingSample))
                            continue;

                        if (parent1 == null)
                            parent1 = trioSample;
                        else if (parent2 == null)
                            parent2 = trioSample;
                        else
                            throw new ReviewedStingException("Cannot be more than 2 parents in TRIO!");
                    }
                    if (parent1 == null || parent2 == null)
                        throw new ReviewedStingException("Must have 2 parents in TRIO!");

                    Genotype samplePrevGtInTrio = prevTrioVc.getGenotype(phasingSample);

                    Genotype parent1PrevGt = prevTrioVc.getGenotype(parent1);
                    Genotype parent1CurGt = curTrioVc.getGenotype(parent1);

                    Genotype parent2PrevGt = prevTrioVc.getGenotype(parent2);
                    Genotype parent2CurGt = curTrioVc.getGenotype(parent2);

                    int prevHomIndex, prevOtherIndex;
                    Allele prevHomAllele;
                    Set<Allele> prevOtherAlleles;
                    if (parent1PrevGt.isHom()) {
                        prevHomIndex = 1;
                        prevOtherIndex = 2;
                        prevHomAllele = parent1PrevGt.getAllele(0);
                        prevOtherAlleles = new TreeSet<Allele>(parent2PrevGt.getAlleles());
                    }
                    else if (parent2PrevGt.isHom()) {
                        prevHomIndex = 2;
                        prevOtherIndex = 1;
                        prevHomAllele = parent2PrevGt.getAllele(0);
                        prevOtherAlleles = new TreeSet<Allele>(parent1PrevGt.getAlleles());
                    }
                    else
                        throw new ReviewedStingException("LOGICAL ERROR: at least one parent is hom!");

                    int curHomIndex, curOtherIndex;
                    Allele curHomAllele;
                    Set<Allele> curOtherAlleles;
                    if (parent1CurGt.isHom()) {
                        curHomIndex = 1;
                        curOtherIndex = 2;
                        curHomAllele = parent1CurGt.getAllele(0);
                        curOtherAlleles = new TreeSet<Allele>(parent2CurGt.getAlleles());
                    }
                    else if (parent2CurGt.isHom()) {
                        curHomIndex = 2;
                        curOtherIndex = 1;
                        curHomAllele = parent2CurGt.getAllele(0);
                        curOtherAlleles = new TreeSet<Allele>(parent1CurGt.getAlleles());
                    }
                    else
                        throw new ReviewedStingException("LOGICAL ERROR: at least one parent is hom!");

                    boolean phased = true;

                    Map<Allele, Integer> prevAlleleToParent = new TreeMap<Allele, Integer>();
                    for (Allele prevAllele : samplePrevGtInTrio.getAlleles()) {
                        if (prevAllele.equals(prevHomAllele))
                            prevAlleleToParent.put(prevAllele, prevHomIndex);
                        else if (prevOtherAlleles.contains(prevAllele))
                            prevAlleleToParent.put(prevAllele, prevOtherIndex);
                        else {
                            logger.warn("CANNOT phase, due to inconsistent inheritance of alleles!");
                            phased = false;
                            break;
                        }
                    }

                    Map<Integer, Allele> parentToCurAllele = new HashMap<Integer, Allele>();
                    for (Allele curAllele : sampleCurGtInTrio.getAlleles()) {
                        if (curAllele.equals(curHomAllele))
                            parentToCurAllele.put(curHomIndex, curAllele);
                        else if (curOtherAlleles.contains(curAllele))
                            parentToCurAllele.put(curOtherIndex, curAllele);
                        else {
                            logger.warn("CANNOT phase, due to inconsistent inheritance of alleles!");
                            phased = false;
                            break;
                        }
                    }

                    if (phased) {
                        List<Allele> phasedCurAlleles = new LinkedList<Allele>();
                        for (Allele prevAllele : prevPhasingGt.getAlleles()) {
                            Integer prevIndex = prevAlleleToParent.get(prevAllele);
                            if (prevIndex == null)
                                throw new ReviewedStingException("LOGICAL error: expecting to find prev allele in trio parents");
                            Allele curAllele = parentToCurAllele.get(prevIndex);
                            if (curAllele == null)
                                throw new ReviewedStingException("LOGICAL error: expecting to find cur allele in trio parents");
                            phasedCurAlleles.add(curAllele);
                        }

                        boolean useTrioPhase = true;
                        Genotype phasedGt = new Genotype(phasingSample, phasedCurAlleles, curPhasingGt.getNegLog10PError(), curPhasingGt.getFilters(), curPhasingGt.getAttributes(), phased);

                        if (curPhasingGt.isPhased()) {
                            stats.bothCanPhase++;
                            useTrioPhase = false;

                            if (!phasedCurAlleles.equals(curPhasingGt.getAlleles())) {
                                String contradictMessage = "Phase from " + PHASING_ROD_NAME + " track at " + curLoc + " contradicts the trio-based phasing.";
                                stats.contradictoryPhaseSites++;
                                addToOutput += "\tcontradictory";

                                if (phasingVc.hasAttribute(ReadBackedPhasingWalker.PHASING_INCONSISTENT_KEY)) {
                                    stats.contradictoryPhaseSitesWithPhaseInconsistency++;
                                    addToOutput += "\tphaseInconsistent";                                    
                                    useTrioPhase = true;
                                    contradictMessage += " Ignoring " + PHASING_ROD_NAME + " phase due to phase-inconsistency.";
                                }
                                else {
                                    contradictMessage += " Maintaining phase from " + PHASING_ROD_NAME + ".";
                                }
                                logger.warn(contradictMessage);
                            }
                        }

                        if (useTrioPhase) { // trio phasing adds PREVIOUSLY UNKNOWN phase information:
                            Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
                            genotypes.put(phasingSample, phasedGt);

                            phasingVc = VariantContext.modifyGenotypes(phasingVc, genotypes);
                            result.phasedVc = phasingVc;
                        }
                    }
                }
            }
            out.println(prevLoc + "\t" + curLoc + "\t" + trioPhaseStatus + "\t" + curPhasingGt.isPhased() + addToOutput);
            
            if (diffTrioAndPhasingTracks != null && prevTrioStatus != TrioStatus.MISSING && currentTrioStatus != TrioStatus.MISSING && sampleCurGtInTrio.isPhased() && curPhasingGt.isPhased()) {
                AllelePair prevTrioAll = new AllelePair(prevTrioVc.getGenotype(phasingSample));
                AllelePair curTrioAll = new AllelePair(sampleCurGtInTrio);
                
                AllelePair prevPhasingAll = new AllelePair(prevPhasingGt);
                AllelePair curPhasingAll = new AllelePair(curPhasingGt);
                
                boolean topsMatch = (GenotypePhasingEvaluator.topMatchesTop(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.topMatchesTop(curTrioAll, curPhasingAll));
                boolean bottomsMatch = (GenotypePhasingEvaluator.bottomMatchesBottom(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.bottomMatchesBottom(curTrioAll, curPhasingAll));

                boolean topMatchesBottom = (GenotypePhasingEvaluator.topMatchesBottom(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.topMatchesBottom(curTrioAll, curPhasingAll));
                boolean bottomMatchesTop = (GenotypePhasingEvaluator.bottomMatchesTop(prevTrioAll, prevPhasingAll) && GenotypePhasingEvaluator.bottomMatchesTop(curTrioAll, curPhasingAll));

                boolean phasesAgree = ((topsMatch && bottomsMatch) || (topMatchesBottom && bottomMatchesTop));

                diffTrioAndPhasingTracks.println(prevLoc + "\t" + curLoc + "\t" + trioPhaseStatus + "\t" + phasesAgree);
                diffTrioAndPhasingCounts.addComparison(trioPhaseStatus, phasesAgree);
            }
        }

        prevLoc = curLoc;
        prevTrioVc = curTrioVc;
        prevTrioStatus = currentTrioStatus;
        prevPhasingGt = curPhasingGt;

        return result;
    }

    public CompareToTrioPhasingStats reduce(CompareResult addIn, CompareToTrioPhasingStats runningCount) {
        if (addIn == null)
            addIn = new CompareResult();

        if (writer != null && addIn.phasedVc != null)
            WriteVCF.writeVCF(addIn.phasedVc, writer, logger);

        return runningCount.addIn(addIn.stats);
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(CompareToTrioPhasingStats result) {
        System.out.println("Compared " + result.comparedSites + " sites.");
        System.out.println("Trio can phase " + result.trioPhaseableSites + " sites.");
        System.out.println("Trio and " + PHASING_ROD_NAME + " track can both phase " + result.bothCanPhase + " sites.");
        System.out.println("Contradiction between phase inferred from " + TRIO_ROD_NAME + " and phase present in " + PHASING_ROD_NAME + " tracks at " + result.contradictoryPhaseSites + " sites.");
        System.out.println("Of those, " + PHASING_ROD_NAME + " track is phase-inconsistent at " + result.contradictoryPhaseSitesWithPhaseInconsistency + " sites.");

        if (diffTrioAndPhasingCounts != null) {
            System.out.println("");
            diffTrioAndPhasingCounts.printSummary(System.out);
        }
    }
}

class CompareToTrioPhasingStats {
    public int comparedSites;
    public int trioPhaseableSites;
    public int contradictoryPhaseSites;
    public int contradictoryPhaseSitesWithPhaseInconsistency;
    public int bothCanPhase;

    public CompareToTrioPhasingStats() {
        this.comparedSites = 0;
        this.trioPhaseableSites = 0;
        this.contradictoryPhaseSites = 0;
        this.contradictoryPhaseSitesWithPhaseInconsistency = 0;
        this.bothCanPhase = 0;
    }

    public CompareToTrioPhasingStats addIn(CompareToTrioPhasingStats other) {
        this.comparedSites += other.comparedSites;
        this.trioPhaseableSites += other.trioPhaseableSites;
        this.contradictoryPhaseSites += other.contradictoryPhaseSites;
        this.contradictoryPhaseSitesWithPhaseInconsistency += other.contradictoryPhaseSitesWithPhaseInconsistency;
        this.bothCanPhase += other.bothCanPhase;

        return this;
    }
}

class CompareResult {
    public VariantContext phasedVc;
    public CompareToTrioPhasingStats stats;

    public CompareResult() {
        this.phasedVc = null;
        this.stats = new CompareToTrioPhasingStats();
    }

    public CompareResult(VariantContext phasedVc, CompareToTrioPhasingStats stats) {
        this.phasedVc = phasedVc;
        this.stats = stats;
    }
}

class CompareTrioAndPhasingTracks {
    private Map<String, AgreeDisagreeCounts> trioStatusToAgreement;

    public CompareTrioAndPhasingTracks() {
        this.trioStatusToAgreement = new HashMap<String, AgreeDisagreeCounts>();
    }

    public void addComparison(String trioStatus, boolean agree) {
        AgreeDisagreeCounts counts = trioStatusToAgreement.get(trioStatus);
        if (counts == null) {
            counts = new AgreeDisagreeCounts();
            trioStatusToAgreement.put(trioStatus, counts);
        }

        if (agree)
            counts.incrementAgree();
        else
            counts.incrementDisagree();
    }

    public void printSummary(PrintStream out) {
        out.println("--------------------------------------------");
        out.println("Summary of trio vs. phasing tracks' phasing:");
        out.println("--------------------------------------------");        

        int globalAgree = 0;
        int globalDisagree = 0;
        for (AgreeDisagreeCounts counts : trioStatusToAgreement.values()) {
            globalAgree += counts.agree;
            globalDisagree += counts.disagree;
        }
        int globalTotal = globalAgree + globalDisagree;

        out.println("Concordant phase:\t" + percentString(globalAgree, globalTotal));
        out.println("Discordant phase:\t" + percentString(globalDisagree, globalTotal));

        for (Map.Entry<String, AgreeDisagreeCounts> statusCounts : trioStatusToAgreement.entrySet()) {
            String status = statusCounts.getKey();
            AgreeDisagreeCounts counts = statusCounts.getValue();

            out.println("");
            out.println("'" + status + "'" + " Concordant phase:\t" + percentString(counts.agree, counts.total()));
            out.println("'" + status + "'" + " Discordant phase:\t" + percentString(counts.disagree, counts.total()));
        }
        out.println("--------------------------------------------");
        out.println("");
    }

    private static String percentString(int numerator, int denominator) {
        int NUM_DECIMAL_PLACES = 1;
        String percent = new Formatter().format("%." + NUM_DECIMAL_PLACES + "f", MathUtils.percentage(numerator, denominator)).toString();

        StringBuilder sb = new StringBuilder();
        sb.append(numerator).append(" (").append(percent).append("%)");

        return sb.toString();
    }
}

class AgreeDisagreeCounts {
    protected int agree;
    protected int disagree;

    public AgreeDisagreeCounts() {
        this.agree = 0;
        this.disagree = 0;
    }

    public void incrementAgree() {
        agree++;
    }

    public void incrementDisagree() {
        disagree++;
    }

    public int total() {
        return agree + disagree;
    }
}