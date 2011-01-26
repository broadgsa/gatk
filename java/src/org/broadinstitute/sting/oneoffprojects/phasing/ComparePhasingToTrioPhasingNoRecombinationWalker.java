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

package org.broadinstitute.sting.oneoffprojects.phasing;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeSet;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = ComparePhasingToTrioPhasingNoRecombinationWalker.TRIO_ROD_NAME, type = ReferenceOrderedDatum.class), @RMD(name = ComparePhasingToTrioPhasingNoRecombinationWalker.PHASING_ROD_NAME, type = ReferenceOrderedDatum.class)})

@ReadFilters({ZeroMappingQualityReadFilter.class})
// Filter out all reads with zero mapping quality

public class ComparePhasingToTrioPhasingNoRecombinationWalker extends RodWalker<Integer, Integer> {
    public final static String TRIO_ROD_NAME = "trio";
    public final static String PHASING_ROD_NAME = "phasing";

    private final static int NUM_IN_TRIO = 3;

    @Output
    protected PrintStream out;

    private String phasingSample = null;

    private enum TrioStatus {
        PRESENT, MISSING, TRIPLE_HET
    }

    private GenomeLoc prevLoc = null;
    private TrioStatus prevTrioStatus = TrioStatus.MISSING;


    public void initialize() {
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLoc = ref.getLocus();
        VariantContext phasingVc = tracker.getVariantContext(ref, PHASING_ROD_NAME, curLoc);
        if (phasingVc == null || phasingVc.isFiltered())
            return null;

        Map<String, Genotype> phasingSampleToGt = phasingVc.getGenotypes();
        if (phasingSampleToGt.size() != 1)
            throw new UserException("Must provide EXACTLY one sample in " + PHASING_ROD_NAME + " track!");
        Map.Entry<String, Genotype> phasingSampGt = phasingSampleToGt.entrySet().iterator().next();
        String sample = phasingSampGt.getKey();
        if (phasingSample == null)
            phasingSample = sample;
        if (!sample.equals(phasingSample))
            throw new UserException("Must provide EXACTLY one sample!");
        Genotype phasingGt = phasingSampGt.getValue();
        if (!phasingGt.isHet())
            return null;

        VariantContext trioVc = tracker.getVariantContext(ref, TRIO_ROD_NAME, curLoc);
        boolean useTrioVc = (trioVc != null && !trioVc.isFiltered());

        Genotype sampleGtInTrio = null;
        if (useTrioVc) {
            sampleGtInTrio = trioVc.getGenotype(phasingSample);

            if (trioVc.getNSamples() > NUM_IN_TRIO || sampleGtInTrio == null)
                throw new UserException("Must provide trio data for sample: " + phasingSample);

            if (!new TreeSet<Allele>(phasingGt.getAlleles()).equals(new TreeSet<Allele>(sampleGtInTrio.getAlleles()))) {
                logger.warn("Locus " + curLoc + " breaks phase, since " + PHASING_ROD_NAME + " and " + TRIO_ROD_NAME + " tracks have different genotypes for " + phasingSample + "!");
                prevLoc = null;
                return null;
            }
        }

        // Now, we have a [trio-consistent] het genotype that may be phased or not [and we want to know if it could be phased based on trio information]:
        int processed = 1;

        TrioStatus currentTrioStatus = TrioStatus.MISSING;
        if (useTrioVc && trioVc.getNSamples() == NUM_IN_TRIO) {
            boolean allHet = true;
            for (int i = 0; i < NUM_IN_TRIO; i++) {
                if (!trioVc.getGenotype(i).isHet()) {
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
            }

            out.println(prevLoc + "\t" + curLoc + "\t" + trioPhaseStatus + "\t" + phasingGt.isPhased());
        }

        prevLoc = curLoc;
        prevTrioStatus = currentTrioStatus;

        return processed;
    }

    public Integer reduce(Integer addIn, Integer runningCount) {
        if (addIn == null)
            addIn = 0;

        return runningCount + addIn;
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(Integer result) {
        System.out.println("Processed " + result + " sites.");
    }
}