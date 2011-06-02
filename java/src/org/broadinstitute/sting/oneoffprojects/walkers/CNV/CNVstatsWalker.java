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

package org.broadinstitute.sting.oneoffprojects.walkers.CNV;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.PrintStream;
import java.util.*;


/**
 * Walks along all variant ROD loci, and tabulates the statistics of the CNVs detected.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = @RMD(name = "variant", type = ReferenceOrderedDatum.class))
@By(DataSource.REFERENCE_ORDERED_DATA)

public class CNVstatsWalker extends RodWalker<CNVstatistics, CNVstatistics> {

    @Output(doc = "File to which copy number counts should be written", required = true)
    protected PrintStream out;

    @Argument(fullName = "alleleCountsCopyCounts", shortName = "AC_CC", doc = "File to which discovered allele copy and copy number counts should be written", required = false)
    private PrintStream alleleCountsCopyCounts = null;

    @Argument(fullName = "allowFilteredGenotypes", shortName = "allowFilteredGenotypes", doc = "Include filtered genotypes in the output/calculations", required = false)
    private boolean ALLOW_FILTERED_GENOTYPES = false;

    private LinkedList<String> rodNames = null;

    public static String CNV_TAG = "<CNV>";
    public static String CN_FIELD = "CN";

    public static String SVLEN_FIELD = "SVLEN";
    public static String AC_FIELD = "AC";

    public static int DIPLOID = 2;

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add("variant");
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public CNVstatistics reduceInit() {
        return new CNVstatistics();
    }

    /**
     * For each site, calculate the CNV stats.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return dummy Integer
     */
    public CNVstatistics map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        logger.debug("REF:" + ref.getLocus());
        CNVstatistics stats = new CNVstatistics();

        boolean requireStartHere = true; // only see each VariantContext once
        boolean takeFirstOnly = false; // take as many entries as the VCF file has
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            if (vc.isMixed() && vc.isBiallelic()) {
                Allele altAll = vc.getAlternateAllele(0);
                if (altAll.isSymbolic() && altAll.getDisplayString().equals(CNV_TAG)) {
                    logger.debug("Found CNV at locus...");
                    stats.cnvDeclaredLoci++;

                    CopyNumberCounts cnc = new CopyNumberCounts();

                    boolean hasDiploidGt = false;
                    boolean hasNonDiploidGt = false;
                    for (Map.Entry<String, Genotype> gtEntry : vc.getGenotypes().entrySet()) {
                        Genotype gt = gtEntry.getValue();
                        Integer copyNum = gt.getAttributeAsIntegerNoException(CN_FIELD);
                        if (copyNum != null && (ALLOW_FILTERED_GENOTYPES || gt.isNotFiltered())) {
                            cnc.incrementCopyNumber(copyNum);

                            if (copyNum == DIPLOID)
                                hasDiploidGt = true;
                            else
                                hasNonDiploidGt = true;
                        }
                    }

                    if (hasDiploidGt && hasNonDiploidGt) {
                        stats.diploidAndNonDiploidLoci++;
                    }
                    else {
                        if (hasDiploidGt)
                            stats.diploidOnlyLoci++;
                        if (hasNonDiploidGt)
                            stats.nonDiploidOnlyLoci++;
                    }

                    int cnvEnd = vc.getEnd();
                    Integer cnvLength = vc.getAttributeAsIntegerNoException(SVLEN_FIELD);
                    if (cnvLength != null)
                        cnvEnd = vc.getStart() + cnvLength - 1;
                    GenomeLoc vcLoc = getToolkit().getGenomeLocParser().createGenomeLoc(vc.getChr(), vc.getStart(), cnvEnd, true);
                    out.print(vcLoc);

                    for (Map.Entry<Integer, Integer> copyNumEntry : cnc.entrySet()) {
                        out.print("\t" + copyNumEntry.getKey() + ":" + copyNumEntry.getValue());
                    }
                    out.println();

                    if (alleleCountsCopyCounts != null) {
                        Integer ac = vc.getAttributeAsIntegerNoException(AC_FIELD);
                        CopyNumberCounts.DeletionDuplicationCounts counts = cnc.deletionDuplicationCounts();
                        double cnvCount = counts.deletionCount + counts.duplicationCount;

                        alleleCountsCopyCounts.println(vcLoc + "\t" + ac + "\t" + counts.deletionCount + "\t" + counts.duplicationCount + "\t" + cnvCount);
                    }
                }
            }
        }

        return stats;
    }

    public CNVstatistics reduce(CNVstatistics result, CNVstatistics total) {
        if (result == null)
            return total;

        return total.addIn(result);
    }

    /**
     * @param result statistics of CNV sites
     */
    public void onTraversalDone(CNVstatistics result) {
        System.out.println();
        System.out.println("--------------------------------------");
        System.out.println("CNV summary:");
        System.out.println("--------------------------------------");

        System.out.println("cnvDeclaredLoci: " + result.cnvDeclaredLoci);

        System.out.println();
        System.out.println("noGenotypesLoci: " + result.noGenotypesLoci());

        System.out.println();
        System.out.println("nonDiploidOnlyLoci: " + result.nonDiploidOnlyLoci);

        System.out.println();
        System.out.println("lociWithDiploid: " + result.lociWithDiploid());
        System.out.println("diploidOnlyLoci: " + result.diploidOnlyLoci);
        System.out.println("diploidAndNonDiploidLoci: " + result.diploidAndNonDiploidLoci);
        String onlyDiploidRateStr = percentageString(result.diploidOnlyLoci, result.lociWithDiploid());
        System.out.println("onlyDiploidRate = " + onlyDiploidRateStr + "%");

        System.out.println();
        int noDiploidGenotypes = result.noGenotypesLoci() + result.nonDiploidOnlyLoci;
        System.out.println("loci with no diploid genotypes: " + noDiploidGenotypes);
        String noDiploidGtRateStr = percentageString(noDiploidGenotypes, result.cnvDeclaredLoci);
        System.out.println("noDiploidGtRate = " + noDiploidGtRateStr + "%");
    }

    private static String percentageString(int numerator, int denominator) {
        int NUM_DECIMAL_PLACES = 2;

        return new Formatter().format("%." + NUM_DECIMAL_PLACES + "f", MathUtils.percentage(numerator, denominator)).toString();
    }
}


class CNVstatistics {
    protected int cnvDeclaredLoci, diploidOnlyLoci, nonDiploidOnlyLoci, diploidAndNonDiploidLoci;

    public CNVstatistics() {
        this.cnvDeclaredLoci = 0;
        this.diploidOnlyLoci = 0;
        this.nonDiploidOnlyLoci = 0;
        this.diploidAndNonDiploidLoci = 0;
    }

    public CNVstatistics addIn(CNVstatistics other) {
        this.cnvDeclaredLoci += other.cnvDeclaredLoci;
        this.diploidOnlyLoci += other.diploidOnlyLoci;
        this.nonDiploidOnlyLoci += other.nonDiploidOnlyLoci;
        this.diploidAndNonDiploidLoci += other.diploidAndNonDiploidLoci;

        return this;
    }

    public int noGenotypesLoci() {
        return cnvDeclaredLoci - (diploidOnlyLoci + nonDiploidOnlyLoci + diploidAndNonDiploidLoci);
    }

    public int lociWithDiploid() {
        return diploidOnlyLoci + diploidAndNonDiploidLoci;
    }
}

class CopyNumberCounts {
    private Map<Integer, Integer> copyNumToCountsMap;

    public CopyNumberCounts() {
        this.copyNumToCountsMap = new TreeMap<Integer, Integer>();
    }

    public void incrementCopyNumber(int copyNum) {
        Integer count = copyNumToCountsMap.get(copyNum);
        if (count == null)
            count = 0;

        copyNumToCountsMap.put(copyNum, count + 1);
    }
    
    public Set<Map.Entry<Integer, Integer>> entrySet() {
        return copyNumToCountsMap.entrySet();
    }

    class DeletionDuplicationCounts {
        public double deletionCount;
        public double duplicationCount;

        public DeletionDuplicationCounts() {
            this.deletionCount = 0;
            this.duplicationCount = 0;
        }
    }

    public DeletionDuplicationCounts deletionDuplicationCounts() {
        int total = 0;
        DeletionDuplicationCounts counts = new DeletionDuplicationCounts();

        for (Map.Entry<Integer, Integer> copyNumEntry : this.entrySet()) {
            int copyNum = copyNumEntry.getKey();
            int count = copyNumEntry.getValue();

            if (copyNum < CNVstatsWalker.DIPLOID) {
                counts.deletionCount += count;
            }
            else if (copyNum > CNVstatsWalker.DIPLOID) {
                counts.duplicationCount += count;
            }

            total += count;
        }

        //counts.deletionCount /= total;
        //counts.duplicationCount /= total;

        return counts;
    }
}