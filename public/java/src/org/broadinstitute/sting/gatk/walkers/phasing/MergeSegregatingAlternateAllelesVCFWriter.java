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

package org.broadinstitute.sting.gatk.walkers.phasing;

import net.sf.picard.reference.ReferenceSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

// Streams in VariantContext objects and streams out VariantContexts produced by merging phased segregating polymorphisms into MNP VariantContexts

public class MergeSegregatingAlternateAllelesVCFWriter implements VCFWriter {
    private VCFWriter innerWriter;

    private GenomeLocParser genomeLocParser;

    private ReferenceSequenceFile referenceFileForMNPmerging;

    private VariantContextMergeRule vcMergeRule;
    private VariantContextUtils.AlleleMergeRule alleleMergeRule;

    private String useSingleSample = null;

    private boolean emitOnlyMergedRecords;

    private VCFRecord vcfrWaitingToMerge;
    private List<VCFRecord> filteredVcfrList;

    private int numRecordsAttemptToMerge;
    private int numRecordsSatisfyingMergeRule;
    private int numMergedRecords;
    private AltAlleleStatsForSamples altAlleleStats = null;

    private Logger logger;

    // Should we call innerWriter.close() in close()
    private boolean takeOwnershipOfInner;

    public MergeSegregatingAlternateAllelesVCFWriter(VCFWriter innerWriter, GenomeLocParser genomeLocParser, File referenceFile, VariantContextMergeRule vcMergeRule, VariantContextUtils.AlleleMergeRule alleleMergeRule, String singleSample, boolean emitOnlyMergedRecords, Logger logger, boolean takeOwnershipOfInner, boolean trackAltAlleleStats) {
        this.innerWriter = innerWriter;
        this.genomeLocParser = genomeLocParser;
        try {
            this.referenceFileForMNPmerging = new CachingIndexedFastaSequenceFile(referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }        
        this.vcMergeRule = vcMergeRule;
        this.alleleMergeRule = alleleMergeRule;
        this.useSingleSample = singleSample;
        this.emitOnlyMergedRecords = emitOnlyMergedRecords;

        this.vcfrWaitingToMerge = null;
        this.filteredVcfrList = new LinkedList<VCFRecord>();
        this.numRecordsSatisfyingMergeRule = 0;
        this.numMergedRecords = 0;

        if (trackAltAlleleStats)
            this.altAlleleStats = new AltAlleleStatsForSamples();

        this.logger = logger;
        this.takeOwnershipOfInner = takeOwnershipOfInner;
    }

    public MergeSegregatingAlternateAllelesVCFWriter(VCFWriter innerWriter, GenomeLocParser genomeLocParser, File referenceFile, int maxGenomicDistanceForMNP, Logger logger, boolean takeOwnershipOfInner) {
        this(innerWriter, genomeLocParser, referenceFile, new DistanceMergeRule(maxGenomicDistanceForMNP, genomeLocParser), new SegregatingMNPmergeAllelesRule(), null, false, logger, takeOwnershipOfInner, true); // by default: consider all samples, emit all records, keep track of alt allele statistics
    }

    public void writeHeader(VCFHeader header) {
        if (useSingleSample != null) { // only want to output context for one sample
            Set<String> singSampSet = new TreeSet<String>();
            singSampSet.add(useSingleSample);
            header = new VCFHeader(header.getMetaData(), singSampSet);
        }

        innerWriter.writeHeader(header);
    }

    public void close() {
        stopWaitingToMerge();

        if (takeOwnershipOfInner)
            innerWriter.close();
    }

    public void add(VariantContext vc) {
        if (useSingleSample != null) { // only want to output context for one sample
            Genotype sampGt = vc.getGenotype(useSingleSample);
            if (sampGt != null) // TODO: subContextFromGenotypes() does not handle any INFO fields [AB, HaplotypeScore, MQ, etc.].  Note that even SelectVariants.subsetRecord() only handles AC,AN,AF, and DP!
                vc = vc.subContextFromGenotypes(sampGt);
            else // asked for a sample that this vc does not contain, so ignore this vc:
                return;
        }

        logger.debug("Next VC input = " + VariantContextUtils.getLocation(genomeLocParser, vc));
        boolean curVcIsNotFiltered = vc.isNotFiltered();

        if (vcfrWaitingToMerge == null) {
            logger.debug("NOT Waiting to merge...");

            if (!filteredVcfrList.isEmpty())
                throw new ReviewedStingException("filteredVcfrList should be empty if not waiting to merge a vc!");

            if (curVcIsNotFiltered) { // still need to wait before can release vc
                logger.debug("Waiting for new variant " + VariantContextUtils.getLocation(genomeLocParser, vc));
                vcfrWaitingToMerge = new VCFRecord(vc, false);
            }
            else if (!emitOnlyMergedRecords) { // filtered records are never merged
                logger.debug("DIRECTLY output " + VariantContextUtils.getLocation(genomeLocParser, vc));
                innerWriter.add(vc);
            }
        }
        else { // waiting to merge vcfrWaitingToMerge
            logger.debug("Waiting to merge " + VariantContextUtils.getLocation(genomeLocParser, vcfrWaitingToMerge.vc));

            if (!curVcIsNotFiltered) {
                if (!emitOnlyMergedRecords) { // filtered records are never merged
                    logger.debug("Caching unprocessed output " + VariantContextUtils.getLocation(genomeLocParser, vc));
                    filteredVcfrList.add(new VCFRecord(vc, false));
                }
            }
            else { // waiting to merge vcfrWaitingToMerge, and curVcIsNotFiltered. So, attempt to merge them:
                numRecordsAttemptToMerge++;
                boolean shouldAttemptToMerge = vcMergeRule.shouldAttemptToMerge(vcfrWaitingToMerge.vc, vc);
                logger.debug("shouldAttemptToMerge? = " + shouldAttemptToMerge);

                /*
                TODO: -- CONSIDER THE FOLLOWING EXAMPLE: WHAT DO WE WANT HERE??? --
                If the following 3 genotypes originally exist for a sample [at sites 1, 2, and 3]:
                1/1
                0|1
                0|1

                Then, after merging the first two, we have [at sites 1 and 3]:
                1/2
                0|1

                Then, not having merged would consider sites 2 and 3 as a MNP (since it's a diploid het site with haplotypes: REF-REF and ALT-ALT)
                But, since we merged sites 1 and 2, we get that sites 1-2 and 3 are counted as two haplotypes of: ALT-REF and ALT-ALT
                 */
                if (altAlleleStats != null)
                    altAlleleStats.updateSampleStats(vcfrWaitingToMerge.vc, vc, shouldAttemptToMerge);

                boolean mergedRecords = false;
                if (shouldAttemptToMerge) {
                    numRecordsSatisfyingMergeRule++;
                    VariantContext mergedVc = VariantContextUtils.mergeIntoMNP(genomeLocParser, vcfrWaitingToMerge.vc, vc, referenceFileForMNPmerging, alleleMergeRule);

                    if (mergedVc != null) {
                        mergedRecords = true;

                        Map<String, Object> addedAttribs = vcMergeRule.addToMergedAttributes(vcfrWaitingToMerge.vc, vc);
                        addedAttribs.putAll(mergedVc.getAttributes());
                        mergedVc = VariantContext.modifyAttributes(mergedVc, addedAttribs);

                        vcfrWaitingToMerge = new VCFRecord(mergedVc, true);
                        numMergedRecords++;
                    }
                }

                if (!mergedRecords) {
                    stopWaitingToMerge();
                    vcfrWaitingToMerge = new VCFRecord(vc, false);
                }
                logger.debug("Merged? = " + mergedRecords);
            }
        }
    }

    private void stopWaitingToMerge() {
        if (vcfrWaitingToMerge == null) {
            if (!filteredVcfrList.isEmpty())
                throw new ReviewedStingException("filteredVcfrList should be empty if not waiting to merge a vc!");
            return;
        }

        if (!emitOnlyMergedRecords || vcfrWaitingToMerge.resultedFromMerge)
            innerWriter.add(vcfrWaitingToMerge.vc);
        vcfrWaitingToMerge = null;

        for (VCFRecord vcfr : filteredVcfrList)
            innerWriter.add(vcfr.vc);
        filteredVcfrList.clear();
    }

    public int getNumRecordsAttemptToMerge() {
        return numRecordsAttemptToMerge;
    }

    public int getNumRecordsSatisfyingMergeRule() {
        return numRecordsSatisfyingMergeRule;
    }

    public int getNumMergedRecords() {
        return numMergedRecords;
    }

    public VariantContextMergeRule getVcMergeRule() {
        return vcMergeRule;
    }

    public VariantContextUtils.AlleleMergeRule getAlleleMergeRule() {
        return alleleMergeRule;
    }

    /**
     * Gets a string representation of this object.
     *
     * @return
     */
    @Override
    public String toString() {
        return getClass().getName();
    }

    public String getAltAlleleStats() {
        if (altAlleleStats == null)
            return "";

        return "\n" + altAlleleStats.toString();
    }

    private static class VCFRecord {
        public VariantContext vc;
        public boolean resultedFromMerge;

        public VCFRecord(VariantContext vc, boolean resultedFromMerge) {
            this.vc = vc;
            this.resultedFromMerge = resultedFromMerge;
        }
    }

    private class AltAlleleStats {
        public int numSuccessiveGenotypes;
        public int numSuccessiveGenotypesAttemptedToBeMerged;

        public int oneSampleMissing;
        public int atLeastOneSampleNotCalledOrFiltered;
        public int segregationUnknown;
        public int eitherNotVariant;

        public int bothInPairHaveVariant;

        public int ref_ref_pair;
        public int ref_alt_pair;
        public int alt_ref_pair;
        public int alt_alt_pair;

        public int MNPsites;
        public int CHetSites;

        public AltAlleleStats() {
            this.numSuccessiveGenotypes = 0;
            this.numSuccessiveGenotypesAttemptedToBeMerged = 0;

            this.oneSampleMissing = 0;
            this.atLeastOneSampleNotCalledOrFiltered = 0;
            this.segregationUnknown = 0;
            this.eitherNotVariant = 0;

            this.bothInPairHaveVariant = 0;

            this.ref_ref_pair = 0;
            this.ref_alt_pair = 0;
            this.alt_ref_pair = 0;
            this.alt_alt_pair = 0;

            this.MNPsites = 0;
            this.CHetSites = 0;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append("Sample missing:\t" + oneSampleMissing + "\n");
            sb.append("Not called or filtered:\t" + atLeastOneSampleNotCalledOrFiltered + "\n");

            sb.append("* Number of successive pairs of genotypes:\t" + numSuccessiveGenotypes + "\n");
            sb.append("Number of successive pairs of genotypes with " + vcMergeRule + ":\t" + numSuccessiveGenotypesAttemptedToBeMerged + "\n");

            sb.append("Unknown segregation, " + vcMergeRule + ":\t" + segregationUnknown + "\n");
            sb.append("Not variant at least one of pair, segregation known, " + vcMergeRule + ":\t" + eitherNotVariant + "\n");
            sb.append("* Variant at both, segregation known, " + vcMergeRule + ":\t" + percentageString(bothInPairHaveVariant, numSuccessiveGenotypes) + "\n");

            sb.append("[Total haplotypes at pairs:\t" + (ref_ref_pair + ref_alt_pair + alt_ref_pair + alt_alt_pair) + "\n");
            sb.append("REF-REF:\t" + ref_ref_pair + "\n");
            sb.append("REF-ALT:\t" + ref_alt_pair + "\n");
            sb.append("ALT-REF:\t" + alt_ref_pair + "\n");
            sb.append("ALT-ALT:\t" + alt_alt_pair + "]\n");

            int hetAfterHetSites = MNPsites + CHetSites;
            sb.append("* Het-Het sites (with REF allele present at each):\t" + percentageString(hetAfterHetSites, bothInPairHaveVariant) + "\n");
            sb.append("* MNPs:\t" + percentageString(MNPsites, hetAfterHetSites) + "\n");
            sb.append("Compound Hets:\t" + CHetSites + "\n");

            return sb.toString();
        }

        private String percentageString(int count, int baseCount) {
            int NUM_DECIMAL_PLACES = 1;
            String percent = new Formatter().format("%." + NUM_DECIMAL_PLACES + "f", MathUtils.percentage(count, baseCount)).toString();
            return count + " (" + percent + "%)";
        }
    }

    private class AltAlleleStatsForSamples {
        private Map<String, AltAlleleStats> sampleStats;

        public AltAlleleStatsForSamples() {
            this.sampleStats = new HashMap<String, AltAlleleStats>();
        }

        public void updateSampleStats(VariantContext vc1, VariantContext vc2, boolean shouldAttemptToMerge) {
            if (vc1.isFiltered() || vc2.isFiltered())
                return;

            Set<String> allSamples = new TreeSet<String>(vc1.getSampleNames());
            allSamples.addAll(vc2.getSampleNames());

            for (String samp : allSamples) {
                AltAlleleStats aas = sampleStats.get(samp);
                if (aas == null) {
                    aas = new AltAlleleStats();
                    sampleStats.put(samp, aas);
                }

                Genotype gt1 = vc1.getGenotype(samp);
                Genotype gt2 = vc2.getGenotype(samp);
                if (gt1 == null || gt2 == null) {
                    aas.oneSampleMissing++;
                }
                else if (gt1.isNoCall() || gt1.isFiltered() || gt2.isNoCall() || gt2.isFiltered()) {
                    aas.atLeastOneSampleNotCalledOrFiltered++;
                }
                else {
                    aas.numSuccessiveGenotypes++;

                    if (shouldAttemptToMerge) {
                        aas.numSuccessiveGenotypesAttemptedToBeMerged++;

                        if (!VariantContextUtils.alleleSegregationIsKnown(gt1, gt2)) {
                            aas.segregationUnknown++;
                            logger.debug("Unknown segregation of alleles [not phased] for " + samp + " at " + VariantContextUtils.getLocation(genomeLocParser, vc1) + ", " + VariantContextUtils.getLocation(genomeLocParser, vc2));
                        }
                        else if (gt1.isHomRef() || gt2.isHomRef()) {
                            logger.debug("gt1.isHomRef() || gt2.isHomRef() for " + samp + " at " + VariantContextUtils.getLocation(genomeLocParser, vc1) + ", " + VariantContextUtils.getLocation(genomeLocParser, vc2));
                            aas.eitherNotVariant++;
                        }
                        else { // BOTH gt1 and gt2 have at least one variant allele (so either hets, or homozygous variant):
                            aas.bothInPairHaveVariant++;

                            List<Allele> site1Alleles = gt1.getAlleles();
                            List<Allele> site2Alleles = gt2.getAlleles();

                            Iterator<Allele> all2It = site2Alleles.iterator();
                            for (Allele all1 : site1Alleles) {
                                Allele all2 = all2It.next(); // this is OK, since alleleSegregationIsKnown(gt1, gt2)

                                if (all1.isReference()) {
                                    if (all2.isReference())
                                        aas.ref_ref_pair++;
                                    else
                                        aas.ref_alt_pair++;
                                }
                                else { // all1.isNonReference()
                                    if (all2.isReference())
                                        aas.alt_ref_pair++;
                                    else
                                        aas.alt_alt_pair++;
                                }
                            }

                            // Check MNPs vs. CHets:
                            if (containsRefAllele(site1Alleles) && containsRefAllele(site2Alleles)) {
                                logger.debug("HET-HET for " + samp + " at " + VariantContextUtils.getLocation(genomeLocParser, vc1) + ", " + VariantContextUtils.getLocation(genomeLocParser, vc2));
                                if (logger.isDebugEnabled() && !(gt1.isHet() && gt2.isHet()))
                                    throw new ReviewedStingException("Since !gt1.isHomRef() && !gt2.isHomRef(), yet both have ref alleles, they BOTH must be hets!");

                                // There's the potential to only have REF-ALT, ALT-REF (CHet), or possibly ALT-ALT together (MNP)
                                boolean hasMNP = false;

                                all2It = site2Alleles.iterator();
                                for (Allele all1 : site1Alleles) {
                                    Allele all2 = all2It.next(); // this is OK, since alleleSegregationIsKnown(gt1, gt2)

                                    if (all1.isNonReference() && all2.isNonReference()) {
                                        hasMNP = true; // has at least one haplotype of ALT-ALT that segregates!
                                        break;
                                    }
                                }

                                if (hasMNP)
                                    aas.MNPsites++;
                                else
                                    aas.CHetSites++;
                            }
                        }
                    }
                }
            }
        }

        private boolean containsRefAllele(List<Allele> siteAlleles) {
            for (Allele all : siteAlleles) {
                if (all.isReference())
                    return true;
            }

            return false;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("-------------------------------------------------------------------------\n");
            sb.append("Per-sample alternate allele statistics [" + vcMergeRule + "]\n");
            sb.append("-------------------------------------------------------------------------");

            for (Map.Entry<String, AltAlleleStats> sampAltAllStatsEntry : sampleStats.entrySet()) {
                String samp = sampAltAllStatsEntry.getKey();
                AltAlleleStats stats = sampAltAllStatsEntry.getValue();
                sb.append("\n* Sample:\t" + samp + "\n" + stats);
            }

            return sb.toString();
        }
    }
}



/*
 External classes:
 */

abstract class VariantContextMergeRule {
    abstract public boolean shouldAttemptToMerge(VariantContext vc1, VariantContext vc2);

    public Map<String, Object> addToMergedAttributes(VariantContext vc1, VariantContext vc2) {
        return new HashMap<String, Object>();
    }
}

class DistanceMergeRule extends VariantContextMergeRule {
    private int maxGenomicDistanceForMNP;
    private GenomeLocParser genomeLocParser;

    public DistanceMergeRule(int maxGenomicDistanceForMNP, GenomeLocParser genomeLocParser) {
        this.maxGenomicDistanceForMNP = maxGenomicDistanceForMNP;
        this.genomeLocParser = genomeLocParser;
    }

    public boolean shouldAttemptToMerge(VariantContext vc1, VariantContext vc2) {
        return minDistance(vc1, vc2) <= maxGenomicDistanceForMNP;
    }

    public String toString() {
        return "Merge distance <= " + maxGenomicDistanceForMNP;
    }

    public int minDistance(VariantContext vc1, VariantContext vc2) {
        return VariantContextUtils.getLocation(genomeLocParser, vc1).minDistance(VariantContextUtils.getLocation(genomeLocParser, vc2));
    }
}


class ExistsDoubleAltAlleleMergeRule extends VariantContextUtils.AlleleMergeRule {
    public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2) {
        return VariantContextUtils.someSampleHasDoubleNonReferenceAllele(vc1, vc2);
    }

    public String toString() {
        return super.toString() + ", some sample has a MNP of ALT alleles";
    }
}

class SegregatingMNPmergeAllelesRule extends ExistsDoubleAltAlleleMergeRule {
    public SegregatingMNPmergeAllelesRule() {
        super();
    }

    public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2) {
        // Must be interesting AND consistent:
        return super.allelesShouldBeMerged(vc1, vc2) && VariantContextUtils.doubleAllelesSegregatePerfectlyAmongSamples(vc1, vc2);
    }

    public String toString() {
        return super.toString() + ", all alleles segregate consistently";
    }
}