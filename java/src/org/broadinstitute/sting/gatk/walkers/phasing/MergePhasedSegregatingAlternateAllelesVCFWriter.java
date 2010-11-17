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

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.util.*;

// Streams in VariantContext objects and streams out VariantContexts produced by merging phased segregating polymorphisms into MNP VariantContexts

public class MergePhasedSegregatingAlternateAllelesVCFWriter implements VCFWriter {
    private VCFWriter innerWriter;

    private GenomeLocParser genomeLocParser;

    private ReferenceSequenceFile referenceFileForMNPmerging;
    private MergeRule mergeRule;

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

    public MergePhasedSegregatingAlternateAllelesVCFWriter(VCFWriter innerWriter, GenomeLocParser genomeLocParser, File referenceFile, MergeRule mergeRule, String singleSample, boolean emitOnlyMergedRecords, Logger logger, boolean takeOwnershipOfInner, boolean trackAltAlleleStats) {
        this.innerWriter = innerWriter;
        this.genomeLocParser = genomeLocParser;
        this.referenceFileForMNPmerging = new IndexedFastaSequenceFile(referenceFile);
        this.mergeRule = mergeRule;
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

    public MergePhasedSegregatingAlternateAllelesVCFWriter(VCFWriter innerWriter, GenomeLocParser genomeLocParser, File referenceFile, int maxGenomicDistanceForMNP, Logger logger) {
        this(innerWriter, genomeLocParser, referenceFile, new DistanceMergeRule(maxGenomicDistanceForMNP, genomeLocParser), null, false, logger, false, false); // by default: consider all samples, emit all records, don't own inner, don't keep track of alt allele statistics
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

    public void add(VariantContext vc, byte refBase) {
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
                vcfrWaitingToMerge = new VCFRecord(vc, refBase, false);
            }
            else if (!emitOnlyMergedRecords) { // filtered records are never merged
                logger.debug("DIRECTLY output " + VariantContextUtils.getLocation(genomeLocParser, vc));
                innerWriter.add(vc, refBase);
            }
        }
        else { // waiting to merge vcfrWaitingToMerge
            logger.debug("Waiting to merge " + VariantContextUtils.getLocation(genomeLocParser, vcfrWaitingToMerge.vc));

            if (!curVcIsNotFiltered) {
                if (!emitOnlyMergedRecords) { // filtered records are never merged
                    logger.debug("Caching unprocessed output " + VariantContextUtils.getLocation(genomeLocParser, vc));
                    filteredVcfrList.add(new VCFRecord(vc, refBase, false));
                }
            }
            else { // waiting to merge vcfrWaitingToMerge, and curVcIsNotFiltered. So, attempt to merge them:
                numRecordsAttemptToMerge++;
                boolean shouldMerge = mergeRule.shouldMerge(vcfrWaitingToMerge.vc, vc);

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
                    altAlleleStats.updateSampleStats(vcfrWaitingToMerge.vc, vc, shouldMerge);

                boolean mergedRecords = false;
                if (shouldMerge) {
                    numRecordsSatisfyingMergeRule++;
                    VariantContext mergedVc = VariantContextUtils.mergeIntoMNP(genomeLocParser, vcfrWaitingToMerge.vc, vc, referenceFileForMNPmerging);

                    if (mergedVc != null) {
                        mergedRecords = true;

                        Map<String, Object> updatedAttribs = RefSeqData.getMergedRefSeqAttributes(vcfrWaitingToMerge.vc, vc);
                        updatedAttribs.putAll(mergedVc.getAttributes());
                        mergedVc = VariantContext.modifyAttributes(mergedVc, updatedAttribs);

                        vcfrWaitingToMerge = new VCFRecord(mergedVc, vcfrWaitingToMerge.refBase, true);
                        numMergedRecords++;
                    }
                }

                if (!mergedRecords) {
                    stopWaitingToMerge();
                    vcfrWaitingToMerge = new VCFRecord(vc, refBase, false);
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
            innerWriter.add(vcfrWaitingToMerge.vc, vcfrWaitingToMerge.refBase);
        vcfrWaitingToMerge = null;

        for (VCFRecord vcfr : filteredVcfrList)
            innerWriter.add(vcfr.vc, vcfr.refBase);
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

    public MergeRule getMergeRule() {
        return mergeRule;
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
        public byte refBase;
        public boolean resultedFromMerge;

        public VCFRecord(VariantContext vc, byte refBase, boolean resultedFromMerge) {
            this.vc = vc;
            this.refBase = refBase;
            this.resultedFromMerge = resultedFromMerge;
        }
    }

    private class AltAlleleStats {
        public int numSuccessiveGenotypes;
        public int numSuccessiveGenotypesThatCouldBeMerged;

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
            this.numSuccessiveGenotypesThatCouldBeMerged = 0;

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
            sb.append("Number of successive pairs of genotypes with " + mergeRule + ":\t" + numSuccessiveGenotypesThatCouldBeMerged + "\n");

            sb.append("Unknown segregation, " + mergeRule + ":\t" + segregationUnknown + "\n");
            sb.append("Not variant at least one of pair, segregation known, " + mergeRule + ":\t" + eitherNotVariant + "\n");
            sb.append("* Variant at both, segregation known, " + mergeRule + ":\t" + percentageString(bothInPairHaveVariant, numSuccessiveGenotypes) + "\n");

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

        public void updateSampleStats(VariantContext vc1, VariantContext vc2, boolean shouldMerge) {
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

                    if (shouldMerge) {
                        aas.numSuccessiveGenotypesThatCouldBeMerged++;

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
            sb.append("Per-sample alternate allele statistics [" + mergeRule + "]\n");
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


/* Some methods for extracting and merging RefSeq-related data from annotated VCF INFO fields:
 */

class RefSeqData {
    private static String REFSEQ_PREFIX = "refseq.";

    private static String NUM_RECORDS_KEY = REFSEQ_PREFIX + "numMatchingRecords";
    private static String NAME_KEY = REFSEQ_PREFIX + "name";
    private static String NAME2_KEY = REFSEQ_PREFIX + "name2";

    private static String CODON_KEY = REFSEQ_PREFIX + "codonCoord";

    private static Map<String, String> getRefSeqEntriesToNames(VariantContext vc, boolean getName2) {
        Map<String, String> entriesToNames = new HashMap<String, String>();

        Integer numRecords = vc.getAttributeAsIntegerNoException(NUM_RECORDS_KEY);
        if (numRecords != null) {
            for (int i = 1; i <= numRecords; i++) {
                String key = NAME_KEY + "_" + i;
                String name = vc.getAttributeAsStringNoException(key);
                if (name != null)
                    entriesToNames.put(key, name);
            }
        }
        else {
            String name = vc.getAttributeAsStringNoException(NAME_KEY);
            if (name != null) {
                entriesToNames.put(NAME_KEY, name);
            }
            else { // Check all INFO fields for a match:
                for (Map.Entry<String, Object> entry : vc.getAttributes().entrySet()) {
                    String key = entry.getKey();
                    if (getName2 && key.startsWith(NAME2_KEY))
                        entriesToNames.put(key, entry.getValue().toString());
                    else if (key.startsWith(NAME_KEY) && !key.startsWith(NAME2_KEY))
                        entriesToNames.put(key, entry.getValue().toString());
                }
            }
        }

        return entriesToNames;
    }

    private static Map<String, String> getRefSeqEntriesToNames(VariantContext vc) {
        return getRefSeqEntriesToNames(vc, false);
    }

    public static Set<String> getRefSeqNames(VariantContext vc, boolean getName2) {
        return new TreeSet<String>(getRefSeqEntriesToNames(vc, getName2).values());
    }

    public static Set<String> getRefSeqNames(VariantContext vc) {
        return getRefSeqNames(vc, false);
    }

    public static Map<String, Object> getMergedRefSeqAttributes(VariantContext vc1, VariantContext vc2) {
        Map<String, Object> refSeqAttribs = new HashMap<String, Object>();

        List<RefSeqEntry> list1 = getAllRefSeqEntries(vc1);
        List<RefSeqEntry> list2 = getAllRefSeqEntries(vc2);
        boolean addSuffix = list1.size() > 1 || list2.size() > 1;
        int count = 1;

        for (RefSeqEntry refseq1 : list1) {
            for (RefSeqEntry refseq2 : list2) {
                Set<String> keys = new HashSet<String>();
                keys.addAll(refseq1.info.keySet());
                keys.addAll(refseq2.info.keySet());

                String keySuffix = "";
                if (addSuffix)
                    keySuffix = "_" + count++;

                Object name1 = refseq1.info.get(NAME_KEY);
                Object name2 = refseq2.info.get(NAME_KEY);
                boolean sameGene = name1 != null && name2 != null && name1.equals(name2);

                for (String key : keys) {
                    Object obj1 = refseq1.info.get(key);
                    Object obj2 = refseq2.info.get(key);
                    if (obj1 == null)
                        obj1 = "";
                    if (obj2 == null)
                        obj2 = "";

                    if (sameGene && key.equals(CODON_KEY) && obj1.equals(obj2)) // vc1 and vc2 have variants in the same codon in the same gene
                        System.out.println(vc1.getChr() + ":" + vc1.getStart() + " --> CODON: obj1 = " + obj1);

                    String useKey = key + keySuffix;
                    String mergedVal = obj1 + "\\" + obj2;
                    refSeqAttribs.put(useKey, mergedVal);
                }
            }
        }

        return refSeqAttribs;
    }

    private static List<RefSeqEntry> getAllRefSeqEntries(VariantContext vc) {
        List<RefSeqEntry> allRefSeq = new LinkedList<RefSeqEntry>();

        for (Map.Entry<String, String> entryToName : getRefSeqEntriesToNames(vc).entrySet()) {
            String entry = entryToName.getKey();
            String entrySuffix = entry.replaceFirst(NAME_KEY, "");
            allRefSeq.add(new RefSeqEntry(vc, entrySuffix));
        }

        return allRefSeq;
    }

    private static class RefSeqEntry {
        public Map<String, Object> info;

        public RefSeqEntry(VariantContext vc, String entrySuffix) {
            this.info = new HashMap<String, Object>();

            for (Map.Entry<String, Object> attribEntry : vc.getAttributes().entrySet()) {
                String key = attribEntry.getKey();
                if (key.startsWith(REFSEQ_PREFIX) && key.endsWith(entrySuffix)) {
                    String genericKey = key.replaceAll(entrySuffix, "");
                    this.info.put(genericKey, attribEntry.getValue());
                }
            }
        }
    }
}