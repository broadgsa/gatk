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

package org.broadinstitute.sting.playground.gatk.walkers.phasing;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.phasing.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.*;
import java.util.*;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.READS, DataSource.REFERENCE})
@Requires(value = {DataSource.READS, DataSource.REFERENCE}, referenceMetaData = @RMD(name = "variant", type = ReferenceOrderedDatum.class))
@By(DataSource.READS)

@ReadFilters({ZeroMappingQualityReadFilter.class})
// Filter out all reads with zero mapping quality

public class ReadBasedPhasingValidationWalker extends RodWalker<Integer, Integer> {
    private LinkedList<String> rodNames = null;

    @Argument(fullName = "sitePairsFile", shortName = "sitePairsFile", doc = "File of pairs of variants for which phasing in ROD should be assessed using input reads", required = true)
    protected File sitePairsFile = null;

    @Output
    protected PrintStream out;

    private Set<SitePair> sitePairs = null;
    private String sampleName = null;

    SiteGenotypeAndReads prevSiteAndReads = null;

    private final static int NUM_IN_PAIR = 2; // trivial

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add("variant");

        sitePairs = new TreeSet<SitePair>();
        GenomeLocParser locParser = getToolkit().getGenomeLocParser();

        InputStream sitePairsStream = null;
        try {
            sitePairsStream = new FileInputStream(sitePairsFile);
        }
        catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
            throw new UserException("Problem opening file: " + sitePairsFile);
        }

        AsciiLineReader sitePairsReader = new AsciiLineReader(sitePairsStream);
        while (true) {
            String line = null;
            try {
                line = sitePairsReader.readLine();
            }
            catch (IOException ioe) {
                ioe.printStackTrace();
                throw new UserException("Problem reading file: " + sitePairsFile);
            }
            if (line == null)
                break; // reached end of file

            String[] twoSites = line.split("\t");
            if (twoSites.length != 2)
                throw new UserException("Must have PAIRS of sites in line " + line + " of " + sitePairsFile);

            SitePair sp = new SitePair(locParser.parseGenomeLoc(twoSites[0]), locParser.parseGenomeLoc(twoSites[1]));
            sitePairs.add(sp);
        }
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

        boolean relevantSitePair = false;
        SitePair sp = null;
        if (prevSiteAndReads != null) {
            // all vc's below start at ref.getLocus() [due to requireStartHere = true]:
            sp = new SitePair(prevSiteAndReads.site, ref.getLocus());
            relevantSitePair = sitePairs.contains(sp);
        }

        if (context == null || !context.hasBasePileup())
            return null;
        ReadBackedPileup pileup = context.getBasePileup();
        String nextName = null;

        //
        // TODO: ASSUMES THAT ALL READS COME FROM A SINGLE SAMPLE:
        // TODO: CODE SHOULD BE AS FOLLOWS:
        /*
        Collection<String> sampNames = pileup.getSampleNames();
        if (sampNames.size() != 1)
            throw new UserException("Reads must be for exactly one sample [not multi-sample]");
        nextName = sampNames.iterator().next();
        if (nextName == null)
            throw new UserException("Reads must be for exactly one sample");

        if (sampleName == null)
            sampleName = nextName;
        else if (!nextName.equals(sampleName))
            throw new UserException("Reads must have a single consistent sample name");

        pileup = pileup.getPileupForSampleName(sampleName);
        */
        //

        ReadBasesAtPosition readBases = new ReadBasesAtPosition();
        for (PileupElement p : pileup)
            readBases.putReadBase(p);

        ReadCounter rdCounts = null;
        if (relevantSitePair) { // otherwise, processed the reads for their possible use in the future:
            PhasingReadList buildReads = new PhasingReadList(NUM_IN_PAIR);
            buildReads.updateBases(0, prevSiteAndReads.readBases);
            buildReads.updateBases(1, readBases);

            List<PhasingRead> reads = new LinkedList<PhasingRead>();
            for (PhasingRead rd : buildReads) {
                if (rd.getNonNullIndices().length == NUM_IN_PAIR) // only want reads with BOTH bases called [possibly as deleted ("D")]
                    reads.add(rd);
            }

            // Count the occurence of each "haplotype":
            rdCounts = new ReadCounter();
            for (PhasingRead rd : reads)
                rdCounts.incrementCount(rd);
        }


        // Now, read the ROD and note the genotypes and their phase to be validated:
        Set<Haplotype> calledHaplotypes = null;
        List<Haplotype> allPossibleHaplotypes = null;

        boolean requireStartHere = true; // only see each VariantContext once
        boolean takeFirstOnly = true; // take only the first entry from the ROD file
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            if (vc.isFiltered() || !vc.isSNP())
                continue;

            if (vc.getNSamples() != 1)
                throw new UserException("ROD file must have exactly one sample [not multi-sample]");
            nextName = vc.getSampleNames().iterator().next();
            if (sampleName == null)
                sampleName = nextName;
            else if (!nextName.equals(sampleName))
                throw new UserException("ROD must have a single consistent sample name");

            Genotype gt = vc.getGenotype(sampleName);

            if (relevantSitePair) {
                Genotype prevGt = prevSiteAndReads.gt;
                List<Allele> prevAlleles = prevGt.getAlleles();
                List<Allele> curAlleles = gt.getAlleles();

                calledHaplotypes = new TreeSet<Haplotype>(); // implemented Haplotype.compareTo()
                if (gt.isPhased()) {
                    if (gt.getPloidy() != prevGt.getPloidy())
                        throw new UserException("Invalid ROD file: cannot be phased AND have different ploidys!");

                    // Consider only the haplotypes called to be phased
                    Iterator<Allele> curAllIt = curAlleles.iterator();
                    for (Allele prevAll : prevAlleles) {
                        Allele curAll = curAllIt.next();
                        calledHaplotypes.add(successiveAllelesToHaplotype(prevAll, curAll));
                    }
                }

                // Consider EVERY combination of alleles as haplotypes [IF PHASED, this will give the contingency table in the CORRECT order]:
                allPossibleHaplotypes = new LinkedList<Haplotype>();
                for (Allele prevAll : prevAlleles) {
                    for (Allele curAll : curAlleles) {
                        allPossibleHaplotypes.add(successiveAllelesToHaplotype(prevAll, curAll));
                    }
                }
            }

            prevSiteAndReads = new SiteGenotypeAndReads(ref.getLocus(), gt, readBases);
        }

        int processedPairs = 0;
        if (relevantSitePair) {
            Map<Haplotype, Integer> haplotypeCounts = new TreeMap<Haplotype, Integer>(); // implemented Haplotype.compareTo()

            processedPairs = 1;
            int totalCount = rdCounts.totalCount();
            System.out.println("\nPair: " + sp + " [# reads = " + totalCount + "]");

            int matchCount = 0;
            for (Map.Entry<PhasingRead, Integer> rdEntry : rdCounts.entrySet()) {
                PhasingRead read = rdEntry.getKey();
                int count = rdEntry.getValue();

                Haplotype readsHaplotype = new Haplotype(read);
                haplotypeCounts.put(readsHaplotype, count);

                boolean readMatchesCalledHaplotype = calledHaplotypes != null && calledHaplotypes.contains(readsHaplotype);
                if (readMatchesCalledHaplotype)
                    matchCount += count;

                System.out.println("read" + ": " + read + (readMatchesCalledHaplotype ? "*" : "") + "\tcount: " + count);
            }

            double percentMatchingReads = 100 * (matchCount / (double) totalCount);
            System.out.println("% MATCHING reads: " + percentMatchingReads + " [of " + totalCount + " TOTAL reads]");

            out.print(sp);
            for (Haplotype hap : allPossibleHaplotypes) {
                Integer count = haplotypeCounts.get(hap);
                if (count == null) // haplotype may not have been observed in ANY reads
                    count = 0;

                out.print("\t" + count);
            }
            out.println();
        }

        return processedPairs;
    }

    private Haplotype successiveAllelesToHaplotype(Allele prevAll, Allele curAll) {
        byte prevBase = SNPallelePair.getSingleBase(prevAll);
        byte curBase = SNPallelePair.getSingleBase(curAll);

        byte[] hapBases = new byte[NUM_IN_PAIR];
        hapBases[0] = prevBase;
        hapBases[1] = curBase;
        return new Haplotype(hapBases);
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
        System.out.println("Validated " + result + " pairs of sites.");
    }
}

class SitePair implements Comparable<SitePair> {
    public GenomeLoc site1;
    public GenomeLoc site2;

    public SitePair(GenomeLoc site1, GenomeLoc site2) {
        if (site1.size() > 1 || site2.size() > 1)
            throw new UserException("Must give pairs of SINGLE-LOCUS record start sites");

        this.site1 = site1;
        this.site2 = site2;
    }

    public String toString() {
        return site1.toString() + "\t" + site2.toString();
    }

    public int compareTo(SitePair other) {
        int comp1 = site1.compareTo(other.site1);
        if (comp1 != 0)
            return comp1;

        return site2.compareTo(other.site2);
    }
}

class SiteGenotypeAndReads {
    public GenomeLoc site;
    public Genotype gt;
    public ReadBasesAtPosition readBases;

    public SiteGenotypeAndReads(GenomeLoc site, Genotype gt, ReadBasesAtPosition readBases) {
        this.site = site;
        this.gt = gt;
        this.readBases = readBases;
    }
}

class PhasingReadList implements Iterable<PhasingRead> {
    private Map<String, PhasingRead> readsAtSites = null;
    private int numSites;

    public PhasingReadList(int numSites) {
        this.readsAtSites = new HashMap<String, PhasingRead>();
        this.numSites = numSites;
    }

    public void updateBases(int index, ReadBasesAtPosition readBases) {
        if (readBases == null)
            return;

        for (ReadBase rb : readBases) {
            String readName = rb.readName;

            PhasingRead rd = readsAtSites.get(readName);
            if (rd == null) {
                rd = new PhasingRead(numSites, rb.mappingQual);
                readsAtSites.put(readName, rd);
            }

            // Arbitrarily updates to the last base observed for this sample and read (rb.base):
            rd.updateBaseAndQuality(index, rb.base, rb.baseQual);
        }
    }

    public Iterator<PhasingRead> iterator() {
        return readsAtSites.values().iterator();
    }

    public int size() {
        return readsAtSites.size();
    }
}

class ReadCounter {
    private Map<PhasingRead, Integer> counts;
    private int totalCount;

    public ReadCounter() {
        this.counts = new TreeMap<PhasingRead, Integer>(); // implemented PhasingRead.compareTo()
    }

    public void incrementCount(PhasingRead rd) {
        Integer cnt = counts.get(rd);
        if (cnt == null)
            cnt = 0;

        counts.put(rd, cnt + 1);
        totalCount++;
    }

    public Set<Map.Entry<PhasingRead, Integer>> entrySet() {
        return counts.entrySet();
    }

    public int totalCount() {
        return totalCount;
    }
}