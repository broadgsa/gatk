/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.gatk.tools.walkers.rnaseq;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.engine.walkers.DisabledReadFilters;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.tools.walkers.coverage.CoverageUtils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.io.PrintStream;
import java.util.List;

/**
 * Calculated allele counts in specific loci, to allow allele specific expression analysis.
 *
 * <p>
 * Calculated allele counts in specific loci after filtering based on mapping quality, base quality, depth of coverage, overlapping paired reads and deletion on the position.
 * All thresholds and option can be set as command line arguments.
 * The output of this tool is a tab delimited txt file, compatible with Mamba, a downstream tool (http://www.well.ox.ac.uk/~rivas/mamba/).
 * A user can use -drf DuplicateRead to allow the tool to count also duplicated reads  (a useful option for RNAseq data tools)
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * BAM files (with proper headers) to be analyzed for ASE * </p>
 * <p>
 * a VCF file with specific sites to process
 * </p>
 * </p></p>
 * <h3>Output</h3>
 * <p>
 * Table (tab delimited file) with the specific allele counts
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ASEReadCounter \
 *   -o file_name \
 *   -I input.bam \
 *   -sites sites.vcf
 *   -U ALLOW_N_CIGAR_READS
 *   [-minDepth 10] \
 *   [--minMappingQuality 10] \
 *   [--minBaseQuality 2]
 *   [-drf DuplicateRead]
 * </pre>
 */
@Downsample(by = DownsampleType.BY_SAMPLE, toCoverage = 10000)
//@DisabledReadFilters({DuplicateReadFilter.class})  //currently can be disabled using the command line argument -drf DuplicateRead
public class ASEReadCounter extends LocusWalker<String, Integer> {

    @Output
    public PrintStream out;

    @Input (fullName = "sitesVCFFile",shortName = "sites")
    public RodBinding<VariantContext> sites;

    /**
     * Loci with total depth lower then this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minDepthOfNonFilteredBase", shortName = "minDepth", doc = "Minimum number of filtered-pass bases that we need to see for reporting this site", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    public int minDepthOfNonFilteredBases = -1;

    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    public int minMappingQuality = 0;

    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    public byte minBaseQuality = 0;

    /**
     * There are few option to deal with an overlapping read pair.
     * COUNT_READS -  Count all reads independently (even if from the same fragment).
     * COUNT_FRAGMENTS - Count all fragments (even if the reads that compose the fragment are not consistent at that base).
     * COUNT_FRAGMENTS_REQUIRE_SAME_BASEQUIRE_SAME_BASE - Count all fragments (but only if the reads that compose the fragment are consistent at that base).
     * defaults to COUNT_FRAGMENTS_REQUIRE_SAME_BASEQUIRE_SAME_BASE
     */
    @Argument(fullName = "countOverlapReadsType", shortName = "overlap", doc = "How should overlapping reads from the same fragment be handled. Current options: COUNT_READS, COUNT_FRAGMENTS and COUNT_FRAGMENTS_REQUIRE_SAME_BASEQUIRE_SAME_BASE (default)", required = false)
    public CoverageUtils.CountPileupType countType = CoverageUtils.CountPileupType.COUNT_FRAGMENTS_REQUIRE_SAME_BASE;

    /**
     * Output file format (e.g. csv, table, rtable); defaults to r-readable table.
     */
    @Argument(fullName = "outputFormat", doc = "The format of the output file, can be csv, table, rtable", required = false)
    public String outputFormat = "rtable";

    /**
     * Consider a spanning deletion as contributing to coverage. Also enables deletion counts in per-base output.
     */
    @Advanced
    @Argument(fullName = "includeDeletions", shortName = "dels", doc = "Include information on deletions", required = false)
    public boolean includeDeletions = false;

    @Advanced
    @Argument(fullName = "ignoreDeletionSites", doc = "Ignore sites consisting only of deletions", required = false)
    public boolean ignoreDeletionSites = false;

    final String[] OUTPUT_FORMATS = {"table","rtable","csv"};
    public String separator = "\t";

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public boolean includeReadsWithDeletionAtLoci() { return includeDeletions && ! ignoreDeletionSites; }

    public void initialize() {

        // Check the output format
        boolean goodOutputFormat = false;
        for ( final String f : OUTPUT_FORMATS ) {
            goodOutputFormat = goodOutputFormat || f.equals(outputFormat);
        }

        if ( ! goodOutputFormat ) {
            throw new IllegalArgumentException("Improper output format. Can be one of table,rtable,csv. Was "+outputFormat);
        }

        if ( outputFormat.equals("csv") ) {
            separator = ",";
        }
        final String header = "contig"+separator+"position"+separator+"variantID"+separator+"refAllele"+separator+"altAllele"+separator+"refCount"+separator+"altCount"+separator+"totalCount"+separator+"lowMAPQDepth"+separator+"lowBaseQDepth"+separator+"rawDepth"+separator+"otherBases"+separator+"improperPairs";
        out.println(header);

    }


    @Override
    public String map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null )
            return null;
        final String contig = context.getLocation().getContig();
        final long position = context.getPosition();

        final char refAllele = (char)ref.getBase();

        final List<VariantContext> VCs =  tracker.getValues(sites, context.getLocation());
        if(VCs != null && VCs.size() > 1)
            throw new UserException("more then one variant context in position: "+contig+":"+position);
        if(VCs == null || VCs.size() == 0)
            return null;

        final VariantContext vc = VCs.get(0);
        if(!vc.isBiallelic()) {
            logger.warn("ignore site: cannot run ASE on non-biallelic sites: " + vc.toString());
            return null;
        }

        final char altAllele = (char)vc.getAlternateAllele(0).getBases()[0];
        final String siteID = vc.getID();
        final ReadBackedPileup pileup = filterPileup(context.getBasePileup(), countType, includeReadsWithDeletionAtLoci());

        // count up the depths of all and QC+ bases
        return calculateLineForSite(pileup, contig, position, siteID, refAllele, altAllele);

    }

    protected ReadBackedPileup filterPileup(final ReadBackedPileup originalPileup, final CoverageUtils.CountPileupType countType, final boolean includeDeletions){

        ReadBackedPileup pileupWithDeletions;
        if(countType.equals(CoverageUtils.CountPileupType.COUNT_FRAGMENTS_REQUIRE_SAME_BASE))
            pileupWithDeletions = originalPileup.getOverlappingFragmentFilteredPileup(true,true);
        else if(countType.equals(CoverageUtils.CountPileupType.COUNT_READS))
            pileupWithDeletions = originalPileup;
        else if(countType.equals(CoverageUtils.CountPileupType.COUNT_FRAGMENTS))
            pileupWithDeletions = originalPileup.getOverlappingFragmentFilteredPileup(false,true);
        else
            throw new UserException("Must use valid CountPileupType");

        return includeDeletions ? pileupWithDeletions: pileupWithDeletions.getPileupWithoutDeletions();

    }

    protected String calculateLineForSite(final ReadBackedPileup pileup, final String contig, final long position, final String siteID, final char refAllele, final char altAllele){

        int rawDepth = 0, lowBaseQDepth = 0, lowMAPQDepth = 0, refCount = 0, altCount = 0, totalNonFilteredCount = 0, otherBasesCount = 0, improperPairsCount = 0 ;

        for (final PileupElement base : pileup) {
            rawDepth++;

            if (!base.getRead().getProperPairFlag()){
                improperPairsCount++;
                continue;
            }
            if (base.getMappingQual() < minMappingQuality) {
                lowMAPQDepth++;
                continue;
            }

            if (base.getQual() < minBaseQuality) {
                lowBaseQDepth++;
                continue;
            }

            if(base.getBase() == refAllele)
                refCount++;
            else if(base.getBase() == altAllele)
                altCount++;
            else {
                otherBasesCount++;
                continue;
            }
            totalNonFilteredCount++;
        }

        if(totalNonFilteredCount < minDepthOfNonFilteredBases)
            return null;

        return contig +separator+
                position +separator+
                siteID +separator+
                refAllele +separator+
                altAllele +separator+
                refCount +separator+
                altCount +separator+
                totalNonFilteredCount +separator+
                lowMAPQDepth +separator+
                lowBaseQDepth +separator+
                rawDepth +separator+
                otherBasesCount +separator+
                improperPairsCount;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(String results, Integer sum) {
        if(results!= null)
            out.println(results);
        return ++sum;
    }

    @Override
    public void onTraversalDone(Integer sum) {
        logger.info("Done processing "+sum+" loci");
        out.close();
    }





}
