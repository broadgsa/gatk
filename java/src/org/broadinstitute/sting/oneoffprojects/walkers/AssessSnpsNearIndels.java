/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.iterators.DownsampleIterator;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.io.Writer;
import java.security.PrivilegedActionException;
import java.util.*;

/**
 * Assesses distance of SNP calls to nearest indel call in Exomes.
 * Use -B:snps,vcf and -B:indels,vcf
 */
public class AssessSnpsNearIndels extends RodWalker<Integer, Integer> {

    private class TiTv {
        int TiCount = 0, TvCount = 0;

        public TiTv() {};
    }

    private static final int Bin0to5 = 0;
    private static final int Bin6to10 = 1;
    private static final int Bin11to15 = 2;
    private static final int Bin16to20 = 3;
    private static final int Bin21to25 = 4;
    private static final int Bin26to30 = 5;
    private static final int BinMoreThan30 = 6;

    private GenomeLoc previousIndel = null;
    private ArrayList<VariantContext> snpQueue = new ArrayList<VariantContext>();
    private TiTv[] counts = new TiTv[7];
    private GenomeLocParser GLparser = null;

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out = null;

    public void initialize() {
        GLparser = getToolkit().getGenomeLocParser();

        for (int i = 0; i < 7; i++)
            counts[i] = new TiTv();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        VariantContext snp = tracker.getVariantContext(ref, "snps", null, ref.getLocus(), true);
        VariantContext indel = tracker.getVariantContext(ref, "indels", null, ref.getLocus(), true);

        // first add the snp if available
        if ( snp != null && !snp.isFiltered() ) {

            // flush the queue on a new contig
            if ( !snpQueue.isEmpty() && !snpQueue.get(0).getChr().equals(snp.getChr()) ) {
                for ( VariantContext vc : snpQueue )
                    calculateDistance(vc, previousIndel, null);
                snpQueue.clear();
            }

            snpQueue.add(snp);
        }

        // then look for the indel
        if ( indel != null && !indel.isFiltered() ) {

            GenomeLoc loc = GLparser.createGenomeLoc(indel.getChr(), indel.getStart(), indel.getEnd());
            boolean sameContig = !snpQueue.isEmpty() && snpQueue.get(0).getChr().equals(indel.getChr());

            // flush the queue
            for ( VariantContext vc : snpQueue )
                calculateDistance(vc, previousIndel, sameContig ? loc : null);
            snpQueue.clear();

            previousIndel = loc;
        }

        return 1;
    }

    private void calculateDistance(VariantContext snp, GenomeLoc previousIndel, GenomeLoc nextIndel) {

        GenomeLoc loc = GLparser.createGenomeLoc(snp.getChr(), snp.getStart(), snp.getEnd());

        int previousDistance = -1, nextDistance = -1;

        if ( previousIndel != null ) {
            // watch out for spanning deletions
            if ( previousIndel.getStop() > snp.getStart() )
                previousDistance = 0;
            else
                previousDistance = snp.getStart() - previousIndel.getStop();
        }

        if ( nextIndel != null )
            nextDistance = nextIndel.getStart() - loc.getStart();

        if ( previousDistance == -1 && nextDistance == -1 )
            return;

        int distance = -1;
        if ( previousDistance == -1 )
            distance = nextDistance;
        else if ( nextDistance == -1 )
            distance = previousDistance;
        else
            distance = Math.min(previousDistance, nextDistance);

        TiTv obj;
        if ( distance < 0 )
             throw new IllegalStateException("Found a negative distance at " + loc);
        else if ( distance < 6 )
            obj = counts[Bin0to5];
        else if ( distance < 11 )
            obj = counts[Bin6to10];
        else if ( distance < 16 )
            obj = counts[Bin11to15];
        else if ( distance < 21 )
            obj = counts[Bin16to20];
        else if ( distance < 26 )
            obj = counts[Bin21to25];
        else if ( distance < 31 )
            obj = counts[Bin26to30];
        else
            obj = counts[BinMoreThan30];

        if ( BaseUtils.isTransition(snp.getReference().getBases()[0], snp.getAlternateAllele(0).getBases()[0]) )
            obj.TiCount++;
        else
            obj.TvCount++;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        // flush the queue
        for ( VariantContext vc : snpQueue )
            calculateDistance(vc, previousIndel, null);

        out.println("Bin\tnumTi\tnumTv\tTi/Tv");
        printLine(counts[Bin0to5], "0to5");
        printLine(counts[Bin6to10], "6to10");
        printLine(counts[Bin11to15], "11to15");
        printLine(counts[Bin16to20], "16to20");
        printLine(counts[Bin21to25], "21to25");
        printLine(counts[Bin26to30], "26to30");
        printLine(counts[BinMoreThan30], ">30");
    }

    private void printLine(TiTv obj, String s) {
        out.println(String.format("%s\t%d\t%d\t%.2f", s, obj.TiCount, obj.TvCount, ((double)obj.TiCount/(double)obj.TvCount)));
    }
}
