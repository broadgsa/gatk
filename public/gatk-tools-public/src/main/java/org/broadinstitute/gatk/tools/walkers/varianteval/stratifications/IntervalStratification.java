/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import htsjdk.samtools.util.IntervalTree;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Stratifies the variants by whether they overlap an interval in the set provided on the command line.
 *
 * The primary use of this stratification is to provide a mechanism to divide asssessment of a call set up
 * by whether a variant overlaps an interval or not.  I use this to differentiate between variants occurring
 * in CCDS exons vs. those in non-coding regions, in the 1000G call set, using a command line that looks like:
 *
 * -T VariantEval -R human_g1k_v37.fasta -eval 1000G.vcf -stratIntervals:BED ccds.bed -ST IntervalStratification
 *
 * Note that the overlap algorithm properly handles symbolic alleles with an INFO field END value.  In order to
 * safely use this module you should provide entire contigs worth of variants, and let the interval strat decide
 * overlap, as opposed to using -L which will not properly work with symbolic variants.
 */
public class IntervalStratification extends VariantStratifier {
    final protected static Logger logger = Logger.getLogger(IntervalStratification.class);
    Map<String, IntervalTree<GenomeLoc>> intervalTreeByContig = null;

    final List<Object> OVERLAPPING = Arrays.asList((Object)"all", (Object)"overlaps.intervals");
    final List<Object> NOT_OVERLAPPING = Arrays.asList((Object)"all", (Object)"outside.intervals");

    
    @Override
    public void initialize() {
        if ( getVariantEvalWalker().intervalsFile == null )
            throw new UserException.MissingArgument("stratIntervals", "Must be provided when IntervalStratification is enabled");

        final List<GenomeLoc> locs = getVariantEvalWalker().intervalsFile.getIntervals(getVariantEvalWalker().getToolkit().getGenomeLocParser());

        if ( locs.isEmpty() )
            throw new UserException.BadArgumentValue("stratIntervals", "Contains no intervals.  Perhaps the file is malformed or empty?");

        intervalTreeByContig = getVariantEvalWalker().createIntervalTreeByContig(getVariantEvalWalker().intervalsFile);

        logger.info(String.format("Creating IntervalStratification %s containing %d intervals covering %d bp",
                getVariantEvalWalker().intervalsFile.getSource(), locs.size(), IntervalUtils.intervalSize(locs)));

        states.addAll(Arrays.asList("all", "overlaps.intervals", "outside.intervals"));
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null) {
            final GenomeLoc loc = getVariantEvalWalker().getToolkit().getGenomeLocParser().createGenomeLoc(eval);
            IntervalTree<GenomeLoc> intervalTree = intervalTreeByContig.get(loc.getContig());
            IntervalTree.Node<GenomeLoc> node = intervalTree.minOverlapper(loc.getStart(), loc.getStop());
            //logger.info(String.format("Overlap %s found %s", loc, node));
            if ( node != null )
                return OVERLAPPING;
            else
                return NOT_OVERLAPPING;
        }

        return Collections.emptyList();
    }
}
