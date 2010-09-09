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

package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMReadGroupRecord;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.BatchedCallsMerger;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;


/**
 * A walker that merges multiple batches of calls, by calling into the Genotyper to fill in sites that
 * were called in one batch but not another.
 */
//@Reference(window=@Window(start=-20,stop=20))
@By(DataSource.REFERENCE)
@Allows(value={DataSource.READS, DataSource.REFERENCE})
//@Requires(value={},referenceMetaData=@RMD(name="genotypes", type=VariantContext.class))
public class QSampleWalker extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    @Argument(shortName = "E", doc="If true, all sites will have pHet emitted into output file")
    boolean SHOW_EACH_SITE;

    @Output(doc = "QSample Log File", required = true)
    public PrintStream out;

    Map<String, ArrayList<Double>> qSampleScores = new HashMap<String, ArrayList<Double>>();

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    Set<String> samples;

    public void initialize() {
        Map<String, VCFHeader> headers = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList("genotypes"));
        for ( String sample : SampleUtils.getSampleList(headers) ) {
            qSampleScores.put(sample, new ArrayList<Double>());
        }

        // update the engine
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, null);

        Set<String> bamSamples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        Set<String> rodSamples = qSampleScores.keySet();
        samples = new LinkedHashSet<String>(bamSamples);
        samples.retainAll(rodSamples);

        logger.info("# bam samples    : " + bamSamples);
        logger.info("# rod samples    : " + rodSamples);
        logger.info("# common samples : " + samples);

        for ( String sample : samples ) {
            out.print(sample);
            out.print("\t");
        }
        out.println();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        // get the calls at this locus
        VariantContext genotypesVC = tracker.getVariantContext(ref, "genotypes", null, context.getLocation(), true);

        if ( genotypesVC == null )
            return 0; // protect ourselves for runs with BTI

        // set of het samples
        Set<String> hetSamples = new HashSet<String>();
        for ( String sample : samples ) {
            if ( genotypesVC.getGenotype(sample).isHet() )
                hetSamples.add(sample);
        }

        AlignmentContext hetsContext = BatchedCallsMerger.filterForSamples(context.getBasePileup(), hetSamples);
        VariantCallContext vcc = UG_engine.runGenotyper(tracker, ref, hetsContext);

        logger.info(String.format("Visiting site %s", context.getLocation()));
        for ( String sample : samples ) {
            Genotype g = genotypesVC.getGenotype(sample);
            double pHet = -1;

            if ( g.isHet() ) {
                pHet = scoreSample(sample, vcc);
                qSampleScores.get(sample).add(pHet);
                //logger.info(String.format("  het sample %s with score %.3f", sample, pHet));
            }

            if ( SHOW_EACH_SITE ) out.printf("%s\t", pHet == -1 ? "NA" : String.format("%.3f", pHet));
        }

        if ( SHOW_EACH_SITE ) out.println();


        return 1;
    }

    private static int HET_INDEX = 1;
    private double scoreSample(String sample, VariantCallContext vcc) {
        if ( vcc == null || vcc.vc == null || vcc.vc.getGenotype(sample) == null )
            return 0.33;
        else {
            // get GLs
            GenotypeLikelihoods gl = vcc.vc.getGenotype(sample).getLikelihoods();
            double[] log10gl = gl.getAsVector();
            return Math.pow(10, log10gl[HET_INDEX]) / MathUtils.sumLog10(log10gl);
        }
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public Integer treeReduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer sum) {
        // print out summary statistics
        for ( String sample : samples ) {
            ArrayList<Double> scores = qSampleScores.get(sample);
            double qSample = MathUtils.averageDouble(scores);
            out.printf("%.3f\t", qSample);
        }
        out.println();
    }
}