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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.StandardVCFWriter;

import java.io.File;
import java.util.*;

import net.sf.picard.liftover.LiftOver;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

/**
 * Lifts a VCF file over from one build to another.  Note that the resulting VCF could be mis-sorted.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
public class LiftoverVariants extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected File file = null;
    protected StandardVCFWriter writer = null;

    @Argument(fullName="chain", shortName="chain", doc="Chain file", required=true)
    protected File CHAIN = null;

    @Argument(fullName="newSequenceDictionary", shortName="dict", doc="Sequence .dict file for the new build", required=true)
    protected File NEW_SEQ_DICT = null;

    private LiftOver liftOver;

    private long successfulIntervals = 0, failedIntervals = 0;

    public void initialize() {
        liftOver = new LiftOver(CHAIN);
        liftOver.setLiftOverMinMatch(LiftOver.DEFAULT_LIFTOVER_MINMATCH);

        final SAMFileHeader toHeader = new SAMFileReader(NEW_SEQ_DICT).getFileHeader();
        liftOver.validateToSequences(toHeader.getSequenceDictionary());

        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList("variant"));
        Map<String, VCFHeader> vcfHeaders = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList("variant"));

        final VCFHeader vcfHeader = new VCFHeader(vcfHeaders.containsKey("variant") ? vcfHeaders.get("variant").getMetaData() : null, samples);
        writer = new StandardVCFWriter(file, false);
        writer.writeHeader(vcfHeader);
    }

    private void convertAndWrite(VariantContext vc, ReferenceContext ref) {

        final Interval fromInterval = new Interval(vc.getChr(), vc.getStart(), vc.getStart(), false, String.format("%s:%d", vc.getChr(), vc.getStart()));
        final int length = vc.getEnd() - vc.getStart();
        final Interval toInterval = liftOver.liftOver(fromInterval);
        VariantContext originalVC = vc;

        if ( toInterval != null ) {
            // check whether the strand flips, and if so reverse complement everything
            // TODO -- make this work for indels (difficult because the 'previous base' context needed will be changing based on indel type/size)
            if ( fromInterval.isPositiveStrand() != toInterval.isPositiveStrand() && (vc.isSNP() || !vc.isVariant()) ) {
                vc = VariantContextUtils.reverseComplement(vc);
            }

            vc = VariantContextUtils.modifyLocation(vc, GenomeLocParser.createPotentiallyInvalidGenomeLoc(toInterval.getSequence(), toInterval.getStart(), toInterval.getStart() + length));
            VariantContext newVC = VariantContext.createVariantContextWithPaddedAlleles(vc, ref.getBase(), false);

            if ( originalVC.isVariant() && VariantContextUtils.getSNPSubstitutionType(originalVC) != VariantContextUtils.getSNPSubstitutionType(newVC) ) {
                logger.warn(String.format("VCF at %s / %d => %s / %d is switching substitution type %s/%s to %s/%s",
                        originalVC.getChr(), originalVC.getStart(), newVC.getChr(), newVC.getStart(),
                        originalVC.getReference(), originalVC.getAlternateAllele(0), newVC.getReference(), newVC.getAlternateAllele(0)));
            }

            writer.add(vc, ref.getBase());
            successfulIntervals++;
        } else {
            failedIntervals++;
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, false);
        for ( VariantContext vc : VCs )
            convertAndWrite(vc, ref);

        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public void onTraversalDone(Integer result) {
        System.out.println("Converted " + successfulIntervals + " records; failed to convert " + failedIntervals + " records.");
        writer.close();
    }
}
