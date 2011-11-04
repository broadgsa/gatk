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

import net.sf.picard.PicardException;
import net.sf.picard.liftover.LiftOver;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Lifts a VCF file over from one build to another.  Note that the resulting VCF could be mis-sorted.
 */
public class LiftoverVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written",required=true)
    protected File file = null;
    protected StandardVCFWriter writer = null;

    @Argument(fullName="chain", shortName="chain", doc="Chain file", required=true)
    protected File CHAIN = null;

    @Argument(fullName="newSequenceDictionary", shortName="dict", doc="Sequence .dict file for the new build", required=true)
    protected File NEW_SEQ_DICT = null;

    @Argument(fullName="recordOriginalLocation", shortName="recordOriginalLocation", doc="Should we record what the original location was in the INFO field?", required=false)
    protected Boolean RECORD_ORIGINAL_LOCATION = false;

    private LiftOver liftOver;

    private long successfulIntervals = 0, failedIntervals = 0;

    public void initialize() {
        try {
            liftOver = new LiftOver(CHAIN);
        } catch (PicardException e) {
            throw new UserException.BadInput("there is a problem with the chain file you are using: " + e.getMessage());
        }

        liftOver.setLiftOverMinMatch(LiftOver.DEFAULT_LIFTOVER_MINMATCH);

        try {
            final SAMFileHeader toHeader = new SAMFileReader(NEW_SEQ_DICT).getFileHeader();
            liftOver.validateToSequences(toHeader.getSequenceDictionary());
        } catch (PicardException e) {
            throw new UserException.BadInput("the chain file you are using is not compatible with the reference you are trying to lift over to; please use the appropriate chain file for the given reference");    
        }

        String trackName = variantCollection.variants.getName();
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(trackName));
        Map<String, VCFHeader> vcfHeaders = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));

        Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
        if ( vcfHeaders.containsKey(trackName) )
            metaData.addAll(vcfHeaders.get(trackName).getMetaData());
        if ( RECORD_ORIGINAL_LOCATION ) {
            metaData.add(new VCFInfoHeaderLine("OriginalChr", 1, VCFHeaderLineType.String, "Original contig name for the record"));
            metaData.add(new VCFInfoHeaderLine("OriginalStart", 1, VCFHeaderLineType.Integer, "Original start position for the record"));
        }


        final VCFHeader vcfHeader = new VCFHeader(metaData, samples);
        writer = new StandardVCFWriter(file, getMasterSequenceDictionary(), false);
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
            if ( fromInterval.isPositiveStrand() != toInterval.isPositiveStrand() && vc.isPointEvent() ) {
                vc = VariantContextUtils.reverseComplement(vc);
            }

            vc = VariantContext.modifyLocation(vc, toInterval.getSequence(), toInterval.getStart(), toInterval.getStart() + length);

            if ( RECORD_ORIGINAL_LOCATION ) {
                HashMap<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
                attrs.put("OriginalChr", fromInterval.getSequence());
                attrs.put("OriginalStart", fromInterval.getStart());
                vc = VariantContext.modifyAttributes(vc, attrs);
            }

            VariantContext newVC = VariantContext.createVariantContextWithPaddedAlleles(vc, false);
            if ( originalVC.isSNP() && originalVC.isBiallelic() && VariantContextUtils.getSNPSubstitutionType(originalVC) != VariantContextUtils.getSNPSubstitutionType(newVC) ) {
                logger.warn(String.format("VCF at %s / %d => %s / %d is switching substitution type %s/%s to %s/%s",
                        originalVC.getChr(), originalVC.getStart(), newVC.getChr(), newVC.getStart(),
                        originalVC.getReference(), originalVC.getAlternateAllele(0), newVC.getReference(), newVC.getAlternateAllele(0)));
            }

            writer.add(vc);
            successfulIntervals++;
        } else {
            failedIntervals++;
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
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
