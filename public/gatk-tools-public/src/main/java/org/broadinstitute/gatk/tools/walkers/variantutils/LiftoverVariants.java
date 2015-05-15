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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.File;
import java.util.*;

/**
 * Lifts a VCF file over from one build to another
 *
 * <p>"Lifting over" variants means adjusting variant calls from one reference to another. Specifically, the process
 * adjusts the position of the call to match the corresponding position on the target reference. For example, if you
 * have variants called from reads aligned to the hg19 reference, and you want to compare them to calls made based on
 * the b37 reference, you need to liftover one of the callsets to the other reference.</p>
 *
 * <p>LiftoverVariants is intended to be the first of two processing steps for the liftover process.
 * The second step is to run FilterLiftedVariants on the output of LiftoverVariants. This will produce valid
 * well-behaved VCF files, where you'll see that the contig names in the header have all been correctly replaced.</p>
 *
 * <h3>Caveat</h3>
 * <p>To be clear, the VCF resulting from the LiftoverVariants run is not guaranteed to be valid according to the official specification.  The file could
 * possibly be mis-sorted and the header may not be complete. That is why you need to run FilterLiftedVariants on it.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant call set to lift over, the sequence  dictionary of the new reference build and the appropriate liftover
 * chain file.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * The lifted-over call set.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T LiftoverVariants \
 *   -R reference_hg19.fasta \
 *   -V input_hg19.vcf \
 *   -chain liftover_hg19_to_b37.txt \
 *   -dict reference_b37.dict \
 *   -o liftedover_output_b37.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class LiftoverVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written", required=true, defaultToStdout=false)
    protected File file = null;
    protected VariantContextWriter writer = null;

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
        } catch (RuntimeException e) {
            throw new UserException.BadInput("there is a problem with the chain file you are using: " + e.getMessage());
        }

        liftOver.setLiftOverMinMatch(LiftOver.DEFAULT_LIFTOVER_MINMATCH);

        try {
            final SAMFileHeader toHeader = new SAMFileReader(NEW_SEQ_DICT).getFileHeader();
            liftOver.validateToSequences(toHeader.getSequenceDictionary());
        } catch (RuntimeException e) {
            throw new UserException.BadInput("the chain file you are using is not compatible with the reference you are trying to lift over to; please use the appropriate chain file for the given reference");    
        }

        String trackName = variantCollection.variants.getName();
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(trackName));
        Map<String, VCFHeader> vcfHeaders = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));

        Set<VCFHeaderLine> metaData = new HashSet<>();
        if ( vcfHeaders.containsKey(trackName) )
            metaData.addAll(vcfHeaders.get(trackName).getMetaDataInSortedOrder());
        if ( RECORD_ORIGINAL_LOCATION ) {
            metaData.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_CONTIG_KEY));
            metaData.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_START_KEY));
        }


        final VCFHeader vcfHeader = new VCFHeader(metaData, samples);
        writer = VariantContextWriterFactory.create(file, getMasterSequenceDictionary(), EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER));
        writer.writeHeader(vcfHeader);
    }

    private void convertAndWrite(VariantContext vc, ReferenceContext ref) {

        final Interval fromInterval = new Interval(vc.getChr(), vc.getStart(), vc.getStart(), false, String.format("%s:%d", vc.getChr(), vc.getStart()));
        final int length = vc.getEnd() - vc.getStart();
        final Interval toInterval = liftOver.liftOver(fromInterval);
        VariantContext originalVC = vc;

        if ( toInterval != null ) {
            // check whether the strand flips, and if so reverse complement everything
            if ( fromInterval.isPositiveStrand() != toInterval.isPositiveStrand() && vc.isPointEvent() ) {
                vc = GATKVariantContextUtils.reverseComplement(vc);
            }

            vc = new VariantContextBuilder(vc).loc(toInterval.getSequence(), toInterval.getStart(), toInterval.getStart() + length).make();

            if ( RECORD_ORIGINAL_LOCATION ) {
                vc = new VariantContextBuilder(vc)
                        .attribute(GATKVCFConstants.ORIGINAL_CONTIG_KEY, fromInterval.getSequence())
                        .attribute(GATKVCFConstants.ORIGINAL_START_KEY, fromInterval.getStart()).make();
            }

            if ( originalVC.isSNP() && originalVC.isBiallelic() && GATKVariantContextUtils.getSNPSubstitutionType(originalVC) != GATKVariantContextUtils.getSNPSubstitutionType(vc) ) {
                logger.warn(String.format("VCF at %s / %d => %s / %d is switching substitution type %s/%s to %s/%s",
                        originalVC.getChr(), originalVC.getStart(), vc.getChr(), vc.getStart(),
                        originalVC.getReference(), originalVC.getAlternateAllele(0), vc.getReference(), vc.getAlternateAllele(0)));
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
