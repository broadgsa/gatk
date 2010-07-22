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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFHeader;

import java.util.*;

import net.sf.samtools.SAMReadGroupRecord;


/**
 * A walker that merges multiple batches of calls, by calling into the Genotyper to fill in sites that
 * were called in one batch but not another.
 */
@Reference(window=@Window(start=-20,stop=20))
@By(DataSource.REFERENCE)
@Requires(value={},referenceMetaData=@RMD(name="trigger", type=ReferenceOrderedData.class))
public class BatchedCallsMerger extends LocusWalker<VariantContext, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    @Argument(doc = "VCF file to which variants should be written", required = false)
    public GenotypeWriter writer = null;

    @Argument(fullName="rod_list", shortName="rods", doc="A comma-separated string describing the rod names representing individual call batches", required=true)
    protected String ROD_STRING = null;
    private HashSet<String> targetRods;

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

    // all samples to be used
    private Set<String> samples = new HashSet<String>();

    // mapping from rod name to set of samples coming from it
    private Map<String, Set<String>> rodsToSamples = new HashMap<String, Set<String>>();

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable extended events for indels
    public boolean generateExtendedEvents() { return UAC.genotypeModel == GenotypeCalculationModel.Model.INDELS; }

    public void initialize() {

        targetRods = new HashSet<String>(Arrays.asList(ROD_STRING.split(",")));

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            RMDTrack rod = source.getReferenceOrderedData();
            // if it's one of our target rods
            if ( targetRods.contains(rod.getName()) ) {
                // read the samples and store them
                VCFReader reader = new VCFReader(rod.getFile());
                HashSet<String> mySamples = new HashSet<String>(reader.getHeader().getGenotypeSamples());
                rodsToSamples.put(rod.getName(), mySamples);
                samples.addAll(mySamples);

                // while we're here, pull out the header lines
                headerLines.addAll(reader.getHeader().getMetaData());

                reader.close();
            }
        }

        // update the engine
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, writer, null, null);
        UG_engine.samples = samples;        

        // initialize the header
        GenotypeWriterFactory.writeHeader(writer, GenomeAnalysisEngine.instance.getSAMFileHeader(), new VCFHeader(headerLines, samples));
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        // get the calls at this locus
        Collection<VariantContext> VCs = tracker.getAllVariantContexts(ref, null, context.getLocation(), true, false);

        Set<VariantContext> calls = new HashSet<VariantContext>();
        Set<String> seenRods = new HashSet<String>();
        for ( VariantContext vc : VCs ) {
            if ( targetRods.contains(vc.getName()) ) {
                calls.add(vc);
                seenRods.add(vc.getName());
            }
        }

        // if there are no calls, ignore this site
        if ( seenRods.size() == 0 )
            return null;

        // figure out which samples still need to be called
        Set<String> missedSamples = new HashSet<String>();
        for ( String rod : targetRods ) {
            if ( !seenRods.contains(rod) )
                missedSamples.addAll(rodsToSamples.get(rod));
        }

        // we are missing samples, call them
        if ( missedSamples.size() > 0 ) {
            AlignmentContext prunedContext = filterForSamples(context.getBasePileup(), missedSamples);
            VariantCallContext vcc = UG_engine.runGenotyper(tracker, ref, prunedContext);
            if ( vcc != null && vcc.vc != null )
                calls.add(vcc.vc);
        }

        // merge the variant contexts
        return VariantContextUtils.simpleMerge(calls, ref.getBases());
    }

    private AlignmentContext filterForSamples(ReadBackedPileup pileup, Set<String> samples) {

        ArrayList<PileupElement> newPileup = new ArrayList<PileupElement>();

        for (PileupElement p : pileup ) {
            SAMReadGroupRecord readGroup = p.getRead().getReadGroup();
            if ( readGroup != null && samples.contains(readGroup.getSample()) )
                newPileup.add(p);
        }
        return new AlignmentContext(pileup.getLocation(), new ReadBackedPileupImpl(pileup.getLocation(), newPileup));

    }

    public Integer reduceInit() { return 0; }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public Integer reduce(VariantContext value, Integer sum) {
        // can't call the locus because of no coverage or no confidence
        if ( value == null )
            return sum;

        try {
            writer.addCall(value, value.getReference().getBases()[0]);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage() + "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
        }

        return sum;
    }

    public void onTraversalDone(Integer sum) {}
}