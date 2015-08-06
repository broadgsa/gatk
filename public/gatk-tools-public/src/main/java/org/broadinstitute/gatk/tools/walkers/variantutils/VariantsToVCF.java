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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.VariantContextAdaptors;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantOverlapAnnotator;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.codecs.hapmap.RawHapMapFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.File;
import java.util.*;

/**
 * Convert variants from other file formats to VCF format
 *
 * <p>
 * Note that there must be a Tribble feature/codec available for the file format as well as an adaptor.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant file to convert.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A VCF file.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T VariantsToVCF \
 *   -R reference.fasta \
 *   -o output.vcf \
 *   --variant:RawHapMap input.hapmap
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-40,stop=40))
public class VariantsToVCF extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter baseWriter = null;
    private VariantContextWriter vcfwriter; // needed because hapmap/dbsnp indel records move

    /**
     * Variants from this input file are used by this tool as input.
     */
    @Input(fullName="variant", shortName = "V", doc="Input variant file", required=true)
    public RodBinding<Feature> variants;

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * This argument is used for data (like GELI) with genotypes but no sample names encoded within.
     */
    @Argument(fullName="sample", shortName="sample", doc="The sample name represented by the variant rod", required=false)
    protected String sampleName = null;

    private Set<String> allowedGenotypeFormatStrings = new HashSet<String>();
    private boolean wroteHeader = false;
    private Set<String> samples;

    // for dealing with indels in hapmap
    CloseableIterator<GATKFeature> dbsnpIterator = null;
    VariantOverlapAnnotator variantOverlapAnnotator = null;

    public void initialize() {
        vcfwriter = VariantContextWriterFactory.sortOnTheFly(baseWriter, 40, false);
        variantOverlapAnnotator = new VariantOverlapAnnotator(dbsnp.dbsnp, getToolkit().getGenomeLocParser());
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || !BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        Collection<VariantContext> contexts = getVariantContexts(tracker, ref);

        for ( VariantContext vc : contexts ) {
            VariantContextBuilder builder = new VariantContextBuilder(vc);

            // set the appropriate sample name if necessary
            if ( sampleName != null && vc.hasGenotypes() && vc.hasGenotype(variants.getName()) ) {
                Genotype g = new GenotypeBuilder(vc.getGenotype(variants.getName())).name(sampleName).make();
                builder.genotypes(g);
            }

            final VariantContext withID = variantOverlapAnnotator.annotateRsID(tracker, builder.make());
            writeRecord(withID, tracker, ref.getLocus());
        }

        return 1;
    }

    private Collection<VariantContext> getVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref) {

        List<Feature> features = tracker.getValues(variants, ref.getLocus());
        List<VariantContext> VCs = new ArrayList<VariantContext>(features.size());

        for ( Feature record : features ) {
            if ( VariantContextAdaptors.canBeConvertedToVariantContext(record) ) {
                // we need to special case the HapMap format because indels aren't handled correctly
                if ( record instanceof RawHapMapFeature) {

                    // is it an indel?
                    RawHapMapFeature hapmap = (RawHapMapFeature)record;
                    if ( hapmap.getAlleles()[0].equals(RawHapMapFeature.NULL_ALLELE_STRING) || hapmap.getAlleles()[1].equals(RawHapMapFeature.NULL_ALLELE_STRING) ) {
                        // get the dbsnp object corresponding to this record (needed to help us distinguish between insertions and deletions)
                        VariantContext dbsnpVC = getDbsnp(hapmap.getName());
                        if ( dbsnpVC == null || dbsnpVC.isMixed() )
                            continue;

                        Map<String, Allele> alleleMap = new HashMap<String, Allele>(2);
                        alleleMap.put(RawHapMapFeature.DELETION, Allele.create(ref.getBase(), dbsnpVC.isSimpleInsertion()));
                        alleleMap.put(RawHapMapFeature.INSERTION, Allele.create((char)ref.getBase() + ((RawHapMapFeature)record).getAlleles()[1], !dbsnpVC.isSimpleInsertion()));
                        hapmap.setActualAlleles(alleleMap);

                        // also, use the correct positioning for insertions
                        hapmap.updatePosition(dbsnpVC.getStart());

                        if ( hapmap.getStart() < ref.getWindow().getStart() ) {
                            logger.warn("Hapmap record at " + ref.getLocus() + " represents an indel too large to be converted; skipping...");
                            continue;
                        }
                    }
                }

                // ok, we might actually be able to turn this record in a variant context
                VariantContext vc = VariantContextAdaptors.toVariantContext(variants.getName(), record, ref);

                if ( vc != null ) // sometimes the track has odd stuff in it that can't be converted
                    VCs.add(vc);
            }
        }

        return VCs;
    }

    private VariantContext getDbsnp(String rsID) {
        if ( dbsnpIterator == null ) {

            if ( dbsnp == null )
                throw new UserException.BadInput("No dbSNP rod was provided, but one is needed to decipher the correct indel alleles from the HapMap records");

            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                                                          getToolkit().getGenomeLocParser(),
                                                          getToolkit().getArguments().unsafe,
                                                          getToolkit().getArguments().disableAutoIndexCreationAndLockingWhenReadingRods,
                                                          null);
            dbsnpIterator = builder.createInstanceOfTrack(VCFCodec.class, new File(dbsnp.dbsnp.getSource())).getIterator();
            // Note that we should really use some sort of seekable iterator here so that the search doesn't take forever
            // (but it's complicated because the hapmap location doesn't match the dbsnp location, so we don't know where to seek to)
        }

        while ( dbsnpIterator.hasNext() ) {
            GATKFeature feature = dbsnpIterator.next();
            VariantContext vc = (VariantContext)feature.getUnderlyingObject();
            if ( vc.getID().equals(rsID) )
                return vc;
        }

        return null;
    }

    private void writeRecord(VariantContext vc, RefMetaDataTracker tracker, GenomeLoc loc) {
        if ( !wroteHeader ) {
            wroteHeader = true;

            // setup the header fields
            Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
            hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), Arrays.asList(variants.getName())));
            hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));

            allowedGenotypeFormatStrings.add(VCFConstants.GENOTYPE_KEY);
            for ( VCFHeaderLine field : hInfo ) {
                if ( field instanceof VCFFormatHeaderLine) {
                    allowedGenotypeFormatStrings.add(((VCFFormatHeaderLine)field).getID());
                }
            }

            samples = new LinkedHashSet<String>();
            if ( sampleName != null ) {
                samples.add(sampleName);
            } else {
                // try VCF first
                samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(variants.getName()));

                if ( samples.isEmpty() ) {
                    List<Feature> features = tracker.getValues(variants, loc);
                    if ( features.size() == 0 )
                        throw new IllegalStateException("No rod data is present, but we just created a VariantContext");

                    Feature f = features.get(0);
                    if ( f instanceof RawHapMapFeature )
                        samples.addAll(Arrays.asList(((RawHapMapFeature)f).getSampleIDs()));
                    else
                        samples.addAll(vc.getSampleNames());
                }
            }

            vcfwriter.writeHeader(new VCFHeader(hInfo, samples));
        }

        vc = GATKVariantContextUtils.purgeUnallowedGenotypeAttributes(vc, allowedGenotypeFormatStrings);
        vcfwriter.add(vc);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer sum) {
        if ( dbsnpIterator != null )
            dbsnpIterator.close();
        vcfwriter.close();
    }
}
