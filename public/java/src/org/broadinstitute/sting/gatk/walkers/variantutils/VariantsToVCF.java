/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.hapmap.RawHapMapFeature;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Converts variants from other file formats to VCF format.
 *
 * <p>
 * Note that there must be a Tribble feature/codec for the file format as well as an adaptor.
 *
 * <h2>Input</h2>
 * <p>
 * A variant file to filter.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A VCF file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantsToVCF \
 *   -o output.vcf \
 *   --variant:RawHapMap input.hapmap \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@Reference(window=@Window(start=-40,stop=40))
public class VariantsToVCF extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter baseWriter = null;
    private SortingVCFWriter vcfwriter; // needed because hapmap/dbsnp indel records move

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

    /**
     * This argument is useful for fixing input VCFs with bad reference bases (the output will be a fixed version of the VCF).
     */
    @Argument(fullName="fixRef", shortName="fixRef", doc="Fix common reference base in case there's an indel without padding", required=false)
    protected boolean fixReferenceBase = false;

    private Set<String> allowedGenotypeFormatStrings = new HashSet<String>();
    private boolean wroteHeader = false;

    // for dealing with indels in hapmap
    CloseableIterator<GATKFeature> dbsnpIterator = null;

    public void initialize() {
        vcfwriter = new SortingVCFWriter(baseWriter, 40, false);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || !BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        String rsID = dbsnp == null ? null : VCFUtils.rsIDOfFirstRealVariant(tracker.getValues(dbsnp.dbsnp, context.getLocation()), VariantContext.Type.SNP);

        Collection<VariantContext> contexts = getVariantContexts(tracker, ref);

        for ( VariantContext vc : contexts ) {
            Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
            if ( rsID != null && !vc.hasID() ) {
                attrs.put(VariantContext.ID_KEY, rsID);
                vc = VariantContext.modifyAttributes(vc, attrs);
            }

            // set the appropriate sample name if necessary
            if ( sampleName != null && vc.hasGenotypes() && vc.hasGenotype(variants.getName()) ) {
                Genotype g = Genotype.modifyName(vc.getGenotype(variants.getName()), sampleName);
                Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
                genotypes.put(sampleName, g);
                vc = VariantContext.modifyGenotypes(vc, genotypes);
            }

            if ( fixReferenceBase ) {
                vc = VariantContext.modifyReferencePadding(vc, ref.getBase());
            }

            writeRecord(vc, tracker, ref.getLocus());
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
                        alleleMap.put(RawHapMapFeature.DELETION, Allele.create(Allele.NULL_ALLELE_STRING, dbsnpVC.isSimpleInsertion()));
                        alleleMap.put(RawHapMapFeature.INSERTION, Allele.create(((RawHapMapFeature)record).getAlleles()[1], !dbsnpVC.isSimpleInsertion()));
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

            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),getToolkit().getGenomeLocParser(),getToolkit().getArguments().unsafe);
            dbsnpIterator = builder.createInstanceOfTrack(VCFCodec.class, new File(dbsnp.dbsnp.getSource())).getIterator();
            // Note that we should really use some sort of seekable iterator here so that the search doesn't take forever
            // (but it's complicated because the hapmap location doesn't match the dbsnp location, so we don't know where to seek to)
        }

        while ( dbsnpIterator.hasNext() ) {
            GATKFeature feature = dbsnpIterator.next();
            VariantContext vc = (VariantContext)feature.getUnderlyingObject();
            if ( vc.hasID() && vc.getID().equals(rsID) )
                return vc;
        }

        return null;
    }

    private void writeRecord(VariantContext vc, RefMetaDataTracker tracker, GenomeLoc loc) {
        if ( !wroteHeader ) {
            wroteHeader = true;

            // setup the header fields
            Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
            hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), Arrays.asList(variants.getName())));
            //hInfo.add(new VCFHeaderLine("source", "VariantsToVCF"));
            //hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

            allowedGenotypeFormatStrings.add(VCFConstants.GENOTYPE_KEY);
            for ( VCFHeaderLine field : hInfo ) {
                if ( field instanceof VCFFormatHeaderLine) {
                    allowedGenotypeFormatStrings.add(((VCFFormatHeaderLine)field).getName());
                }
            }

            Set<String> samples = new LinkedHashSet<String>();
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

        vc = VariantContextUtils.purgeUnallowedGenotypeAttributes(vc, allowedGenotypeFormatStrings);
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
