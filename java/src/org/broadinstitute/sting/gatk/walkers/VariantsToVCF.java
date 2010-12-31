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

package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.dbsnp.DbSNPCodec;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.hapmap.HapMapFeature;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.SampleUtils;

import java.util.*;

/**
 * Converts variants from other file formats to VCF format.
 */
@Requires(value={},referenceMetaData=@RMD(name=VariantsToVCF.INPUT_ROD_NAME, type=VariantContext.class))
@Reference(window=@Window(start=-40,stop=40))
public class VariantsToVCF extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter baseWriter = null;
    private SortingVCFWriter vcfwriter; // needed because hapmap indel records move

    public static final String INPUT_ROD_NAME = "variant";

    @Argument(fullName="sample", shortName="sample", doc="The sample name represented by the variant rod (for data like GELI with genotypes)", required=false)
    protected String sampleName = null;

    private Set<String> allowedGenotypeFormatStrings = new HashSet<String>();
    private boolean wroteHeader = false;

    // Don't allow mixed types for now
    private EnumSet<VariantContext.Type> ALLOWED_VARIANT_CONTEXT_TYPES = EnumSet.of(VariantContext.Type.SNP,
            VariantContext.Type.NO_VARIATION, VariantContext.Type.INDEL, VariantContext.Type.MNP);

    // for dealing with indels in hapmap
    CloseableIterator<GATKFeature> dbsnpIterator = null;

    public void initialize() {
        vcfwriter = new SortingVCFWriter(baseWriter, 40, false);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || !BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));

        Collection<VariantContext> contexts = getVariantContexts(tracker, ref);

        for ( VariantContext vc : contexts ) {
            Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
            if ( rsID != null && !vc.hasID() ) {
                attrs.put(VariantContext.ID_KEY, rsID);
                vc = VariantContext.modifyAttributes(vc, attrs);
            }

            // set the appropriate sample name if necessary
            if ( sampleName != null && vc.hasGenotypes() && vc.hasGenotype(INPUT_ROD_NAME) ) {
                Genotype g = Genotype.modifyName(vc.getGenotype(INPUT_ROD_NAME), sampleName);
                Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
                genotypes.put(sampleName, g);
                vc = VariantContext.modifyGenotypes(vc, genotypes);
            }

            writeRecord(vc, tracker, ref.getBase());
        }

        return 1;
    }

    private Collection<VariantContext> getVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref) {
        // we need to special case the HapMap format because indels aren't handled correctly
        List<Object> features = tracker.getReferenceMetaData(INPUT_ROD_NAME, true);
        if ( features.size() > 0 && features.get(0) instanceof HapMapFeature ) {
            ArrayList<VariantContext> hapmapVCs = new ArrayList<VariantContext>(features.size());
            for ( Object feature : features ) {
                HapMapFeature hapmap = (HapMapFeature)feature;
                Byte refBase = null;

                // if it's an indel, we need to figure out the alleles
                if ( hapmap.getAlleles()[0].equals("-") ) {
                    Map<String, Allele> alleleMap = new HashMap<String, Allele>(2);

                    // get the dbsnp object corresponding to this record, so we can learn whether this is an insertion or deletion
                    DbSNPFeature dbsnp = getDbsnpFeature(hapmap.getName());
                    if ( dbsnp == null || dbsnp.getVariantType().equalsIgnoreCase("mixed") )
                        continue;

                    boolean isInsertion = dbsnp.getVariantType().equalsIgnoreCase("insertion");

                    alleleMap.put(HapMapFeature.DELETION, Allele.create(Allele.NULL_ALLELE_STRING, isInsertion));
                    alleleMap.put(HapMapFeature.INSERTION, Allele.create(hapmap.getAlleles()[1], !isInsertion));
                    hapmap.setActualAlleles(alleleMap);

                    // also, use the correct positioning for insertions
                    if ( isInsertion )
                        hapmap.updatePosition(dbsnp.getStart());                        
                    else
                        hapmap.updatePosition(dbsnp.getStart() - 1);

                    if ( hapmap.getStart() < ref.getWindow().getStart() ) {
                        logger.warn("Hapmap record at " + ref.getLocus() + " represents an indel too large to be converted; skipping...");
                        continue;
                    }
                    refBase = ref.getBases()[hapmap.getStart() - ref.getWindow().getStart()];
                }
                VariantContext vc = VariantContextAdaptors.toVariantContext(INPUT_ROD_NAME, hapmap, ref);
                if ( vc != null ) {
                    if ( refBase != null ) {
                        Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
                        attrs.put(VariantContext.REFERENCE_BASE_FOR_INDEL_KEY, refBase);
                        vc = VariantContext.modifyAttributes(vc, attrs);
                    }
                    hapmapVCs.add(vc);
                }
            }
            return hapmapVCs;
        }

        // for everything else, we can just convert to VariantContext
        return tracker.getVariantContexts(ref, INPUT_ROD_NAME, ALLOWED_VARIANT_CONTEXT_TYPES, ref.getLocus(), true, false);
    }

    private DbSNPFeature getDbsnpFeature(String rsID) {
        if ( dbsnpIterator == null ) {
            ReferenceOrderedDataSource dbsnpDataSource = null;
            for ( ReferenceOrderedDataSource ds : getToolkit().getRodDataSources() ) {
                if ( ds.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                    dbsnpDataSource = ds;
                    break;
                }
            }

            if ( dbsnpDataSource == null )
                throw new UserException.BadInput("No dbSNP rod was provided, but one is needed to decipher the correct indel alleles from the HapMap records");

            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),getToolkit().getGenomeLocParser(),getToolkit().getArguments().unsafe);
            dbsnpIterator = builder.createInstanceOfTrack(DbSNPCodec.class, dbsnpDataSource.getFile()).getIterator();
            // Note that we should really use some sort of seekable iterator here so that the search doesn't take forever
            // (but it's complicated because the hapmap location doesn't match the dbsnp location, so we don't know where to seek to)
        }

        while ( dbsnpIterator.hasNext() ) {
            GATKFeature feature = dbsnpIterator.next();
            DbSNPFeature dbsnp = (DbSNPFeature)feature.getUnderlyingObject();
            if ( dbsnp.getRsID().equals(rsID) )
                return dbsnp;
        }

        return null;
    }

    private void writeRecord(VariantContext vc, RefMetaDataTracker tracker, byte ref) {
        if ( !wroteHeader ) {
            wroteHeader = true;

            // setup the header fields
            Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
            hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
            hInfo.add(new VCFHeaderLine("source", "VariantsToVCF"));
            hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

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
                samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(INPUT_ROD_NAME));

                if ( samples.isEmpty() ) {
                    List<Object> rods = tracker.getReferenceMetaData(INPUT_ROD_NAME);
                    if ( rods.size() == 0 )
                        throw new IllegalStateException("No rod data is present");

                    Object rod = rods.get(0);
                    if ( rod instanceof HapMapFeature)
                        samples.addAll(Arrays.asList(((HapMapFeature)rod).getSampleIDs()));
                    else
                        samples.addAll(vc.getSampleNames());
                }
            }

            vcfwriter.writeHeader(new VCFHeader(hInfo, samples));
        }

        vc = VariantContextUtils.purgeUnallowedGenotypeAttributes(vc, allowedGenotypeFormatStrings);
        vcfwriter.add(vc, ref);
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
