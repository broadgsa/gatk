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

import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;

/**
 * Converts variants from other file formats to VCF format.
 */
@Requires(value={},referenceMetaData=@RMD(name=VariantsToVCF.INPUT_ROD_NAME,type= ReferenceOrderedDatum.class))
@Reference(window=@Window(start=0,stop=40))
public class VariantsToVCF extends RodWalker<Integer, Integer> {

    public static final String INPUT_ROD_NAME = "variant";

    @Argument(fullName="sample", shortName="sample", doc="The sample name represented by the variant rod (for data like GELI with genotypes)", required=false)
    protected String sampleName = null;

    private VCFWriter vcfwriter = null;

    // Don't allow mixed types for now
    private EnumSet<VariantContext.Type> ALLOWED_VARIANT_CONTEXT_TYPES = EnumSet.of(VariantContext.Type.SNP, VariantContext.Type.NO_VARIATION, VariantContext.Type.INDEL);

    private String[] ALLOWED_FORMAT_FIELDS = {VCFConstants.GENOTYPE_KEY, VCFConstants.GENOTYPE_QUALITY_KEY, VCFConstants.DEPTH_KEY, VCFConstants.GENOTYPE_LIKELIHOODS_KEY };

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || !BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        DbSNPFeature dbsnp = DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));

        Collection<VariantContext> contexts = tracker.getVariantContexts(ref, INPUT_ROD_NAME, ALLOWED_VARIANT_CONTEXT_TYPES, context.getLocation(), true, false);

        for ( VariantContext vc : contexts ) {
            VCFRecord vcf = VariantContextAdaptors.toVCF(vc, ref.getBase(), Arrays.asList(ALLOWED_FORMAT_FIELDS), false, false);
            if ( dbsnp != null )
                vcf.setID(dbsnp.getRsID());
            // set the appropriate sample name if necessary
            if ( sampleName != null && vcf.hasGenotypeData() && vcf.getGenotype(INPUT_ROD_NAME) != null )
                 vcf.getGenotype(INPUT_ROD_NAME).setSampleName(sampleName);
            writeRecord(vcf, tracker);
        }

        return 1;
    }

    private void writeRecord(VCFRecord rec, RefMetaDataTracker tracker) {
        if ( vcfwriter == null ) {
            // setup the header fields
            Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
            hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
            hInfo.add(new VCFHeaderLine("source", "VariantsToVCF"));
            hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

            TreeSet<String> samples = new TreeSet<String>();
            if ( sampleName != null ) {
                samples.add(sampleName);
            } else {

                List<Object> rods = tracker.getReferenceMetaData(INPUT_ROD_NAME);
                if ( rods.size() == 0 )
                    throw new IllegalStateException("VCF record was created, but no rod data is present");

                Object rod = rods.get(0);
                if ( rod instanceof VCFRecord )
                    samples.addAll(Arrays.asList(((VCFRecord)rod).getSampleNames()));
                else if ( rod instanceof HapMapROD )
                    samples.addAll(Arrays.asList(((HapMapROD)rod).getSampleIDs()));
                else
                    samples.addAll(Arrays.asList(rec.getSampleNames()));
            }

            vcfwriter = new VCFWriter(out);
            vcfwriter.writeHeader(new VCFHeader(hInfo, samples));
        }
        vcfwriter.addRecord(rec);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer sum) {
        if ( vcfwriter != null )
            vcfwriter.close();
    }
}
