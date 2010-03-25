package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

/**
 * Converts variants from other file formats to VCF format.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant",type= ReferenceOrderedDatum.class))
public class VariantsToVCF extends RodWalker<Integer, Integer> {

    @Argument(fullName="sample", shortName="sample", doc="The sample name represented by the variant rod (for data like GELI with genotypes)", required=false)
    protected String sampleName = null;

    private VCFWriter vcfwriter = null;

    // Don't allow mixed types for now
    private EnumSet<VariantContext.Type> ALLOWED_VARIANT_CONTEXT_TYPES = EnumSet.of(VariantContext.Type.SNP, VariantContext.Type.NO_VARIATION, VariantContext.Type.INDEL);

    private String[] ALLOWED_FORMAT_FIELDS = {"GT"};

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));        

        Collection<VariantContext> contexts = tracker.getVariantContexts("variant", ALLOWED_VARIANT_CONTEXT_TYPES, context.getLocation(), true, false);
        for ( VariantContext vc : contexts ) {
            VCFRecord vcf = VariantContextAdaptors.toVCF(vc, ref.getBase(), Arrays.asList(ALLOWED_FORMAT_FIELDS), false, false);
            if ( dbsnp != null )
                vcf.setID(dbsnp.getRS_ID());
            if ( sampleName != null && vcf.hasGenotypeData() && vcf.getGenotype("variant") != null )
                 vcf.getGenotype("variant").setSampleName(sampleName);
            writeRecord(vcf);
        }

        return 1;
    }

    private void writeRecord(VCFRecord rec) {
        if ( vcfwriter == null ) {
            // setup the header fields
            Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
            hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
            hInfo.add(new VCFHeaderLine("source", "VariantFiltration"));
            hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

            vcfwriter = new VCFWriter(out);

            TreeSet<String> samples = new TreeSet<String>();
            if ( sampleName != null )
                samples.add(sampleName);
            else
                samples.addAll(Arrays.asList(rec.getSampleNames()));

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
        vcfwriter.close();
    }
}
