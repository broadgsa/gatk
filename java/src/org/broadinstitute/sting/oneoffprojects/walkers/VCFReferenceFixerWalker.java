package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.util.*;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Apr 13, 2010
 */
public class VCFReferenceFixerWalker extends RodWalker<VCFRecord,Long> {

    private VCFWriter vcfWriter;

    public void initialize() {
        TreeSet<String> samples = new TreeSet<String>();
        SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, new HashMap<Pair<String, String>, String>());
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantAnnotator"));
        vcfWriter = new VCFWriter(out);
        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext context, AlignmentContext alicon) {
        if ( tracker == null ) {
            return null;
        }
        List<Object> rods = tracker.getReferenceMetaData("fixme");
        Object rod = rods.get(0);
        RodVCF vcfrod = null;
        if ( rod instanceof RodVCF ) {
            vcfrod = (RodVCF) rod;
        }

        VCFRecord rec = vcfrod.getRecord();
        rec.setReferenceBase(new String(BaseUtils.charSeq2byteSeq(context.getBases())));
        return rec;

        /*
        VariantContext vcon = null;
        if ( rod instanceof RodVCF) {
            vcon = VariantContextAdaptors.toVariantContext("fixme", (RodVCF) rod, new Allele(BaseUtils.charSeq2byteSeq(context.getBases()),true));
        }
        if ( vcon == null ) {
            return null;
        }
        Set<Allele> otherAlleles = vcon.getAlternateAlleles();
        VariantContext fixedContext = new VariantContext(vcon.getName(),context.getLocus(),otherAlleles,vcon.getGenotypes(),vcon.getNegLog10PError(),vcon.getFilters(),vcon.getAttributes());
        return VariantContextAdaptors.toVCF(fixedContext,context.getBase());*/
    }

    public Long reduce(VCFRecord con, Long num) {
        if ( con == null ) {
            return num;
        }

        vcfWriter.addRecord(con);
        return 1 + num;
    }

    public Long reduceInit() {
        return 0l;
    }

    public void onTraversalDone(Long num){
        return;
    }
}
