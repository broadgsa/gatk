package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Counts of bases from SLX, 454, and SOLiD at this site
 */
@Hidden
public class TechnologyComposition extends InfoFieldAnnotation implements ExperimentalAnnotation {
    private String nSLX = "NumSLX";
    private String n454 ="Num454";
    private String nSolid = "NumSOLiD";
    private String nOther = "NumOther";
    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int readsIllumina = 0;
        int readsSolid = 0;
        int reads454 = 0;
        int readsOther = 0;

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            AlignmentContext context = sample.getValue();

            ReadBackedPileup pileup = null;
            if (context.hasExtendedEventPileup())
                pileup = context.getExtendedEventPileup();
            else if (context.hasBasePileup())
                pileup = context.getBasePileup();

            if (pileup != null) {
                for (PileupElement p : pileup ) {
                    if(ReadUtils.is454Read(p.getRead()))
                        reads454++;
                    else if (ReadUtils.isSOLiDRead(p.getRead()))
                        readsSolid++;
                    else if (ReadUtils.isSLXRead(p.getRead()))
                        readsIllumina++;
                    else
                        readsOther++;
                }
            }
        }

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(nSLX, String.format("%d", readsIllumina));
        map.put(n454, String.format("%d", reads454));
        map.put(nSolid, String.format("%d", readsSolid));
        map.put(nOther, String.format("%d", readsOther));
            return map;
    }

     public List<String> getKeyNames() { return Arrays.asList(nSLX,n454,nSolid,nOther); }

     public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(nSLX, 1, VCFHeaderLineType.Integer, "Number of SLX reads"),
             new VCFInfoHeaderLine(n454, 1, VCFHeaderLineType.Integer, "Number of 454 reads"),
             new VCFInfoHeaderLine(nSolid, 1, VCFHeaderLineType.Integer, "Number of SOLiD reads"),
             new VCFInfoHeaderLine(nOther, 1, VCFHeaderLineType.Integer, "Number of Other technology reads")); }

}
