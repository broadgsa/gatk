package org.broadinstitute.sting.oneoffprojects.walkers.vcftools;

import com.google.common.collect.ArrayListMultimap;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import javax.management.NotificationBroadcaster;
import java.util.*;

/**
 * @doc Merges N callsets that have been made on the same set of samples, and averages specific annotations.
 * @Author chartl
 *
 */
public class BootstrapCallsMerger extends RodWalker<BootstrapCallsMerger.VCHolder,VCFWriter> implements TreeReducible<VCFWriter>{

    @Output
    VCFWriter vcfWriter = null;

    // todo -- remove this, can be done just by looking at the type and iterating over ANNOTS_TO_AVG
    // todo -- multi-allelic sites (e.g. what happens here:)
    // todo --      Set 1:    A   G,T    AC=2,4   ??
    // todo --      Set 1:    A   G      AC=2      Set 2:    A    T    AC=4
    // todo -- fix above behavior
    final static private Set<String> INTEGER_ANNOTS_CAN_MEAN = new HashSet<String>(Arrays.asList("AC","AN"));
    final static private Set<String> ANNOTS_TO_AVG =  new HashSet<String>(Arrays.asList(
            "QD","SB","HaplotypeScore","Dels","MQ","MQ0","sumGLByD","AC","AF","AN"));

    public void initialize() {
        // grab the samples
        Set<String> samples = new HashSet<String>();
        // setup the header fields
        // note that if any of the definitions conflict with our new ones, then we want to overwrite the old ones
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        for ( Map.Entry<String,VCFHeader> headers : VCFUtils.getVCFHeadersFromRodPrefix(getToolkit(), "bootstrap").entrySet() ) {
            samples.addAll(headers.getValue().getGenotypeSamples());
            for ( VCFHeaderLine line : headers.getValue().getMetaData() ) {
                logger.debug(line);
                VCFHeaderLine altered = alterHeaderLine(line);
                if ( VariantAnnotator.isUniqueHeaderLine(altered, hInfo) )
                    hInfo.add(altered);
            }
        }
        hInfo.add(new VCFInfoHeaderLine("NB",1,VCFHeaderLineType.Integer,"Number of bootsrap sets site was seen in"));
        hInfo.add(new VCFFormatHeaderLine("BC",4,VCFHeaderLineType.Integer,"Genotype counts across bootsraps: ref,het,var,nocall"));
        HashSet<String> rodName = new HashSet<String>();
        rodName.add("variant");
        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    /**
     * Note: integer annotations will need to become floats, others will not
     * (e.g. HRun won't change, but counts will)
     * @param line
     * @return line with type changed
     */
    private VCFHeaderLine alterHeaderLine(VCFHeaderLine line) {
        if ( line instanceof VCFInfoHeaderLine ) {
            if(INTEGER_ANNOTS_CAN_MEAN.contains(((VCFInfoHeaderLine) line).getName())) {
                return new VCFInfoHeaderLine(((VCFInfoHeaderLine) line).getName(),
                        ((VCFInfoHeaderLine) line).getCount(),
                        VCFHeaderLineType.Float,
                        ((VCFInfoHeaderLine) line).getDescription());
            }
        }
        return line;
    }

    public VCFWriter reduceInit() {
        return vcfWriter;
    }

    public VCHolder map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext con) {
        if ( tracker == null ) { return null; }
        Collection<VariantContext> bootstraps = tracker.getVariantContextsByPrefix(ref,Arrays.asList("bootstrap"),null,ref.getLocus(),true,false);
        int num_bootstraps = bootstraps.size();
        if ( num_bootstraps == 0 ) { return null; }
        Map<String,Double> avgInfo = new HashMap<String,Double>(ANNOTS_TO_AVG.size());
        Map<String,Integer[]> genotypeCountsBySample = new HashMap<String,Integer[]>();
        for ( VariantContext vc : bootstraps ) {
            // update genotype counts
            for ( Map.Entry<String,Genotype> genotype : vc.getGenotypes().entrySet() ) {
                if ( ! genotypeCountsBySample.containsKey(genotype.getKey())) {
                    genotypeCountsBySample.put(genotype.getKey(),new Integer[]{0,0,0,0});
                }
                genotypeCountsBySample.get(genotype.getKey())[genotype2offset(genotype.getValue())]++;
            }
            // update info field annotations
            for ( String anno : ANNOTS_TO_AVG ) {
                if ( ! avgInfo.containsKey(anno) ) {
                    avgInfo.put(anno,0.0);
                }
                Object value = vc.getAttribute(anno);
                if ( value instanceof Number ) {
                    //logger.debug(value);
                    avgInfo.put(anno,avgInfo.get(anno) + ((Number)value).doubleValue()/num_bootstraps);
                }
                if ( value instanceof String ) {
                    //logger.debug("string: "+value.toString());
                    avgInfo.put(anno,avgInfo.get(anno) + Double.valueOf((String)value)/num_bootstraps);
                }
            }
        }
        VariantContext first = bootstraps.iterator().next();
        Map<String,Object> finalInfo = new HashMap<String,Object>(first.getAttributes().size());
        for ( Map.Entry<String,Object> attrib : first.getAttributes().entrySet() ) {
            if ( ANNOTS_TO_AVG.contains(attrib.getKey()) ) {
                finalInfo.put(attrib.getKey(),avgInfo.get(attrib.getKey()));
            } else {
                finalInfo.put(attrib.getKey(),attrib.getValue());
            }
        }
        Map<String,Genotype> finalGenotypes = new HashMap<String,Genotype>(first.getSampleNames().size());
        for ( Map.Entry<String,Genotype> g : first.getGenotypes().entrySet() ) {
            Map<String,Object> att = new HashMap<String,Object>(g.getValue().getAttributes());
            att.put("BC",countsToString(genotypeCountsBySample.get(g.getKey())));
            //logger.debug(g.getValue());
            finalGenotypes.put(g.getKey(),Genotype.modifyAttributes(g.getValue(),att));
            //logger.debug("final:");
            //logger.debug(finalGenotypes.get(g.getKey()));
        }

        finalInfo.put("NB",String.format("%d",num_bootstraps));

        VariantContext attributeModified = VariantContext.modifyAttributes(first,finalInfo);
        logger.debug(attributeModified.hasGenotypes() ? "attributes have genotypes" : "VERY BAD");
        VariantContext genotypeModified = VariantContext.modifyGenotypes(attributeModified,finalGenotypes);
        logger.debug(genotypeModified.hasGenotypes() ? "modified genotypes have genotypes" : "NOT SO BAD");

        return new VCHolder(genotypeModified,ref.getBase());

        //return new VCHolder(VariantContext.modifyGenotypes(VariantContext.modifyAttributes(first, finalInfo), finalGenotypes),
               // ref.getBase());
    }

    private static String countsToString(Integer[] c) {
        return String.format("%d,%d,%d,%d",c[0],c[1],c[2],c[3]);
    }

    public VCFWriter treeReduce(VCFWriter l, VCFWriter r) {
        return l;
    }

    public VCFWriter reduce(VCHolder h, VCFWriter w) {
        if ( h != null ) {
            w.add(h.v,h.b);
        }

        return w;
    }

    private static int genotype2offset(Genotype g) {
        if ( g.isHomRef() ) { return 0; }
        if ( g.isHet() ) { return 1; }
        if ( g.isHomVar() ) { return 2; }
        return 3;
    }

    class VCHolder {
        public VariantContext v;
        public Byte b;
        public VCHolder(VariantContext v, Byte b) {
            this.v = v;
            this.b = b;
        }
    }
}
