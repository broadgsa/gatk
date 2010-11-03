package org.broadinstitute.sting.oneoffprojects.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 20, 2010
 * Time: 3:08:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class HammingDistance implements ExperimentalAnnotation, InfoFieldAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( tracker == null ) {
            return null;
        }
        VariantContext hamCon = tracker.getVariantContext(ref,"hamming",null,ref.getLocus(),true);
        if ( hamCon == null ) {
            return null;
        }

        Set<String> interSamples =  new HashSet<String>(vc.getSampleNames());
        interSamples.retainAll(hamCon.getSampleNames());

        int dist  = 0;
        int nrd_num = 0;
        int hamCalls = 0;
        int vcCallsAtHamCalls = 0;
        int num_variant = 0;

        for ( String s : interSamples ) {
            dist += dist(vc.getGenotype(s),hamCon.getGenotype(s),true);
            nrd_num += dist(vc.getGenotype(s),hamCon.getGenotype(s),false);
            if ( vc.getGenotype(s).isHet() || vc.getGenotype(s).isHomVar() || hamCon.getGenotype(s).isHet() || hamCon.getGenotype(s).isHomVar() ) {
                num_variant ++;
            }
            if ( hamCon.getGenotype(s).isCalled() ) {
                hamCalls++;
                if ( vc.getGenotype(s).isCalled() ) {
                    vcCallsAtHamCalls++;
                }
            }

        }

        HashMap<String,Object> map = new HashMap<String,Object>(1);
        map.put("HMD",dist);
        map.put("HCR",(0.0+vcCallsAtHamCalls)/(0.0+hamCalls));
        map.put("NRD",(0.0+nrd_num)/(0.0+num_variant));
        map.put("OGC",(0.0+nrd_num)/(0.0+interSamples.size()));
        return map;
        
    }

    public int dist(Genotype a, Genotype b, boolean weightByAC) {
        if ( a.isNoCall() || b.isNoCall() ) {
            return 0;
        }
        if ( weightByAC ) {
            if ( a.isHomRef() ) {
                if ( b.isHomVar() ) {
                    return 2;
                } else if ( b.isHet() ) {
                    return 1;
                } else {
                    return 0;
                }
            } else if ( a.isHet() ) {
                if ( b.isHom() ) {
                    return 1;
                } else {
                    return 0;
                }
            } else {
                if ( b.isHomRef() ) {
                    return 2;
                } else if ( b.isHet() ) {
                    return 1;
                } else {
                    return 0;
                }
            }
        } else {
            if ( ! a.equals(b) ) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    public List<String> getKeyNames() { return Arrays.asList("HMD","HCR","NRD","OGC"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("HMD",1, VCFHeaderLineType.Integer,"The hamming distance between record in Hamming ROD and this record"),
            new VCFInfoHeaderLine("HCR",1,VCFHeaderLineType.Float,"The differential call rate between record in Hamming ROD and this record"),
            new VCFInfoHeaderLine("NRD",1,VCFHeaderLineType.Float,"The Non-reference discrepancy between Hamming ROD and this record"),
            new VCFInfoHeaderLine("OGC",1,VCFHeaderLineType.Float,"The Overall Genotype Concordance between Hamming ROD and this one")); }


}