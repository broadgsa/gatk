package org.broadinstitute.sting.oneoffprojects.walkers.vcftools;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 16, 2010
 * Time: 8:20:41 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type= VariantContext.class))
public class FixRefBases extends RodWalker<Integer,Integer> {
    @Output(doc="output file to write to",required=true)
    VCFWriter out;
    @Hidden
    @Argument(doc="",fullName="bypassException",required=false)
    boolean bypass = false;

    public void initialize() {
        if ( ! bypass ) {
            throw new ReviewedStingException("This walker is currently broken. This exception is being thrown because you will not get out what you think you will.");
        }
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList("variant"));
	Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        headerLines.add(new VCFHeaderLine("source", "FixRefBases"));
        headerLines.add(new VCFInfoHeaderLine("FRI",1,VCFHeaderLineType.String,"The Fix-Ref-Info (if present, either \"Flipped\" or \"Fixed\")"));
        out.writeHeader(new VCFHeader(headerLines,vcfSamples));
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker != null && tracker.hasROD("variant") ) {
            VariantContext vc = null;
	    try {
		vc = tracker.getVariantContext(ref,"variant",null,context.getLocation(),true);
	    } catch ( ReviewedStingException e ) {
		logger.warn("Multiple variants found, catching exception ",e);
		return 0;
	    }
            VariantContext newContext = null;
            if ( vc.isSNP() && ref.getBase() != vc.getReference().getBases()[0] && vc.getReference().length() == 1 ) {
                if ( basesAreFlipped(vc,ref) ) {
                    logger.warn(String.format("Variant context at %s has ref and alt bases flipped according to reference",ref.getLocus().toString()));
                    newContext = flipBases(vc,ref);
                } else {
                    HashSet<Allele> newAlleles = new HashSet<Allele>(vc.getAlternateAlleles());
                    Allele refAllele = Allele.create(ref.getBase(),true);
                    newAlleles.add(refAllele);
                    HashMap<String,Object> newAttributes = new HashMap<String,Object>(vc.getAttributes());
                    newAttributes.put("FRI",String.format("Fixed%s-%s",vc.getReference().toString(),refAllele));
                    newContext = new VariantContext("FixRefBasesVC", ref.getLocus().getContig(),
                            ref.getLocus().getStart(), ref.getLocus().getStop(), newAlleles, fixGenotypes(vc.getGenotypes(),refAllele),
                            vc.hasNegLog10PError() ? 10.0*vc.getNegLog10PError() : VCFConstants.MISSING_QUALITY_v3_DOUBLE,
                            vc.isFiltered() ? vc.getFilters() : null, newAttributes);
                }

                if ( ! newContext.hasAttribute("FRI") ) {
                    throw new StingException("FRI for fixed base not propagated. vc="+vc.toString());
                }
                out.add(newContext,ref.getBase());
                return 1;

            } else {
                out.add(vc,ref.getBase());
            }
        }

        return 0;
    }

    public Integer reduce(Integer map, Integer reduce) {
        return map + reduce;
    }

    public void onTraversalDone(Integer fReduce) {
        logger.info(String.format("Fixed %d records",fReduce));
    }

    private boolean basesAreFlipped(VariantContext vc, ReferenceContext reference) {
        for ( Allele a : vc.getAlternateAlleles() ) {
            if ( a.getBases().length == 1 && a.getBases()[0] == reference.getBase() ) {
                return true;
            }
        }
        return false;
    }

    private VariantContext flipBases(VariantContext vc, ReferenceContext ref) {
        logger.info(String.format("Flipping bases at variant position %s:%d",vc.getChr(),vc.getStart()));
        HashSet<Allele> newAlleles = new HashSet<Allele>(vc.getAlleles().size());
        newAlleles.add(Allele.create(ref.getBase(),true));
        newAlleles.add(Allele.create(vc.getReference().getBases()[0],false));
        Map<String,Object> attribs = new HashMap<String,Object>(vc.getAttributes());
        attribs.put("FRI",String.format("Flipped%s-%s",vc.getReference().toString(),Allele.create(ref.getBase(),true).toString()));
        for ( Allele a : vc.getAlternateAlleles() ) {
            if ( a.getBases()[0] != ref.getBase() ) {
                newAlleles.add(a);
            }
        }

        VariantContext newVC = new VariantContext("FixRefBasesVC", ref.getLocus().getContig(),
                ref.getLocus().getStart(), ref.getLocus().getStop(), newAlleles, flipGenotypes(vc.getGenotypes(),newAlleles),
                vc.hasNegLog10PError() ? 10.0*vc.getNegLog10PError() : VCFConstants.MISSING_QUALITY_v3_DOUBLE,
                vc.isFiltered() ? vc.getFilters() : null, attribs);

        VariantContextUtils.calculateChromosomeCounts(newVC,attribs,false);
        VariantContext.modifyAttributes(newVC,attribs);
        return newVC;
    }

    private Map<String, Genotype> fixGenotypes(Map<String,Genotype> old, Allele newRef) {
        HashMap<String,Genotype> newGTs = new HashMap<String,Genotype>(old.size());
        for ( Map.Entry<String,Genotype> e : old.entrySet() ) {
            newGTs.put(e.getKey(),fixGenotype(e.getValue(),newRef));
        }

        return newGTs;
    }

    private Genotype fixGenotype(Genotype g, Allele newRef) {
        List<Allele> newAlleles = new ArrayList<Allele>(g.getAlleles().size());
        for ( Allele a : g.getAlleles() ) {
            if ( a.isReference() ) {
                newAlleles.add(newRef);
            } else {
                newAlleles.add(a);
            }
        }
        return new Genotype(g.getSampleName(),
                newAlleles,g.hasNegLog10PError() ? g.getNegLog10PError() : VCFConstants.MISSING_QUALITY_v3_DOUBLE,
                g.getAttributes().keySet(), g.getAttributes(),g.genotypesArePhased());
    }

    private Map<String,Genotype> flipGenotypes(Map<String,Genotype> old, Set<Allele> newAlleles) {
        HashMap<String,Genotype> newGTs = new HashMap<String,Genotype>(old.size());
        for ( Map.Entry<String,Genotype> e : old.entrySet() ) {
            newGTs.put(e.getKey(),flipGenotype(e.getValue(),newAlleles));
        }

        return newGTs;
    }

    private Genotype flipGenotype(Genotype old, Set<Allele> newAlleles) {
        Allele ref = null;
        for ( Allele a : newAlleles ) {
            if ( a.isReference() ) {
                ref = a;
            }
        }
        if ( ref == null ) {
            throw new StingException("No reference allele in variant context with which to flip genotype alleles");
        }

        List<Allele> newGTAlleles = new ArrayList<Allele>(old.getAlleles().size());

        for ( Allele a : old.getAlleles() ) {
            if ( ! a.isNoCall() && a.getBases()[0] == ref.getBases()[0] ) {
                newGTAlleles.add(ref);
            } else {
                if ( a.isReference() ) {
                    newGTAlleles.add(Allele.create(a.getBases(),false));
                } else {
                    newGTAlleles.add(a);
                }
            }
        }

        return new Genotype(old.getSampleName(),
                newGTAlleles,old.hasNegLog10PError() ? old.getNegLog10PError() : VCFConstants.MISSING_QUALITY_v3_DOUBLE,
                old.getAttributes().keySet(), old.getAttributes(),old.genotypesArePhased());
    }

}
