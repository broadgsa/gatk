package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.varianteval.MendelianViolationEvaluator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;

import java.util.*;


/**
 * Takes in a VCF file for a trio (and optionally, bam files for some or all members) and classifies mendelian violations
 * as deNovo events, or opposite homozygotes, and includes additional information pertaining to event QC (size of
 * homozygous region in child and/or parents, whether a site is likely tri-allelic or not, etc)
 *
 * @Author chartl
 * @Date Jun 8, 2010
 */
public class MendelianViolationClassifier extends LocusWalker<MendelianViolationClassifier.MendelianViolation, VCFWriter> {
    @Argument(shortName="f",fullName="familyPattern",required=true,doc="Pattern for the family structure (usage: mom+dad=child)")
    String familyStr = null;

    /*
     ****** PRIVATE CLASSES
     */

    public class ExtendedTrioStructure extends MendelianViolationEvaluator.TrioStructure {
        public HashMap<String,HomozygosityRegion> homozygousRegions;
        public HashMap<String,Integer> homozygousRegionCounts;

        public ExtendedTrioStructure(String family) {
            MendelianViolationEvaluator.TrioStructure struct = MendelianViolationEvaluator.parseTrioDescription(family);
            this.child = struct.child;
            this.mom = struct.mom;
            this.dad = struct.dad;
            homozygousRegions = new HashMap<String,HomozygosityRegion>(3);
            homozygousRegionCounts = new HashMap<String,Integer>(3);
            homozygousRegions.put(child,null);
            homozygousRegions.put(mom,null);
            homozygousRegions.put(dad,null);
            homozygousRegionCounts.put(child,0);
            homozygousRegionCounts.put(mom,0);
            homozygousRegionCounts.put(dad,0);
        }
    }

    public class HomozygosityRegion {
        public GenomeLoc regionStart;
        public int callsWithinRegion;

        public HomozygosityRegion(GenomeLoc start) {
            regionStart = start;
            callsWithinRegion = 0;
        }
    }

    public class MendelianViolation {
        private VariantContext trio;
        public MendelianViolationType type;
        private HashMap<String,Object> newAttributes;

        public MendelianViolation(VariantContext context, MendelianViolationType violationType) {
            trio = context;
            type = violationType;
            newAttributes = new HashMap<String,Object>();
            newAttributes.put(MendelianInfoKey.ViolationType.getKey(),type);
        }

        public void addAttribute(String key, Object val) {
            newAttributes.put(key,val);
        }

        public VariantContext toVariantContext() {
            newAttributes.putAll(trio.getAttributes());
            return new VariantContext(trio.getName(),trio.getLocation(),trio.getAlleles(),trio.getGenotypes(),trio.getNegLog10PError(),trio.getFilters(),newAttributes);
        }

        public boolean siteIsFiltered() {
            return trio.isFiltered();
        }

    }

    /*
     ********** PRIVATE ENUMS
     */

    public enum MendelianViolationType {
        OPPOSITE_HOMOZYGOTE("oppositeHomozygote"),
        DE_NOVO_SNP("deNovoSNP"),
        NONE("none");

        private String infoString;

        MendelianViolationType(String typeName) {
            infoString=typeName;
        }

        public String toString() {
            return infoString;
        }
    }

    public enum MendelianInfoKey {
        ViolationType("MVT","String","\"The Mendelian violation type\""),
        ParentalDeletion("deletedParent","String","\"The parent from whom the child (putatively) inherited a deletion at opposite homozygous sites\""),
        ChildHomozygosityRegion("CHR","Integer","\"An identifier for the region of homozygosity in the child in which the opposite homozygote is located\""),
        ChildHomozygosityRegionSize("CHRS","Integer","\"The size of the region of homozygosity in the child in which the opposite homozygote is located\""),
        OppositeHomozygotesInRegion("CHROH","Integer","\"The number of opposite-homozygotes located in the region of homozygosity\""),
        TriAllelicQuality("TriAlQ","Integer","\"The variant quality of the third allele at putative tri-allelic sites\""),
        TriAllelicBase("TriAlB","String","\"The third allele at putative tri-allelic sites\""),
        NumCallsInRegion("REGION_NCALLS","Integer","\"The number of SNP calls found in the homozygosity region\""),
        ParentHomozygosityRegion("PHR","String","\"An identifier for the parent, and the region number, for the homozygosity region where the deNovo SNP is located\""),
        ParentHomozygosityRegionSize("PHRS","Integer","\"The size of the parental homozygosity region where the deNovo SNP is located\""),
        DeNovoSNPsInRegion("PHRDN","Integer","\"The number of deNovo SNP events located in the same region of parental homozygosity where the deNovo SNP is located\""),
        ProportionOfParentAllele("PropParent","Double","\"The proportion of bases in the child that were the parent allele at deNovo SNP sites\"");


        String keyName;
        String valueType;
        String valueDescription;

        MendelianInfoKey(String keyStr,String infoType,String description) {
            keyName = keyStr;
            valueType = infoType;
            valueDescription = description;
        }

        public String toString() {
            return String.format("%s,%s,%s",keyName,valueType,valueDescription);
        }

        public String getKey() {
            return keyName;
        }
    }

    /*
     *********** PRIVATE DATA
     */
    private ExtendedTrioStructure trioStructure;
    private UnifiedGenotyperEngine engine;
    private ArrayList<VariantContext> contextBuffer;

    /*
     *********** INITIALIZE
     */
    public void initialize() {
        trioStructure = new ExtendedTrioStructure(familyStr);
        contextBuffer = new ArrayList<VariantContext>(5000);
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.MIN_BASE_QUALTY_SCORE = 10;
        uac.MIN_MAPPING_QUALTY_SCORE = 10;
        uac.STANDARD_CONFIDENCE_FOR_CALLING = 30;
        engine = new UnifiedGenotyperEngine(getToolkit(),uac);
    }

    /*
     *********** REDUCE INIT
     */
    public VCFWriter reduceInit() {
        VCFWriter writer = new VCFWriter(out);
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "OppositeHomozygoteClassifier"));
        for ( MendelianInfoKey key : EnumSet.allOf(MendelianInfoKey.class) ) {
            hInfo.add( new VCFHeaderLine("INFO",key.toString()));
        }
        VCFHeader vcfHeader = new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit()));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /*
     *********** MAP
     */
    public MendelianViolation map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        return assessViolation(tracker.getVariantContext(ref,"trio", EnumSet.of(VariantContext.Type.SNP),ref.getLocus(),true),tracker,ref,context);
    }

    private MendelianViolation assessViolation(VariantContext varContext, RefMetaDataTracker tracker, ReferenceContext reference, AlignmentContext context) {
        MendelianViolation violation;
        if ( MendelianViolationEvaluator.isViolation(varContext,trioStructure) ) {
            if ( isDeNovo(varContext) ) {
                violation = assessDeNovo(varContext,tracker,reference,context);
            } else if ( isOppositeHomozygote(varContext) ) {
                violation = assessOppositeHomozygote(varContext,tracker,reference,context);
            } else {
                throw new StingException("Mendelian violation that is neither deNovo nor opposite homozygote. Should never see this.");
            }
        } else {
            violation = new MendelianViolation(varContext,MendelianViolationType.NONE);
        }

        return violation;
    }

    // Note, mendelian violation is guaranteed at this point; deNovo only happens at child het sites
    private boolean isDeNovo(VariantContext trio) {
        return trio.getGenotype(trioStructure.child).isHet();
    }

    // Note, mendelian violation is guaranteed at this point; technically we do not have to check this...but we are
    private boolean isOppositeHomozygote(VariantContext trio) {
        if ( trio.getGenotype(trioStructure.child).isHet() ) { // not valid at child het sites
            return false;
        } else if ( trio.getHetCount() > 1 ) { // child is not het, so if this is 2, mom and dad are both het, invalid
            return false;
        } else if ( trio.getGenotype(trioStructure.dad) == null || trio.getGenotype(trioStructure.mom) == null ) {
            return false;
        }

        return true;
    }

    private MendelianViolation assessDeNovo(VariantContext trio, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // todo -- implement me
        return new MendelianViolation(trio,MendelianViolationType.DE_NOVO_SNP);
    }

    private MendelianViolation assessOppositeHomozygote(VariantContext trio, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // todo -- implement me
        return new MendelianViolation(trio,MendelianViolationType.OPPOSITE_HOMOZYGOTE);
    }

    /*
     *********** REDUCE
     */
    public VCFWriter reduce(MendelianViolation variant, VCFWriter writer) {
        if ( variant != null ) {
            if ( ! variant.siteIsFiltered() ) {
                // todo -- fill me in
            } else { /* just add filtered variants to output if they're wanted */
                contextBuffer.add(variant.toVariantContext());
            }
        }

        return writer;
    }

    /*
     ********** ON TRAVERSAL DONE
     */
    public void onTraversalDone(VCFWriter writer) {
        while ( contextBuffer.size() > 0 ) {
            VariantContext vc = contextBuffer.remove(0);
            writer.addRecord(VariantContextAdaptors.toVCF(vc,vc.getReference().getBases()[0]));
        }
    }
}
