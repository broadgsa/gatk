package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.StandardVCFWriter;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.varianteval.MendelianViolationEvaluator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.io.PrintStream;
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
    @Output
    PrintStream out;
    @Argument(shortName="f",fullName="familyPattern",required=true,doc="Pattern for the family structure (usage: mom+dad=child)")
    String familyStr = null;
    @Argument(shortName="ob",fullName="outputBed",required=true,doc="Output file to write the homozygous region information to")
    PrintStream bedOutput = null;
    @Argument(fullName="deNovoTriAllelicQ",required=false,doc="Cutoff for quality scores of 3rd allele at denovo mendelian violations to remove it from the denovo set; e.g. Q40 = need Q40 evidence for 3rd allele to toss out deNovo")
    int deNovoTriQ = 20;
    @Argument(fullName="deNovoParentalAllele",required=false,doc="Range for the parental allele at denovo sites to be kept in denovo set, e.g. 0.4-0.6 will toss out denovo sites with parental allele proportions of <0.4 and >0.6")
    String deNovoParentalAllele = "-0.1-1.1";
    @Argument(fullName="oppositeHomozygoteTriAllelicQ",required=false,doc="Cutoff for quality scores of 3rd allele at opposite homozygote sites to remove it from the violation set")
    int opHomTriQ = 20;
    @Argument(fullName="oppositeHomozygoteAlleleProportion",required=false,doc="Range for the parental allele in the parents at opposite homozygote sites for it to be kept in violation set")
    String opHomAlleleProp = "-0.1-1.1";


    /*
     *********** PRIVATE CLASSES
     */

    public class ExtendedTrioStructure extends MendelianViolationEvaluator.TrioStructure {
        public HashMap<String,HomozygosityRegion> homozygousRegions;
        public HashMap<String,Integer> homozygousRegionCounts;
        public HashMap<String,MendelianInfoKey> regionKeys;

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
            regionKeys = new HashMap<String,MendelianInfoKey>(3);
            regionKeys.put(child,MendelianInfoKey.ChildHomozygosityRegion);
            regionKeys.put(mom,MendelianInfoKey.MotherHomozygosityRegion);
            regionKeys.put(dad,MendelianInfoKey.FatherHomozygosityRegion);
        }

        public void updateHomozygosityRegions(MendelianViolation v, PrintStream output) {
            if ( ! v.siteIsFiltered() ) {
                ArrayList<HomozygosityRegion> brokenRegions = new ArrayList<HomozygosityRegion>(3);
                // can only enter or break regions at unfiltered calls
                for( Map.Entry<String, Genotype> memberGenotype : v.getUnderlyingGenotypes().entrySet() ) {
                    // for each family member
                    if ( homozygousRegions.get(memberGenotype.getKey()) == null ) {
                        // currently in a heterozygous region, update if possible
                        if ( memberGenotype.getValue().isHom() ) {
                            homozygousRegionCounts.put(memberGenotype.getKey(),homozygousRegionCounts.get(memberGenotype.getKey())+1);
                            homozygousRegions.put(memberGenotype.getKey(),new HomozygosityRegion(v.getLocus()));
                            if ( v.type != MendelianViolationType.NONE ) {
                                v.addAttribute(regionKeys.get(memberGenotype.getKey()).getKey(),homozygousRegionCounts.get(memberGenotype.getKey()));
                            }
                        }
                    } else {
                        // potentially breaking a homozygous region
                        if ( memberGenotype.getValue().isHom() ) {
                            // no break, update the region
                            HomozygosityRegion r = homozygousRegions.get(memberGenotype.getKey());
                            r.lastSeen = v.getLocus();
                            r.callsWithinRegion++;
                            if ( v.type != MendelianViolationType.NONE && ! v.violationIsFiltered() ) {
                                v.addAttribute(regionKeys.get(memberGenotype.getKey()).getKey(),homozygousRegionCounts.get(memberGenotype.getKey()));
                                if ( v.type == MendelianViolationType.DE_NOVO_SNP ) {
                                    r.deNovoSNPsInRegion++;
                                } else if ( v.type == MendelianViolationType.OPPOSITE_HOMOZYGOTE ) {
                                    r.oppositeHomsInRegion++;
                                }
                            }
                        } else if ( memberGenotype.getValue().isHet() ) {
                            // explicitly check for hets -- no calls are not counted -- this breaks a region so we print it
                            homozygousRegions.get(memberGenotype.getKey()).finalize(v.getLocus(),memberGenotype.getKey(),homozygousRegionCounts.get(memberGenotype.getKey()));
                            brokenRegions.add(homozygousRegions.get(memberGenotype.getKey()));
                            homozygousRegions.put(memberGenotype.getKey(),null);
                        }
                    }
                }

                if ( brokenRegions.size() > 0 ) {
                    Collections.sort(brokenRegions);
                }

                for ( HomozygosityRegion r : brokenRegions ) {
                    output.printf("%s%n",r);
                }
            }
        }
    }

    public class HomozygosityRegion implements Comparable{
        public GenomeLoc regionStart;
        public GenomeLoc lastSeen;
        public GenomeLoc endedBy;
        public int callsWithinRegion;
        public int oppositeHomsInRegion;
        public int deNovoSNPsInRegion;
        private String parent;
        private int id;

        public HomozygosityRegion(GenomeLoc start) {
            regionStart = start;
            lastSeen = start;
            endedBy = null;
            callsWithinRegion = 0;
            oppositeHomsInRegion = 0;
            deNovoSNPsInRegion = 0;
        }

        public void finalize(GenomeLoc regionEnd,String parent, int id) {
            endedBy = regionEnd;
            this.parent = parent;
            this.id = id;
        }

        private int getSizeLowerBound() {
            return lastSeen.distance(regionStart);
        }

        private int getSizeOfFirstHomToFirstHet() {
            return endedBy.distance(regionStart);
        }

        public String toString() {
            return String.format("%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d",regionStart.getContig(),regionStart.getStart(),
                    lastSeen.getStart(),endedBy.getStart(),parent,id,getSizeLowerBound(),getSizeOfFirstHomToFirstHet(),
                    callsWithinRegion,oppositeHomsInRegion,deNovoSNPsInRegion);
        }

        public int compareTo(Object o) {
            if ( ! ( o instanceof HomozygosityRegion) ) {
                return Integer.MIN_VALUE;
            }

            return this.regionStart.compareTo(((HomozygosityRegion) o).regionStart);
        }

        public String getContigStr() {
            return regionStart.getContig();
        }
    }

    public class MendelianViolation {
        private VariantContext trio;
        public MendelianViolationType type;
        private HashMap<String,Object> newAttributes;
        private HashMap<String,Integer> homozygosityRegions;
        private boolean filtered = false;

        public MendelianViolation(VariantContext context, MendelianViolationType violationType) {
            trio = context;
            type = violationType;
            newAttributes = new HashMap<String,Object>();
            newAttributes.put(MendelianInfoKey.ViolationType.getKey(),type);
            homozygosityRegions = new HashMap<String,Integer>(3);
        }

        public void addAttribute(String key, Object val) {
            newAttributes.put(key,val);
        }

        public VariantContext toVariantContext() {
            newAttributes.putAll(trio.getAttributes());
            return new VariantContext(trio.getName(), trio.getChr(), trio.getStart(), trio.getEnd(),trio.getAlleles(),trio.getGenotypes(),trio.getNegLog10PError(),trio.filtersWereApplied()?trio.getFilters():null,newAttributes);
        }

        public boolean siteIsFiltered() {
            return trio.isFiltered();
        }

        public void addRegions(HashMap<String,Integer> regionIDsByName ) {
            for ( Map.Entry<String,Integer> e : regionIDsByName.entrySet() ) {
                setRegion(e.getKey(),e.getValue());
            }
        }

        public void setRegion(String parent, int regionID) {
            homozygosityRegions.put(parent,regionID);
        }

        public boolean isInPreviousRegions(Map<String,Integer> otherRegions) {
            for ( String s : otherRegions.keySet() ) {
                if ( homozygosityRegions.get(s) >= otherRegions.get(s) ) {
                    return false;
                }
            }

            return true;
        }

        public Map<String,Genotype> getUnderlyingGenotypes() {
            return trio.getGenotypes();
        }

        public GenomeLoc getLocus() {
            return VariantContextUtils.getLocation(trio);
        }

        public byte getRefBase() {
            return trio.getReference().getBases()[0];
        }

        public Object getAttribute(String key) {
            if ( newAttributes.keySet().contains(key) ) {
                return newAttributes.get(key);
            } else {
                return trio.getAttribute(key);
            }
        }

        public void filter() {
            filtered = true;
            newAttributes.put(MendelianInfoKey.ViolationType.getKey(),"Filtered_"+newAttributes.get(MendelianInfoKey.ViolationType.getKey()));
        }

        public String getType() {
            return filtered ? "Filtered_"+type.toString() : type.toString();
        }

        public boolean violationIsFiltered() {
            return filtered;
        }
    }

    public class Range {
        private double upper;
        private double lower;
        private double epsilon = 10e-3;

        Range(String rangeStr) {
            String rs = rangeStr.substring(0); // don't clobber original string
            boolean startIsNegative = rangeStr.startsWith("-");
            if ( startIsNegative ) {
                rs = rs.substring(1);
            }
            String[] lu = rs.split("-");
            lower = startIsNegative ? -1*Double.parseDouble(lu[0]) : Double.parseDouble(lu[0]);
            upper = Double.parseDouble(lu[1]);
            //System.out.printf("Lower: %.2f, Upper: %.2f",lower,upper);
        }

        public boolean contains(double p) {
            return p > lower-epsilon && p < upper+epsilon;
        }
    }

    /*
     *************** PRIVATE ENUMS
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
        ViolationType("MVT","String",1,"\"The Mendelian violation type\""),
        ParentalDeletion("deletedParent","String",1,"\"The parent from whom the child (putatively) inherited a deletion at opposite homozygous sites\""),
        TriAllelicQuality("TriAlQ","Integer",1,"\"The variant quality of the third allele at putative tri-allelic sites\""),
        TriAllelicBase("TriAlB","String",1,"\"The third allele at putative tri-allelic sites\""),
        MotherHomozygosityRegion("MHR","String",1,"\"An identifier for the mother's homozygosity region where the violation is located\""),
        FatherHomozygosityRegion("FHR","String",1,"\"An identifier for the father's homozygosity region where the violation is located\""),
        ChildHomozygosityRegion("CHR","Integer",1,"\"An identifier for the child's homozygosity region where the violation is located\""),
        ProportionOfParentAllele("PropParent","Float",1,"\"The proportion of bases in the child that were the parent allele at deNovo SNP sites\""),
        /* ***************************************** UNUSED ************************************************ */
        NumCallsInRegion("REGION_NCALLS","Integer",1,"\"The number of unfiltered SNP calls found in the homozygosity region\""),
        ChildHomozygosityRegionSize("CHRS","Integer",1,"\"The size of the region of homozygosity in the child in which the opposite homozygote is located\""),
        OppositeHomozygotesInRegion("CHROH","Integer",1,"\"The number of opposite-homozygotes located in the region of homozygosity\""),
        ParentHomozygosityRegionSize("PHRS","Integer",1,"\"The size of the parental homozygosity region where the deNovo SNP is located\""),
        DeNovoSNPsInRegion("PHRDN","Integer",1,"\"The number of deNovo SNP events located in the same region of parental homozygosity where the deNovo SNP is located\"");


        String keyName;
        String valueType;
        String valueDescription;
        int numFields;

        MendelianInfoKey(String keyStr,String infoType, int fields, String description) {
            keyName = keyStr;
            valueType = infoType;
            valueDescription = description;
            numFields = fields;
        }

        public String toString() {
            return String.format("%s,%s,%d,%s",keyName,valueType,numFields,valueDescription);
        }

        public String getKey() {
            return keyName;
        }
    }

    /*
     **************** PRIVATE DATA
     */
    private ExtendedTrioStructure trioStructure;
    private UnifiedGenotyperEngine engine;
    private Range deNovoRange;
    private Range opHomRange;

    /*
     ***************** INITIALIZE
     */
    public void initialize() {
        trioStructure = new ExtendedTrioStructure(familyStr);
        deNovoRange = new Range(deNovoParentalAllele);
        opHomRange = new Range(opHomAlleleProp);
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.MIN_BASE_QUALTY_SCORE = 10;
        uac.MIN_MAPPING_QUALTY_SCORE = 10;
        uac.STANDARD_CONFIDENCE_FOR_CALLING = Math.min(deNovoTriQ,opHomTriQ);
        engine = new UnifiedGenotyperEngine(getToolkit(),uac);
        logger.info("Mom: "+trioStructure.mom+" Dad: "+trioStructure.dad+" Child: "+trioStructure.child);
        bedOutput.printf("%s%n",getBedFileHeader());
    }

    /*
     *********** REDUCE INIT
     */
    public VCFWriter reduceInit() {
        VCFWriter writer = new StandardVCFWriter(out);
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "MendelianViolationClassifier"));
        for ( MendelianInfoKey key : EnumSet.allOf(MendelianInfoKey.class) ) {
            hInfo.add( new VCFHeaderLine("INFO",key.toString()));
        }
        VCFHeader vcfHeader = new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit()));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /*
     ***************** FILTER
     */
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return tracker != null;
    }

    /*
     *************** MAP
     */
    public MendelianViolation map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return assessViolation(tracker.getVariantContext(ref,"trio", EnumSet.of(VariantContext.Type.SNP),ref.getLocus(),true),tracker,ref,context);
    }

    private boolean isComplete(VariantContext vc) {
        for ( Genotype g : vc.getGenotypes().values() ) {
            if ( g.isNoCall() || g.isFiltered() ) {
                return false;
            }
        }

        return true;
    }

    private MendelianViolation assessViolation(VariantContext varContext, RefMetaDataTracker tracker, ReferenceContext reference, AlignmentContext context) {
        MendelianViolation violation;
        if ( varContext != null ) {
            if ( isComplete(varContext) && MendelianViolationEvaluator.isViolation(varContext,trioStructure) ) {
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

        } else {
            violation = null;
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

    /*
     ************ ASSESS DE NOVO AND OPPOSITE HOMOZYGOTES
     */

    private MendelianViolation assessDeNovo(VariantContext trio, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        MendelianViolation violation = new MendelianViolation(trio,MendelianViolationType.DE_NOVO_SNP);

        // look for mis-genotyped sites by examining the proportion of the parent-allele bases in the child --
        // as a spoeedup we do this only at non-filtered sites

        if ( ! trio.isFiltered() ) {
            Allele parental = trio.getGenotype(trioStructure.mom).getAllele(0); // guaranteed homozygous
            if ( parental.getBases().length < 1 ) {
                throw new StingException("Parental bases have length zero at "+trio.toString());
            }

            Map<String,StratifiedAlignmentContext> splitContext = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
            Double proportion = getAlleleProportion(parental,splitContext.get(trioStructure.child));
            if ( proportion != null ) {
                violation.addAttribute(MendelianInfoKey.ProportionOfParentAllele.getKey(), proportion);
                if ( ! deNovoRange.contains(proportion) ) {
                    //System.out.println("Filtering deNovo by proportion: is "+proportion+" should be in range "+deNovoRange.lower+"-"+deNovoRange.upper);
                    violation.filter();
                }
            }

            Pair<Allele,Integer> triAl = getTriAllelicQuality(tracker,ref,trio,splitContext);
            if ( triAl != null ) {
                violation.addAttribute(MendelianInfoKey.TriAllelicBase.getKey(),triAl.first.toString());
                violation.addAttribute(MendelianInfoKey.TriAllelicQuality.getKey(),triAl.second);
                if ( triAl.second >= deNovoTriQ ) {
                    violation.filter();
                }
            }

        } else {
            violation.filter();
        }

        return violation;
    }

    private MendelianViolation assessOppositeHomozygote(VariantContext trio, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        MendelianViolation violation = new MendelianViolation(trio,MendelianViolationType.OPPOSITE_HOMOZYGOTE);
        logger.debug(getParentalDeletion(trio));
        violation.addAttribute(MendelianInfoKey.ParentalDeletion.getKey(),getParentalDeletion(trio));

        // look for tri-allelic sites mis-called as hom -- as a speedup we do this only at non-filtered, non genotype error sites

        if ( ! trio.isFiltered()  ) {
            Map<String,StratifiedAlignmentContext> splitCon = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
            Pair<Allele,Integer> triAl = getTriAllelicQuality(tracker, ref, trio, splitCon);
            if ( triAl != null ) {
                violation.addAttribute(MendelianInfoKey.TriAllelicBase.getKey(),triAl.first.toString());
                violation.addAttribute(MendelianInfoKey.TriAllelicQuality.getKey(),triAl.second);
                if ( triAl.second >= opHomTriQ ) {
                    violation.filter();
                }
            }

            Double childProp = getAlleleProportion(trio.getGenotype(trioStructure.child).getAllele(0),splitCon.get(trioStructure.child));
            Double motherProp = getAlleleProportion(trio.getGenotype(trioStructure.mom).getAllele(0),splitCon.get(trioStructure.mom));
            Double fatherProp = getAlleleProportion(trio.getGenotype(trioStructure.dad).getAllele(0),splitCon.get(trioStructure.dad));
            if ( childProp != null ) {
                violation.addAttribute(MendelianInfoKey.ProportionOfParentAllele.getKey(),childProp);
                if ( ! opHomRange.contains(childProp) ) {
                    violation.filter();
                }
            }

            if ( motherProp != null && ! opHomRange.contains(motherProp) ) {
                violation.filter();
            }

            if ( fatherProp != null && ! opHomRange.contains(fatherProp) ) {
                violation.filter();
            }
        } else {
            violation.filter();
        }

        return violation;
    }

    private Double getAlleleProportion(Allele a, StratifiedAlignmentContext context) {
        int numParental = 0;
        int total = 0;
        if ( context != null ) {
            for ( PileupElement e : context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup()) {
                if ( e.getQual() >= 10 && e.getMappingQual() >= 10 ) {
                    total++;
                    if ( e.getBase() == a.getBases()[0]) {
                        numParental ++;
                    }
                }
            }
            return ( (double) numParental )/total;
        } else {
            return null;
        }


    }

    private Pair<Allele,Integer> getTriAllelicQuality(RefMetaDataTracker tracker, ReferenceContext ref, VariantContext var, Map<String,StratifiedAlignmentContext> strat) {
        int conf = 0;
        Allele alt = null;
        for ( Map.Entry<String,StratifiedAlignmentContext> sEntry : strat.entrySet() ) {
            VariantCallContext call = engine.runGenotyper(tracker, ref, sEntry.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE));
            if ( call != null && call.confidentlyCalled && call.vc != null ) {
                if ( call.vc.isSNP() ) {
                    if ( ! call.vc.getAlternateAllele(0).basesMatch(var.getAlternateAllele(0))) {
                        if ( alt == null ) {
                            alt = call.vc.getAlternateAllele(0);
                            conf = (int) Math.floor(10*call.vc.getNegLog10PError());
                        } else {
                            conf += (int) Math.floor(10*call.vc.getNegLog10PError());
                        }
                    }
                }
            }
        }

        if ( alt == null ) {
            return null;
        } else {
            return new Pair<Allele,Integer>(alt,conf);
        }
    }

    private String getParentalDeletion(VariantContext trio) {
        // case 1, mom and dad are both hom, so missing is the one whose alleles don't match the child
        if ( trio.getGenotype(trioStructure.mom).isHom() && trio.getGenotype(trioStructure.dad).isHom() ) {
            if ( trio.getGenotype(trioStructure.mom).isHomRef() == trio.getGenotype(trioStructure.child).isHomRef() ) {
                // mom and child both hom ref or hom var
                return trioStructure.dad;
            } else if ( trio.getGenotype(trioStructure.dad).isHomRef() == trio.getGenotype(trioStructure.child).isHomRef() ) {
                return trioStructure.mom;
            } else {
                // child matches neither parent
                return "genotypeError";
            }
        } else { // case 2, either mom or dad is het - the hom must be the missing allele
            return trio.getGenotype(trioStructure.mom).isHet() ? trioStructure.dad : trioStructure.mom;  
        }
    }

    /*
     *************** REDUCE
     */
    public VCFWriter reduce(MendelianViolation variant, VCFWriter writer) {
        if ( variant != null ) {
            trioStructure.updateHomozygosityRegions(variant,bedOutput);
            writer.add(variant.toVariantContext(),variant.getRefBase());
        }

        return writer;
    }

    /*
     ********** ON TRAVERSAL DONE
     */
    public void onTraversalDone(VCFWriter writer) {
        Map<String,HomozygosityRegion> regions = trioStructure.homozygousRegions;
        Map<String,Integer> counts = trioStructure.homozygousRegionCounts;
        List<HomozygosityRegion> to_print = new ArrayList<HomozygosityRegion>(3);
        for ( Map.Entry<String,HomozygosityRegion> entryRegion : regions.entrySet() ) {
            if ( entryRegion.getValue() != null ) {
                logger.info("---------------- REGION NOT FINALIZED -----------------");
                logger.info(String.format("%s,%s,%s,%d,%d",entryRegion.getKey(),entryRegion.getValue().regionStart,entryRegion.getValue().lastSeen,
                        entryRegion.getValue().deNovoSNPsInRegion,entryRegion.getValue().oppositeHomsInRegion));
                int chr_end = getToolkit().getSAMFileHeader().getSequenceDictionary().getSequence(entryRegion.getValue().getContigStr()).getSequenceLength();
                entryRegion.getValue().endedBy = GenomeLocParser.createGenomeLoc(entryRegion.getValue().getContigStr(),chr_end,chr_end);
                to_print.add(entryRegion.getValue());
            }
        }

        Collections.sort(to_print);
        for ( HomozygosityRegion hr : to_print ) {
            bedOutput.printf("%s%n",hr);
        }
    }

    /*
     ***************** STATIC METHODS
     */

    public static String getBedFileHeader() {
        return String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s","Chrom","first_seen_hom","last_seen_hom","first_seen_het",
                "sample_name","region_id","homozygous_region_size", "size_to_first_het","calls_within_region",
                "opposite_homozygotes_in_region","deNovo_SNPs_in_region");
    }
}
