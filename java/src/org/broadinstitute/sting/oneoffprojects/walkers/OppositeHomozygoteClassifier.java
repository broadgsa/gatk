package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Jun 7, 2010
 */
public class OppositeHomozygoteClassifier extends RodWalker<VariantContext,VCFWriter> {
    @Argument(shortName="f",fullName="familyPattern",required=true,doc="Pattern for the family structure (usage: mom+dad=child)")
    String familyStr = null;

    private static Pattern FAMILY_PATTERN = Pattern.compile("(.*)\\+(.*)=(.*)");
    private TrioStructure trioStructure; /* holds sample names of mom dad child */
    private ArrayList<VariantContext> contextBuffer; /* holds contexts until they're ready to be printed */
    private GenomeLoc homozygousRegionStartChild; /* start of the homozygous region for child */
    private int callsWithinHomozygousRegion; /* number of calls in the current homozygous child region */
    private int childHomozygousRegionCounter; /* holds number of child homozygous regions */

    public void initialize() {
        trioStructure = parseTrioDescription(familyStr);
        homozygousRegionStartChild = null;
        childHomozygousRegionCounter = 0;
        callsWithinHomozygousRegion = 0;
        contextBuffer = new ArrayList<VariantContext>(500);
    }

    public static class TrioStructure {
        public String mom, dad, child;
    }

    public static TrioStructure parseTrioDescription(String family) {
        Matcher m = FAMILY_PATTERN.matcher(family);
        if (m.matches()) {
            TrioStructure trio = new TrioStructure();
            //System.out.printf("Found a family pattern: %s%n", parent.FAMILY_STRUCTURE);
            trio.mom = m.group(1);
            trio.dad = m.group(2);
            trio.child = m.group(3);
            return trio;
        } else {
            throw new IllegalArgumentException("Malformatted family structure string: " + family + " required format is mom+dad=child");
        }
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        VariantContext trioVariants = tracker.getVariantContext(ref,"trio", EnumSet.allOf(VariantContext.Type.class), ref.getLocus(), true);
        // for this to work we need mismatching parents, one a homozyote, and a child homozygote

        if ( trioVariants == null ) {
            return null;
        }

        //System.out.println(": "+trioVariants.getGenotype(trioStructure.mom)+" dad: "+trioVariants.getGenotype(trioStructure.dad)+" child: "+trioVariants.getGenotype(trioStructure.child));
        if ( isOppositeHomozygoteSite(trioVariants) && ! trioVariants.isFiltered()) {
            // find out who the homozygote is in the parents
            if ( trioVariants.getGenotype(trioStructure.mom).isHom() ) {
                return assessVariant(trioStructure.mom,trioStructure.dad,trioStructure.child,trioVariants,ref,context);
            } else if ( trioVariants.getGenotype(trioStructure.dad).isHom() ) {
                return assessVariant(trioStructure.dad,trioStructure.mom,trioStructure.child,trioVariants,ref,context);
            }
        }

        return trioVariants;
    }

    private boolean isOppositeHomozygoteSite(VariantContext trio) {
        if ( trio.getGenotype(trioStructure.child).isHet() ) { // not valid at child het sites
            return false;
        } else if ( trio.getHetCount() > 1 ) { // child is not het, so if this is 2, mom and dad are both het, invalid
            return false;
        } else if ( trio.getGenotype(trioStructure.dad) == null || trio.getGenotype(trioStructure.mom) == null ) {
            return false;
        }

        return true;
    }

    private VariantContext assessVariant(String homParent, String otherParent, String child, VariantContext variCon, ReferenceContext refCon, AlignmentContext aliCon) {
        // see if the child matches the hom parent
        HashMap<String,Object> attributes = new HashMap<String,Object>(variCon.getAttributes());
        //System.out.println(refCon.getLocus()+" homParent: "+variCon.getGenotype(homParent).getGenotypeString()+" otherParent: "+variCon.getGenotype(otherParent).getGenotypeString()+" child: "+variCon.getGenotype(child).getGenotypeString());
        //do child and hom parent NOT match genotypes?
        if ( variCon.getGenotype(child).isHomRef() && variCon.getGenotype(homParent).isHomVar() ||
             variCon.getGenotype(child).isHomVar() && variCon.getGenotype(homParent).isHomRef()  ) {
            // check for genotyping error (other must be het, or opposite of first parent)
            if ( variCon.getGenotype(otherParent).isHet() || variCon.getGenotype(otherParent).isHomRef() != variCon.getGenotype(homParent).isHomRef() ) {
                attributes.put("opHom",homParent);
            } else {
                attributes.put("opHom","genotypeError");
            }
        } else if ( variCon.getGenotype(otherParent).isHom() && variCon.getGenotype(otherParent).isHomRef() != variCon.getGenotype(homParent).isHomRef() ) {
            // is other parent both homozygous and different?
            attributes.put("opHom",otherParent);
        }
        // todo -- assessment of site based on alignment contest (tri allelic? etc)
        return new VariantContext(variCon.getName(), variCon.getLocation(), variCon.getAlleles(), variCon.getGenotypes(), variCon.getNegLog10PError(), variCon.getFilters(), attributes);
    }

    public VCFWriter reduceInit() {
        VCFWriter writer = new VCFWriter(out);
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "OppositeHomozygoteClassifier"));
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        //hInfo.add(new VCFHeaderLine("opHom","Child tentatively inheritied the NULL ALLELE from this parent"));
        // todo -- add info field annotation lines: "opHom", "CHR", "CHRS"
        VCFHeader vcfHeader = new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit()));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    public VCFWriter reduce(VariantContext variCon, VCFWriter writer) {
        if ( variCon == null ) {
            return writer;
        }

        if ( homozygosityRegionIsBroken(variCon) && contextBuffer.size() > 0 ) {
            outputBufferedRecords(contextBuffer,variCon,writer);
        } else if ( ! variCon.isFiltered() && ! variCon.getGenotype(trioStructure.child).isNoCall() ) {
            callsWithinHomozygousRegion++;
        }

        if ( variCon.hasAttribute("opHom") ) {
            //writer.addRecord(VariantContextAdaptors.toVCF(variCon, variCon.getReference().getBases()[0]));
            contextBuffer.add(variCon);
        }

        if ( homozygousRegionStartChild == null ) {
            if ( variCon.getGenotype(trioStructure.child).isHom() && ! variCon.isFiltered() ) {
                homozygousRegionStartChild = variCon.getLocation();
            }
        }

        return writer;

    }

    private boolean homozygosityRegionIsBroken(VariantContext context) {
        // check to see if either the parent or child homozygosity regions have been broken
        if ( homozygousRegionStartChild != null && context.getGenotype(trioStructure.child).isHet() && ! context.isFiltered() ) {
            // NOTE: NO CALLS DO NOT BREAK REGIONS OF HOMOZYGOSITY
            return true;
        }

        return false;
    }

    private void outputBufferedRecords(List<VariantContext> bufCon, VariantContext varCon, VCFWriter writer) {
        // the buffered contexts all share one feature -- come from the same child homozygosity region
        String regionSize;
        if ( varCon != null ) {
            regionSize = Integer.toString(varCon.getLocation().distance(homozygousRegionStartChild));
        } else {
            regionSize = "unknown";
        }
        for ( VariantContext vc : bufCon ) {
            HashMap<String,Object> attributes = new HashMap<String,Object>(vc.getAttributes());
            attributes.put("CHR",childHomozygousRegionCounter);
            attributes.put("CHRS",regionSize);
            attributes.put("CHRNCALL",callsWithinHomozygousRegion);
            attributes.put("CHRNOPHOM",bufCon.size());
            VariantContext newVC = new VariantContext(vc.getName(),vc.getLocation(),vc.getAlleles(),vc.getGenotypes(),vc.getNegLog10PError(),vc.getFilters(),attributes);
            writer.addRecord(VariantContextAdaptors.toVCF(newVC,vc.getReference().getBases()[0]));
        }
        childHomozygousRegionCounter++;
        homozygousRegionStartChild = null;
        callsWithinHomozygousRegion = 0;
        bufCon.clear();
    }

    public void onTraversalDone(VCFWriter w) {
        outputBufferedRecords(contextBuffer,null,w);
    }
}
