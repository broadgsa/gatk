package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Merges read-back-phased and phase-by-transmission files.
 */
public class MergeAndMatchHaplotypes extends RodWalker<Integer, Integer> {
    @Output
    protected VCFWriter vcfWriter = null;

    @Input(fullName="pbt", shortName = "pbt", doc="Input VCF truth file", required=true)
    public RodBinding<VariantContext> pbtTrack;

    @Input(fullName="rbp", shortName = "rbp", doc="Input VCF truth file", required=true)
    public RodBinding<VariantContext> rbpTrack;

    private Map<String, Genotype> pbtCache = new HashMap<String, Genotype>();
    private Map<String, Genotype> rbpCache = new HashMap<String, Genotype>();

    private final String SOURCE_NAME = "MergeReadBackedAndTransmissionPhasedVariants";

    public void initialize() {
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(pbtTrack.getName());

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(VCFUtils.getHeaderFields(this.getToolkit()));

        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            Collection<VariantContext> pbts = tracker.getValues(pbtTrack, ref.getLocus());
            Collection<VariantContext> rbps = tracker.getValues(rbpTrack, ref.getLocus());

            VariantContext pbt = pbts.iterator().hasNext() ? pbts.iterator().next() : null;
            VariantContext rbp = rbps.iterator().hasNext() ? rbps.iterator().next() : null;

            if (pbt != null && rbp != null) {
                Map<String, Genotype> genotypes = pbt.getGenotypes();

                if (!rbp.isFiltered()) {
                    for (String sample : rbp.getSampleNames()) {
                        Genotype rbpg = rbp.getGenotype(sample);
                        Genotype pbtg = pbt.getGenotype(sample);

                        // Propagate read-backed phasing information to genotypes unphased by transmission
                        //if (!pbtg.isPhased() && rbpCache.containsKey(sample)) {
                        if (!pbtg.isPhased() && rbpg.isPhased() && rbpCache.containsKey(sample)) {
                            boolean orientationMatches = rbpCache.get(sample).sameGenotype(pbtCache.get(sample), false);

                            if (orientationMatches) {
                                pbtg = rbpg;
                            } else {
                                List<Allele> fwdAlleles = rbpg.getAlleles();
                                List<Allele> revAlleles = new ArrayList<Allele>();

                                for (int i = fwdAlleles.size() - 1; i >= 0; i--) {
                                    revAlleles.add(fwdAlleles.get(i));
                                }

                                pbtg = new Genotype(sample, revAlleles, rbpg.getNegLog10PError(), rbpg.getFilters(), rbpg.getAttributes(), rbpg.isPhased());
                            }
                        }

                        genotypes.put(sample, pbtg);

                        // Update the cache
                        if (/*rbpg.isPhased() &&*/ rbpg.isHet()) {
                            rbpCache.put(sample, rbpg);
                            pbtCache.put(sample, pbtg);
                        } else if (!rbpg.isPhased()) {
                            rbpCache.remove(sample);
                            pbtCache.remove(sample);
                        }
                    }
                }

                VariantContext newvc = new VariantContext(SOURCE_NAME, pbt.getChr(), pbt.getStart(), pbt.getStart(), pbt.getAlleles(), genotypes, pbt.getNegLog10PError(), pbt.getFilters(), pbt.getAttributes());
                vcfWriter.add(newvc);
            }
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
