package org.broadinstitute.sting.playground.gatk.walkers.secondaryBases;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.playground.utils.NamedTable;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.HashMap;

/**
 * Given a secondary base annotated .bam file and a reference, this walker generates a table of secondary base counts
 * for all called loci in the .bam. Each base call made is an instance included in the table. Specifically, the walker
 * maps the following vector to a count of secondary bases:
 * <HomRef/Het/HomVar, Read Group, Called Genotype, AlleleBalance, Reference Base, Primary Base, Previous Read Base, Secondary Base>.
 */

@Reference(window=@Window(start=-1,stop=1))
public class SecondaryBaseTransitionTableWalker extends LocusWalker<Integer, Integer> {

    HashMap<String,Long> counts = new HashMap<String,Long>();
    private UnifiedGenotyperEngine ug;
    private NamedTable altTable;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.STANDARD_CONFIDENCE_FOR_CALLING = uac.STANDARD_CONFIDENCE_FOR_EMITTING = 50.0;
        uac.ALL_BASES_MODE = true;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        altTable = new NamedTable();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        char refBase = Character.toUpperCase(ref.getBaseAsChar());
        ReadBackedPileup pileup = context.getBasePileup();
        int[] baseCounts = pileup.getBaseCounts();
        int length = 0;
        for (int i : baseCounts) {length += i;}
        byte[] contextBases = ref.getBases();
        byte prevBase = (byte)Character.toUpperCase(contextBases[0]);
        byte nextBase = (byte)Character.toUpperCase(contextBases[contextBases.length - 1]);

        if (contextBases.length == 3 && refBase != 'N' && pileup.getBases() != null && pileup.getSecondaryBases() != null) {
            VariantCallContext ugResult = ug.runGenotyper(tracker,ref,context);
            if (ugResult != null && ugResult.vc != null) {
                Genotype res = ugResult.vc.getGenotype(0);
                String call = res.getGenotypeString();
                String type;
                String alleleBalance = "N/A";
                if (res.isHomRef()) {
                    type = "homref";
                }
                else if (!res.isHet()) {type = "homvar";}
                else if (call.contains(Character.toString(refBase))) {
                    type = "het";
                    char alt;
                    if (call.charAt(0) == refBase) {alt = call.charAt(1);}
                    else {alt = call.charAt(0);}
                    double refCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(refBase)];
                    double altCount = baseCounts[BaseUtils.simpleBaseToBaseIndex(alt)];
                    alleleBalance = Double.toString(Math.round(100.0*refCount/(refCount + altCount))/100.0);
                }
                else {type = "bad";}
                if (!type.equals("bad")) {
                    for (PileupElement element : pileup) {
                        char primaryBase = Character.toUpperCase((char)element.getBase());
                        char secondaryBase = Character.toUpperCase((char)element.getSecondBase());
                        String RG = element.getRead().getReadGroup().getReadGroupId();
                        if (secondaryBase != 'N' && secondaryBase != '.' && primaryBase != 'N') {
                            String strandRef;
                            String strandPrimary;
                            String strandPrev;
                            String strandSecondary;
                            String strandCall;
                            if (!element.getRead().getReadNegativeStrandFlag()) {
                                strandRef = Character.toString(refBase);
                                strandPrimary = Character.toString(primaryBase);
                                strandPrev = Character.toString((char)prevBase);
                                strandSecondary = Character.toString(secondaryBase);
                                strandCall = call;
                            }
                            else {
                                strandRef = Character.toString(BaseUtils.simpleComplement(refBase));
                                strandPrimary = Character.toString(BaseUtils.simpleComplement(primaryBase));
                                strandPrev = Character.toString(BaseUtils.simpleComplement((char)nextBase));
                                strandSecondary = Character.toString(BaseUtils.simpleComplement(secondaryBase));
                                strandCall = BaseUtils.simpleReverseComplement(call);
                            }
                            if (strandPrev.charAt(0) != 'N') {
                                String key = type+' '+RG+' '+strandCall+' '+alleleBalance+' '+strandRef+' '+strandPrimary+' '+strandPrev+' '+strandSecondary;
                                if (counts.containsKey(key)) {
                                    counts.put(key, counts.get(key) + Long.valueOf(1));
                                }
                                else {
                                    counts.put(key, Long.valueOf(1));
                                }
                            }
                        }
                    }
                }
            }
        }
    return 1;
    }

    public Integer reduceInit() {return 0;}

    public Integer reduce(Integer value, Integer sum) {return sum + value;}

    public void onTraversalDone(Integer result) {
        out.println(">>>");
        out.println("Type ReadGroup CalledGenotype AlleleBalance ReferenceBase PrimaryBase PreviousBase SecondaryBase Count");
        for (String key : counts.keySet()) {
            out.println(key + ' ' + counts.get(key).toString());
        }
        out.println("Processed " + result.toString() + " loci.");
    }
}