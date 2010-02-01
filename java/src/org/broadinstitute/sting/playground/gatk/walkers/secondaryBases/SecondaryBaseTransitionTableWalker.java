package org.broadinstitute.sting.playground.gatk.walkers.secondaryBases;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.playground.utils.NamedTable;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import java.util.HashMap;

/**
 * Created By User: Michael Melgar
 * Date: Nov 2, 2009
 * Time: 8:45:32 PM
 * Given a secondary base annotated .bam file and a reference, this walker generates a table of secondary base counts
 * for all called loci in the .bam. Each base call made is an instance included in the table. Specifically, the walker
 * maps the following vector to a count of secondary bases:
 * <Called Genotype, Reference Base, Primary Base, Previous Genome Base, Read Group, Secondary Base>.
 */

@Reference(window=@Window(start=-1,stop=1))
public class SecondaryBaseTransitionTableWalker extends LocusWalker<Integer, Integer> {

    HashMap<String,Long> counts = new HashMap<String,Long>();
    private UnifiedGenotyperEngine ug;
    private NamedTable altTable;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.CONFIDENCE_THRESHOLD = 50;
        uac.ALL_BASES = true;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        altTable = new NamedTable();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        char refBase = Character.toUpperCase(ref.getBase());
        ReadBackedPileup pileup = context.getBasePileup();
        char[] contextBases = ref.getBases();
        char prevBase = Character.toUpperCase(contextBases[0]);
        char nextBase = Character.toUpperCase(contextBases[contextBases.length - 1]);

        if (contextBases.length == 3 && refBase != 'N' && pileup.getBases() != null && pileup.getSecondaryBases() != null) {
            VariantCallContext ugResult = ug.runGenotyper(tracker,ref,context);
            if (ugResult != null && ugResult.variation != null) {
                Genotype res = ugResult.genotypes.get(0);
                String call = res.getBases();
                String type;
                if (!res.isVariant(refBase)) {type = "homref";}
                else if (!res.isHet()) {type = "homvar";}
                else if (call.contains(Character.toString(refBase))) {type = "het";}
                else {type = "bad";}
                if (type != "bad") {
                    for (PileupElement element : pileup) {
                        char primaryBase = Character.toUpperCase((char)element.getBase());
                        char secondaryBase = Character.toUpperCase((char)element.getSecondBase());
                        String RG = element.getRead().getReadGroup().getReadGroupId();
                        if (secondaryBase != 'N' && secondaryBase != '.' && primaryBase != 'N') {
                            String strandRef;
                            String strandPrimary;
                            String strandPrev;
                            String strandSecondary;
                            if (!element.getRead().getReadNegativeStrandFlag()) {
                                strandRef = Character.toString(refBase);
                                strandPrimary = Character.toString(primaryBase);
                                strandPrev = Character.toString(prevBase);
                                strandSecondary = Character.toString(secondaryBase);
                            }
                            else {
                                strandRef = Character.toString(BaseUtils.simpleComplement(refBase));
                                strandPrimary = Character.toString(BaseUtils.simpleComplement(primaryBase));
                                strandPrev = Character.toString(BaseUtils.simpleComplement(nextBase));
                                strandSecondary = Character.toString(BaseUtils.simpleComplement(secondaryBase));
                            }
                            if (strandPrev.charAt(0) != 'N') {
                                String key = RG+' '+type+' '+call+' '+strandRef+' '+strandPrimary+' '+strandPrev+' '+strandSecondary;
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
        out.println("ReadGroup CallType CalledGenotype ReferenceBase PrimaryBase PreviousBase SecondaryBase Count");
        for (String key : counts.keySet()) {
            out.println(key + ' ' + counts.get(key).toString());
        }
        out.println("Processed " + result.toString() + " loci.");
    }
}
