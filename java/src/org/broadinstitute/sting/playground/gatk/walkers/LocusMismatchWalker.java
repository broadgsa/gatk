package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

@By(DataSource.REFERENCE)
public class LocusMismatchWalker extends LocusWalker<String,Integer> implements TreeReducible<Integer> {
    @Argument(fullName="confidentRefThreshold",doc="Set the lod score that defines confidence in ref, defaults to 4", required=false)
    int confidentRefThreshold = 5;
    @Argument(fullName="maxNumMismatches",doc="Set the maximum number of mismatches at a locus before choosing not to use it in calculation. Defaults to 1.", required=false)
    int maxNumMismatches = 100;
    @Argument(fullName="minMappingQuality", doc ="Set the alignment quality below which to ignore reads; defaults to 30", required = false)
    int minMappingQuality = 1;
    @Argument(fullName="minDepth",doc="Set the minimum number of reads at a locus before choosing to use it in calculation. Defaults to 20.", required=false)
    int minDepth = 10;
    @Argument(fullName="maxDepth",doc="Set the minimum number of reads at a locus before choosing to use it in calculation. Defaults to 20.", required=false)
    int maxDepth = 100;
    @Argument(fullName="minBaseQuality", doc = "Set the base quality score below which to ignore bases in the pileup, defaults to 20", required = false)
    int minQualityScore = 1;
    @Argument(fullName="minMismatches", doc = "Minimum number of mismatches at a locus before a site is displayed", required = false)
    int minMismatches = 1;

    private UnifiedGenotyper ug;

    public void initialize() {
        ug = new UnifiedGenotyper();
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        ug.initialize();
        uac.baseModel = BaseMismatchModel.THREE_STATE;
        uac.ALL_BASES = true;
        ug.setUnifiedArgumentCollection(uac);
    }

    public String map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        String result = null;

        ReadBackedPileup pileup = context.getPileup();
        if ( locusIsUsable(tracker, ref, pileup, context) ) {
            result = errorCounts( ref, pileup );
        }

        return result;
    }

    public Integer reduce( String map, Integer reduce  ) {
        if ( map != null )
            out.println(map);
        return reduce;
    }

    public Integer treeReduce( Integer reduce1, Integer reduce2 ) {
        return reduce1 + reduce2;
    }

    public Integer reduceInit() {
        out.printf("loc ref depth nMM qSumMM A C G T%n");
        return null;
    }

    private String errorCounts( ReferenceContext ref, ReadBackedPileup pileup ) {
        int[] baseCounts = { 0, 0, 0, 0 };
        int usableDepth = 0;
        int nMismatches = 0;
        int qSumMismatches = 0;

        for ( PileupElement e : pileup ) {
            if ( useRead(e) ) {
                //System.out.printf("Using %s%n", e.getRead().getReadName());
                baseCounts[e.getBaseIndex()] += 1;
                usableDepth++;
                if ( ! BaseUtils.basesAreEqual(e.getBase(), (byte)ref.getBase()) ) {
                    nMismatches++;
                    qSumMismatches += e.getQual();
                }
            }
        }

        if ( nMismatches < maxNumMismatches && nMismatches >= minMismatches && usableDepth >= minDepth ) {
            String baseCountString = "";
            for ( char b : BaseUtils.BASES ) {
                baseCountString += baseCounts[BaseUtils.simpleBaseToBaseIndex(b)] + " ";
            }
            return String.format("%s %c %d %d %d %s", pileup.getLocation(), ref.getBase(), usableDepth, nMismatches, qSumMismatches, baseCountString);
        }

        return null;
    }

    public boolean useRead( PileupElement e ) {
        if ( e.getRead().getMappingQuality() <= minMappingQuality ) {
            return false;
        } else if ( ! BaseUtils.isRegularBase( e.getBase() ) ) {
            return false;
        } else if ( e.getQual() <= minQualityScore ) {
            return false;
        } else {
            return true;
        }
    }

    private boolean locusIsUsable( RefMetaDataTracker tracker, ReferenceContext ref, ReadBackedPileup pileup, AlignmentContext context ) {
        return BaseUtils.isRegularBase(ref.getBase()) &&
                pileup.size() >= minDepth && pileup.size() < maxDepth &&
                notCoveredByVariations(tracker) &&
                pileupContainsNoNs(pileup) &&
                baseIsConfidentRef(tracker,ref,context);
    }

    private boolean notCoveredByVariations( RefMetaDataTracker tracker ) {
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if ( datum instanceof Variation || datum instanceof Genotype ) {
                //System.out.printf("Ignoring site because of %s%n", datum);
                return false;
            }
        }

        return true;
    }

    private boolean pileupContainsNoNs(ReadBackedPileup pileup) {
        for ( byte c : pileup.getBases() ) {
            if ( c == 'N' ) {
                return false;
            }
        }

        return true;
    }

    private boolean baseIsConfidentRef( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        Pair<VariationCall, List<Genotype>> calls = ug.map(tracker,ref,context);
        if ( calls == null || calls.first == null)
            return false;
        else {
            VariationCall var = calls.getFirst();
            return var.isReference() && var.getNegLog10PError() > confidentRefThreshold;
            //return  ( var.isReference() > 0 && !calls.second.get(0).isVariant(ref.getBase()) && calls.second.get(0).getNegLog10PError() > confidentRefThreshold );
        }
    }

    public void onTraversalDone(Integer result) {
        ;
    }
}
