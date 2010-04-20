/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.genotyper.BaseMismatchModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;

/**
 * Walker to calculate the number of mismatches, their base counts, and their quality sums at confidence ref sites" 
 */
@By(DataSource.REFERENCE)
public class LocusMismatchWalker extends LocusWalker<String,Integer> implements TreeReducible<Integer> {
    //@Argument(fullName="confidentRefThreshold",doc="Set the lod score that defines confidence in ref, defaults to 4", required=false)
    //int confidentRefThreshold = 5;
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
    @Argument(fullName="maxBaseQuality", doc = "Set the base quality score below which to ignore bases in the pileup, defaults to no restriction", required = false)
    int maxQualityScore = 99;
    @Argument(fullName="minMismatches", doc = "Minimum number of mismatches at a locus before a site is displayed", required = false)
    int minMismatches = 1;

    @Argument(fullName="skip", doc = "Only display every skip eligable sites.  Defaults to all sites", required = false)
    int skip = 1;

    private UnifiedGenotyperEngine ug;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.baseModel = BaseMismatchModel.THREE_STATE;
        uac.ALL_BASES = true;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        // print the header
        out.printf("loc ref genotype genotypeQ depth nMM qSumMM A C G T%n");
    }

    public String map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        String result = null;

        ReadBackedPileup pileup = context.getBasePileup();
        if ( locusIsUsable(tracker, ref, pileup, context) ) {
            Genotype g = getGenotype(tracker, ref, context);
            if ( g != null )
                result = errorCounts( ref, pileup, g );
        }

        return result;
    }

    public Integer reduce( String map, Integer reduce  ) {
        if ( map != null && (reduce % skip == 0) )
            out.println(map);

        //if (reduce % skip == 0) System.out.printf("Keeping %d%n", reduce);

        return reduce + (map != null ? 1 : 0);
    }

    public Integer treeReduce( Integer reduce1, Integer reduce2 ) {
        return reduce1 + reduce2;
    }

    public Integer reduceInit() {
        return 1;
    }

    private String errorCounts( ReferenceContext ref, ReadBackedPileup pileup, Genotype g ) {
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
            return String.format("%s %c %10s %5.2f %d %d %d %s",
                    pileup.getLocation(), ref.getBase(),
                    getGenotypeClass(g), 10 * g.getNegLog10PError(),
                    usableDepth, nMismatches, qSumMismatches, baseCountString);
        }

        return null;
    }

    private String getGenotypeClass(Genotype g) {
        if ( g.isHomRef() ) return "HOM-REF";
        else if ( g.isHet() ) return "HET";
        else if ( g.isHom() ) return "HOM-NONREF";
        else throw new StingException("Unexpected genotype in getGenotypeClass " + g);
    }

    public boolean useRead( PileupElement e ) {
        if ( e.getRead().getMappingQuality() <= minMappingQuality ) {
            return false;
        } else if ( ! BaseUtils.isRegularBase( e.getBase() ) ) {
            return false;
        } else if ( e.getQual() <= minQualityScore || e.getQual() > maxQualityScore ) {
            return false;
        } else {
            return true;
        }
    }

    private boolean locusIsUsable( RefMetaDataTracker tracker, ReferenceContext ref, ReadBackedPileup pileup, AlignmentContext context ) {
        return BaseUtils.isRegularBase(ref.getBase()) &&
                pileup.size() >= minDepth && pileup.size() < maxDepth &&
                notCoveredByVariations(tracker, ref) &&
                pileupContainsNoNs(pileup);
//        pileupContainsNoNs(pileup) &&
//        baseIsConfidentRef(tracker,ref,context);
    }

    private boolean notCoveredByVariations( RefMetaDataTracker tracker, ReferenceContext ref ) {
        Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref);
        // TODO: check this logic. I think it's the best approximation of what was here before, but it's a different system
        if (vcs != null && vcs.size() > 0 ) {
                return false;
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

    private Genotype getGenotype( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        VariantCallContext calls = ug.runGenotyper(tracker,ref,context);
        if ( calls == null || calls.vc == null || calls.vc.getNSamples() == 0 || !calls.vc.isSNP() )
            return null;
        else {
            return calls.vc.getGenotype(0);
        }
    }

//    private boolean baseIsConfidentRef( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
//        Pair<VariationCall, List<Genotype>> calls = ug.map(tracker,ref,context);
//        if ( calls == null || calls.first == null)
//            return false;
//        else {
//            VariationCall var = calls.getFirst();
//            return var.isReference() && var.getNegLog10PError() > confidentRefThreshold;
//            //return  ( var.isReference() > 0 && !calls.second.get(0).isVariant(ref.getBase()) && calls.second.get(0).getNegLog10PError() > confidentRefThreshold );
//        }
//    }

    public void onTraversalDone(Integer result) {
        ;
    }
}
