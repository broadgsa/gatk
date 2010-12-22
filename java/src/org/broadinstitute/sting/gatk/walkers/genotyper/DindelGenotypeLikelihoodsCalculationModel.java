/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.indels.HaplotypeIndelErrorModel;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.Haplotype;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broad.tribble.util.variantcontext.Allele;

import java.util.*;

public class DindelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private final int maxReadDeletionLength = 3;
    private final double insertionStartProbability = 1e-3;
    private final double insertionEndProbability = 0.5;
    private final double alphaDeletionProbability = 1e-3;
    private final int HAPLOTYPE_SIZE = 80;

    private final int minIndelCountForGenotyping;
    private final boolean getAlleleListFromVCF;
    // todo - the following  need to be exposed for command line argument control
    private final double indelHeterozygosity = 1.0/8000;
    boolean useFlatPriors = true;
    private static final boolean DEBUGOUT = false;
    private static final boolean DEBUG = false;
    private HaplotypeIndelErrorModel model;
    private GenomeLoc lastSiteVisited;
    private List<Allele> alleleList;
    protected DindelGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        model = new HaplotypeIndelErrorModel(maxReadDeletionLength, insertionStartProbability,
                insertionEndProbability, alphaDeletionProbability, HAPLOTYPE_SIZE, false, DEBUGOUT);
        alleleList = new ArrayList<Allele>();
        getAlleleListFromVCF = UAC.GET_ALLELES_FROM_VCF;
        minIndelCountForGenotyping = UAC.MIN_INDEL_COUNT_FOR_GENOTYPING;
    }


    private ArrayList<Allele> computeConsensusAlleles(ReferenceContext ref,
                                                      Map<String, StratifiedAlignmentContext> contexts,
                                                      StratifiedAlignmentContext.StratifiedContextType contextType) {
        Allele refAllele, altAllele;
        GenomeLoc loc = ref.getLocus();
        ArrayList<Allele> aList = new ArrayList<Allele>();

        if (DEBUG) {
            System.out.println("'''''''''''''''''''''");
            System.out.println("Loc:"+loc.toString());
        }
        HashMap<String,Integer> consensusIndelStrings = new HashMap<String,Integer>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = sample.getValue().getContext(contextType);

            final ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();
            insCount += indelPileup.getNumberOfInsertions();
            delCount += indelPileup.getNumberOfDeletions();
        }

        if (insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping)
            return aList;
        
        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = sample.getValue().getContext(contextType);

            final ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();




            for ( ExtendedEventPileupElement p : indelPileup.toExtendedIterable() ) {
                SAMRecord read = p.getRead();

                if (DEBUG && p.isIndel()) {
                    System.out.format("Read: %s, cigar: %s, aln start: %d, aln end: %d, p.len:%d, Type:%s, EventBases:%s\n",
                            read.getReadName(),read.getCigar().toString(),read.getAlignmentStart(),read.getAlignmentEnd(),
                            p.getEventLength(),p.getType().toString(), p.getEventBases());
                }


                String indelString = p.getEventBases();
                if (p.isInsertion()) {
                    boolean foundKey = false;
                    if (read.getAlignmentEnd() == loc.getStart()) {
                        // first corner condition: a read has an insertion at the end, and we're right at the insertion.
                        // In this case, the read could have any of the inserted bases and we need to build a consensus
                        for (String s : consensusIndelStrings.keySet()) {
                            int cnt = consensusIndelStrings.get(s);
                            if (s.startsWith(indelString)){
                                // case 1: current insertion is prefix of indel in hash map
                                consensusIndelStrings.put(s,cnt+1);
                                foundKey = true;
                                break;
                            }
                            else if (indelString.startsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                consensusIndelStrings.remove(s);
                                consensusIndelStrings.put(indelString,cnt+1);
                                foundKey = true;
                                break;
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            consensusIndelStrings.put(indelString,1);

                    }
                    else if (read.getAlignmentStart() == loc.getStart()+1) {
                        // opposite corner condition: read will start at current locus with an insertion
                        for (String s : consensusIndelStrings.keySet()) {
                            int cnt = consensusIndelStrings.get(s);
                            if (s.endsWith(indelString)){
                                // case 1: current insertion is suffix of indel in hash map
                                consensusIndelStrings.put(s,cnt+1);
                                foundKey = true;
                                break;
                            }
                            else if (indelString.endsWith(s)) {
                                // case 2: indel stored in hash table is suffix of current insertion
                                // In this case, new bases are new key.

                                consensusIndelStrings.remove(s);
                                consensusIndelStrings.put(indelString,cnt+1);
                                foundKey = true;
                                break;
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            consensusIndelStrings.put(indelString,1);

                    }
                    else {
                        // normal case: insertion somewhere in the middle of a read: add count to hash map
                        int cnt = consensusIndelStrings.containsKey(indelString)? consensusIndelStrings.get(indelString):0;
                        consensusIndelStrings.put(indelString,cnt+1);
                    }

                }
                else if (p.isDeletion()) {
                    indelString = String.format("D%d",p.getEventLength());
                    int cnt = consensusIndelStrings.containsKey(indelString)? consensusIndelStrings.get(indelString):0;
                    consensusIndelStrings.put(indelString,cnt+1);

                }
            }

            if (DEBUG) {
                int icount = indelPileup.getNumberOfInsertions();
                 int dcount = indelPileup.getNumberOfDeletions();
                if (icount + dcount > 0)
                {
                    List<Pair<String,Integer>> eventStrings = indelPileup.getEventStringsWithCounts(ref.getBases());
                    System.out.format("#ins: %d, #del:%d\n", insCount, delCount);

                    for (int i=0 ; i < eventStrings.size() ; i++ ) {
                        System.out.format("%s:%d,",eventStrings.get(i).first,eventStrings.get(i).second);
                        //                int k=0;
                    }
                    System.out.println();
                }
            }
        }

        int maxAlleleCnt = 0;
        String bestAltAllele = "";
        for (String s : consensusIndelStrings.keySet()) {
            int curCnt = consensusIndelStrings.get(s);
            if (curCnt > maxAlleleCnt) {
                maxAlleleCnt = curCnt;
                bestAltAllele = s;
            }
            if (DEBUG)
                System.out.format("Key:%s, number: %d\n",s,consensusIndelStrings.get(s)  );
        }         //gdebug-

        if (maxAlleleCnt <  minIndelCountForGenotyping)
            return aList;

        if (bestAltAllele.startsWith("D")) {
            // get deletion length
            int dLen = Integer.valueOf(bestAltAllele.substring(1));
            // get ref bases of accurate deletion
            int startIdxInReference = (int)(1+loc.getStart()-ref.getWindow().getStart());

            byte[] refBases = Arrays.copyOfRange(ref.getBases(),startIdxInReference,startIdxInReference+dLen);
            refAllele = Allele.create(refBases,true);
            altAllele = Allele.create(Allele.NULL_ALLELE_STRING, false);
        }
        else {
            // insertion case
            refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
            altAllele = Allele.create(bestAltAllele, false);
        }

        aList.add(0,refAllele);
        aList.add(1,altAllele);

        return aList;

    }
    public Allele getLikelihoods(RefMetaDataTracker tracker,
                                 ReferenceContext ref,
                                 Map<String, StratifiedAlignmentContext> contexts,
                                 StratifiedAlignmentContext.StratifiedContextType contextType,
                                 GenotypePriors priors,
                                 Map<String, BiallelicGenotypeLikelihoods> GLs,
                                 Allele alternateAlleleToUse) {

        if ( tracker == null )
            return null;


        GenomeLoc loc = ref.getLocus();
        Allele refAllele, altAllele;
        VariantContext vc = null;

        if (!ref.getLocus().equals(lastSiteVisited)) {
            // starting a new site: clear allele list
            alleleList.clear();
            lastSiteVisited = ref.getLocus().clone();


            if (getAlleleListFromVCF) {

                 for( final VariantContext vc_input : tracker.getVariantContexts(ref, "indels", null, ref.getLocus(), false, false) ) {
                     if( vc_input != null && vc_input.isIndel() && ref.getLocus().getStart() == vc_input.getStart()) {
                         vc = vc_input;
                         break;
                     }
                 }
                 // ignore places where we don't have a variant
                 if ( vc == null )
                     return null;


                 if (!vc.isIndel())
                     return null;

                alleleList.clear();
                for (Allele a : vc.getAlleles())
                    alleleList.add(a);

            }
            else {
                alleleList = computeConsensusAlleles(ref,contexts, contextType);
                if (alleleList.isEmpty())
                    return null;
            }
        }
        // protect against having an indel too close to the edge of a contig
        if (loc.getStart() <= HAPLOTYPE_SIZE)
            return null;

        if ( !(priors instanceof DiploidIndelGenotypePriors) )
            throw new StingException("Only diploid-based Indel priors are supported in the DINDEL GL model");

        refAllele = alleleList.get(0);
        altAllele = alleleList.get(1);
        int eventLength = refAllele.getBaseString().length() - altAllele.getBaseString().length();
        // assume only one alt allele for now
        if (eventLength<0)
            eventLength = - eventLength;

        int currentHaplotypeSize = HAPLOTYPE_SIZE;

        // int numSamples = getNSamples(contexts);
        List<Haplotype> haplotypesInVC = Haplotype.makeHaplotypeListFromAlleles( alleleList, loc.getStart(),
                ref, currentHaplotypeSize);
        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.
        // initialize the GenotypeLikelihoods
        GLs.clear();

        double[][] haplotypeLikehoodMatrix;

        if (useFlatPriors) {
            priors = new DiploidIndelGenotypePriors();
        }
        else
            priors = new DiploidIndelGenotypePriors(indelHeterozygosity,eventLength,currentHaplotypeSize);

        //double[] priorLikelihoods = priors.getPriors();

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = sample.getValue().getContext(contextType);

            ReadBackedPileup pileup = null;
            if (context.hasExtendedEventPileup())
                pileup = context.getExtendedEventPileup();
            else if (context.hasBasePileup())
                pileup = context.getBasePileup();

            if (pileup != null ) {
                haplotypeLikehoodMatrix = model.computeReadHaplotypeLikelihoods( pileup, haplotypesInVC);


                double[] genotypeLikelihoods = HaplotypeIndelErrorModel.getHaplotypeLikelihoods( haplotypeLikehoodMatrix);

                GLs.put(sample.getKey(), new BiallelicGenotypeLikelihoods(sample.getKey(),
                        refAllele,
                        altAllele,
                        genotypeLikelihoods[0],
                        genotypeLikelihoods[1],
                        genotypeLikelihoods[2],
                        getFilteredDepth(pileup)));
                //System.out.format("%4.2f %4.2f %4.2f\n",genotypeLikelihoods[0],genotypeLikelihoods[1], genotypeLikelihoods[2]);
            }
        }

        return refAllele;
    }
}