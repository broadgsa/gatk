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

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.walkers.indels.HaplotypeIndelErrorModel;
import org.broadinstitute.sting.gatk.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.Haplotype;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

public class DindelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private final int maxReadDeletionLength = 3;
    private final int HAPLOTYPE_SIZE;

    private final int minIndelCountForGenotyping;
    private final boolean getAlleleListFromVCF;
    // todo - the following  need to be exposed for command line argument control
    private final double indelHeterozygosity = 1.0/8000;
    boolean useFlatPriors = true;

    private boolean DEBUG = false;

    // todo -cleanup
    private HaplotypeIndelErrorModel model;
    private PairHMMIndelErrorModel pairModel;

    private GenomeLoc lastSiteVisited;
    private List<Allele> alleleList;

    protected DindelGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
            pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY,UAC.INDEL_GAP_CONTINUATION_PENALTY,
                    UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.DO_CONTEXT_DEPENDENT_PENALTIES, UAC.dovit);
        alleleList = new ArrayList<Allele>();
        getAlleleListFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        minIndelCountForGenotyping = UAC.MIN_INDEL_COUNT_FOR_GENOTYPING;
        HAPLOTYPE_SIZE = UAC.INDEL_HAPLOTYPE_SIZE;
        DEBUG = UAC.OUTPUT_DEBUG_INDEL_INFO;
    }


    private ArrayList<Allele> computeConsensusAlleles(ReferenceContext ref,
                                                      Map<String, AlignmentContext> contexts,
                                                      AlignmentContextUtils.ReadOrientation contextType) {
        Allele refAllele=null, altAllele=null;
        GenomeLoc loc = ref.getLocus();
        ArrayList<Allele> aList = new ArrayList<Allele>();

        if (DEBUG) {
            System.out.println("'''''''''''''''''''''");
            System.out.println("Loc:"+loc.toString());
        }
        HashMap<String,Integer> consensusIndelStrings = new HashMap<String,Integer>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();
            insCount += indelPileup.getNumberOfInsertions();
            delCount += indelPileup.getNumberOfDeletions();
        }

        if (insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping)
            return aList;
        
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            // todo -- warning, can be duplicating expensive partition here
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();




            for ( ExtendedEventPileupElement p : indelPileup.toExtendedIterable() ) {
                SAMRecord read = p.getRead();
                if(ReadUtils.is454Read(read)) {
                    continue;
                }

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

            //System.out.println(new String(ref.getBases()));
            byte[] refBases = Arrays.copyOfRange(ref.getBases(),startIdxInReference,startIdxInReference+dLen);
            boolean ok = true;
            for (int i=0; i < refBases.length; i++)
                if (!BaseUtils.isRegularBase(refBases[i]))
                    ok = false;

            if (ok) {
                refAllele = Allele.create(refBases,true);
                altAllele = Allele.create(Allele.NULL_ALLELE_STRING, false);
            }
        }
        else {
            // insertion case
            boolean ok = true;
            for (int i=0; i < bestAltAllele.length(); i++)
                if (!BaseUtils.isRegularBase(bestAltAllele.getBytes()[i]))
                    ok = false;
            if (ok)  {
                refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
                altAllele = Allele.create(bestAltAllele, false);
            }
        }
        if (refAllele != null && altAllele != null) {
            aList.add(0,refAllele);
            aList.add(1,altAllele);
        }
        return aList;

    }
    public Allele getLikelihoods(RefMetaDataTracker tracker,
                                 ReferenceContext ref,
                                 Map<String, AlignmentContext> contexts,
                                 AlignmentContextUtils.ReadOrientation contextType,
                                 GenotypePriors priors,
                                 Map<String, BiallelicGenotypeLikelihoods> GLs,
                                 Allele alternateAlleleToUse,
                                 boolean useBAQedPileup) {

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

                 for( final VariantContext vc_input : tracker.getVariantContexts(ref, "alleles", null, ref.getLocus(), false, false) ) {
                     if( vc_input != null && ! vc_input.isFiltered() && vc_input.isIndel() && ref.getLocus().getStart() == vc_input.getStart()) {
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

        // check if there is enough reference window to create haplotypes (can be an issue at end of contigs)
        if (ref.getWindow().getStop() < loc.getStop()+HAPLOTYPE_SIZE)
            return null;
        if ( !(priors instanceof DiploidIndelGenotypePriors) )
            throw new StingException("Only diploid-based Indel priors are supported in the DINDEL GL model");

        if (alleleList.isEmpty())
            return null;
        
        refAllele = alleleList.get(0);
        altAllele = alleleList.get(1);
        int eventLength = altAllele.getBaseString().length() - refAllele.getBaseString().length();
        // assume only one alt allele for now

        List<Haplotype> haplotypesInVC;

        int hsize = (int)ref.getWindow().size()-Math.abs(eventLength)-1;
        int numPrefBases = ref.getLocus().getStart()-ref.getWindow().getStart()+1;
        if (DEBUG)
            System.out.format("hsize: %d eventLength: %d refSize: %d, locStart: %d numpr: %d\n",hsize,eventLength,
                    (int)ref.getWindow().size(), loc.getStart(), numPrefBases);

        haplotypesInVC = Haplotype.makeHaplotypeListFromAlleles( alleleList, loc.getStart(),
            ref, hsize, numPrefBases);

        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.
        // initialize the GenotypeLikelihoods
        GLs.clear();

        double[][] haplotypeLikehoodMatrix;


        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            ReadBackedPileup pileup = null;
            if (context.hasExtendedEventPileup())
                pileup = context.getExtendedEventPileup();
            else if (context.hasBasePileup())
                pileup = context.getBasePileup();

            if (pileup != null ) {

                haplotypeLikehoodMatrix = pairModel.computeReadHaplotypeLikelihoods( pileup, haplotypesInVC, ref, HAPLOTYPE_SIZE, eventLength);


                double[] genotypeLikelihoods = HaplotypeIndelErrorModel.getHaplotypeLikelihoods( haplotypeLikehoodMatrix);

                GLs.put(sample.getKey(), new BiallelicGenotypeLikelihoods(sample.getKey(),
                        refAllele,
                        altAllele,
                        genotypeLikelihoods[0],
                        genotypeLikelihoods[1],
                        genotypeLikelihoods[2],
                        getFilteredDepth(pileup)));
                if (DEBUG)
                    System.out.format("Sample:%s GL:%4.2f %4.2f %4.2f\n",sample.getKey(), genotypeLikelihoods[0],genotypeLikelihoods[1], genotypeLikelihoods[2]);
            }
        }

        return refAllele;
    }
}