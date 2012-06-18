/*
 * Copyright (c) 2012, The Broad Institute 
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Code for determining which indels are segregating among the samples.
 *
 * This code is just a refactor of the original code from Guillermo in the UG.
 *
 * @author Mark DePristo
 * @since 3/26/12
 */
public class ConsensusAlleleCounter {
    final protected static Logger logger = Logger.getLogger(ConsensusAlleleCounter.class);
    private final int minIndelCountForGenotyping;
    private final boolean doMultiAllelicCalls;
    private final double minFractionInOneSample;
    private final GenomeLocParser locParser;

    public ConsensusAlleleCounter(final GenomeLocParser locParser,
                                  final boolean doMultiAllelicCalls,
                                  final int minIndelCountForGenotyping,
                                  final double minFractionInOneSample) {
        this.minIndelCountForGenotyping = minIndelCountForGenotyping;
        this.doMultiAllelicCalls = doMultiAllelicCalls;
        this.minFractionInOneSample = minFractionInOneSample;
        this.locParser = locParser;
    }

    /**
     * Returns a list of Alleles at this locus that may be segregating
     * 
     * @param ref
     * @param contexts
     * @param contextType
     * @return
     */
    public List<Allele> computeConsensusAlleles(ReferenceContext ref,
                                                Map<String, AlignmentContext> contexts,
                                                AlignmentContextUtils.ReadOrientation contextType) {
        final Map<String, Integer> consensusIndelStrings = countConsensusAlleles(ref, contexts, contextType);
//        logger.info("Alleles at " + ref.getLocus());
//        for ( Map.Entry<String, Integer> elt : consensusIndelStrings.entrySet() ) {
//            logger.info("    " + elt.getValue() + " => " + elt.getKey());
//        }
        return consensusCountsToAlleles(ref, consensusIndelStrings);
    }

    //
    // TODO -- WARNING DOESN'T WORK WITH REDUCED READS
    //
    private Map<String, Integer> countConsensusAlleles(ReferenceContext ref,
                                                       Map<String, AlignmentContext> contexts,
                                                       AlignmentContextUtils.ReadOrientation contextType) {
        final GenomeLoc loc = ref.getLocus();
        HashMap<String, Integer> consensusIndelStrings = new HashMap<String, Integer>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            final AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedPileup indelPileup = context.getBasePileup();
            insCount += indelPileup.getNumberOfInsertionsAfterThisElement();
            delCount += indelPileup.getNumberOfDeletionsAfterThisElement();
        }

        if ( insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping )
            return Collections.emptyMap();

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            // todo -- warning, can be duplicating expensive partition here
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedPileup indelPileup = context.getBasePileup();

            final int nIndelReads = indelPileup.getNumberOfInsertionsAfterThisElement() + indelPileup.getNumberOfDeletionsAfterThisElement();
            final int nReadsOverall = indelPileup.getNumberOfElements();

            if ( nIndelReads == 0 || (nIndelReads / (1.0 * nReadsOverall)) < minFractionInOneSample ) {
//                if ( nIndelReads > 0 )
//                    logger.info("Skipping sample " + sample.getKey() + " with nIndelReads " + nIndelReads + " nReads " + nReadsOverall);
                continue;
//            } else {
//                logger.info("### Keeping sample " + sample.getKey() + " with nIndelReads " + nIndelReads + " nReads " + nReadsOverall);
            }


            for (PileupElement p : indelPileup) {
                final GATKSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());
                if (read == null)
                    continue;
                if (ReadUtils.is454Read(read)) {
                    continue;
                }

/*                if (DEBUG && p.isIndel()) {
                    System.out.format("Read: %s, cigar: %s, aln start: %d, aln end: %d, p.len:%d, Type:%s, EventBases:%s\n",
                            read.getReadName(),read.getCigar().toString(),read.getAlignmentStart(),read.getAlignmentEnd(),
                            p.getEventLength(),p.getType().toString(), p.getEventBases());
                }
   */
                String indelString = p.getEventBases();

                if ( p.isBeforeInsertion() ) {
                    // edge case: ignore a deletion immediately preceding an insertion as p.getEventBases() returns null [EB]
                    if ( indelString == null )
                        continue;

                    boolean foundKey = false;
                    // copy of hashmap into temp arrayList
                    ArrayList<Pair<String,Integer>> cList = new ArrayList<Pair<String,Integer>>();
                    for (String s : consensusIndelStrings.keySet()) {
                        cList.add(new Pair<String, Integer>(s,consensusIndelStrings.get(s)));
                    }

                    if (read.getAlignmentEnd() == loc.getStart()) {
                        // first corner condition: a read has an insertion at the end, and we're right at the insertion.
                        // In this case, the read could have any of the inserted bases and we need to build a consensus

                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            // case 1: current insertion is prefix of indel in hash map
                            if (s.startsWith(indelString)) {
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (indelString.startsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(indelString,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(indelString,1));

                    }
                    else if (read.getAlignmentStart() == loc.getStart()+1) {
                        // opposite corner condition: read will start at current locus with an insertion
                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            if (s.endsWith(indelString)) {
                                // case 1: current insertion (indelString) is suffix of indel in hash map (s)
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (indelString.endsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(indelString,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(indelString,1));


                    }
                    else {
                        // normal case: insertion somewhere in the middle of a read: add count to arrayList
                        int cnt = consensusIndelStrings.containsKey(indelString)? consensusIndelStrings.get(indelString):0;
                        cList.add(new Pair<String, Integer>(indelString,cnt+1));
                    }

                    // copy back arrayList into hashMap
                    consensusIndelStrings.clear();
                    for (Pair<String,Integer> pair : cList) {
                        consensusIndelStrings.put(pair.getFirst(),pair.getSecond());
                    }

                }
                else if ( p.isBeforeDeletedBase() ) {
                    indelString = String.format("D%d",p.getEventLength());
                    int cnt = consensusIndelStrings.containsKey(indelString)? consensusIndelStrings.get(indelString):0;
                    consensusIndelStrings.put(indelString,cnt+1);

                }
            }
        }

        return consensusIndelStrings;
    }

    private List<Allele> consensusCountsToAlleles(final ReferenceContext ref,
                                                  final Map<String, Integer> consensusIndelStrings) {
        final GenomeLoc loc = ref.getLocus();
        final Collection<VariantContext> vcs = new ArrayList<VariantContext>();
        int maxAlleleCnt = 0;
        Allele refAllele, altAllele;

        for (final Map.Entry<String, Integer> elt : consensusIndelStrings.entrySet()) {
            final String s = elt.getKey();
            final int curCnt = elt.getValue();
            int stop = 0;

            // if observed count if above minimum threshold, we will genotype this allele
            if (curCnt < minIndelCountForGenotyping)
                continue;

            if (s.startsWith("D")) {
                // get deletion length
                final int dLen = Integer.valueOf(s.substring(1));
                // get ref bases of accurate deletion
                final int startIdxInReference = 1 + loc.getStart() - ref.getWindow().getStart();
                stop = loc.getStart() + dLen;
                final byte[] refBases = Arrays.copyOfRange(ref.getBases(), startIdxInReference, startIdxInReference + dLen);

                if (Allele.acceptableAlleleBases(refBases, false)) {
                    refAllele = Allele.create(refBases, true);
                    altAllele = Allele.create(Allele.NULL_ALLELE_STRING, false);
                }
                else continue; // don't go on with this allele if refBases are non-standard
            } else {
                // insertion case
                if (Allele.acceptableAlleleBases(s, false)) { // don't allow N's in insertions
                    refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
                    altAllele = Allele.create(s, false);
                    stop = loc.getStart();
                }
                else continue; // go on to next allele if consensus insertion has any non-standard base.
            }


            final VariantContextBuilder builder = new VariantContextBuilder().source("");
            builder.loc(loc.getContig(), loc.getStart(), stop);
            builder.alleles(Arrays.asList(refAllele, altAllele));
            builder.referenceBaseForIndel(ref.getBase());
            builder.noGenotypes();
            if (doMultiAllelicCalls) {
                vcs.add(builder.make());
                if (vcs.size() >= GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
                    break;
            } else if (curCnt > maxAlleleCnt) {
                maxAlleleCnt = curCnt;
                vcs.clear();
                vcs.add(builder.make());
            }
        }

        if (vcs.isEmpty())
            return Collections.emptyList(); // nothing else to do, no alleles passed minimum count criterion

        final VariantContext mergedVC = VariantContextUtils.simpleMerge(locParser, vcs, null, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);
        return mergedVC.getAlleles();
    }
}
