/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.variant.variantcontext.Allele;

import java.io.PrintStream;
import java.util.*;

public class AlleleBiasedDownsamplingUtils {

    /**
     * Computes an allele biased version of the given pileup
     *
     * @param pileup                    the original pileup
     * @param downsamplingFraction      the fraction of total reads to remove per allele
     * @param log                       logging output
     * @return allele biased pileup
     */
    public static ReadBackedPileup createAlleleBiasedBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction, final PrintStream log) {
        // special case removal of all or no reads
        if ( downsamplingFraction <= 0.0 )
            return pileup;
        if ( downsamplingFraction >= 1.0 )
            return new ReadBackedPileupImpl(pileup.getLocation(), new ArrayList<PileupElement>());

        final ArrayList<PileupElement>[] alleleStratifiedElements = new ArrayList[4];
        for ( int i = 0; i < 4; i++ )
            alleleStratifiedElements[i] = new ArrayList<PileupElement>();

        // start by stratifying the reads by the alleles they represent at this position
        for ( final PileupElement pe : pileup ) {
            // we do not want to remove a reduced read
            if ( !pe.getRead().isReducedRead() ) {
                final int baseIndex = BaseUtils.simpleBaseToBaseIndex(pe.getBase());
                if ( baseIndex != -1 )
                    alleleStratifiedElements[baseIndex].add(pe);
            }
        }

        // make a listing of allele counts
        final int[] alleleCounts = new int[4];
        for ( int i = 0; i < 4; i++ )
            alleleCounts[i] = alleleStratifiedElements[i].size();

        // do smart down-sampling
        int numReadsToRemove = (int)(pileup.getNumberOfElements() * downsamplingFraction); // floor
        final int[] targetAlleleCounts = runSmartDownsampling(alleleCounts, numReadsToRemove);

        final HashSet<PileupElement> readsToRemove = new HashSet<PileupElement>(numReadsToRemove);
        for ( int i = 0; i < 4; i++ ) {
            final ArrayList<PileupElement> alleleList = alleleStratifiedElements[i];
            // if we don't need to remove any reads, then don't
            if ( alleleList.size() > targetAlleleCounts[i] )
                readsToRemove.addAll(downsampleElements(alleleList, alleleList.size() - targetAlleleCounts[i], log));
        }

        // clean up pointers so memory can be garbage collected if needed
        for ( int i = 0; i < 4; i++ )
            alleleStratifiedElements[i].clear();

        // we need to keep the reads sorted because the FragmentUtils code will expect them in coordinate order and will fail otherwise
        final List<PileupElement> readsToKeep = new ArrayList<PileupElement>(pileup.getNumberOfElements() - numReadsToRemove);
        for ( final PileupElement pe : pileup ) {
            if ( !readsToRemove.contains(pe) ) {
                readsToKeep.add(pe);
            }
        }

        return new ReadBackedPileupImpl(pileup.getLocation(), new ArrayList<PileupElement>(readsToKeep));
    }

    private static int scoreAlleleCounts(final int[] alleleCounts) {
        if ( alleleCounts.length < 2 )
            return 0;

        // sort the counts (in ascending order)
        final int[] alleleCountsCopy = alleleCounts.clone();
        Arrays.sort(alleleCountsCopy);

        final int maxCount = alleleCountsCopy[alleleCounts.length - 1];
        final int nextBestCount = alleleCountsCopy[alleleCounts.length - 2];

        int remainderCount = 0;
        for ( int i = 0; i < alleleCounts.length - 2; i++ )
            remainderCount += alleleCountsCopy[i];

        // try to get the best score:
        //    - in the het case the counts should be equal with nothing else
        //    - in the hom case the non-max should be zero
        return Math.min(maxCount - nextBestCount + remainderCount, Math.abs(nextBestCount + remainderCount));
    }

    /**
     * Computes an allele biased version of the given pileup
     *
     * @param alleleCounts              the original pileup
     * @param numReadsToRemove          fraction of total reads to remove per allele
     * @return allele biased pileup
     */
    protected static int[] runSmartDownsampling(final int[] alleleCounts, final int numReadsToRemove) {
        final int numAlleles = alleleCounts.length;

        int maxScore = scoreAlleleCounts(alleleCounts);
        int[] alleleCountsOfMax = alleleCounts;

        final int numReadsToRemovePerAllele = numReadsToRemove / 2;

        for ( int i = 0; i < numAlleles; i++ ) {
            for ( int j = i; j < numAlleles; j++ ) {
                final int[] newCounts = alleleCounts.clone();

                // split these cases so we don't lose on the floor (since we divided by 2)
                if ( i == j ) {
                    newCounts[i] = Math.max(0, newCounts[i] - numReadsToRemove);
                } else {
                    newCounts[i] = Math.max(0, newCounts[i] - numReadsToRemovePerAllele);
                    newCounts[j] = Math.max(0, newCounts[j] - numReadsToRemovePerAllele);
                }

                final int score = scoreAlleleCounts(newCounts);

                if ( score < maxScore ) {
                    maxScore = score;
                    alleleCountsOfMax = newCounts;
                }
            }
        }

        return alleleCountsOfMax;
    }

    /**
     * Performs allele biased down-sampling on a pileup and computes the list of elements to remove
     *
     * @param elements                  original list of records
     * @param numElementsToRemove       the number of records to remove
     * @param log                       logging output
     * @return the list of pileup elements TO REMOVE
     */
    private static <T> List<T> downsampleElements(final List<T> elements, final int numElementsToRemove, final PrintStream log) {
        ArrayList<T> elementsToRemove = new ArrayList<T>(numElementsToRemove);

        // are there no elements to remove?
        if ( numElementsToRemove == 0 )
            return elementsToRemove;

        // should we remove all of the elements?
        final int pileupSize = elements.size();
        if ( numElementsToRemove == pileupSize ) {
            logAllElements(elements, log);
            elementsToRemove.addAll(elements);
            return elementsToRemove;
        }

        // create a bitset describing which elements to remove
        final BitSet itemsToRemove = new BitSet(pileupSize);
        for ( Integer selectedIndex : MathUtils.sampleIndicesWithoutReplacement(pileupSize, numElementsToRemove) ) {
            itemsToRemove.set(selectedIndex);
        }

        for ( int i = 0; i < pileupSize; i++ ) {
            if ( itemsToRemove.get(i) ) {
                final T element = elements.get(i);
                logElement(element, log);
                elementsToRemove.add(element);
            }
        }

        return elementsToRemove;
    }

    /**
     * Computes reads to remove based on an allele biased down-sampling
     *
     * @param alleleReadMap             original list of records per allele
     * @param downsamplingFraction      the fraction of total reads to remove per allele
     * @param log                       logging output
     * @return list of reads TO REMOVE from allele biased down-sampling
     */
    public static List<GATKSAMRecord> selectAlleleBiasedReads(final Map<Allele, List<GATKSAMRecord>> alleleReadMap, final double downsamplingFraction, final PrintStream log) {
        int totalReads = 0;
        for ( final List<GATKSAMRecord> reads : alleleReadMap.values() )
            totalReads += reads.size();

        int numReadsToRemove = (int)(totalReads * downsamplingFraction);

        // make a listing of allele counts
        final List<Allele> alleles = new ArrayList<Allele>(alleleReadMap.keySet());
        alleles.remove(Allele.NO_CALL);    // ignore the no-call bin
        final int numAlleles = alleles.size();
        final int[] alleleCounts = new int[numAlleles];
        for ( int i = 0; i < numAlleles; i++ )
            alleleCounts[i] = alleleReadMap.get(alleles.get(i)).size();

        // do smart down-sampling
        final int[] targetAlleleCounts = runSmartDownsampling(alleleCounts, numReadsToRemove);

        final List<GATKSAMRecord> readsToRemove = new ArrayList<GATKSAMRecord>(numReadsToRemove);
        for ( int i = 0; i < numAlleles; i++ ) {
            final List<GATKSAMRecord> alleleBin = alleleReadMap.get(alleles.get(i));

            if ( alleleBin.size() > targetAlleleCounts[i] ) {
                readsToRemove.addAll(downsampleElements(alleleBin, alleleBin.size() - targetAlleleCounts[i], log));
            }
        }

        return readsToRemove;
    }

    private static <T> void logAllElements(final List<T> elements, final PrintStream log) {
        if ( log != null ) {
            for ( final T obj : elements ) {
                logElement(obj, log);
            }
        }
    }

    private static <T> void logElement(final T obj, final PrintStream log) {
        if ( log != null ) {

            final GATKSAMRecord read;
            if ( obj instanceof PileupElement )
                read = ((PileupElement)obj).getRead();
            else
                read = (GATKSAMRecord)obj;

            final SAMReadGroupRecord readGroup = read.getReadGroup();
            log.println(String.format("%s\t%s\t%s\t%s", read.getReadName(), readGroup.getSample(), readGroup.getLibrary(), readGroup.getPlatformUnit()));
        }
    }
}
