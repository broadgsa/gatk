/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.utils.BaseUtils;
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
