/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.downsampling;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.collections.DefaultHashMap;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.text.XReadLines;
import htsjdk.variant.variantcontext.Allele;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class AlleleBiasedDownsamplingUtils {

    // define this class so that we can use Java generics below
    private final static class PileupElementList extends ArrayList<PileupElement> {}

    /**
     * Computes an allele biased version of the given pileup
     *
     * @param pileup                    the original pileup
     * @param downsamplingFraction      the fraction of total reads to remove per allele
     * @return allele biased pileup
     */
    public static ReadBackedPileup createAlleleBiasedBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction) {
        // special case removal of all or no reads
        if ( downsamplingFraction <= 0.0 )
            return pileup;
        if ( downsamplingFraction >= 1.0 )
            return new ReadBackedPileupImpl(pileup.getLocation(), new ArrayList<PileupElement>());

        final PileupElementList[] alleleStratifiedElements = new PileupElementList[4];
        for ( int i = 0; i < 4; i++ )
            alleleStratifiedElements[i] = new PileupElementList();

        // start by stratifying the reads by the alleles they represent at this position
        for ( final PileupElement pe : pileup ) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(pe.getBase());
            if ( baseIndex != -1 )
                alleleStratifiedElements[baseIndex].add(pe);
        }

        // make a listing of allele counts and calculate the total count
        final int[] alleleCounts = calculateAlleleCounts(alleleStratifiedElements);
        final int totalAlleleCount = (int)MathUtils.sum(alleleCounts);

        // do smart down-sampling
        final int numReadsToRemove = (int)(totalAlleleCount * downsamplingFraction); // floor
        final int[] targetAlleleCounts = runSmartDownsampling(alleleCounts, numReadsToRemove);

        final HashSet<PileupElement> readsToRemove = new HashSet<PileupElement>(numReadsToRemove);
        for ( int i = 0; i < 4; i++ ) {
            final PileupElementList alleleList = alleleStratifiedElements[i];
            // if we don't need to remove any reads, then don't
            if ( alleleCounts[i] > targetAlleleCounts[i] )
                readsToRemove.addAll(downsampleElements(alleleList, alleleCounts[i], alleleCounts[i] - targetAlleleCounts[i]));
        }

        // we need to keep the reads sorted because the FragmentUtils code will expect them in coordinate order and will fail otherwise
        final List<PileupElement> readsToKeep = new ArrayList<PileupElement>(totalAlleleCount - numReadsToRemove);
        for ( final PileupElement pe : pileup ) {
            if ( !readsToRemove.contains(pe) ) {
                readsToKeep.add(pe);
            }
        }

        return new ReadBackedPileupImpl(pileup.getLocation(), new ArrayList<PileupElement>(readsToKeep));
    }

    /**
     * Calculates actual allele counts for each allele (which can be different than the list size when reduced reads are present)
     *
     * @param alleleStratifiedElements       pileup elements stratified by allele
     * @return non-null int array representing allele counts
     */
    private static int[] calculateAlleleCounts(final PileupElementList[] alleleStratifiedElements) {
        final int[] alleleCounts = new int[alleleStratifiedElements.length];
        for ( int i = 0; i < alleleStratifiedElements.length; i++ ) {
            alleleCounts[i] = alleleStratifiedElements[i].size();
        }
        return alleleCounts;
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
     * Computes an allele biased version of the allele counts for a given pileup
     *
     * @param alleleCounts              the allele counts for the original pileup
     * @param numReadsToRemove          number of total reads to remove per allele
     * @return non-null array of new counts needed per allele
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
     * @param elements                  original list of pileup elements
     * @param originalElementCount      original count of elements (taking reduced reads into account)
     * @param numElementsToRemove       the number of records to remove
     * @return the list of pileup elements TO REMOVE
     */
    protected static List<PileupElement> downsampleElements(final List<PileupElement> elements, final int originalElementCount, final int numElementsToRemove) {
        // are there no elements to remove?
        if ( numElementsToRemove == 0 )
            return Collections.<PileupElement>emptyList();

        final ArrayList<PileupElement> elementsToRemove = new ArrayList<PileupElement>(numElementsToRemove);

        // should we remove all of the elements?
        if ( numElementsToRemove >= originalElementCount ) {
            elementsToRemove.addAll(elements);
            return elementsToRemove;
        }

        // create a bitset describing which elements to remove
        final BitSet itemsToRemove = new BitSet(originalElementCount);
        for ( final Integer selectedIndex : MathUtils.sampleIndicesWithoutReplacement(originalElementCount, numElementsToRemove) ) {
            itemsToRemove.set(selectedIndex);
        }

        int currentBitSetIndex = 0;
        for ( final PileupElement element : elements ) {
            if ( itemsToRemove.get(currentBitSetIndex++) ) {
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
     * @return list of reads TO REMOVE from allele biased down-sampling
     */
    public static <A extends Allele> List<GATKSAMRecord> selectAlleleBiasedReads(final Map<A, List<GATKSAMRecord>> alleleReadMap, final double downsamplingFraction) {
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
            if ( alleleCounts[i] > targetAlleleCounts[i] ) {
                readsToRemove.addAll(downsampleElements(alleleReadMap.get(alleles.get(i)), alleleCounts[i] - targetAlleleCounts[i]));
            }
        }

        return readsToRemove;
    }

    /**
     * Performs allele biased down-sampling on a pileup and computes the list of elements to remove
     *
     * @param reads                     original list of records
     * @param numElementsToRemove       the number of records to remove
     * @return the list of pileup elements TO REMOVE
     */
    protected static List<GATKSAMRecord> downsampleElements(final List<GATKSAMRecord> reads, final int numElementsToRemove) {
        // are there no elements to remove?
        if ( numElementsToRemove == 0 )
            return Collections.<GATKSAMRecord>emptyList();

        final ArrayList<GATKSAMRecord> elementsToRemove = new ArrayList<GATKSAMRecord>(numElementsToRemove);
        final int originalElementCount = reads.size();

        // should we remove all of the elements?
        if ( numElementsToRemove >= originalElementCount ) {
            elementsToRemove.addAll(reads);
            return elementsToRemove;
        }

        // create a bitset describing which elements to remove
        final BitSet itemsToRemove = new BitSet(originalElementCount);
        for ( final Integer selectedIndex : MathUtils.sampleIndicesWithoutReplacement(originalElementCount, numElementsToRemove) ) {
            itemsToRemove.set(selectedIndex);
        }

        int currentBitSetIndex = 0;
        for ( final GATKSAMRecord read : reads ) {
            if ( itemsToRemove.get(currentBitSetIndex++) )
                elementsToRemove.add(read);
        }

        return elementsToRemove;
    }

    /**
     * Create sample-contamination maps from file
     *
     * @param ContaminationFractionFile   Filename containing two columns: SampleID and Contamination
     * @param AvailableSampleIDs          Set of Samples of interest (no reason to include every sample in file) or null to turn off checking
     * @param logger                      for logging output
     * @return sample-contamination Map
     */

    public static DefaultHashMap<String, Double> loadContaminationFile(File ContaminationFractionFile, final Double defaultContaminationFraction, final Set<String> AvailableSampleIDs, Logger logger) throws GATKException {
        DefaultHashMap<String, Double> sampleContamination = new DefaultHashMap<String, Double>(defaultContaminationFraction);
        Set<String> nonSamplesInContaminationFile = new HashSet<String>(sampleContamination.keySet());
        try {

            XReadLines reader = new XReadLines(ContaminationFractionFile, true);
            for (String line : reader) {

                if (line.length() == 0) {
                    continue;
                }

                StringTokenizer st = new StringTokenizer(line,"\t");

                String fields[] = new String[2];
                try {
                    fields[0] = st.nextToken();
                    fields[1] = st.nextToken();
                } catch(NoSuchElementException e){
                    throw new UserException.MalformedFile("Contamination file must have exactly two, tab-delimited columns. Offending line:\n" + line);
                }
                if(st.hasMoreTokens()) {
                    throw new UserException.MalformedFile("Contamination file must have exactly two, tab-delimited columns. Offending line:\n" + line);
                }

                if (fields[0].length() == 0 || fields[1].length() == 0) {
                    throw new UserException.MalformedFile("Contamination file can not have empty strings in either column. Offending line:\n" + line);
                }

                if (sampleContamination.containsKey(fields[0])) {
                    throw new UserException.MalformedFile("Contamination file contains duplicate entries for input name " + fields[0]);
                }

                try {
                    final Double contamination = Double.valueOf(fields[1]);
                    if (contamination < 0 || contamination > 1){
                        throw new UserException.MalformedFile("Contamination file contains unacceptable contamination value (must be 0<=x<=1): " + line);
                    }
                    if (AvailableSampleIDs==null || AvailableSampleIDs.contains(fields[0])) {// only add samples if they are in the sampleSet (or if it is null)
                        sampleContamination.put(fields[0], contamination);
                    }
                    else {
                        nonSamplesInContaminationFile.add(fields[0]);
                    }
                } catch (NumberFormatException e) {
                    throw new UserException.MalformedFile("Contamination file contains unparsable double in the second field. Offending line: " + line);
                }
            }


            //output to the user info lines telling which samples are in the Contamination File
            if (sampleContamination.size() > 0) {
                logger.info(String.format("The following samples were found in the Contamination file and will be processed at the contamination level therein: %s", sampleContamination.keySet().toString()));

                //output to the user info lines telling which samples are NOT in the Contamination File
                if(AvailableSampleIDs!=null){
                    Set<String> samplesNotInContaminationFile = new HashSet<String>(AvailableSampleIDs);
                    samplesNotInContaminationFile.removeAll(sampleContamination.keySet());
                    if (samplesNotInContaminationFile.size() > 0)
                        logger.info(String.format("The following samples were NOT found in the Contamination file and will be processed at the default contamination level: %s", samplesNotInContaminationFile.toString()));
                }
            }

            //output to the user Samples that do not have lines in the Contamination File
            if (nonSamplesInContaminationFile.size() > 0) {
                logger.info(String.format("The following entries were found in the Contamination file but were not SAMPLEIDs. They will be ignored: %s", nonSamplesInContaminationFile.toString()));
            }

            return sampleContamination;

        } catch (IOException e) {
            throw new GATKException("I/O Error while reading sample-contamination file " + ContaminationFractionFile.getName() + ": " + e.getMessage());
        }

    }
}
