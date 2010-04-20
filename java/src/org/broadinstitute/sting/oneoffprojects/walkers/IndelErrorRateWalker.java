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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.CircularArray;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.commandline.Argument;

import java.util.List;
import java.util.LinkedList;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 5, 2010
 * Time: 5:25:02 PM
 * To change this template use File | Settings | File Templates.
 */
@Reference(window=@Window(start=-10,stop=10))
public class IndelErrorRateWalker extends LocusWalker<Integer,Integer> {
    @Argument(fullName="minCoverage",shortName="minC",doc="Assess only sites with coverage at or above the specified value.",required=true)
    int MIN_COVERAGE = 0;
    @Argument(fullName="maxCoverage",shortName="maxC",doc="Assess only sites with coverage at or below the specified value.",required=false)
    int MAX_COVERAGE = 1000000000;
    @Argument(fullName="maxIndels",shortName="maxI",doc="Assess only sites with no more indels than the specified value.",required=false)
    int MAX_INDELS = 1;
    @Argument(fullName="minSeparation",shortName="minS", doc="Ignore reference sites within that many bases of a NON-countable indel sites ( > MAX_INDELS ).",
    required=true)
    int MIN_DISTANCE = 10;
    private GenomeLoc skipToLoc = null;
    private List<ReadBackedExtendedEventPileup> countableIndelBuffer = new LinkedList<ReadBackedExtendedEventPileup>();

    private long totalObservationsMade = 0; // total number of observations at each reference base (regradless of the outcome);
                                            // in other words, it is Sum_{all assessed ref positions R} coverage(R)

    private CircularArray.Int coverageBuffer = new CircularArray.Int(MIN_DISTANCE);

    private int MAX_LENGTH = 40;
    private int[] delCounts = new int[MAX_LENGTH];
    private int[] insCounts = new int[MAX_LENGTH];


    @Override
    public boolean generateExtendedEvents() { return true; }

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    @Override
    public void initialize() {
        
    }

    private void countIndels(ReadBackedExtendedEventPileup p) {
        for ( ExtendedEventPileupElement pe : p ) {
            if ( ! pe.isIndel() ) continue;
            if ( pe.getEventLength() > MAX_LENGTH ) continue;
            if ( pe.isInsertion() ) insCounts[pe.getEventLength()-1]++;
            else delCounts[pe.getEventLength()-1]++;
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        // if we got to ignore a stretch of reference bases (because we have seen a non-countable indel recently), do it now:
        if ( skipToLoc != null ) {
            if ( ref.getLocus().compareTo(skipToLoc) < 0 ) return 0;
            else skipToLoc = null; // reached it, no need to do extra checks anymore
        }

        // if we previously cached indels that looked like something we'd like to count:
        if ( countableIndelBuffer.size() != 0 ) {
            // when we are at least MIN_DISTANCE bases away, we are guaranteed that we are not going
            // to run into a NON-countable indel anymore that is so close that it will render last countable indel useless.

            Iterator<ReadBackedExtendedEventPileup> iter = countableIndelBuffer.iterator();
            while ( iter.hasNext() ) {

                ReadBackedExtendedEventPileup p = iter.next();

                if ( ref.getLocus().distance(p.getLocation()) >= MIN_DISTANCE ) {
                    countIndels(p);
                    iter.remove();
                } else {
                    break;
                }
            }
        }

        // at this point we have counted (and discarded from the buffer) all indels that are sufficiently far behind
        // the current position. Now it's time to examine the pileup at the current position in more details:

        if ( context.hasExtendedEventPileup() ) {
            // if we got indels at current position:

            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup().getPileupWithoutMappingQualityZeroReads();
            if ( pileup.size() < MIN_COVERAGE ) return 0;

            if ( pileup.getNumberOfDeletions() + pileup.getNumberOfInsertions() > MAX_INDELS ) {
                // we got too many indel events. Maybe it's even a true event, and what we are looking for are
                // errors rather than true calls. Hence, we do not need these indels. We have to 1) discard
                // all remaining indels from the buffer: if they are still in the buffer, they are too close
                // to the current position; and 2) make sure that the next position at which we attempt to count again is
                // sufficiently far *after* the current position.
      //                      System.out.println("Non countable indel event at "+pileup.getLocation());
                countableIndelBuffer.clear();
                coverageBuffer.clear(); // we do not want to count observations (read bases) around non-countable indel as well
                skipToLoc = GenomeLocParser.createGenomeLoc(pileup.getLocation().getContigIndex(),pileup.getLocation().getStop()+pileup.getMaxDeletionLength()+MIN_DISTANCE+1);
 //                       System.out.println("Skip to "+skipToLoc);
            } else {
                // pileup does not contain too many indels, we need to store them in the buffer and count them later,
                // if a non-countable indel event(s) do not show up too soon:
                countableIndelBuffer.add(pileup);
            }
            return 0;
        }

        // we are here only if we have a "regular" base pileup; let's count coverage:

        ReadBackedPileup pileup = context.getBasePileup().getPileupWithoutMappingQualityZeroReads();

        int coverage = pileup.size() - pileup.getNumberOfDeletions(); // do not count bases that we did not sequence (deletions)

        if ( coverage < MIN_COVERAGE ) return 0;
 //                 System.out.println("at "+ref.getLocus()+"; adding "+coverageBuffer.get(0));

        if ( MIN_DISTANCE > 0 ) {

            totalObservationsMade += coverageBuffer.get(0);
            coverageBuffer.shiftData(1);
            coverageBuffer.set(MIN_DISTANCE-1,coverage);
        } else totalObservationsMade += coverage;
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return value+sum;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void onTraversalDone(Integer result) {
        for ( ReadBackedExtendedEventPileup p : countableIndelBuffer ) {
            countIndels(p);
        }
        for ( int i = 0 ; i < MIN_DISTANCE ; i++ ) {
            //System.out.println("done and printing "+coverageBuffer.get(i));
            totalObservationsMade += coverageBuffer.get(i);
        }
        super.onTraversalDone(result);


        out.println("Total observations (bases): "+totalObservationsMade);
        out.println("Indel error events:");
        out.println("len\tins_count\tins_rate\tdel_count\tdel_rate");
        int totalIns = 0;
        int totalDels = 0;
        for ( int i = 0 ; i < MAX_LENGTH ; i++ ) {
            out.printf("%d\t%d\t%.3g\t%d\t%.3g%n",i+1,
                                insCounts[i],((double)insCounts[i])/totalObservationsMade,
                                delCounts[i],((double)delCounts[i])/totalObservationsMade
                    );
            totalIns += insCounts[i];
            totalDels += delCounts[i];
        }
        out.println();
        out.print("Total indel errors found: "+(totalIns+totalDels));
        out.printf(" (rate: %.3g)%n",((double)(totalIns+totalDels))/totalObservationsMade);
        out.print("              insertions: "+totalIns);
        out.printf(" (rate: %.3g)%n",((double)totalIns)/totalObservationsMade);
        out.print("               deletions: "+totalDels);
        out.printf(" (rate: %.3g)%n",((double)totalDels)/totalObservationsMade);
    }

}
