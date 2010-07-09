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

package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;

/**
 * Emits intervals for the Local Indel Realigner to target for cleaning.  Ignores 454 and MQ0 reads.
 */
@ReadFilters({Platform454Filter.class, ZeroMappingQualityReadFilter.class})
@Reference(window=@Window(start=-1,stop=50))
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@By(DataSource.REFERENCE)
public class RealignerTargetCreator extends RodWalker<RealignerTargetCreator.Event, RealignerTargetCreator.Event> {

    // mismatch/entropy/SNP arguments
    @Argument(fullName="windowSize", shortName="window", doc="window size for calculating entropy or SNP clusters", required=false)
    protected int windowSize = 10;

    @Argument(fullName="mismatchFraction", shortName="mismatch", doc="fraction of base qualities needing to mismatch for a position to have high entropy; to disable set to <= 0 or > 1", required=false)
    protected double mismatchThreshold = 0.15;

    @Argument(fullName="minReadsAtLocus", shortName="minReads", doc="minimum reads at a locus to enable using the entropy calculation", required=false)
    protected int minReadsAtLocus = 4;

    // interval merging arguments
    @Argument(fullName="maxIntervalSize", shortName="maxInterval", doc="maximum interval size", required=false)
    protected int maxIntervalSize = 500;

    @Argument(fullName="realignReadsWithBadMates", required=false, doc="Should we try to realign paired-end reads whose mates map to other chromosomes?")
    protected boolean REALIGN_BADLY_MATED_READS = false;

    @Override
    public boolean generateExtendedEvents() { return true; }

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }


    public void initialize() {
        if ( windowSize < 2 )
            throw new StingException("Window Size must be an integer greater than 1");
    }

    public Event map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        boolean hasIndel = false;
        boolean hasInsertion = false;
        boolean hasPointEvent = false;

        long furthestStopPos = -1;

        // look for insertions in the extended context (we'll get deletions from the normal context)
        if ( context.hasExtendedEventPileup() ) {
            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup();
            if ( pileup.getNumberOfInsertions() > 0 ) {
                hasIndel = hasInsertion = true;
                // check the ends of the reads to see how far they extend
                for (ExtendedEventPileupElement p : pileup.toExtendedIterable() )
                    furthestStopPos = Math.max(furthestStopPos, p.getRead().getAlignmentEnd());
            }
        }

        // look at the rods for indels or SNPs
        if ( tracker != null ) {
            for ( VariantContext vc : tracker.getAllVariantContexts(ref) ) {
                switch ( vc.getType() ) {
                    case INDEL:
                        hasIndel = true;
                        if ( vc.isInsertion() )
                            hasInsertion = true;
                        break;
                    case SNP:
                        hasPointEvent = true;
                        break;
                    case MIXED:
                        hasPointEvent = true;
                        hasIndel = true;
                        if ( vc.isInsertion() )
                            hasInsertion = true;
                        break;
                    default:
                        break;
                }
                if ( hasIndel )
                    furthestStopPos = vc.getLocation().getStop();
            }
        }

        // look at the normal context to get deletions and positions with high entropy
        ReadBackedPileup pileup = context.getBasePileup();
        if ( pileup != null ) {

            int mismatchQualities = 0, totalQualities = 0;
            byte refBase = ref.getBase();
            for (PileupElement p : pileup ) {
                if ( !REALIGN_BADLY_MATED_READS && BadMateFilter.hasBadMate(p.getRead()) )
                    continue;

                // check the ends of the reads to see how far they extend
                furthestStopPos = Math.max(furthestStopPos, p.getRead().getAlignmentEnd());

                // is it a deletion? (sanity check in case extended event missed it)
                if ( p.isDeletion() ) {
                    hasIndel = true;
                }

                // look for mismatches
                else {
                    if ( p.getBase() != refBase )
                        mismatchQualities += p.getQual();
                    totalQualities += p.getQual();
                }
            }

            // make sure we're supposed to look for high entropy
            if ( mismatchThreshold > 0.0 &&
                    mismatchThreshold <= 1.0 &&
                    pileup.size() >= minReadsAtLocus &&
                    (double)mismatchQualities / (double)totalQualities >= mismatchThreshold )
                hasPointEvent = true;
        }

        // return null if no event occurred
        if ( !hasIndel && !hasPointEvent )
            return null;

        // return null if we didn't find any usable reads/rods associated with the event
        if ( furthestStopPos == -1 )
            return null;

        GenomeLoc eventLoc = context.getLocation();
        if ( hasInsertion )
            eventLoc =  GenomeLocParser.createGenomeLoc(eventLoc.getContigIndex(), eventLoc.getStart(), eventLoc.getStart()+1);
        else if ( hasIndel && (context.getBasePileup() == null || context.getBasePileup().size() == 0) )
            eventLoc =  GenomeLocParser.createGenomeLoc(eventLoc.getContigIndex(), eventLoc.getStart(), furthestStopPos);        

        EVENT_TYPE eventType = (hasIndel ? (hasPointEvent ? EVENT_TYPE.BOTH : EVENT_TYPE.INDEL_EVENT) : EVENT_TYPE.POINT_EVENT);

        return new Event(eventLoc, furthestStopPos, eventType);
    }

    public void onTraversalDone(Event sum) {
        if ( sum != null && sum.isReportableEvent() )
            out.println(sum.toString());
    }

    public Event reduceInit() {
        return null;
    }

    public Event reduce(Event value, Event sum) {
        // ignore no new events
        if ( value == null )
            return sum;

        // if it's the first good value, use it
        if ( sum == null )
            return value;

        // if we hit a new contig or they have no overlapping reads, then they are separate events - so clear sum
        if ( sum.loc.getContigIndex() != value.loc.getContigIndex() || sum.furthestStopPos < value.loc.getStart() ) {
            if ( sum.isReportableEvent() )
                out.println(sum.toString());
            return value;
        }

        // otherwise, merge the two events
        sum.merge(value);
        return sum;
    }

    private enum EVENT_TYPE { POINT_EVENT, INDEL_EVENT, BOTH }

    class Event {
        public long furthestStopPos;

        public GenomeLoc loc;
        public long eventStartPos;
        private long eventStopPos;
        private EVENT_TYPE type;
        private ArrayList<Long> pointEvents = new ArrayList<Long>();

        public Event(GenomeLoc loc, long furthestStopPos, EVENT_TYPE type) {
            this.loc = loc;
            this.furthestStopPos = furthestStopPos;
            this.type = type;

            if ( type == EVENT_TYPE.INDEL_EVENT || type == EVENT_TYPE.BOTH ) {
                eventStartPos = loc.getStart();
                eventStopPos = loc.getStop();
            } else {
                eventStartPos = -1;
                eventStopPos = -1;
            }

            if ( type == EVENT_TYPE.POINT_EVENT || type == EVENT_TYPE.BOTH ) {
                pointEvents.add(loc.getStart());
            }
        }

        public void merge(Event e) {

            // merges only get called for events with certain types
            if ( e.type == EVENT_TYPE.INDEL_EVENT || e.type == EVENT_TYPE.BOTH ) {
                if ( eventStartPos == -1 )
                    eventStartPos = e.eventStartPos;
                eventStopPos = e.eventStopPos;
                furthestStopPos = e.furthestStopPos;
            }

            if ( e.type == EVENT_TYPE.POINT_EVENT || e.type == EVENT_TYPE.BOTH ) {
                long newPosition = e.pointEvents.get(0);
                if ( pointEvents.size() > 0 ) {
                    long lastPosition = pointEvents.get(pointEvents.size()-1);
                    if ( newPosition - lastPosition < windowSize ) {
                        eventStopPos = Math.max(eventStopPos, newPosition);
                        furthestStopPos = e.furthestStopPos;

                        if ( eventStartPos == -1 )
                            eventStartPos = lastPosition;
                        else
                            eventStartPos = Math.min(eventStartPos, lastPosition);
                    }
                }
                pointEvents.add(newPosition);
            }
        }

        public boolean isReportableEvent() {
            return GenomeLocParser.validGenomeLoc(loc.getContig(), eventStartPos, eventStopPos) && eventStopPos >= 0 && eventStopPos - eventStartPos < maxIntervalSize;
        }

        public String toString() {
            return String.format("%s:%d-%d", loc.getContig(), eventStartPos, eventStopPos);
        }
    }
}