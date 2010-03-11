package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

import java.util.List;
/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Provides access to the portion of the reference covering a single locus.
 */
public class LocusReferenceView extends ReferenceView {
    /**
     * Bound the reference view to make sure all accesses are within the shard.
     */
    private GenomeLoc bounds;

    /**
     * Start of the expanded window for which the reference context should be provided,
     * relative to the locus in question.
     */
    private final int windowStart;


    /**
     * Start of the expanded window for which the reference context should be provided,
     * relative to the locus in question.
     */
    private final int windowStop;

    /**
     * Track the reference sequence and the last point accessed.  Used to
     * track state when traversing over the reference.
     */
    private ReferenceSequence referenceSequence;

    /**
     * Create a LocusReferenceView given no other contextual information about
     * the walkers, etc.
     * @param provider  source for locus data.
     */
    public LocusReferenceView( LocusShardDataProvider provider ) {
        super(provider);
        initializeBounds(provider);
        windowStart = windowStop = 0;
        initializeReferenceSequence(bounds);
    }

    /**
     * Create a new locus reference view.
     * @param provider source for locus data.
     */
    public LocusReferenceView( Walker walker, LocusShardDataProvider provider ) {
        super( provider );
        initializeBounds(provider);

        // Retrieve information about the window being accessed.
        if( walker.getClass().isAnnotationPresent(Reference.class) ) {
            Window window = walker.getClass().getAnnotation(Reference.class).window();

            if( window.start() > 0 ) throw new StingException( "Reference window starts after current locus" );
            if( window.stop() < 0 ) throw new StingException( "Reference window ends before current locus" );

            windowStart = window.start();
            windowStop = window.stop();
        }
        else {
            windowStart = 0;
            windowStop = 0;
        }

        if(bounds != null) {
            long expandedStart = getWindowStart( bounds );
            long expandedStop  = getWindowStop( bounds );
            initializeReferenceSequence(GenomeLocParser.createGenomeLoc(bounds.getContig(), expandedStart, expandedStop));
        }
    }

    /** Returns true if the specified location is fully within the bounds of the reference window loaded into
     *  this LocusReferenceView object.
     */
    public boolean isLocationWithinBounds(GenomeLoc loc) {
        return bounds.containsP(loc);
    }

    /** Ensures that specified location is within the bounds of the reference window loaded into this
     * LocusReferenceView object. If the location loc is within the current bounds (or if it is null), then nothing is done.
     * Otherwise, the bounds are expanded on either side, as needed, to accomodate the location, and the reference seuqence for the
     * new bounds is reloaded (can be costly!).
     * @param loc
     */
    public void expandBoundsToAccomodateLoc(GenomeLoc loc) {
        if ( bounds==null || loc==null) return; // can bounds be null actually???
        if ( isLocationWithinBounds(loc) ) return;
        if ( loc.getContigIndex() != bounds.getContigIndex() )
            throw new StingException("Illegal attempt to expand reference view bounds to accomodate location on a different contig.");

        bounds = GenomeLocParser.createGenomeLoc(bounds.getContigIndex(),
                                                 Math.min(bounds.getStart(),loc.getStart()),
                                                 Math.max(bounds.getStop(),loc.getStop()));
        long expandedStart = getWindowStart( bounds );
        long expandedStop  = getWindowStop( bounds );
        initializeReferenceSequence(GenomeLocParser.createGenomeLoc(bounds.getContig(), expandedStart, expandedStop));
    }

    /**
     * Initialize the bounds of this shard, trimming the bounds so that they match the reference.
     * @param provider Provider covering the appropriate locus.
     */
    private void initializeBounds(LocusShardDataProvider provider) {
        if(provider.getLocus() != null) {
            long sequenceLength = reference.getSequenceDictionary().getSequence(provider.getLocus().getContig()).getSequenceLength();
            bounds = GenomeLocParser.createGenomeLoc(provider.getLocus().getContig(),
                    Math.max(provider.getLocus().getStart(),1),
                    Math.min(provider.getLocus().getStop(),sequenceLength));
        }
        else
            bounds = null;
    }

    /**
     * Initialize reference sequence data using the given locus.
     * @param locus
     */
    private void initializeReferenceSequence( GenomeLoc locus ) {
        this.referenceSequence = reference.getSubsequenceAt( locus.getContig(), locus.getStart(), locus.getStop() );
    }

    /**
     * Gets the reference context associated with this particular point on the genome.
     * @param genomeLoc Region for which to retrieve the base.  GenomeLoc must represent a 1-base region.
     * @return The base at the position represented by this genomeLoc.
     */
    public ReferenceContext getReferenceContext( GenomeLoc genomeLoc ) {
        validateLocation( genomeLoc );

        GenomeLoc window = GenomeLocParser.createGenomeLoc( genomeLoc.getContig(), getWindowStart(genomeLoc), getWindowStop(genomeLoc) );
        char[] bases = null;

        if(bounds != null) {
            bases = StringUtil.bytesToString( referenceSequence.getBases(), (int)(window.getStart() - getWindowStart(bounds)), (int)window.size() ).toCharArray();
        }
        else {
            if(referenceSequence == null || referenceSequence.getContigIndex() != genomeLoc.getContigIndex())
                referenceSequence = reference.getSequence(genomeLoc.getContig());
            bases = StringUtil.bytesToString( referenceSequence.getBases(), (int)window.getStart()-1, (int)window.size()).toCharArray();
        }
        return new ReferenceContext( genomeLoc, window, bases );
    }

    /**
     * Allow the user to pull reference info from any arbitrary region of the reference.
     * @param genomeLoc The locus.
     * @return A list of the bases starting at the start of the locus (inclusive) and ending
     *         at the end of the locus (inclusive).
     */
    public char[] getReferenceBases( GenomeLoc genomeLoc ) {
        return super.getReferenceBases(genomeLoc);
    }

    /**
     * Validates that the genomeLoc is one base wide and is in the reference sequence.
     * @param genomeLoc location to verify.
     */
    private void validateLocation( GenomeLoc genomeLoc ) throws InvalidPositionException {
        if( bounds != null && !bounds.containsP(genomeLoc) )
            throw new InvalidPositionException(
                    String.format("Requested position %s not within interval %s", genomeLoc, bounds));
    }

    /**
     * Gets the start of the expanded window, bounded if necessary by the contig.
     * @param locus The locus to expand.
     * @return The expanded window.
     */
    private long getWindowStart( GenomeLoc locus ) {
        // If the locus is not within the bounds of the contig it allegedly maps to, don't expand the locus at all. 
        if(locus.getStart() < 1) return locus.getStart();
        return Math.max( locus.getStart() + windowStart, 1 );
    }

    /**
     * Gets the stop of the expanded window, bounded if necessary by the contig.
     * @param locus The locus to expand.
     * @return The expanded window.
     */    
    private long getWindowStop( GenomeLoc locus ) {
        // If the locus is not within the bounds of the contig it allegedly maps to, don't expand the locus at all.
        long sequenceLength = reference.getSequenceDictionary().getSequence(locus.getContig()).getSequenceLength();
        if(locus.getStop() > sequenceLength) return locus.getStop();
        return Math.min( locus.getStop() + windowStop, sequenceLength );
    }
}
