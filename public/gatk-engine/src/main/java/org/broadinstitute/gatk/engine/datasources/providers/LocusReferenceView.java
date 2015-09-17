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

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
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
     * Stop of the expanded window for which the reference context should be provided,
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
     * @param walker input walker
     * @param provider source for locus data.
     */
    public LocusReferenceView( Walker walker, LocusShardDataProvider provider ) {
        super(provider);
        initializeBounds(provider);

        // Retrieve information about the window being accessed.
        if( walker.getClass().isAnnotationPresent(Reference.class) ) {
            Window window = walker.getClass().getAnnotation(Reference.class).window();

            if( window.start() > 0 ) throw new ReviewedGATKException( "Reference window starts after current locus" );
            if( window.stop() < 0 ) throw new ReviewedGATKException( "Reference window ends before current locus" );

            windowStart = window.start();

            if ( walker.getArguments() == null ){
                windowStop = window.stop();
            } else {
                // Use reference arguments if set, otherwise use the annotation
                windowStop = walker.getArguments().reference_window_stop != GATKArgumentCollection.DEFAULT_REFERENCE_WINDOW_STOP ?
                        walker.getArguments().reference_window_stop : window.stop();
            }
        }
        else {
            windowStart = 0;
            if ( walker.getArguments() == null ){
                windowStop = 0;
            } else {
                windowStop = walker.getArguments().reference_window_stop;
            }
        }

        if(bounds != null) {
            int expandedStart = getWindowStart( bounds );
            int expandedStop  = getWindowStop( bounds );
            initializeReferenceSequence(genomeLocParser.createGenomeLoc(bounds.getContig(), bounds.getContigIndex(), expandedStart, expandedStop));
        }
    }

    /**
     * Initialize the bounds of this shard, trimming the bounds so that they match the reference.
     * @param provider Provider covering the appropriate locus.
     */
    private void initializeBounds(LocusShardDataProvider provider) {
        if(provider.getLocus() != null) {
            int sequenceLength = reference.getSequenceDictionary().getSequence(provider.getLocus().getContig()).getSequenceLength();
            bounds = genomeLocParser.createGenomeLoc(provider.getLocus().getContig(),
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

    protected GenomeLoc trimToBounds(GenomeLoc l) {
        int expandedStart = getWindowStart( bounds );
        int expandedStop  = getWindowStop( bounds );
        if ( l.getStart() < expandedStart ) l = genomeLocParser.setStart(l, expandedStart);
        if ( l.getStop() > expandedStop  ) l = genomeLocParser.setStop(l, expandedStop);
        return l;
    }

    public class Provider implements ReferenceContext.ReferenceContextRefProvider {
        int refStart, len;

        public Provider( int refStart, int len ) {
            this.refStart = refStart;
            this.len = len;
        }

        public byte[] getBases() {
            //System.out.printf("Getting bases for location%n");
            byte[] bases = new byte[len];
            System.arraycopy(referenceSequence.getBases(), refStart, bases, 0, len);
            return bases;
        }
    }

    /**
     * Gets the reference context associated with this particular point or extended interval on the genome.
     * @param genomeLoc Region for which to retrieve the base(s). If region spans beyond contig end or beyond current bounds, it will be trimmed down.
     * @return The base at the position represented by this genomeLoc.
     */
    public ReferenceContext getReferenceContext( GenomeLoc genomeLoc ) {
        //validateLocation( genomeLoc );

        GenomeLoc window = genomeLocParser.createGenomeLoc( genomeLoc.getContig(), genomeLoc.getContigIndex(),
                getWindowStart(genomeLoc), getWindowStop(genomeLoc) );

        int refStart = -1;
        if (bounds != null) {
            window = trimToBounds(window);
            refStart = (int)(window.getStart() - getWindowStart(bounds));
        }
        else {
            if(referenceSequence == null || referenceSequence.getContigIndex() != genomeLoc.getContigIndex())
                referenceSequence = reference.getSequence(genomeLoc.getContig());
            refStart = (int)window.getStart()-1;
        }

        int len = (int)window.size();
        return new ReferenceContext( genomeLocParser, genomeLoc, window, new Provider(refStart, len));
    }

    /**
     * Allow the user to pull reference info from any arbitrary region of the reference.
     * @param genomeLoc The locus.
     * @return A list of the bases starting at the start of the locus (inclusive) and ending
     *         at the end of the locus (inclusive).
     */
    public byte[] getReferenceBases( GenomeLoc genomeLoc ) {
        return super.getReferenceBases(genomeLoc);
    }

    /**
     * Gets the start of the expanded window, bounded if necessary by the contig.
     * @param locus The locus to expand.
     * @return The expanded window.
     */
    private int getWindowStart( GenomeLoc locus ) {
        // If the locus is not within the bounds of the contig it allegedly maps to, expand only as much as we can.
        if(locus.getStart() < 1) return 1;
//        if(locus.getStart() < 1) return locus.getStart();
        return Math.max( locus.getStart() + windowStart, 1 );
    }

    /**
     * Gets the stop of the expanded window, bounded if necessary by the contig.
     * @param locus The locus to expand.
     * @return The expanded window.
     */    
    private int getWindowStop( GenomeLoc locus ) {
        // If the locus is not within the bounds of the contig it allegedly maps to, expand only as much as we can.
        int sequenceLength = reference.getSequenceDictionary().getSequence(locus.getContig()).getSequenceLength();
        if(locus.getStop() > sequenceLength) return sequenceLength;
        return Math.min( locus.getStop() + windowStop, sequenceLength );
    }
}
