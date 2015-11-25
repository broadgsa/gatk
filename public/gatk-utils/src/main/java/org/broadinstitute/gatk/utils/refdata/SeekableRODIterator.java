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

package org.broadinstitute.gatk.utils.refdata;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.iterators.PushbackIterator;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Wrapper class for iterators over ROD objects. It is assumed that the underlying iterator can only
 * perform standard next() operation, which advances it to the next ROD in the stream (i.e. reads the data file
 * line by line). This iterator 1) shifts the focus from record-based traversal to position-based traversal,
 * and 2) adds querying seekForward() method.
 *
 * Namely, this iterator's next() method advances not to the next ROD in the underlying stream, but to the next
 * genomic position covered by (at least one) ROD, and returns all RODs overlapping with that position as a RODRecordList
 * collection-like object. Similarly, when seekForward(interval) is called, this iterator skips all the RODs from the
 * underlying stream, until it reaches specified genomic interval, and returns the list of all RODs overlapping with that interval.
 *
 * NOTE: this iterator has a STATE: next() operation is not allowed after a seekForward() to a non-point (extended) interval
 * of length > 1. Such a call would leave the iterator in an inconsistent state. seekForward() can always be called after
 * either seekForward() or next() (as long as usual ordering criteria are satisfied: the query interval location can neither
 * start before the current position, nor end before the previous query end). seekForward to an interval of length 1
 * reenables next() operation. 
 *
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 10, 2009
 * Time: 6:20:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class SeekableRODIterator implements LocationAwareSeekableRODIterator {
    /**
     * Header for the datasource backing this iterator.
     */
    private final Object header;

    /**
     * The parser, used to construct new genome locs.
     */
    private final GenomeLocParser parser;

    private final SAMSequenceDictionary sequenceDictionary;

    private PushbackIterator<GATKFeature> it;
    List<GATKFeature> records = null;  // here we will keep a pile of records overlaping with current position; when we iterate
                               // and step out of record's scope, we purge it from the list
    String name = null; // name of the ROD track wrapped by this iterator. Will be pulled from underlying iterator.

    int curr_position = 0; // where the iterator is currently positioned on the genome
    int max_position = 0;  // the rightmost stop position of currently loaded records
    String curr_contig = null;   // what contig the iterator is currently on
    boolean next_is_allowed = true; // see discussion below. next() is illegal after seek-forward queries of length > 1

    // the stop position of the last query. We can query only in forward direction ("seek forward");
    // it is not only the start position of every successive query that can not be before the start
    // of the previous one (curr_start), but it is also illegal for a query interval to *end* before
    // the end of previous query, otherwise we can end up in an inconsistent state
    int curr_query_end = -1;

    // EXAMPLE of inconsistency curr_query_end guards against:
    //              record 1      record 2
    //             ----------     -----------
    // -------------------------------------------------- REF
    //         ------------------------- query 1 (interval 1)
    //               ----------  query 2 (interval 2)
    //                     --------------- query 3
    //
    // If we query first for interval 1, both record 1 and record 2 will be loaded.
    // Query for interval 2, on the other hand, should return only record 1, but after
    // query 1 was performed, record 2 is already loaded from the file. If, on the other hand,
    // we try to un-load it from memory, we won't be able to read it again. Hence query 2 is not
    // allowed after query 1. Note also, that curr_query_end is not equivalent to max_position:
    // the latter only tracks where currently loaded records end (and hence helps to re-load records);
    // after query 1 is performed, max_position will be the end of record 2, but query 3 is still
    // perfectly legal after query 1.
    //
    // IMPORTANT NOTE: it follows from the above discussion and example that next() is illegal after ANY
    // seek-forward query EXCEPT those that are performed with length-1 intervals (queryInterval.start=queryinteval.stop).
    // Indeed, in the example above, after, e.g., query 1 is performed, the iterator is "located" at the start
    // of interval 1, but record1 and record 2 are already loaded. On the other hand, a subsequent call to next() would
    // need to shift iterator's position by 1 base and return only record 1.
    //
    // This implementation tracks the query history and makes next() illegal after a seekforward query of length > 1,
    // but re-enables next() again after a length-1 query.

    public SeekableRODIterator(Object header,SAMSequenceDictionary rodDictionary,SAMSequenceDictionary referenceDictionary,GenomeLocParser parser,CloseableIterator<GATKFeature> it) {
        this.header = header;
        this.parser = parser;
        this.sequenceDictionary = rodDictionary;
        this.it = new PushbackIterator<GATKFeature>(it);
        records = new LinkedList<GATKFeature>();
        // the following is a trick: we would like the iterator to know the actual name assigned to
        // the ROD implementing object we are working with. But the only way to do that is to
        // get an instance of that ROD and query it for its name. Now, the only generic way we have at this point to instantiate
        // the ROD is to make the underlying stream iterator to do it for us. So we are reading (or rather peeking into)
        // the first line of the track data file just to get the ROD object created.
        GATKFeature r = null;
        if (this.it.hasNext()) r = this.it.element();
        name = (r==null?null:r.getName());

        curr_contig = referenceDictionary.getSequence(0).getSequenceName();
    }

    /**
     * Gets the header associated with the backing input stream.
     * @return the ROD header.
     */
    @Override
    public Object getHeader() {
        return header;
    }

    /**
     * Gets the sequence dictionary associated with the backing input stream.
     * @return sequence dictionary from the ROD header.
     */
    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }


    /**
     * Returns true if the data we iterate over has records associated with (any, not necessarily adjacent)
     * genomic position farther along the reference.
     * @return
     */
    public boolean hasNext() {

        // if we did not walk to the very end of the interval(s) covered by currently loaded
        // annotations (records), then we definitely have data for next genomic location
        if ( curr_position < max_position ) return true;

        // we are past currently loaded stuff; we have next if there are more lines to load:
        return it.hasNext();
    }

    // Returns point location (i.e. genome loc of length 1) on the reference, to which this iterator will advance
    // upon next call to next().
    public GenomeLoc peekNextLocation() {
        if ( curr_position + 1 <= max_position ) return parser.createGenomeLoc(curr_contig,curr_position+1);

        // sorry, next reference position is not covered by the RODs we are currently holding. In this case,
        // the location we will jump to upon next call to next() is the start of the next ROD record that we did
        // not read yet:
        if ( it.hasNext() ) {
            GATKFeature r = it.element(); // peek, do not load!
            return parser.createGenomeLoc(r.getLocation().getContig(),r.getLocation().getStart());
        }
        return null; // underlying iterator has no more records, there is no next location!
    }

    /** Advances iterator to the next genomic position that has ROD record(s) associated with it,
     * and returns all the records overlapping with that position as a RODList. The location of the whole
     * RODList object will be set to the smallest interval subsuming genomic intervals of all returned records.
     * Note that next() is disabled (will throw an exception) after seekForward() operation with query length > 1.
     * @return list of all RODs overlapping with the next "covered" genomic position
     */
     public RODRecordList next() {
         if ( ! next_is_allowed )
             throw new ReviewedGATKException("Illegal use of iterator: Can not advance iterator with next() after seek-forward query of length > 1");

         curr_position++;
 //        curr_query_end = -1;

         if ( curr_position <= max_position ) {

             // we still have bases covered by at least one currently loaded record;
             // we have to purge only subset of records, on which we moved past the end
             purgeOutOfScopeRecords();
         } else {
             // ooops, we are past the end of all loaded records - kill them all at once,
             // load next record and reinitialize by fastforwarding current position to the start of next record
             records.clear();
             GATKFeature r = it.next(); // if hasNext() previously returned true, we are guaranteed that this call to reader.next() is safe
             records.add( r );
             curr_contig = r.getLocation().getContig();
             curr_position = r.getLocation().getStart();
             max_position = r.getLocation().getStop();
         }

         // current position is ste and at this point 'records' only keeps those annotations, on which we did not reach the end yet
         // (we might have reloaded records completely if it was necessary); but we are not guaranteed yet that we
         // hold ALL the records overlapping with the current position. Time to check if we just walked into the interval(s)
         // covered by new records, so we need to load them too:

         while ( it.hasNext() ) {
             GATKFeature r = it.element();
             if ( r == null ) {
                 it.next();
                 continue;
             }

             GenomeLoc currentContig = parser.createOverEntireContig(curr_contig);
             GenomeLoc thatContig = r.getLocation();

             if ( currentContig.isPast(thatContig) )
                 throw new UserException("LocationAwareSeekableRODIterator: contig " +r.getLocation().getContig() +
                         " occurs out of order in track " + r.getName() );
             if ( currentContig.isBefore(thatContig) ) break; // next record is on a higher contig, we do not need it yet...

             if ( r.getLocation().getStart() < curr_position )
                 throw new UserException("LocationAwareSeekableRODIterator: track "+r.getName() +
                         " is out of coordinate order on contig "+r.getLocation() + " compared to " + curr_contig + ":" + curr_position);

             if ( r.getLocation().getStart() > curr_position ) break; // next record starts after the current position; we do not need it yet

             r = it.next(); // we got here only if we do need next record, time to load it for real

             int stop = r.getLocation().getStop();
             if ( stop < curr_position ) throw new ReviewedGATKException("DEBUG: encountered contig that should have been loaded earlier"); // this should never happen
             if ( stop > max_position ) max_position = stop; // max_position keeps the rightmost stop position across all loaded records
             records.add(r);
         }

         // 'records' and current position are fully updated. Last, we need to set the location of the whole track
        // (collection of ROD records) to the genomic site we are currently looking at, and return the list

         return new RODRecordListImpl(name,records, parser.createGenomeLoc(curr_contig,curr_position));
     }

    /**
     * Removes from the underlying collection the last element returned by the
     * iterator (optional operation).  This method can be called only once per
     * call to <tt>next</tt>.  The behavior of an iterator is unspecified if
     * the underlying collection is modified while the iteration is in
     * progress in any way other than by calling this method.
     *
     * @throws UnsupportedOperationException if the <tt>remove</tt>
     *                                       operation is not supported by this Iterator.
     * @throws IllegalStateException         if the <tt>next</tt> method has not
     *                                       yet been called, or the <tt>remove</tt> method has already
     *                                       been called after the last call to the <tt>next</tt>
     *                                       method.
     */
    public void remove() {
        throw new UnsupportedOperationException("LocationAwareSeekableRODIterator does not implement remove() operation");
    }


    /**
     * Returns the current "position" (not location!! ;) ) of this iterator. This method is used by the sharding
     * system when it searches for available iterators in the pool that can be reused to resume traversal.
     * When iterator is advanced using next(), current position
     * is the same as 'location'. However, after a seekForward() query with extended interval, returned position
     * will be set to the last position of the query interval, to disable (illegal) attempts to roll the iterator
     * back and re-start traversal from current location.
     * @return Current ending position of the iterator, or null if no position exists.
     */
    public GenomeLoc position() {
        if ( curr_contig == null ) return null;
        if ( curr_query_end > curr_position )  {
            // do not attempt to reuse this iterator if the position we need it for lies before the end of last query performed
            return parser.createGenomeLoc(curr_contig,curr_query_end,curr_query_end);
        }
        else {
            return parser.createGenomeLoc(curr_contig,curr_position);
        }
    }

    /**
     * Seeks forward through the file until the specified interval is reached.
     * The location object <code>interval</code> can be either a single point or an extended interval. All
     * ROD records overlapping with the whole interval will be returned, or null if no such records exist.
     *
     * Query interval must start at or after the iterator's current location, or exception will be thrown.
     *
     * Query interval must end at or after the stop position of the previous query, if any, or an exception will
     * be thrown: subsequent queries that end before the stop of previous ones are illegal.
     *
     * If seekForward() is performed to an extended (length > 1 i.e. start != stop) interval, next() operation becomes
     * illegal (the iterator changes state). Only seekForward() calls are allowed thereafter, until a seekForward() call
     * to a length-1 interval is performed, which re-enables next(). seekForward() queries with length-1 intervals can
     * always be safely intermixed with next() (as long as ordering is respected and query intervals are at or after the
     * current position).
     *
     * Note that in contrast to
     * next() (which always advances current position of the iterator on the reference), this method scrolls
     * forward ONLY if the specified interval is ahead of the current location of
     * the iterator. However, if called again with the same 'interval' argument as before, seekForward will NOT
     * advance, but will simply return the same ROD list as before.
     *
     *
     * @param interval point-like genomic location to fastforward to.
     * @return ROD object at (or overlapping with) the specified position, or null if no such ROD exists.
     */
    public RODRecordList seekForward(GenomeLoc interval) {

        if ( interval.isBefore(parser.createOverEntireContig(curr_contig)) &&
             !(interval.getStart() == 0 && interval.getStop() == 0 && interval.getContig().equals(curr_contig)) ) // This criteria is syntactic sugar for 'seek to right before curr_contig'
            throw new ReviewedGATKException("Out of order query: query contig "+interval.getContig()+" is located before "+
                                     "the iterator's current contig");
        if ( interval.getContig().equals(curr_contig) ) {
            if ( interval.getStart() < curr_position )
                throw new ReviewedGATKException("Out of order query: query position "+interval +" is located before "+
                        "the iterator's current position "+curr_contig + ":" + curr_position);
            if ( interval.getStop() < curr_query_end )
                throw new ReviewedGATKException("Unsupported querying sequence: current query interval " +
                        interval+" ends before the end of previous query interval ("+curr_query_end+")");
        }

        curr_position = interval.getStart();
        curr_query_end = interval.getStop();

        next_is_allowed = ( curr_position == curr_query_end ); // we can call next() later only if interval length is 1

        if (  interval.getContig().equals(curr_contig) &&  curr_position <= max_position ) {
            // some of the intervals we are currently keeping do overlap with the query interval

            purgeOutOfScopeRecords();
        } else {
            // clean up and get ready for fast-forwarding towards the requested position
            records.clear();
            max_position = -1;
            curr_contig = interval.getContig();
        }

        // curr_contig and curr_position are set to where we asked to scroll to

        while ( it.hasNext() ) {
            GATKFeature r = it.next();
            if ( r == null ) continue;

            GenomeLoc currentContig = parser.createOverEntireContig(curr_contig);
            GenomeLoc thatContig = r.getLocation();

            if ( currentContig.isPast(thatContig) ) continue; // did not reach requested contig yet
            if ( currentContig.isBefore(thatContig) ) {
                it.pushback(r); // next record is on the higher contig, we do not need it yet...
                break;
            }

            // we get here if we are on the requested contig:

            if ( r.getLocation().getStop() < curr_position ) continue; // did not reach the requested interval yet

            if ( r.getLocation().getStart() > curr_query_end ) {
                // past the query interval
                it.pushback(r);
                break;
            }

            // we get here only if interval of the record r overlaps with query interval, so the record should be loaded
            if ( r.getLocation().getStop() > max_position ) max_position = r.getLocation().getStop();
            records.add(r);
        }

        if ( records.size() > 0 ) {
            return new RODRecordListImpl(name,records,interval);
        } else {
            return null;
        }

    }

    /**
     * Removes records that end before the curr_position from the list of currently kept records. This is a
     * convenience (private) shortcut that does not perform extensive checking. In particular, it assumes that
     * curr_position <= max_position, as well as that we are still on the same contig.
     */
    private void purgeOutOfScopeRecords() {
        Iterator<GATKFeature> i = records.iterator();
        while ( i.hasNext() ) {
            GATKFeature r = i.next();
            if ( r.getLocation().getStop() < curr_position ) {
                i.remove(); // we moved past the end of interval the record r is associated with, purge the record forever
            }
        }

    }

    @Override
    public void close() {
        if (this.it != null) ((CloseableIterator)this.it.getUnderlyingIterator()).close();
    }

}
