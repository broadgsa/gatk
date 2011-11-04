package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 *
 * User: aaron
 * Date: May 22, 2009
 * Time: 10:54:40 AM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 *         <p/>
 *         Class GenomeLocCollection
 *         <p/>
 *         a set of genome locations. This collection is self sorting,
 *         and will merge genome locations that are overlapping. The remove function
 *         will also remove a region from the list, if the region to remove is a
 *         partial interval of a region in the collection it will remove the region from
 *         that element.
 */
public class GenomeLocSortedSet extends AbstractSet<GenomeLoc> {
    private static Logger logger = Logger.getLogger(GenomeLocSortedSet.class);

    private GenomeLocParser genomeLocParser;

    // our private storage for the GenomeLoc's
    private List<GenomeLoc> mArray = new ArrayList<GenomeLoc>();

    /** default constructor */
    public GenomeLocSortedSet(GenomeLocParser parser) {
        this.genomeLocParser = parser;
    }

    public GenomeLocSortedSet(GenomeLocParser parser,GenomeLoc e) {
        this(parser);
        add(e);
    }

    public GenomeLocSortedSet(GenomeLocParser parser,Collection<GenomeLoc> l) {
        this(parser);

        for ( GenomeLoc e : l )
            add(e);
    }

    /**
     * Gets the GenomeLocParser used to create this sorted set.
     * @return The parser.  Will never be null.
     */
    public GenomeLocParser getGenomeLocParser() {
        return genomeLocParser;
    }

    /**
     * get an iterator over this collection
     *
     * @return an iterator<GenomeLoc>
     */
    public Iterator<GenomeLoc> iterator() {
        return mArray.iterator();
    }

    /**
     * return the size of the collection
     *
     * @return the size of the collection
     */
    public int size() {
        return mArray.size();
    }

    /**
     * Return the size, in bp, of the genomic regions by all of the regions in this set
     * @return size in bp of the covered regions
     */
    public long coveredSize() {
        long s = 0;
        for ( GenomeLoc e : this )
            s += e.size();
        return s;
    }

    /**
     * Return the number of bps before loc in the sorted set
     *
     * @param loc the location before which we are counting bases
     * @return
     */
    public long sizeBeforeLoc(GenomeLoc loc) {
        long s = 0;

        for ( GenomeLoc e : this ) {
            if ( e.isBefore(loc) )
                s += e.size();
            else if ( e.isPast(loc) )
                ; // don't do anything
            else // loc is inside of s
                s += loc.getStart() - e.getStart();
        }

        return s;
    }

    /**
     * determine if the collection is empty
     *
     * @return true if we have no elements
     */
    public boolean isEmpty() {
        return mArray.isEmpty();
    }

    /**
     * add a genomeLoc to the collection, simply inserting in order into the set
     *
     * @param e the GenomeLoc to add
     *
     * @return true
     */
    public boolean add(GenomeLoc e) {
        // assuming that the intervals coming arrive in order saves us a fair amount of time (and it's most likely true)
        if (mArray.size() > 0 && e.isPast(mArray.get(mArray.size() - 1))) {
            mArray.add(e);
            return true;
        } else {
            int loc = Collections.binarySearch(mArray,e);
            if (loc >= 0) {
                throw new ReviewedStingException("Genome Loc Sorted Set already contains the GenomicLoc " + e.toString());
            } else {
                mArray.add((loc+1) * -1,e);
                return true;
            }
        }
    }

    /**
     * Adds a GenomeLoc to the collection, merging it if it overlaps another region.
     * If it's not overlapping then we add it in sorted order.
     *
     * @param e the GenomeLoc to add to the collection
     *
     * @return true, if the GenomeLoc could be added to the collection
     */
    public boolean addRegion(GenomeLoc e) {
        if (e == null) {
            return false;
        }
        // have we added it to the collection?
        boolean haveAdded = false;

        /**
         * check if the specified element overlaps any current locations, if so
         * we should merge the two.
         */
        for (GenomeLoc g : mArray) {
            if (g.contiguousP(e)) {
                GenomeLoc c = g.merge(e);
                mArray.set(mArray.indexOf(g), c);
                haveAdded = true;
            } else if ((g.getContigIndex() == e.getContigIndex()) &&
                    (e.getStart() < g.getStart()) && !haveAdded) {
                mArray.add(mArray.indexOf(g), e);
                return true;
            } else if (haveAdded && ((e.getContigIndex() > e.getContigIndex()) ||
                    (g.getContigIndex() == e.getContigIndex() && e.getStart() > g.getStart()))) {
                return true;
            }
        }
        /** we're at the end and we haven't found locations that should fall after it,
         * so we'll put it at the end
         */
        if (!haveAdded) {
            mArray.add(e);
        }
        return true;
    }

    public GenomeLocSortedSet subtractRegions(GenomeLocSortedSet toRemoveSet) {
        LinkedList<GenomeLoc> good = new LinkedList<GenomeLoc>();
        Stack<GenomeLoc> toProcess = new Stack<GenomeLoc>();
        Stack<GenomeLoc> toExclude = new Stack<GenomeLoc>();

        // initialize the stacks
        toProcess.addAll(mArray);
        Collections.reverse(toProcess);
        toExclude.addAll(toRemoveSet.mArray);
        Collections.reverse(toExclude);

        int i = 0;
        while ( ! toProcess.empty() ) {    // while there's still stuff to process
            if ( toExclude.empty() ) {
                good.addAll(toProcess);         // no more excludes, all the processing stuff is good
                break;
            }

            GenomeLoc p = toProcess.peek();
            GenomeLoc e = toExclude.peek();

            if ( p.overlapsP(e) ) {
                toProcess.pop();
                for ( GenomeLoc newP : p.subtract(e) )
                    toProcess.push(newP);
            } else if ( p.compareContigs(e) < 0 ) {
                good.add(toProcess.pop());         // p is now good
            } else if ( p.compareContigs(e) > 0 ) {
                toExclude.pop();                 // e can't effect anything
            } else if ( p.getStop() < e.getStart() ) {
                good.add(toProcess.pop());         // p stops before e starts, p is good
            } else if ( e.getStop() < p.getStart() ) {
                toExclude.pop();                 // p starts after e stops, e is done
            } else {
                throw new ReviewedStingException("BUG: unexpected condition: p=" + p + ", e=" + e);
            }

            if ( i++ % 10000 == 0 )
                logger.debug("removeRegions operation: i = " + i);
        }

        return createSetFromList(genomeLocParser,good);
    }


    /**
     * a simple removal of an interval contained in this list.  The interval must be identical to one in the list (no partial locations or overlapping)
     * @param location the GenomeLoc to remove
     */
    public void remove(GenomeLoc location) {
        if (!mArray.contains(location)) throw new IllegalArgumentException("Unable to remove location: " + location + ", not in the list");
        mArray.remove(location);
    }

    /**
     * create a list of genomic locations, given a reference sequence
     *
     * @param dict the sequence dictionary to create a collection from
     *
     * @return the GenomeLocSet of all references sequences as GenomeLoc's
     */
    public static GenomeLocSortedSet createSetFromSequenceDictionary(SAMSequenceDictionary dict) {
        GenomeLocParser parser = new GenomeLocParser(dict);
        GenomeLocSortedSet returnSortedSet = new GenomeLocSortedSet(parser);
        for (SAMSequenceRecord record : dict.getSequences()) {
            returnSortedSet.add(parser.createGenomeLoc(record.getSequenceName(), 1, record.getSequenceLength()));
        }
        return returnSortedSet;
    }

    /**
     * Create a sorted genome location set from a list of GenomeLocs.
     *
     * @param locs the list<GenomeLoc>
     *
     * @return the sorted genome loc list
     */
    public static GenomeLocSortedSet createSetFromList(GenomeLocParser parser,List<GenomeLoc> locs) {
        GenomeLocSortedSet set = new GenomeLocSortedSet(parser);
        set.addAll(locs);
        return set;
    }


    /**
     * return a deep copy of this collection.
     *
     * @return a new GenomeLocSortedSet, identical to the current GenomeLocSortedSet.
     */
    public GenomeLocSortedSet clone() {
        GenomeLocSortedSet ret = new GenomeLocSortedSet(genomeLocParser);
        for (GenomeLoc loc : this.mArray) {
            // ensure a deep copy
            ret.mArray.add(genomeLocParser.createGenomeLoc(loc.getContig(), loc.getStart(), loc.getStop()));
        }
        return ret;
    }

    /**
     * convert this object to a list
     * @return the lists
     */
    public List<GenomeLoc> toList() {
        return this.mArray;
    }

    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append("[");
        for ( GenomeLoc e : this ) {
            s.append(" ");
            s.append(e.toString());
        }
        s.append("]");

        return s.toString();
    }

     /**
     * Check to see whether two genomeLocSortedSets are equal.
     * Note that this implementation ignores the contigInfo object.
     *
     */  /*
    @Override
    public boolean equals(Object other) {
        if(other == null)
            return false;
        if(other instanceof GenomeLocSortedSet) {
            // send to a list, so we can ensure order correct
            List otherList = ((GenomeLocSortedSet)other).toList();
            List thisList = this.toList();
            if (otherList.size() != this.size())
                return false;

            for (Integer i=0;i<thisList.size();i++) {
                if (otherList.get(i).equals(thisList.get(i)))
                    return false;
            }
            return true;
        }
        return false;

    }   */

}
