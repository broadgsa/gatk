package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;

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

    // our private storage for the GenomeLoc's
    private List<GenomeLoc> mArray = new ArrayList<GenomeLoc>();

    /** default constructor */
    public GenomeLocSortedSet() {
    }

    public GenomeLocSortedSet(GenomeLoc e) {
        this();
        add(e);
    }

    public GenomeLocSortedSet(Collection<GenomeLoc> l) {
        this();

        for ( GenomeLoc e : l )
            add(e);
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
                throw new StingException("Genome Loc Sorted Set already contains the GenomicLoc " + e.toString());
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

    public GenomeLocSortedSet substractRegions(GenomeLocSortedSet toRemoveSet) {
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
                for ( GenomeLoc newP : subtractRegion(p, e) )
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
                throw new StingException("BUG: unexpected condition: p=" + p + ", e=" + e);
            }

            if ( i++ % 10000 == 0 )
                logger.debug("removeRegions operation: i = " + i);
        }

        return GenomeLocSortedSet.createSetFromList(good);
    }

    private static final List<GenomeLoc> EMPTY_LIST = new ArrayList<GenomeLoc>();
    private List<GenomeLoc> subtractRegion(GenomeLoc g, GenomeLoc e) {
        if (g.equals(e)) {
            return EMPTY_LIST;
        } else if (g.containsP(e)) {
            List<GenomeLoc> l = new ArrayList<GenomeLoc>();

            /**
             * we have to create two new region, one for the before part, one for the after
             * The old region:
             * |----------------- old region (g) -------------|
             *        |----- to delete (e) ------|
             *
             * product (two new regions):
             * |------|  + |--------|
             *
             */
            GenomeLoc before = GenomeLocParser.createGenomeLoc(g.getContigIndex(), g.getStart(), e.getStart() - 1);
            GenomeLoc after = GenomeLocParser.createGenomeLoc(g.getContigIndex(), e.getStop() + 1, g.getStop());
            if (after.getStop() - after.getStart() >= 0) {
                l.add(after);
            }
            if (before.getStop() - before.getStart() >= 0) {
                l.add(before);
            }

            return l;
        } else if (e.containsP(g)) {
            /**
             * e completely contains g, delete g, but keep looking, there may be more regions
             * i.e.:
             *   |--------------------- e --------------------|
             *       |--- g ---|    |---- others ----|
             */
            return EMPTY_LIST;   // don't need to do anything
        } else {
            /**
             * otherwise e overlaps some part of g
             *
             * figure out which region occurs first on the genome.  I.e., is it:
             * |------------- g ----------|
             *       |------------- e ----------|
             *
             * or:
             *       |------------- g ----------|
             * |------------ e -----------|
             *
             */

            GenomeLoc n;
            if (e.getStart() < g.getStart()) {
                n = GenomeLocParser.createGenomeLoc(g.getContigIndex(), e.getStop() + 1, g.getStop());
            } else {
                n = GenomeLocParser.createGenomeLoc(g.getContigIndex(), g.getStart(), e.getStart() - 1);
            }

            // replace g with the new region
            return Arrays.asList(n);
        }
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
        GenomeLocSortedSet returnSortedSet = new GenomeLocSortedSet();
        for (SAMSequenceRecord record : dict.getSequences()) {
            returnSortedSet.add(GenomeLocParser.createGenomeLoc(record.getSequenceIndex(), 1, record.getSequenceLength()));
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
    public static GenomeLocSortedSet createSetFromList(List<GenomeLoc> locs) {
        GenomeLocSortedSet set = new GenomeLocSortedSet();
        set.addAll(locs);
        return set;
    }


    /**
     * return a deep copy of this collection.
     *
     * @return a new GenomeLocSortedSet, indentical to the current GenomeLocSortedSet.
     */
    public GenomeLocSortedSet clone() {
        GenomeLocSortedSet ret = new GenomeLocSortedSet();
        for (GenomeLoc loc : this.mArray) {
            // ensure a deep copy
            ret.mArray.add(GenomeLocParser.createGenomeLoc(loc.getContigIndex(), loc.getStart(), loc.getStop()));
        }
        return ret;
    }


    public boolean addAllRegions(List<GenomeLoc> locations) {
        this.mArray.addAll(locations);
        Collections.sort(this.mArray);
        this.mArray = GenomeLocSortedSet.mergeOverlappingLocations(this.mArray);
        return true;
    }

/**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     *
     * @return the list of merged locations
     */
    public static List<GenomeLoc> mergeOverlappingLocations(final List<GenomeLoc> raw) {
        logger.debug("  Raw locations are: " + Utils.join(", ", raw));
        if (raw.size() <= 1)
            return raw;
        else {
            ArrayList<GenomeLoc> merged = new ArrayList<GenomeLoc>();
            Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while (it.hasNext()) {
                GenomeLoc curr = it.next();
                if (prev.contiguousP(curr)) {
                    prev = prev.merge(curr);
                } else {
                    merged.add(prev);
                    prev = curr;
                }
            }
            merged.add(prev);
            return merged;
        }
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
}
