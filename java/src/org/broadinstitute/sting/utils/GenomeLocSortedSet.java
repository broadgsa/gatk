package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
    // our private storage for the GenomeLoc's
    private final ArrayList<GenomeLoc> mArray = new ArrayList<GenomeLoc>();

    /** default constructor */
    public GenomeLocSortedSet() {
    }

    /**
     * get an iterator over this collection
     *
     * @return
     */
    public Iterator<GenomeLoc> iterator() {
        return mArray.iterator();
    }

    /**
     * return the size of the collection
     *
     * @return
     */
    public int size() {
        return mArray.size();
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
        if (mArray.contains(e)) {
            throw new IllegalArgumentException("attempting to add a duplicate object to the set");
        }
        int index = 0;
        while (index < mArray.size()) {
            if (!e.isPast(mArray.get(index))) {
                mArray.add(index, e);
                return true;
            }
            ++index;
        }
        this.mArray.add(e);
        return true;
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

    /**
     * remove an element from the set.  Given a specific genome location, this function will
     * remove all regions in the element set that overlap the specified region.
     *
     * @param e the genomic range to remove
     *
     * @return true if a removal action was performed, false if the collection was unchanged.
     */
    public boolean removeRegion(GenomeLoc e) {
        if (e == null) {
            return false;
        }

        // sometimes we can't return right away, this holds the value for those cases
        boolean returnValue = false;
        /**
         * check if the specified element overlaps any current locations, subtract the removed
         * region and reinsert what is left.
         */
        for (GenomeLoc g : mArray) {
            if (g.overlapsP(e)) {
                if (g.equals(e)) {
                    mArray.remove(mArray.indexOf(g));
                    return true;
                } else if (g.containsP(e)) {
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
                    GenomeLoc before = new GenomeLoc(g.getContigIndex(), g.getStart(), e.getStart() - 1);
                    GenomeLoc after = new GenomeLoc(g.getContigIndex(), e.getStop() + 1, g.getStop());
                    int index = mArray.indexOf(g);
                    if (after.getStop() - after.getStart() > 0) {
                        mArray.add(index, after);
                    }
                    if (before.getStop() - before.getStart() > 0) {
                        mArray.add(index, before);
                    }
                    mArray.remove(mArray.indexOf(g));
                    return true;
                } else if (e.containsP(g)) {
                    /**
                     * e completely contains g, delete g, but keep looking, there may be more regions
                     * i.e.:
                     *   |--------------------- e --------------------|
                     *       |--- g ---|    |---- others ----|
                     */
                    mArray.remove(mArray.indexOf(g));
                    returnValue = true;
                } else {

                    /**
                     * otherwise e overlaps some part of g
                     */
                    GenomeLoc l;

                    /**
                     * figure out which region occurs first on the genome.  I.e., is it:
                     * |------------- g ----------|
                     *       |------------- e ----------|
                     *
                     * or:
                     *       |------------- g ----------|
                     * |------------ e -----------|
                     *
                     */

                    if (e.getStart() < g.getStart()) {
                        l = new GenomeLoc(g.getContigIndex(), e.getStop() + 1, g.getStop());
                    } else {
                        l = new GenomeLoc(g.getContigIndex(), g.getStart(), e.getStart() - 1);
                    }
                    // replace g with the new region
                    mArray.set(mArray.indexOf(g), l);
                    returnValue = true;
                }
            }
        }
        return returnValue;
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
            returnSortedSet.add(new GenomeLoc(record.getSequenceIndex(), 1, record.getSequenceLength()));
        }
        return returnSortedSet;
    }

    /**
     * Create a sorted genome location set from a list of GenomeLocs.
     * @param locs the list<GenomeLoc>
     * @return the sorted genome loc list
     */
    public static GenomeLocSortedSet createSetFromList(List<GenomeLoc> locs) {
        GenomeLocSortedSet set = new GenomeLocSortedSet();
        for (GenomeLoc l: locs) {
            set.add(l);
        }
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
            ret.mArray.add(new GenomeLoc(loc.getContigIndex(), loc.getStart(), loc.getStop()));
        }
        return ret;
    }

}
