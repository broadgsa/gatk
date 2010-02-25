package org.broadinstitute.sting.gatk.refdata.utils;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class FlashBackIterator
 *
 * better than acid washed jeans...more like a Delorean that flies through time
 *
 * This iterator buffers a certain amount of ROD data to 'flash back' to.  This
 * is needed for using ROD's in read traversals, because between shards we sometimes
 * (actually often) need to go back to before the current iterators location and
 * get RODs that overlap the current read.
 */
public class FlashBackIterator implements LocationAwareSeekableRODIterator {
    private LocationAwareSeekableRODIterator iterator;
    private LinkedList<ComparableList> list = new LinkedList<ComparableList>();
    private int MAX_QUEUE = 5000;
    private boolean usingQueue = false;

    public FlashBackIterator(LocationAwareSeekableRODIterator iterator) {
        this.iterator = iterator;
    }

    @Override
    public GenomeLoc peekNextLocation() {
        return iterator.peekNextLocation();
    }

    @Override
    public GenomeLoc position() {
        return (usingQueue) ? list.getFirst().getLocation() : iterator.position();
    }

    @Override
    public RODRecordList seekForward(GenomeLoc interval) {
        RODRecordList lt = iterator.seekForward(interval);
        if (lt != null) list.addLast(new ComparableList(lt));
        return lt;
    }

    @Override
    public boolean hasNext() {
        if (usingQueue) return (list.size() > 0 || iterator.hasNext());
        return iterator.hasNext();
    }

    @Override
    public RODRecordList next() {
        RODRecordList ret;
        if (!usingQueue || list.size() < 1) {
            usingQueue = false;
            ret = iterator.next();
            list.addLast(new ComparableList(ret));
            if (list.size() > MAX_QUEUE) list.removeFirst();
        } else {
            ret = list.getFirst().getList();
            list.removeFirst();
        }
        return ret;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("We don't support remove");
    }

    /**
     * can we flash back to the specified location?
     *
     * @param location the location to try and flash back to
     *
     * @return true if we can, false otherwise
     */
    public boolean canFlashBackTo(GenomeLoc location) {
        GenomeLoc farthestBack = (list.size() > 0) ? list.getFirst().getLocation() : iterator.peekNextLocation();
        System.err.println("farthestBack = " + farthestBack + " loc = " + location);
        return (!farthestBack.isPast(location));
    }

    /**
     * flashback! Throws an unsupported operation exception
     *
     * @param location where to flash back to
     */
    public void flashBackTo(GenomeLoc location) {
        if (!canFlashBackTo(location)) throw new UnsupportedOperationException("we can't flash back to " + location);
        if (list.size() > 0 && !list.getLast().getLocation().isBefore(location))
            usingQueue = true;
    }
}

class ComparableList implements Comparator<ComparableList> {
    private RODRecordList list;
    private GenomeLoc location = null;
    public ComparableList(RODRecordList list) {
        this.list = list;
        if (list != null && list.size() != 0) location = list.get(0).getLocation();
        else throw new IllegalStateException("Bad voodoo!");
    }

    @Override
    public int compare(ComparableList list1, ComparableList list2) {
        if (list1.location == null && list2.location == null)
            return 0;
        if (list1.location == null) return 1;
        if (list2.location == null) return -1;
        return (list1.location.compareTo(list2.location));
    }

    public GenomeLoc getLocation() {
        return location;
    }

    public RODRecordList getList() {
        return list;
    }
}