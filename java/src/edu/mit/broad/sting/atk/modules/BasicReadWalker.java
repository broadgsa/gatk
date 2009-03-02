package edu.mit.broad.sting.atk.modules;

import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sting.atk.LocusContext;
import edu.mit.broad.sting.atk.ReadWalker;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class BasicReadWalker<MapType, ReduceType> implements ReadWalker<MapType, ReduceType> {
    public void initialize() { }
    public String walkerType() { return "ByRead"; }

    public boolean filter(LocusContext context, SAMRecord read) {
        // We are keeping all the reads
        return true;
    }

    public void onTraveralDone() {

    }

    // Three basic abstract function that *must* be overridden
    public abstract MapType map(LocusContext context, SAMRecord read);
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);
}