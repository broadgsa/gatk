package edu.mit.broad.sting.atk.modules;

import edu.mit.broad.sting.atk.LocusWalker;
import edu.mit.broad.sting.atk.LocusIterator;
import edu.mit.broad.sam.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class EmptyLocusWalker implements LocusWalker<Integer, Integer> {
    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(char ref, LocusIterator context) {
        return true;    // We are keeping all the reads
    }

    // Map over the edu.mit.broad.sting.atk.LocusContext
    public Integer map(char ref, LocusIterator context) {
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraveralDone() {
    }
}