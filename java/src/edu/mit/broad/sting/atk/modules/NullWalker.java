package edu.mit.broad.sting.atk.modules;

import edu.mit.broad.sting.atk.LocusWalker;
import edu.mit.broad.sting.atk.LocusIterator;
import edu.mit.broad.sting.utils.ReferenceOrderedDatum;
import edu.mit.broad.sting.utils.rodDbSNP;
import edu.mit.broad.sting.utils.Utils;
import net.sf.samtools.SAMRecord;

import java.util.List;

// Null traversal. For ATK performance measuring.
// j.maguire 3-7-2009

public class NullWalker implements LocusWalker<Integer, Integer> {
    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusIterator context) {
        return true;    // We are keeping all the reads
    }

    // Map over the edu.mit.broad.sting.atk.LocusContext
    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusIterator context) 
    {
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() 
    {
        return 0; 
    }
    public Integer reduce(Integer value, Integer sum) 
    {
        return 0;
    }

    public void onTraveralDone() {
    }
}
