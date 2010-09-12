package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.io.IOException;

public class IntervalRod extends BasicReferenceOrderedDatum {
    private GenomeLoc loc = null;

    public IntervalRod(String name) {
		super(name);
	}

    public IntervalRod(String name, GenomeLoc loc) {
        super(name);
        this.loc = loc;
    }

    public GenomeLoc getLocation() { return loc; }

	/** Required by ReferenceOrderedDatum interface; this method does nothing (always returns false),
	 * since this rod provides its own iterator for reading underlying data files.
	 */
	public boolean parseLine(Object header, String[] parts) {
		return false; // this rod has its own iterator
	}

	public String repl() {
		throw new GATKException("repl() is not implemented yet");
	}

	public String toSimpleString() { return loc.toString(); }
    public String toString() { return toSimpleString(); }

	public static IntervalRodIterator createIterator(String trackName, File f) throws IOException {
		return IntervalRodIterator.IntervalRodIteratorFromLocsFile(trackName, f);
	}
}