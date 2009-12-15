package org.broadinstitute.sting.gatk.refdata;

import java.io.File;
import java.util.*;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.GenomeLoc;

public class IntervalRodIterator implements Iterator<IntervalRod> {
    private List<GenomeLoc> locations = null;
    private Iterator<GenomeLoc> iter;
    private String trackName = null;

    //public static RODIterator<IntervalRod> IntervalRodIteratorFromLocs(List<GenomeLoc> locs) {
    //    IntervalRodIterator it = new IntervalRodIterator("interval", locs);
    //    return new RODIterator<IntervalRod>(it);
    //}

    public static IntervalRodIterator IntervalRodIteratorFromLocsFile(final String trackName, final File file) {
        //System.out.printf("Parsing %s for intervals %s%n", file, trackName);
        List<GenomeLoc> locs = GenomeAnalysisEngine.parseIntervalRegion(Collections.singletonList(file.getPath()));
        //System.out.printf(" => got %d entries %n", locs.size());
        return new IntervalRodIterator(trackName, locs);
    }

	public IntervalRodIterator(String trackName, List<GenomeLoc> locs) {
        this.trackName = trackName;
        locations = locs;
        iter = locations.iterator();
	}

	@Override
	public boolean hasNext() {
		return iter.hasNext();
	}

    /**
     * @return the next element in the iteration.
     * @throws NoSuchElementException - iterator has no more elements.
     */
	@Override
	public IntervalRod next() {
        if (!this.hasNext()) throw new NoSuchElementException("IntervalRodIterator next called on iterator with no more elements");
        IntervalRod r = new IntervalRod(trackName, iter.next());
        //System.out.printf("IntervalRod next is %s%n", r);
        return r;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}