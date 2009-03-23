package org.broadinstitute.sting.playground.indels;

import org.broadinstitute.sting.playground.utils.Interval;


public class IntervalComparator implements java.util.Comparator<Interval> {
	public int compare(Interval r1, Interval r2) {
		if ( r1.getStart() < r2.getStart() ) return -1;
		if ( r1.getStart() == r2.getStart() ) {
			if ( r1.getStop() < r2.getStop() ) return -1;
			if ( r1.getStop() == r2.getStop() ) return 0;
		}
		return 1;
	}
	
	@Override
	public boolean equals(Object o) {
		if ( o instanceof IntervalComparator) return true;
		return false;
	}

}
