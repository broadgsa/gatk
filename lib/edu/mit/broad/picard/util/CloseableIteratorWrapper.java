/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.util;

import java.util.Iterator;

import edu.mit.broad.sam.util.CloseableIterator;

public class CloseableIteratorWrapper<T> implements CloseableIterator<T> {
	Iterator<T> wrappedIterator;

	public CloseableIteratorWrapper(Iterator<T> wrappedIterator) {
		this.wrappedIterator = wrappedIterator;
	}

	@Override
	public boolean hasNext() {
		return wrappedIterator.hasNext();
	}

	@Override
	public T next() {
		return wrappedIterator.next();
	}

	@Override
	public void remove() {
		wrappedIterator.remove();
	}

	@Override
	public void close() {
	}
}