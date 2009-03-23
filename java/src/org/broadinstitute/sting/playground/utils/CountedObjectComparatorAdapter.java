package org.broadinstitute.sting.playground.utils;


/** Support class for counted objects. This comparator is an adapter: it is initialized with an arbitrary
 * comparator for objects of type T and can be used to directly compare counted objects of type CountedObject<T>
 * (the underlying Comparator<T> will be used to compare the "object" part of the counted objects, the counter values
 * will be ignored). This comparator also provides additional, non-standard methods that allow direct
 * comparison between a CountedObject<T> and "raw" object of type T (the same underlying Comparator<T> will be used,
 * and the value of the counter in the counted object wil be ignored).
  * @param <T>
 */
public class CountedObjectComparatorAdapter<T>  implements java.util.Comparator<CountedObject<T>> {
	
	private java.util.Comparator<? super T> mComp;

    /** Initializes comparator adapter with a comparator for objects of trype T */
	public CountedObjectComparatorAdapter(java.util.Comparator<? super T> adaptee) {
		mComp = adaptee;
	}

    @Override
	public int compare(CountedObject<T> o1, CountedObject<T> o2) {
		return mComp.compare(o1.getObject(),o2.getObject());
	}
	
	public int compare(T o1, CountedObject<T> o2) {
		return mComp.compare(o1,o2.getObject());
	}

	public int compare(CountedObject<T> o1, T o2) {
		return mComp.compare(o1.getObject(),o2);
	}
	
	@Override
	public boolean equals(Object o) {
		if ( o instanceof CountedObjectComparatorAdapter) {
			if ( ((CountedObjectComparatorAdapter) o).mComp.getClass() == mComp.getClass() ) return true;
		}
		return false;
	}

}
