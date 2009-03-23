package org.broadinstitute.sting.indels;

import java.util.List;
import java.util.ArrayList;

public class MSAColumn {
	private List<Character> mBytes;
	
	public MSAColumn() {
		mBytes = new ArrayList<Character>();
	}
	
	/** Adds specified byte to the end of the column */
	public void add(char b) throws IllegalArgumentException {
		if ( b == 'A' || b == 'C' || b == 'T' || b == 'G' || b == '-' || b==' ' || b == '*' || b=='N') {
			mBytes.add(b);
		} else {
			throw new IllegalArgumentException("Invalid base letter passed to MSAColumn");
		}
	}
	
	/** Removes first element from the column */
	public void removeFirst() throws IndexOutOfBoundsException {
		mBytes.remove(0);
	}
	
	/** Removes value at the specified position from the column */ 
	public void remove (int index) throws IndexOutOfBoundsException {
		mBytes.remove(index);
	}
	
	public int size() { return mBytes.size(); }
	
	public Character charAt(int offset) { return mBytes.get(offset); }
}
