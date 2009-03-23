package org.broadinstitute.sting.indels;

import java.util.List;
import java.util.ArrayList;
import net.sf.samtools.*;

public class SequencePile {
	private List<MSAColumn> mSeqGrid;
	private StringBuffer mRefGrid;
	private int mDepth;
	private List<Boolean> mSeqRC;

	
	public SequencePile(String ref) {
		mRefGrid = new StringBuffer( ref );
		mSeqGrid = new ArrayList<MSAColumn>();
		for ( int i = 0 ; i < mRefGrid.length(); i++ ) {
			mSeqGrid.add(new MSAColumn());
		}
		mDepth = 0;
		mSeqRC = new ArrayList<Boolean>();
	}
	
	/** Adds to the pile nucleotide sequence <seq> that aligns at zero-based position <refpos>
		relative to the original reference stretch the pile is built upon; the detailed alignment
		of the sequence to that reference stretch is specified by the <cigar>. 
		
		@param seq nucleotide sequence
		@param isRC true indicates that RC of the <seq> is being aligned
		@param cigar specification of the alignment of the sequence <seq> to the reference
		@param refpos 0-based position of the alignment with respect to the original stretch of the reference
			that was passed to the pile's constructor. Either <pos> or <pos>+sequence_length can be outside of 
			the pile's boundaries, the SequencePile class will deal with such situations correctly. 
		*/
	public void addAlignedSequence(String seq, boolean isRC, Cigar cigar, int refpos) {

		String alignedSeq = seq ;
//		if ( isRC ) {
//			alignedSeq = ReverseComplement(seq);
//		} else alignedSeq = seq;
		mSeqRC.add(isRC);
		
		int pos = 0; // actual position on the grid (reference can have insertions)
		for ( int i = 0 ; i < refpos ; i++ ) { // i is the position on the original reference
			// if we got some insertions on the reference prior to refpos, we need to take care of them:
			while( mRefGrid.charAt(pos) == '+' ) {
				mSeqGrid.get(pos).add(' ');
				pos++;
			}
			mSeqGrid.get(pos).add(' '); // fill with ' ' to the left of the read
			pos++;
		}
		
		// we reached start position of the alignment

		int readpos = 0; // position on the read
		
		for ( int i = 0 ; i < cigar.numCigarElements() ; i++ ) {

			final CigarElement ce = cigar.getCigarElement(i);
			
			switch(ce.getOperator()) {
    		case I: // read has an insertion
    				for ( int j = 0 ; j < ce.getLength() ; j++ ) {
						if ( pos >= mRefGrid.length() ) break;
						if ( pos >= 0 ) { 
							if ( mRefGrid.charAt(pos) !='+' ) {  // there was no insertion here yet: add it now!
								mRefGrid.insert(pos, '+');
								MSAColumn c = new MSAColumn();
								for ( int k = 0 ; k < mDepth ; k++ ) {
									if ( mSeqGrid.get(pos-1).charAt(k) == ' ') c.add(' ');
									else c.add('*');
								}
								mSeqGrid.add(pos, c);
							}
							try {
								mSeqGrid.get(pos).add(alignedSeq.charAt(readpos));
							} catch (IllegalArgumentException e) {
								throw new IllegalArgumentException(e.getMessage()+": "+seq);
							}
						}
						readpos++;
						pos++;
    				}
    				break;
    		case D: // read has a deletion
    				for ( int j = 0 ; j < ce.getLength() ; j++ ) {
						while( pos < mRefGrid.length() && mRefGrid.charAt(pos) == '+' ) { // skip insertions on the ref
							mSeqGrid.get(pos).add('*');
							pos++;
						}    					
						if ( pos >= mRefGrid.length() ) break;
						mSeqGrid.get(pos).add('-'); // mark deletion
						pos++;
    				}
    				break;
    		case M: 
    				for ( int j = 0 ; j < ce.getLength() ; j++ ) {
						// if ref has an insertion, but the read does not: skip the insertion and continue with "gapless" alignment
						while( pos < mRefGrid.length() && mRefGrid.charAt(pos) == '+' ) {
							mSeqGrid.get(pos).add('*');
							pos++;
						}
						if ( pos >= mRefGrid.length() ) break;
						try {
							mSeqGrid.get(pos).add(alignedSeq.charAt(readpos));
						} catch (IllegalArgumentException e) {
							throw new IllegalArgumentException(e.getMessage()+": "+seq);							
						}
						pos++;
						readpos++;
    				}
    				break; 
    		default : throw new IllegalArgumentException("Unknown cigar element");
			}
		}
		for ( int i = pos ; i < mRefGrid.length() ; i++ ) { // i is the position on the modified reference
			mSeqGrid.get(i).add(' '); // fill with ' ' to the left of the read
		}
		mDepth++;
	}
	
	public String format() {
		StringBuffer b = new StringBuffer();
		b.append("  ");
		b.append(mRefGrid);
		b.append('\n');
		
		try {
		for ( int i = 0 ; i < mDepth; i++ ) {
			if ( mSeqRC.get(i).booleanValue() ) b.append("<-");
			else b.append("->");
			for ( int j = 0 ; j < mRefGrid.length() ; j++) {
				b.append(mSeqGrid.get(j).charAt(i));
			}
			b.append('\n');
		}
		} catch (Exception e) {}
		return b.toString();
	}
	
	private String ReverseComplement(String s) {
		StringBuffer b = new StringBuffer();
		char [] data = s.toCharArray();
		for ( int i = data.length - 1 ; i >= 0 ; i-- ) b.append(BaseComplement(data[i]));
		return b.toString();
	}
	
	private char BaseComplement(char b) {
		switch ( b ) {
		case 'A' : return 'T';
		case 'C': return 'G'; 
		case 'G': return 'C';
		case 'T': return 'A';
		default: throw new IllegalArgumentException(b + " is not a DNA base");
		}
	}

    public void colorprint() { colorprint(false); }

	public void colorprint(boolean printId) {
        if ( printId ) System.out.print("      ");
		else System.out.print("  ");
		System.out.println(mRefGrid);

        StringBuilder sb = new StringBuilder();
        java.util.Formatter frmt = new java.util.Formatter(sb);
        
		try {
		for ( int i = 0 ; i < mDepth; i++ ) {
            if ( printId ) {
                sb.delete(0,sb.length());
                frmt.format("%3d:", i);
                System.out.print(frmt.out().toString());
            }
			if ( mSeqRC.get(i).booleanValue() ) System.out.print("<-");
			else System.out.print("->");
			for ( int j = 0 ; j < mRefGrid.length() ; j++) {
				char seqbase = mSeqGrid.get(j).charAt(i);
				char refbase = mRefGrid.charAt(j);
				if ( isBase(refbase) && isBase(seqbase) && refbase != seqbase ) System.out.print("\033[31m"+seqbase+"\033[30m");
				else System.out.print(seqbase);
			}
			System.out.print('\n');
		}
		} catch (Exception e) {}
	}
	
	private boolean isBase(char b) {
		return ( b=='A' ||b == 'C' || b=='G' || b=='T' || b=='N');
	}
}
