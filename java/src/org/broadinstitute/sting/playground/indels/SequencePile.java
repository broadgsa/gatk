package org.broadinstitute.sting.playground.indels;

import java.util.List;
import java.util.ArrayList;
import net.sf.samtools.*;

public class SequencePile {
	private List<MSAColumn> mSeqGrid;
	private StringBuilder mRefGrid;
    private StringBuilder headerGrid;
	private int mDepth;
	private List<Boolean> mSeqRC;

	
	public SequencePile(String ref) {
		mRefGrid = new StringBuilder( ref );
        headerGrid = new StringBuilder();
        for ( int i = 0; i < ref.length(); i++ ) headerGrid.append(' ');
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

        // will hold actual position on the grid; reference can have insertions on the grid,
        // so position on the grid where we should start placing the read is not refpos!
		int pos = 0;
		for ( int i = 0 ; i < refpos ; i++ ) { // i is the position on the original reference
			// if we got some insertions on the reference prior to refpos, we need to count them in:
			while( mRefGrid.charAt(pos) == '+' ) {
				mSeqGrid.get(pos).add(' '); // add additional spaces in the line that will hold sequence seq
				pos++;
			}
			mSeqGrid.get(pos).add(' '); // fill with ' ' to the left of the read
			pos++;
		}
		
		// we reached start position of the alignment on the reference grid

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
                                headerGrid.insert(pos,'+');
								MSAColumn c = new MSAColumn();
                                // reads up to the previous depth (prior to adding current read) did not
                                // have an insertion here, so we insert '*' into all of them:
								for ( int k = 0 ; k < mDepth ; k++ ) {
									if ( mSeqGrid.get(pos-1).charAt(k) == ' ') c.add(' ');
									else c.add('*');
								}
								mSeqGrid.add(pos, c); // finally, add the base from the current read
							}
							mSeqGrid.get(pos).add(alignedSeq.charAt(readpos));
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
                        headerGrid.setCharAt(pos,'-');
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
						mSeqGrid.get(pos).add(alignedSeq.charAt(readpos));
                        if ( Character.toUpperCase(alignedSeq.charAt(readpos)) !=
                                Character.toUpperCase(mRefGrid.charAt(pos))
                                && headerGrid.charAt(pos)== ' ') headerGrid.setCharAt(pos,'*');
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

    public void dotprint(boolean printId) {

        String skip = null;
        if ( printId ) skip = new String("     ");
        else skip = new String("  ");

        System.out.print(formatHeader(skip));
        System.out.print(skip);
        System.out.println(mRefGrid);

        try {
        for ( int i = 0 ; i < mDepth; i++ ) {
            if ( printId ) System.out.printf("%3d",i);
            if ( mSeqRC.get(i).booleanValue() ) System.out.print("<-");
            else System.out.print("->");
            for ( int j = 0 ; j < mRefGrid.length() ; j++) {
                char seqbase = mSeqGrid.get(j).charAt(i);
                char refbase = mRefGrid.charAt(j);
                if ( isBase(refbase) && isBase(seqbase) &&
                        Character.toUpperCase(refbase) ==
                        Character.toUpperCase(seqbase) ) {
                    if ( mSeqRC.get(i) ) System.out.print(',');
                    else System.out.print('.');
                }
                else System.out.print(seqbase);
            }
            System.out.print('\n');
        }
        } catch (Exception e) {}
    }


	public void colorprint(boolean printId) {

        String skip = null;
        if ( printId ) skip = new String("     ");
		else skip = new String("  ");

        System.out.print(formatHeader(skip));
        System.out.print(skip);
        System.out.println(mRefGrid);
 
		try {
		for ( int i = 0 ; i < mDepth; i++ ) {
            if ( printId ) System.out.printf("%3d",i);
			if ( mSeqRC.get(i).booleanValue() ) System.out.print("<-");
			else System.out.print("->");
			for ( int j = 0 ; j < mRefGrid.length() ; j++) {
				char seqbase = mSeqGrid.get(j).charAt(i);
				char refbase = mRefGrid.charAt(j);
				if ( isBase(refbase) && isBase(seqbase) &&
                        Character.toUpperCase(refbase) !=
                                Character.toUpperCase(seqbase) ) System.out.print("\033[31m"+seqbase+"\033[30m");
				else System.out.print(seqbase);
			}
			System.out.print('\n');
		}
		} catch (Exception e) {}
	}

    private String formatHeader(String leadString) {
        char [][] mm_strings = new char[2][mRefGrid.length()];
        for ( int i = 0 ; i < mRefGrid.length() ; i++ ) {
            int count = 0;
            char refC = mRefGrid.charAt(i);
            MSAColumn col = mSeqGrid.get(i);
            if ( refC == '+' ) {
                 // count number of observations for insertion
                for ( int j = 0 ; j < col.size() ; j++ ) {
                     if ( col.charAt(j) != '*' && col.charAt(j) != ' ') count++;
                }
            } else {
                if ( headerGrid.charAt(i) == '-' ) {
                    // count number of observations for deletion
                    for ( int j = 0 ; j < col.size() ; j++ ) {
                         if ( col.charAt(j) == '-' ) count++;
                    }
                } else {
                    if ( headerGrid.charAt(i) == '*') {
                        for ( int j = 0 ; j < col.size() ; j++ ) {
                             if ( col.charAt(j)!=' ' &&
                                     Character.toUpperCase(col.charAt(j)) !=
                                     Character.toUpperCase(refC) ) count++;
                        }
                    }
                }
            }
            if ( count > 9 ) mm_strings[0][i] = Character.forDigit(count/10,10);
            else mm_strings[0][i] = ' ';
            if ( count > 0 ) mm_strings[1][i] = Character.forDigit(count%10,10);
            else mm_strings[1][i] = ' ';
        }

        StringBuilder b = new StringBuilder();
        b.append(leadString);
        b.append(mm_strings[0]);
        b.append('\n');
        b.append(leadString);
        b.append(mm_strings[1]);
        b.append('\n');
        b.append(leadString);
        b.append(headerGrid);
        b.append('\n');
        return b.toString();
    }

	private boolean isBase(char b) {
        b = Character.toUpperCase(b);
		return ( b=='A' ||b == 'C' || b=='G' || b=='T' || b=='N');
	}
}
