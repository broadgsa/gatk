package org.broadinstitute.sting.indels;

import java.util.Iterator;

public class KmerCodeIterator implements Iterator<Short> {
	private byte[] m_seq;
	private int m_last_offset; // offset of the last base to be added to current code when the next() shift occurs
	private short m_code;
	private short m_mask; // used to mask out all bits in m_code except the lowest 2*K bits used for kmers
	
	public KmerCodeIterator(String s, int K) {
		assert K <= 8 : "Currently only Kmers of length K <=8 are supported";
		assert K > 0 : "Kmer length must be positive";
		m_seq = s.getBytes();
		m_last_offset = K-1;
		m_mask = 0; 
		for ( int i = 0 ; i < K ; i++ ) {
			m_mask <<= 2; 
			m_mask |= 0x03;
		}
		if ( K <= m_seq.length ) m_code = kmerCode(m_seq, 0, m_last_offset );
		// m_code now contains first K-1 bases encoded (the last, K-th base will be added when next() is called
	}
	
	@Override
	public boolean hasNext()  {
		return m_last_offset <  m_seq.length ;
	}
		
	@Override
	public Short next() {
		m_code <<= 2;
		m_code |= toBaseCode(m_seq[m_last_offset]);
		m_code &= m_mask;
		m_last_offset++;
		return m_code;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

	/**
	 * Converts base letter (character passed as byte: 'A', 'C', 'T', or 'G', case-insensitive) into 
	 * numerical code according to A=0, C=1, G=2, T=3 (i.e. two bits per nucleotide)
	 * @param b
	 * @return
	 */
	private byte toBaseCode(byte b) {
		// the following transformation is performed in the lines below:
		// A,a->0, C,c->1, G,g->3, T,t->2; (we rely on the ASCII codes here!)
		b >>= 1;
		b &= 0x03;
		// to get conventional base codes (A=0, C=1, G=2, T=3), we need to flip 3 (G) and 2 (T). 
		// In order to do that we xor the lowest bit
		// with the second one (the latter is 1 for 2 and 3, but 0 for 0 and 1, so A anC will not be affected)
		b  ^= (b >> 1);
		return  b;
	}
	
	/** Returns compact code that uniquely represents nucleotide sequence specified
	 *  as an array of character codes (bytes) such as {'A','C','G','G','T',...}. Case insensitive.
	 *  Sequence can not be longer than 8 nucleotides. 
	 * @param s Nucleotide sequence
	 * @return unique code
	 */
//	private short kmerCode(byte [] s) { return kmerCode(s, 0, s.length); }
	
	/** Returns compact code that uniquely represents nucleotide sequence found in the interval
	 * [start,stop) of the specified  array of character codes (bytes) such as {'A','C','G','G','T',...}. Case insensitive.
	 *  Sequence can not be longer than 8 nucleotides. 
	 * @param s Nucleotide sequence
	 * @param start index of the start position of the sub-sequence to be encoded
	 * @param stop next position after the last element of the sub-sequence to be encoded
	 * @return unique code
	 */
	private short kmerCode(byte [] s, int start, int stop) {
		assert start <= stop : "Start position of the subsequence can not be greater than the stop position";
		assert start >= 0 : "Negative subsequence start positions are not allowed";
		assert stop <= s.length : "Stop position of the subsequence can not extend beyond the sequence end";
		short code = 0;
		for ( int i = start ; i < stop ; i++ ) {
			code  <<= 2;
			code |= toBaseCode(s[i]);
		}
		return code;
	}


}
