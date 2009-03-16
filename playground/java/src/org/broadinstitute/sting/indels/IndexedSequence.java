package org.broadinstitute.sting.indels;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

public class IndexedSequence implements Iterable<Map.Entry<Short,List<Integer>>> {
	private Map<Short,List<Integer> > m_locs;
	private String m_seq;
	private int m_K;
	
	public IndexedSequence(String seq, int K) {
		m_locs = new Hashtable<Short,List<Integer> >();
		m_seq = new String(seq);
		m_K = K;
		Iterator<Short> iter = new KmerCodeIterator(seq, K);
		int offset = 0;
		while ( iter.hasNext() ) {
			Short k = iter.next();
			List<Integer> offset_list = m_locs.get(k);
			if ( offset_list == null ) { 
				offset_list = new ArrayList<Integer>();
				m_locs.put(k,offset_list);
			}
			offset_list.add(offset++);
		}
	}
	
	public Iterator<Entry<Short, List<Integer>>> iterator() {
		return m_locs.entrySet().iterator();
	}

	public int length() { return m_seq.length(); }
	
	public List<Integer> getOffsets(short k) { return m_locs.get(k); }
	
	String getSequence() { return m_seq; }
	
	public int getK() { return m_K; }
	
	public static void testMe() {
		String s = "ACCGTGCGGGCACCTGC";
		int K = 3;
		IndexedSequence is = new IndexedSequence(s,K);
		System.out.println("Sequence: "+ s);
		System.out.print("          ");
		for ( int i= 0 ; i < s.length() ; i++ ) {
			if ( i % 10 == 0 ) System.out.print(i/10);
			else System.out.print(' ');
		}
		System.out.println();
		System.out.print("          ");
		for ( int i= 0 ; i < s.length() ; i++ ) System.out.print(i%10);
		System.out.println();
		System.out.println();
	
		System.out.println("Indexing with K="+K+":");
		
		Set< Map.Entry<Short, List<Integer> > > data  = is.m_locs.entrySet(); 
		for ( Map.Entry<Short, List<Integer> > e : data ) {
//			System.out.print("("+e.getKey()+") ");
			System.out.print(kmerToString(e.getKey().shortValue(),K));
			System.out.print("-->");
			for ( Integer offset : e.getValue() ) {
				System.out.print(" "+offset.toString());
			}
			System.out.println();
		}
	}
	
	private static String kmerToString(short code, int K) {
		StringBuffer b = new StringBuffer(K);
		for ( int i = 0 ; i < K ; i++ ) {
			char c='N';
			switch( code & 0x3 ) {
			case 0 : c = 'A' ; break; 
			case 1 : c = 'C' ; break; 
			case 2 : c = 'G' ; break; 
			case 3 : c = 'T' ; break; 
			}
			b.append(c);
			code >>= 2;
		}
		return b.reverse().toString(); 
	}
	
	public static void main(String argv[]) {
		testMe();
	}

}
