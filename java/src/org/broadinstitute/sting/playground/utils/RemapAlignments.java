package org.broadinstitute.sting.playground.utils;

import java.io.File;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.broadinstitute.sting.utils.GenomeLocParser;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;


public class RemapAlignments extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Investigates indels called in the alignment data\n";
	@Option(shortName="M", 
			doc="Map file: from the reference the reads were aligned to, to the master reference the alignments should be remapped to", 
			optional=false)
	public File MAP_FILE = null;
	@Option(shortName="I", 
			doc="Input file (bam or sam) with alignments to be remapped", 
			optional=false)
	public File IN = null;
	@Option(shortName="O", 
			doc="File to write remapped reads to.", 
			optional=false)
	public File OUT = null;
	@Option(shortName="R", 
			doc="Reference to remap alignments onto.", 
			optional=false)
	public File REFERENCE = null;
	@Option(
			doc="If a read has multiple alignments that are exactly the same after remapping, "+
			"keep only one copy of such alignment in output file. Multiple alignments that are "+
			"not equivalent after remaping are not affected by this flag. "+
			"Multiple alignments for the same query must be grouped on adjacent lines of the input file to be detected," +
			"otherwise REDUCE will have no effect.", 
			optional=true)
	public boolean REDUCE = false;


	private GenomicMap map = null;
	private String lastReadName = null;
	private int totalReads = 0;
	private int totalRecords = 0;
	private int badRecords = 0;
	private int totalUnmappedReads = 0;
	private int writtenRecords = 0;

	private Set<SAMRecord> remappedReads = null;
	private SAMFileWriter writer = null;
	private SAMFileReader reader = null;
	
	private static int [] g_log_n; // copied from bwa
	
	
    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new RemapAlignments().instanceMain(argv));
    }
    
    protected int doWork() {
    			
    	g_log_n = new int[256];
    	for (int i = 1; i < 256; ++i) g_log_n[i] = (int)(4.343 * Math.log(i) + 0.5);
    	
    	reader = new SAMFileReader(IN);
    	reader.setValidationStringency(ValidationStringency.SILENT);
		SAMFileHeader oldHeader = reader.getFileHeader();
		if ( oldHeader == null ) throw new RuntimeException("Failed to retrieve SAM file header from the input bam file");
		
		if ( REDUCE && oldHeader.getSortOrder() != SortOrder.queryname ) 
			System.out.println("WARNING: Input file is not sorted by query name, REDUCE may have no effect. Sort order: "
					+oldHeader.getSortOrder());
		
		remappedReads = new TreeSet<SAMRecord>(new AlignmentComparator());
		
		SAMFileHeader h = new SAMFileHeader();
		
		for ( Entry<String, Object> attr : oldHeader.getAttributes() ) h.setAttribute(attr.getKey(), attr.getValue());
		h.setGroupOrder(oldHeader.getGroupOrder());
		h.setProgramRecords(oldHeader.getProgramRecords());
		h.setReadGroups(oldHeader.getReadGroups());
		
		if ( oldHeader.getSortOrder() == SortOrder.queryname ) {
			h.setSortOrder(SortOrder.queryname);
		} else {
			h.setSortOrder(SortOrder.unsorted);
		}
		
		ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(REFERENCE);

        if ( reference.getSequenceDictionary() == null ) {
        	System.out.println("No reference sequence dictionary found. Aborting.");
        	reader.close();
        	System.exit(1);
        }
		
		h.setSequenceDictionary(reference.getSequenceDictionary());
		GenomeLocParser.setupRefContigOrdering(reference.getSequenceDictionary());
		
		map = new GenomicMap(10000);
		map.read(MAP_FILE);
		System.out.println("Map loaded successfully: "+map.size()+" contigs");
				
		
		writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, OUT);
		
		for ( SAMRecord read : reader ) {
			
			
			if ( map.remapToMasterReference(read,h,true) == null ) {
				badRecords++;
				continue;
			}
			if ( read.getReadUnmappedFlag() ) totalUnmappedReads++;
			else {
				// destroy mate pair mapping information, if any (we will need to reconstitute pairs after remapping both ends):
				read.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
				read.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
//				if ( read.getReadPairedFlag() ) System.out.println("PAIRED READ!!");
			}
			totalRecords++;
			
			if ( totalRecords % 1000000 == 0 ) System.out.println(totalRecords + " valid records processed");
			

			if ( ! read.getReadName().equals(lastReadName) ) {
				totalReads++;
				lastReadName = read.getReadName();
			
						
				if ( REDUCE ) {
					
					updateCountsAndQuals(remappedReads);
					
					for ( SAMRecord r : remappedReads ) {
						writer.addAlignment(r); // emit non-redundant alignments for previous query
						writtenRecords++;
					}
					remappedReads.clear(); 
				}
			} 
			if ( REDUCE ) remappedReads.add(read); 
			else {
				writer.addAlignment(read);
				writtenRecords++;
			}
		}

		// write remaining bunch of reads:
		if ( REDUCE ) {
			updateCountsAndQuals(remappedReads);
			for ( SAMRecord r : remappedReads ) {
				writer.addAlignment(r); // emit non-redundant alignments for previous query
				writtenRecords++;
			}
		}
		
		System.out.println("Total valid records processed: "+totalRecords);
		System.out.println("Incorrect records (alignments across contig boundary) detected: "+badRecords + 
				" (discarded and excluded from any other stats)");
		System.out.println("Total reads processed: "+totalReads);
		System.out.println("Total mapped reads: "+(totalReads-totalUnmappedReads));
		System.out.println("Average hits per mapped read: "+((double)(totalRecords-totalUnmappedReads))/(totalReads-totalUnmappedReads));
		System.out.println("Records written: "+writtenRecords);
		System.out.println("Average hits per mapped read written (after reduction): "
				+((double)(writtenRecords-totalUnmappedReads))/(totalReads-totalUnmappedReads));
		reader.close();
		writer.close();
		return 0;
	}
	
    class AlignmentComparator implements Comparator<SAMRecord> {

    	@Override
    	public int compare(SAMRecord r1, SAMRecord r2) {
    		if ( r1.getReferenceIndex() < r2.getReferenceIndex() ) return -1; 
    		if ( r1.getReferenceIndex() > r2.getReferenceIndex() ) return  1;
    		if ( r1.getAlignmentStart() < r2.getAlignmentStart() ) return -1;
    		if ( r1.getAlignmentStart() > r2.getAlignmentStart() ) return 1;
    		return r1.getCigarString().compareTo(r2.getCigarString());
    	}
    	
    }

    private void updateCountsAndQuals(Set<SAMRecord> reads) {
    	if ( reads.size() == 1 ) {
    		SAMRecord r = reads.iterator().next();
 
        	// technically, if edit distance of the read is equal to max_diff used in alignments, 
        	// we should have set 25... 
    		r.setMappingQuality(37);
    		r.setAttribute("X0", new Integer(1));
    		r.setAttribute("X1", new Integer(0));
    		r.setNotPrimaryAlignmentFlag(false);
    		
    	} else {
    		
    		// we have multiple alignments for the read
    		// need to figure out how many best vs inferior alignments are there:
    		int minNM = 1000000;
    		int cnt = 0; // count of best alignments
    		for ( SAMRecord r : reads ) {
    			int nm = (Integer)r.getAttribute("NM"); 
    			if ( nm < minNM  ) { 
    				minNM = nm;
    				cnt = 1;
    			} else if ( nm == minNM ) cnt++;
    		}
    		// now reset counts and quals:
    		for ( SAMRecord r : reads ) {
    			
    			int cnt2 = reads.size() - cnt; // count of inferior alignments
    			
    	   		r.setAttribute("X0", new Integer(cnt));
        		r.setAttribute("X1", new Integer(cnt2));
    			
        		if ( cnt2 > 255 ) cnt2 = 255; // otherwise we will be out of bounds in g_log_n
        		
    			if ( ((Integer)r.getAttribute("NM")).intValue() == minNM ) { 
    				
    				// one of the best alignments:

    				r.setNotPrimaryAlignmentFlag(false);
    				if ( cnt == 1 ) {    					
    					// single best alignment; additional inferior alignments will only affect mapping qual
    					r.setMappingQuality( 23 < g_log_n[cnt2] ? 0 : 23 - g_log_n[cnt2] ); // this recipe for Q is copied from bwa
    				} else {
    					r.setMappingQuality(0); // multiple best alignments - mapping quality is 0
    				}
    			} else {
    				
    				// secondary alignment ( we know we hold a better one)
    				r.setNotPrimaryAlignmentFlag(true);
    				r.setMappingQuality(0); // ??? should we set 0 for secondary??
    			}
    		}
    	}
    	
    }
    
/*    
    private int bwa_approx_mapQ(SAMRecord r, int max_diff) {
    	int c1 = (Integer)r.getAttribute("X0");
    	int c2 = (Integer)r.getAttribute("X1");
    	int mm = (Integer)r.getAttribute("NM");
    	if ( c1 > 0 ) return 0;
    	if ( c1 == 0 ) return 23;
    	if ( mm == max_diff ) return 25;
    	return 0;
    }
*/
}


