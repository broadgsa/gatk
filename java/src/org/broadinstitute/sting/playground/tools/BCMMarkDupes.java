package org.broadinstitute.sting.playground.tools;

import net.sf.samtools.*;
import java.io.*;
import java.util.*;


/**
 * 
 * Marks duplicate reads in a sam or bam file based on start point and strand of the alignment.
 * Requires: java samtools "picard".  java 1.5 or better.
 * 
 * Note, will not unmark duplicates.  The final number of duplicates marked at the end is not necess. the total number of marked duplicates in the final output file.
 * bam or sam type is inferred by the number of arguments, read order, header, filetype is preserved in the output file.
 * 
 * Usage: <sam or bam file> [index file if using bam] <output file>
 * @author bainbrid
 *
 */
public class BCMMarkDupes
{
	final static String inHT = "ht";
	
	/**
	 *  Usage: <sam or bam file> [index file if using bam] <output file>
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		int cnt = 0;
		int dcnt = 0;

        // begin hack by kiran
        
        /*
		if(args.length != 2 && args.length !=3 )
		{
			throw new Exception("Invalid number of arguments.\n"+usage());
		}
		boolean isBam = args.length == 3;
		String inputfile = args[0];
		String index = null;
		if(isBam) index = args[1];
		String outputfile = args[1];
		if(isBam) index = args[2];
		*/

        String inputfile = args[0];
        String index = inputfile + ".bai";
        String outputfile = args[1];
        boolean isBam = true;
        
        // end hack by kiran
		
		//set the reader and writer, and header
		SAMFileReader sfr = null;
		if(isBam)
			sfr = new SAMFileReader(new File(inputfile), new File(index));
		else
			sfr = new SAMFileReader(new File(inputfile));
		Iterator<SAMRecord> iter = sfr.iterator();
		SAMFileWriter fw = null;
		SAMFileHeader header = sfr.getFileHeader();
		if(isBam)
			fw = (new SAMFileWriterFactory()).makeBAMWriter(header, true, new File(outputfile));
		else
			fw = (new SAMFileWriterFactory()).makeSAMWriter(header, true, new File(outputfile));
		Hashtable<String,String> ht = new Hashtable<String,String>(2000000);
		
		//iterate through file and mark dupes
		while(iter.hasNext())
		{
			
			SAMRecord rec = iter.next();
			int start = 0;
			String strand = "N";
			if(rec.getReadNegativeStrandFlag())
			{
				start = rec.getAlignmentEnd();
				//start = rec.getUnclippedEnd();
			}
			else
			{
				strand = "P";
				start = rec.getAlignmentStart();
				//start = rec.getUnclippedStart();
			}
			if(start == 0) continue;
			cnt++;
			String chromo = rec.getReferenceName();
			String key = chromo+strand+start;
			String s = ht.get(key);
			if(s == null)
			{
				ht.put(key, inHT);
				//rec.setDuplicateReadFlag(false);
			}
			else
			{
				rec.setDuplicateReadFlag(true);
				dcnt++;
			}
			fw.addAlignment(rec);
		}
		fw.close();
		System.err.println("Total records:"+cnt);
		System.err.println("Records marked as duplicates:"+dcnt);
		System.err.println("HT:"+ht.size());
		
	}
	
	
	
	
	public static String usage()
	{
		return "usage:\n\t<sam or bam file> [index file if using bam] <output file>";
		
	}

}
