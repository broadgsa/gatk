package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;


import java.io.*;
import java.util.*;

import net.sf.picard.PicardException;
import net.sf.picard.util.Interval;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

//import org.apache.commons.math.optimization.*;
//import org.apache.commons.math.optimization.direct.*;
//import org.apache.commons.math.analysis.MultivariateRealFunction;

// Program for frequency-specific  VCF-files.


/**
 * @author jmaguire
 */


class VCFMerge extends CommandLineProgram 
{

		@Argument(fullName = "vcf1", shortName = "vcf1", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "vcf2", shortName = "vcf2", doc = "file to open", required = true) public String filename2;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = true) public String output_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print way too much debugging output", required = false) public Boolean verbose = false;

		@Override
		protected int execute() 
		{
			VCFReader reader1;
			VCFReader reader2;

			if (autocorrect) 
			{ 
				reader1 = new VCFReader(VCFHomogenizer.create(filename1)); 
				reader2 = new VCFReader(VCFHomogenizer.create(filename2)); 
			}
			else 
			{ 
				reader1 = new VCFReader(new File(filename1)); 
				reader2 = new VCFReader(new File(filename2)); 
			}
		
			VCFHeader header1 = reader1.getHeader();
			VCFHeader header2 = reader2.getHeader();

			VCFRecord record1 = reader1.next();
			VCFRecord record2 = reader2.next();

			VCFWriter writer = new VCFWriter(new File(output_filename));
			writer.writeHeader(header1);

			while(true)
			{
				if ((record1 == null) && (record2 == null)) { break; }
				else if (record1 == null) { writer.addRecord(record2); record2 = reader2.next(); continue; }
				else if (record2 == null) { writer.addRecord(record1); record1 = reader1.next(); continue; }

				if (verbose)
				{
					System.out.printf("RECORD1: %s\n", record1.toStringEncoding(header1));
					System.out.printf("RECORD2: %s\n", record2.toStringEncoding(header2));
				}

				if (record1.isFiltered()) { record1 = reader1.next(); continue; }
				if (record2.isFiltered()) { record2 = reader2.next(); continue; }

				Interval interval1 = VCFTool.getIntervalFromRecord(record1);
				Interval interval2 = VCFTool.getIntervalFromRecord(record2);

				int comparison = VCFTool.compareIntervals(interval1, interval2);
				
				if (comparison == 0)
				{
					// records match! Emit one.
					writer.addRecord(record1);
					record1 = reader1.next();
					record2 = reader2.next();
				}
				else if (comparison > 0)
				{
					writer.addRecord(record2);
					record2 = reader2.next();
				}
				else if (comparison < 0)
				{
					writer.addRecord(record1);
					record1 = reader1.next();
				}
			}

			writer.close();

			return 0;
		}
}


