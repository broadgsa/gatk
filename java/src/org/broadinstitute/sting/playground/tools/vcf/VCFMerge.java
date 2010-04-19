/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ÓSoftwareÓ), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ÓAS ISÓ, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;


import java.io.*;

import net.sf.picard.util.Interval;

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


