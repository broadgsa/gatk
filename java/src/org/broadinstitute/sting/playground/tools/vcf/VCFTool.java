package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;

import edu.mit.broad.picard.util.Interval;
import org.broadinstitute.sting.utils.GenomeLoc;


import java.io.*;
import java.util.*;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;


// First draft of a program for working with VCF files in various ways.


/**
 * @author jmaguire
 */


class VCFValidate extends CommandLineProgram 
{
		@Argument(fullName = "vcf", shortName = "vcf", doc = "file to open", required = true) public String filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "print", shortName = "print", doc = "print the vcf records to stdout", required = false) public Boolean print = false;
		@Argument(fullName = "profile", shortName = "profile", doc = "print performance information", required = false) public Boolean profile = false;

		@Override
		protected int execute() 
		{
			System.out.println("Validating " + filename + "...");
		
			VCFReader reader = null;

			if (autocorrect) { reader = new VCFReader(VCFHomogenizer.create(filename)); }
			else { reader = new VCFReader(new File(filename)); }

			VCFHeader header = reader.getHeader();

			Date start_time = new Date();
			int n_records_processed = 0;
			while(reader.hasNext())
			{
				VCFRecord record = reader.next();
				if (print) { System.out.println(record.toStringEncoding(header)); }

				if ((profile) && (n_records_processed % 10000 == 0))
				{
					Date current_time = new Date();
					long elapsed = current_time.getTime() - start_time.getTime();
					System.out.printf("RUNTIME: %d records processed in %f seconds; %f seconds per record.\n",
						n_records_processed,
						(double)elapsed/1000.0,
						((double)elapsed/1000.0)/(double)n_records_processed);
				}
				n_records_processed += 1;
			}

			if (autocorrect) { System.out.println(filename + " is VALID (after auto-correction)."); }
			else { System.out.println(filename + " is VALID."); }
			
			return 0;
		}
}

class VCFStats extends CommandLineProgram 
{
		@Argument(fullName = "input", shortName = "input", doc = "file to read", required = true) public String in_filename;
		//@Argument(fullName = "output", shortName = "output", doc = "file to write", required = true) public String out_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "locus", shortName = "locus", doc = "file listing loci to extract", required = true) public String locus_string;


		@Override
		protected int execute() 
		{
			VCFReader reader = null;

			String[] tokens = locus_string.split("\\:|\\-");
			String chr   = tokens[0];
			String start = tokens[1];
			String stop  = tokens[2];
			Interval locus = new Interval(chr, Integer.parseInt(start), Integer.parseInt(stop));

			if (autocorrect) { reader = new VCFReader(VCFHomogenizer.create(in_filename)); }
			else { reader = new VCFReader(new File(in_filename)); }

			VCFHeader header = reader.getHeader();


			//////////////
			// Stats collectors
			int transitions = 0;
			int transversions = 0;
			int total_snps = 0;
			int[] AC_histogram = new int[1000]; int highest_AC = 0;

			boolean before = true;

			while(reader.hasNext())
			{
				VCFRecord record = reader.next();
				Interval this_locus = VCFTool.getIntervalFromRecord(record);
				if (locus.intersects(this_locus)) 
				{ 
					before = false;

					Map<String,String> info = record.getInfoValues();
					int AC = Integer.parseInt(info.get("AC"));
					AC_histogram[AC] += 1;
					if (AC > highest_AC) { highest_AC = AC; }

					if (VCFTool.isTransition(record)) { transitions += 1; }
					else { transversions += 1; }

					total_snps += 1;
					//System.out.printf("%s\n", record.toStringEncoding(header));
				}
				else if ((before == false) && (this_locus.compareTo(locus) > 0)) { break; }
			}

			System.out.printf("Total SNPs : %d\n", total_snps);
			System.out.printf("Ts/Tv      : %.02f\n", (double)transitions / (double)transversions);
			System.out.printf("\n");
			System.out.printf("AAF\tCount\n");
			for (int AC = 1; AC <= highest_AC; AC++)
			{	
				System.out.printf("%d\t%d\n", AC, AC_histogram[AC]);
			}
			System.out.printf("\n");

			return 0;
		}

}

class CheckRefFields extends CommandLineProgram
{
		@Argument(fullName = "vcf", shortName = "vcf", doc = "file to open", required = true) public String filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "fasta", shortName = "fasta", doc = "reference FASTA", required = true) public String fasta_filename;

		@Override
		protected int execute() 
		{
			System.out.println("Checking " + filename + "...");
		
			VCFReader reader = null;

			if (autocorrect) { reader = new VCFReader(VCFHomogenizer.create(filename)); }
			else { reader = new VCFReader(new File(filename)); }

			ReferenceSequenceFileWalker ref = new ReferenceSequenceFileWalker(new File(fasta_filename));
			String ref_seq_name = "";
			byte[] ref_seq = null;
			SAMSequenceDictionary ref_dict = ref.getSequenceDictionary();

			VCFHeader header = reader.getHeader();

			Date start_time = new Date();
			int n_records_processed = 0;
			while(reader.hasNext())
			{
				VCFRecord record = reader.next();

				String chr = record.getLocation().getContig();
				if (! chr.equals(ref_seq_name))
				{
					System.out.println("Loading " + chr);
					ref_seq = ref.get(ref_dict.getSequence(chr).getSequenceIndex()).getBases();
					ref_seq_name = chr;
				}	

				long offset   = record.getLocation().getStart();
				char vcf_ref_base = record.getReferenceBase();
				char fasta_ref_base = (char)ref_seq[(int)offset-1];

				List<VCFGenotypeEncoding> alleles = record.getAlternateAlleles();
				char vcf_alt_base = alleles.get(0).getBases().charAt(0); 

				//System.out.println(chr + " " + offset + " " + fasta_ref_base + " " + vcf_ref_base + " " + vcf_alt_base);

				String ans = null;
				if (vcf_ref_base != fasta_ref_base)
				{
					System.out.println("Error! Ref field does not match fasta. Fasta says " + fasta_ref_base);
					System.out.println(record.toStringEncoding(header));
				}
			}

			System.out.println("All reference fields correct.");
			return 0;
		}
}


class FixRefFields extends CommandLineProgram
{
		@Argument(fullName = "vcf", shortName = "vcf", doc = "file to open", required = true) public String filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "fasta", shortName = "fasta", doc = "reference FASTA", required = true) public String fasta_filename;
		@Argument(fullName = "output", shortName = "output", doc = "output file", required = true) public String output_filename;

		@Override
		protected int execute() 
		{
			System.out.println("Fixing " + filename + "...");
		
			VCFReader reader = null;

			if (autocorrect) { reader = new VCFReader(VCFHomogenizer.create(filename)); }
			else { reader = new VCFReader(new File(filename)); }

			ReferenceSequenceFileWalker ref = new ReferenceSequenceFileWalker(new File(fasta_filename));
			String ref_seq_name = "";
			byte[] ref_seq = null;
			SAMSequenceDictionary ref_dict = ref.getSequenceDictionary();


			VCFHeader header = reader.getHeader();

			PrintStream output;
			try
			{
				VCFWriter writer = new VCFWriter(new File(output_filename));
                writer.writeHeader(header);
				writer.close();
				output = new PrintStream(new FileOutputStream(output_filename, true));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			Date start_time = new Date();
			int n_records_processed = 0;
			while(reader.hasNext())
			{
				VCFRecord record = reader.next();

				String chr = record.getLocation().getContig();
				if (! chr.equals(ref_seq_name))
				{
					System.out.println("Loading " + chr);
					ref_seq = ref.get(ref_dict.getSequence(chr).getSequenceIndex()).getBases();
					ref_seq_name = chr;
				}	

				long offset   = record.getLocation().getStart();
				char vcf_ref_base = record.getReferenceBase();
				char fasta_ref_base = (char)ref_seq[(int)offset-1];

				List<VCFGenotypeEncoding> alleles = record.getAlternateAlleles();
				char vcf_alt_base = alleles.get(0).getBases().charAt(0); 

				//System.out.println(chr + " " + offset + " " + fasta_ref_base + " " + vcf_ref_base + " " + vcf_alt_base);

				String ans = null;
				if ((vcf_ref_base != fasta_ref_base) && ((vcf_alt_base == fasta_ref_base) || (vcf_alt_base == '.')))
				{
					// swap!
					String s = record.toStringEncoding(header);
					String[] tokens = s.split("\\s+");
					tokens[3] = Character.toString(fasta_ref_base);
					tokens[4] = Character.toString(vcf_ref_base);
					for (int i = 9; i < tokens.length; i++)
					{
						tokens[i] = tokens[i].replaceAll("0", "A");
						tokens[i] = tokens[i].replaceAll("1", "B");
						tokens[i] = tokens[i].replaceAll("B", "0");
						tokens[i] = tokens[i].replaceAll("A", "1");
					}

					ans = "";
					for (int i = 0; i < tokens.length; i++)
					{
						ans = ans + tokens[i] + "\t";
					}
					ans.replaceAll("\\s+$", "");

					//System.out.println("from: " + s);
					//System.out.println("to: " + ans);
				}
				else
				{
					ans = record.toStringEncoding(header);
				}

				output.println(ans);
			}

			output.flush();
			output.close();

			System.out.println("Done.");
			return 0;
		}
}


class VCFGrep extends CommandLineProgram 
{
		@Argument(fullName = "input", shortName = "input", doc = "file to read", required = true) public String in_filename;
		@Argument(fullName = "output", shortName = "output", doc = "file to write", required = true) public String out_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "loci", shortName = "loci", doc = "file listing loci to extract", required = true) public String loci_filename;

		@Override
		protected int execute() 
		{
			VCFReader reader = null;
			VCFWriter writer = null;

			HashSet<Interval> loci = new HashSet<Interval>();
			try
			{
				Scanner loci_reader = new Scanner(new File(loci_filename));
				while(loci_reader.hasNextLine())
				{
					String line = loci_reader.nextLine();
					String[] tokens = line.split("\\:");
	
					String chr = tokens[0];
					String off = tokens[1];
					loci.add(new Interval(chr, Integer.parseInt(off), Integer.parseInt(off)));
				}
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			if (autocorrect) { reader = new VCFReader(VCFHomogenizer.create(in_filename)); }
			else { reader = new VCFReader(new File(in_filename)); }


            writer = new VCFWriter(new File(out_filename));
            writer.writeHeader(reader.getHeader());

			while(reader.hasNext())
			{
				VCFRecord record = reader.next();
				Interval locus = VCFTool.getIntervalFromRecord(record);
				if (loci.contains(locus)) { writer.addRecord(record); }
			}
			writer.close();

			return 0;
		}

}

class PrintGQ extends CommandLineProgram 
{
		@Argument(fullName = "vcf", shortName = "vcf", doc = "file to open", required = true) public String filename;

		@Override
		protected int execute() 
		{
			VCFReader reader;
			VCFReader reader2;

				reader = new VCFReader(VCFHomogenizer.create(filename)); 
		
			VCFHeader header = reader.getHeader();
			VCFRecord record = reader.next();

			while(true)
			{
				if (record == null) { break; }

				Interval interval = VCFTool.getIntervalFromRecord(record);

					if (record.isFiltered())
					{
						record = reader.next();
					}
					
					char ref = record.getReferenceBase();
					
					String[] sample_names = record.getSampleNames();

					List<VCFGenotypeRecord> genotypes = record.getVCFGenotypeRecords();

					for (int i = 0; i < sample_names.length; i++)
					{
						VCFGenotypeRecord rec = genotypes.get(i);

						String gq = rec.getFields().get("GQ");

						List<VCFGenotypeEncoding> alleles = rec.getAlleles();

						String g = "";

						for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
						char[] c = g.toCharArray();

						Arrays.sort(c);

						g = new String(c);

						System.out.println(g + " " + gq);
					}

					record = reader.next();
			}

			return 0;
		}
}

class VCFSimpleStats extends CommandLineProgram 
{
		@Argument(fullName = "vcf1", shortName = "vcf1", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = true) public String output_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print extremely detailed stats", required = false) public Boolean verbose = false;
		@Argument(fullName = "min_call_rate", shortName = "min_call_rate", doc = "what fraction of samples must have a call", required = false) public double min_call_rate = 0.9;

		@Override
		protected int execute() 
		{
			//System.out.println("Loading " + filename + "...");
		
			PrintStream output = null;
			try
			{
				output = new PrintStream(new FileOutputStream(output_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			VCFReader reader1;

			if (autocorrect) 
			{ 
				reader1 = new VCFReader(VCFHomogenizer.create(filename1)); 
			}
			else 
			{ 
				reader1 = new VCFReader(new File(filename1)); 
			}
		
			VCFHeader header1 = reader1.getHeader();

			VCFRecord record1 = reader1.next();

			int TP = 0;
			int FP = 0;
			int TN = 0;
			int FN = 0;
			int total = 0;
			int dropped = 0;

			int ts = 0;
			int tv = 0;

			while(true)
			{
				if (record1 == null) { break; }

				Interval interval1 = VCFTool.getIntervalFromRecord(record1);
					
					// (unless it is "filtered")
					if (record1.isFiltered())
					{
						record1 = reader1.next();
					}
					
					char ref = record1.getReferenceBase();

					
					String[] sample_names1 = record1.getSampleNames();

					List<VCFGenotypeRecord> genotypes1 = record1.getVCFGenotypeRecords();

					long n_ref_1 = 0;
					long n_alt_1 = 0;
					long n_total_1 = 0;
					long n_calls_1 = 0;
					long n_dropped_1 = 0;

					for (int i = 0; i < sample_names1.length; i++)
					{
						VCFGenotypeRecord rec1 = genotypes1.get(i);

						//if (rec2 == null) { continue; }

						Long gq1;

						if (rec1.getFields().get("GQ") != null)
						{
							Double gq1_double = Double.parseDouble(rec1.getFields().get("GQ"));
							gq1 = gq1_double.longValue();
						}
						else
						{
							gq1 = 0L;
						}

						List<VCFGenotypeEncoding> alleles1 = rec1.getAlleles();

						String g1 = "";

						for (int j = 0; j < alleles1.size(); j++) { g1 += alleles1.get(j).getBases(); }

						char[] c1 = g1.toCharArray();

						Arrays.sort(c1);

						g1 = new String(c1);

						n_total_1 += 1;

						if (g1.equals(".."))
						{
							n_dropped_1 += 1;
							continue;
						}

						n_calls_1 += 1;

						if (g1.charAt(0) == ref) { n_ref_1 += 1; } else { n_alt_1 += 1; }
						if (g1.charAt(1) == ref) { n_ref_1 += 1; } else { n_alt_1 += 1; }
					}

					if (((double)n_calls_1 / (double)n_total_1) >= min_call_rate)
					{
						if (n_alt_1 == 0) { FP += 1; }
						if (n_alt_1 >  0) { TP += 1; }
						total += 1;

						if (VCFTool.isTransition(record1)) { ts += 1; }
						else { tv += 1; }
					}
					else
					{
						dropped += 1;
					}

					if ((verbose) && (((double)n_calls_1 / (double)n_total_1) >= min_call_rate))
					{
						//output.printf("SNP " 
						//				+ interval1.toString() 
						//				+ " " + n_total_1 + " " + n_calls_1 + " " +  (double)n_calls_1/(double)n_total_1 + " " + n_ref_1 + " " + n_alt_1 + "\n"); 
						if (n_alt_1 == 0) { output.printf("FP: %s\n", interval1.toString()); }
						if (n_alt_1 != 0) { output.printf("TP: %s\n", interval1.toString()); }
					}

					record1 = reader1.next();
			}


			// Now output the statistics.

			output.printf("TP FP dropped ts tv ts/tv\n%d(%f) %d(%f) %d %d %d %f\n",
							TP, (double)TP/(double)total,		
							FP, (double)FP/(double)total,
							dropped,
							ts, tv,
							(double)ts/(double)tv);	

			output.flush();
			output.close();


			return 0;
		}
}
/*
{
		@Argument(fullName = "vcf1", shortName = "vcf1", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "vcf2", shortName = "vcf2", doc = "file to open", required = true) public String filename2;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = true) public String output_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print extremely detailed stats", required = false) public Boolean verbose = false;

		@Override
		protected int execute() 
		{
			//System.out.println("Loading " + filename + "...");
		
			PrintStream output = null;
			try
			{
				output = new PrintStream(new FileOutputStream(output_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

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

			int TP = 0;
			int FP = 0;
			int TN = 0;
			int FN = 0;
			int total = 0;

			while(true)
			{
				if (record1 == null) { break; }
				if (record2 == null) { break; }

				Interval interval1 = VCFTool.getIntervalFromRecord(record1);
				Interval interval2 = VCFTool.getIntervalFromRecord(record2);

				int comparison = interval1.compareTo(interval2);
				
				if (comparison == 0)
				{
					// records match! compute concordance.
					
					// (unless one of them is "filtered")
					if (record1.isFiltered() || record2.isFiltered())
					{
						record1 = reader1.next();
						record2 = reader2.next();
					}
					
					char ref = record1.getReferenceBase();
					
					String[] sample_names1 = record1.getSampleNames();
					String[] sample_names2 = record2.getSampleNames();

					List<VCFGenotypeRecord> genotypes1 = record1.getVCFGenotypeRecords();
					List<VCFGenotypeRecord> genotypes2 = record2.getVCFGenotypeRecords();

					Map<String, VCFGenotypeRecord> map2 = new HashMap<String, VCFGenotypeRecord>();
					for (int i = 0; i < genotypes2.size(); i++)
					{
						map2.put(genotypes2.get(i).getSampleName(), genotypes2.get(i));
					}

					long n_ref_1 = 0;
					long n_alt_1 = 0;

					long n_ref_2 = 0;
					long n_alt_2 = 0;

					for (int i = 0; i < sample_names1.length; i++)
					{
						VCFGenotypeRecord rec1 = genotypes1.get(i);
						VCFGenotypeRecord rec2 = map2.get(sample_names1[i]);

						if (rec2 == null) { continue; }

						Long gq1;

						if (rec1.getFields().get("GQ") != null)
						{
							Double gq1_double = Double.parseDouble(rec1.getFields().get("GQ"));
							gq1 = gq1_double.longValue();
						}
						else
						{
							gq1 = 0L;
						}

						List<VCFGenotypeEncoding> alleles1 = rec1.getAlleles();
						List<VCFGenotypeEncoding> alleles2 = rec2.getAlleles();

						String g1 = "";
						String g2 = "";

						for (int j = 0; j < alleles1.size(); j++) { g1 += alleles1.get(j).getBases(); }
						for (int j = 0; j < alleles2.size(); j++) { g2 += alleles2.get(j).getBases(); }

						char[] c1 = g1.toCharArray();
						char[] c2 = g2.toCharArray();

						Arrays.sort(c1);
						Arrays.sort(c2);

						g1 = new String(c1);
						g2 = new String(c2);

						if ((g1.equals("..")) ||
							(g2.equals("..")))
						{
							continue;
						}

						if (g1.charAt(0) == ref) { n_ref_1 += 1; } else { n_alt_1 += 1; }
						if (g1.charAt(1) == ref) { n_ref_1 += 1; } else { n_alt_1 += 1; }

						if (g2.charAt(0) == ref) { n_ref_2 += 1; } else { n_alt_2 += 1; }
						if (g2.charAt(1) == ref) { n_ref_2 += 1; } else { n_alt_2 += 1; }
					}

					if ((n_alt_1 > 0) && (n_alt_2 > 0)) { TP += 1; }
					if ((n_alt_1 > 0) && (n_alt_2 == 0)) { FP += 1; }
					if ((n_alt_1 == 0) && (n_alt_2 > 0)) { FN += 1; }
					if ((n_alt_1 == 0) && (n_alt_2 == 0)) { TN += 1; }
					total += 1;

					if (verbose)
					{
						output.printf("SNP " 
										+ interval1.toString() 
										+ " " + n_ref_1 + " " + n_alt_1 
										+ " " + n_ref_2 + " " + n_alt_2 + "\n"); 
					}

					record1 = reader1.next();
					record2 = reader2.next();
				}
				else if (comparison > 0)
				{
					// interval1 is later than interval2.
					record2 = reader2.next();
				}
				else if (comparison < 0)
				{
					// interval2 is later than interval1.
					record1 = reader1.next();
				}

			}


			// Now output the statistics.

			output.printf("TP FP TN FN\n%d(%f) %d(%f) %d(%f) %d(%f)\n",
							TP, (double)TP/(double)total,		
							FP, (double)FP/(double)total,		
							TN, (double)TN/(double)total,		
							FN, (double)FN/(double)total);

			output.flush();
			output.close();


			return 0;
		}
}
*/


class VCFConcordance extends CommandLineProgram 
{
		@Argument(fullName = "vcf1", shortName = "vcf1", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "vcf2", shortName = "vcf2", doc = "file to open", required = true) public String filename2;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = true) public String output_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print extremely detailed stats", required = false) public Boolean verbose = false;
		@Argument(fullName = "qual_threshold", shortName = "qual_threshold", doc = "minimum genotype quality to consider", required = false) public long qual_threshold = 1;


		@Override
		protected int execute() 
		{
			//System.out.println("Loading " + filename + "...");
		
			/////////////////////////////////
			// All the various concordance counters
		
			HashMap<String,GenotypeConcordance> individual = new HashMap<String,GenotypeConcordance>();
			HashMap<Long,GenotypeConcordance>   AAF        = new HashMap<Long,GenotypeConcordance>();
			HashMap<Long,GenotypeConcordance>   Qual       = new HashMap<Long,GenotypeConcordance>();

			// 
			/////////////////////////////////
		
			PrintStream output = null;
			try
			{
				output = new PrintStream(new FileOutputStream(output_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

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

			while(true)
			{
				if (record1 == null) { break; }
				if (record2 == null) { break; }


				Interval interval1 = VCFTool.getIntervalFromRecord(record1);
				Interval interval2 = VCFTool.getIntervalFromRecord(record2);

				int comparison = interval1.compareTo(interval2);
				
				if (comparison == 0)
				{
					// records match! compute concordance.
					
					// (unless one of them is "filtered")
					if (record1.isFiltered() || record2.isFiltered())
					{
						record1 = reader1.next();
						record2 = reader2.next();
					}
					
					char ref = record1.getReferenceBase();
					
					String[] sample_names1 = record1.getSampleNames();
					String[] sample_names2 = record2.getSampleNames();

					List<VCFGenotypeRecord> genotypes1 = record1.getVCFGenotypeRecords();
					List<VCFGenotypeRecord> genotypes2 = record2.getVCFGenotypeRecords();

					Map<String, VCFGenotypeRecord> map2 = new HashMap<String, VCFGenotypeRecord>();
					for (int i = 0; i < genotypes2.size(); i++)
					{
						map2.put(genotypes2.get(i).getSampleName(), genotypes2.get(i));
					}

					GenotypeConcordance SNP = new GenotypeConcordance(interval1.toString());

					long n_ref = 0;
					long n_alt = 0;

					for (int i = 0; i < sample_names1.length; i++)
					{
						VCFGenotypeRecord rec1 = genotypes1.get(i);
						VCFGenotypeRecord rec2 = map2.get(sample_names1[i]);

						if (rec2 == null) { continue; }

						Long gq1;

						if (rec1.getFields().get("GQ") != null)
						{
							Double gq1_double = Double.parseDouble(rec1.getFields().get("GQ"));
							gq1 = gq1_double.longValue();
						}
						else
						{
							gq1 = 0L;
						}

						List<VCFGenotypeEncoding> alleles1 = rec1.getAlleles();
						List<VCFGenotypeEncoding> alleles2 = rec2.getAlleles();

						String g1 = "";
						String g2 = "";

						for (int j = 0; j < alleles1.size(); j++) { g1 += alleles1.get(j).getBases(); }
						for (int j = 0; j < alleles2.size(); j++) { g2 += alleles2.get(j).getBases(); }

						char[] c1 = g1.toCharArray();
						char[] c2 = g2.toCharArray();

						Arrays.sort(c1);
						Arrays.sort(c2);

						g1 = new String(c1);
						g2 = new String(c2);

						if ((g1.equals("..")) ||
							(g2.equals("..")))
						{
							continue;
						}

						if (g1.charAt(0) == ref) { n_ref += 1; } else { n_alt += 1; }
						if (g1.charAt(1) == ref) { n_ref += 1; } else { n_alt += 1; }

						if (! individual.containsKey(sample_names1[i])) { individual.put(sample_names1[i], new GenotypeConcordance(sample_names1[i])); }
						if (! Qual.containsKey(gq1)) { Qual.put(gq1, new GenotypeConcordance(Long.toString(gq1))); }

						individual.get(sample_names1[i]).add(ref, g1, g2);
						Qual.get(gq1).add(ref, g1, g2);
						SNP.add(ref, g1, g2);

						/*
						if (verbose) 
						{ 
							String flag = "";
							if (! g1.equals(g2)) { flag = "X"; }
							output.printf("GENOTYPES " + interval1.toString() + " " + sample_names1[i] + " " + g1 + " " + g2 + " " + flag +  "\n"); 
						}
						*/

					}

					if (verbose) { output.printf("SNP " + SNP.toString()); }

					if (! AAF.containsKey(n_alt)) { AAF.put(n_alt, new GenotypeConcordance(Long.toString(n_alt))); }
					AAF.get(n_alt).add(SNP);

					record1 = reader1.next();
					record2 = reader2.next();
				}
				else if (comparison > 0)
				{
					// interval1 is later than interval2.
					record2 = reader2.next();
				}
				else if (comparison < 0)
				{
					// interval2 is later than interval1.
					record1 = reader1.next();
				}

			}


			// Now output the statistics.
			{
				if (verbose)
				{
					output.printf("\n");
					Object[] individuals = individual.keySet().toArray();
					for (int i = 0; i < individuals.length; i++)
					{
						String ind = (String)individuals[i];
						output.print("INDIVIDUAL " + individual.get(ind).toString());
					}

					output.printf("\n");
					Object[] AAFs = AAF.keySet().toArray();
					for (int i = 0; i < AAFs.length; i++)
					{
						Long aaf = (Long)AAFs[i];
						output.print("AAF " + AAF.get(aaf).toString());
					}

					output.printf("\n");
					Object[] quals = Qual.keySet().toArray();
					for (int i = 0; i < quals.length; i++)
					{
						Long qual = (Long)quals[i];
						output.print("QUAL " + Qual.get(qual).toString());
					}
					output.printf("\n");
				}

				output.printf("\n");
				Object[] individuals = individual.keySet().toArray();
				for (int i = 0; i < individuals.length; i++)
				{
					String ind = (String)individuals[i];
					output.printf("INDIVIDUAL %s %f %d %d\n", ind, individual.get(ind).errorRate(), individual.get(ind).total(), individual.get(ind).totalNonHomRef());
				}

				output.printf("\n");
				Object[] AAFs = AAF.keySet().toArray();
				for (int i = 0; i < AAFs.length; i++)
				{
					Long aaf = (Long)AAFs[i];
					output.printf("AAF %d %f %d %d\n", aaf, AAF.get(aaf).errorRate(), AAF.get(aaf).total(), AAF.get(aaf).totalNonHomRef());
				}

				output.printf("\n");
				Object[] quals = Qual.keySet().toArray();
				for (int i = 0; i < quals.length; i++)
				{
					Long qual = (Long)quals[i];
					output.printf("QUAL %d %f %d %d\n", qual, Qual.get(qual).errorRate(), Qual.get(qual).total(), Qual.get(qual).totalNonHomRef());
				}

			}

			output.flush();
			output.close();


			return 0;
		}
}


public class VCFTool 
{
		public static void main(String args[]) 
		{
			String mode = args[0];
			String[] realArgs = Arrays.copyOfRange(args, 1, args.length);

			if (mode.equals("validate"))
			{
				VCFValidate cm = new VCFValidate();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("grep"))
			{
				VCFGrep cm = new VCFGrep();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("concordance"))
			{
				VCFConcordance cm = new VCFConcordance();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("simple_stats"))
			{
				VCFSimpleStats cm = new VCFSimpleStats();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("printGQ"))
			{
				PrintGQ cm = new PrintGQ();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("fix_ref_fields"))
			{
				FixRefFields cm = new FixRefFields();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("check_ref_fields"))
			{
				CheckRefFields cm = new CheckRefFields();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("stats"))
			{
				VCFStats cm = new VCFStats();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("sequenom"))
			{
				VCFSequenomAnalysis cm = new VCFSequenomAnalysis();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("sequenom2"))
			{
				VCFSequenomAnalysis2 cm = new VCFSequenomAnalysis2();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			if (mode.equals("call_rates"))
			{
				VCFCallRates cm = new VCFCallRates();
				CommandLineProgram.start(cm,realArgs);
				System.exit(0);
			}

			System.out.printf("ERROR: mode %s not defined.\n", mode);			
			System.exit(-1);

		}


		/////////////////////////
		// Some helpful utilities.
		
		public static Interval getIntervalFromRecord(VCFRecord record)
		{
			String chr = record.getLocation().getContig();
			long   off = record.getLocation().getStart();
			return new Interval(chr, (int)off, (int)off);
		}

		public static char getAlt(VCFRecord record)
		{
			List<VCFGenotypeEncoding> alleles = record.getAlternateAlleles();
			char alt = alleles.get(0).getBases().charAt(0); 
			return alt;
		}

		public static boolean isTransition(VCFRecord record)
		{
			char ref = record.getReferenceBase();
			List<VCFGenotypeEncoding> alleles = record.getAlternateAlleles();
			char alt = alleles.get(0).getBases().charAt(0); 

			if (((ref == 'A') && (alt == 'G')) ||
				((ref == 'G') && (alt == 'A')) ||
				((ref == 'C') && (alt == 'T')) ||
				((ref == 'T') && (alt == 'C')))
			{
				return true;
			}
			else
			{
				return false; 
			}
		}


		public static int Compute_n_total(VCFRecord record)
		{
			return VCFTool.Compute_n_total(record, (String[])null);
		}

		public static int Compute_n_total(VCFRecord record, String[] sample_names)
		{
			HashSet<String> set = null;
			if (sample_names != null)
			{
				set = new HashSet<String>();
				for (int i = 0; i < sample_names.length; i++) { set.add(sample_names[i]); }
			}
			return VCFTool.Compute_n_total(record, set);
		}

		public static int Compute_n_total(VCFRecord record, Set<String> sample_mask)
		{
			String[] sample_names = record.getSampleNames();
			List<VCFGenotypeRecord> genotypes = record.getVCFGenotypeRecords();
			int n_ref = 0;
			int n_alt = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				if ((sample_mask != null) && (! sample_mask.contains(sample_names[i])))
				{
					continue;
				}

				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
				if (g.charAt(1) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
			}
			return n_alt + n_ref;
		}

		public static int Compute_n_alt(VCFRecord record)
		{
			return VCFTool.Compute_n_alt(record, (String[])null);
		}

		public static int Compute_n_alt(VCFRecord record, String[] sample_names)
		{
			HashSet<String> set = null;
			if (sample_names != null)
			{
				set = new HashSet<String>();
				for (int i = 0; i < sample_names.length; i++) { set.add(sample_names[i]); }
			}
			return VCFTool.Compute_n_alt(record, set);
		}

		public static int Compute_n_alt(VCFRecord record, Set<String> sample_mask)
		{
			String[] sample_names = record.getSampleNames();
			List<VCFGenotypeRecord> genotypes = record.getVCFGenotypeRecords();
			int n_ref = 0;
			int n_alt = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				// Skip samples we should skip.
				if ((sample_mask != null) && (! sample_mask.contains(sample_names[i])))
				{
					continue;
				}

				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
				if (g.charAt(1) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
			}
			return n_alt;
		}


		public static int Compute_n_het(VCFRecord record)
		{
			return VCFTool.Compute_n_het(record, (String[])null);
		}

		public static int Compute_n_het(VCFRecord record, String[] sample_names)
		{
			HashSet<String> set = null;
			if (sample_names != null)
			{
				set = new HashSet<String>();
				for (int i = 0; i < sample_names.length; i++) { set.add(sample_names[i]); }
			}
			return VCFTool.Compute_n_het(record, set);
		}

		public static int Compute_n_het(VCFRecord record, Set<String> sample_mask)
		{
			String[] sample_names = record.getSampleNames();
			List<VCFGenotypeRecord> genotypes = record.getVCFGenotypeRecords();
			int n_het = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				// Skip samples we should skip.
				if ((sample_mask != null) && (! sample_mask.contains(sample_names[i])))
				{
					continue;
				}

				int n_ref = 0;
				int n_alt = 0;

				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
				if (g.charAt(1) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
				if (n_alt == 1) { n_het += 1; }
			}
			return n_het;
		}

		public static double Compute_failure_rate(VCFRecord record)
		{
			String[] sample_names = record.getSampleNames();
			List<VCFGenotypeRecord> genotypes = record.getVCFGenotypeRecords();
			double failure_rate = 0.0;
			for (int i = 0; i < sample_names.length; i++)
			{
				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { failure_rate += 1; continue; }
			}
			return failure_rate / (double)sample_names.length;
		}

		public static double Compute_HWE(VCFRecord record)
		{
			return VCFTool.Compute_HWE(record, (String[])null);
		}

		public static double Compute_HWE(VCFRecord record, String[] sample_names)
		{
			HashSet<String> set = null;
			if (sample_names != null)
			{
				set = new HashSet<String>();
				for (int i = 0; i < sample_names.length; i++) { set.add(sample_names[i]); }
			}
			return VCFTool.Compute_HWE(record, set);
		}

		public static double Compute_HWE(VCFRecord record, Set<String> sample_mask)
		{
			int ref = 0;
			int het = 0;
			int hom = 0;
			int N   = 0;

			String[] sample_names = record.getSampleNames();
			List<VCFGenotypeRecord> genotypes = record.getVCFGenotypeRecords();
			for (int i = 0; i < sample_names.length; i++)
			{
				// Skip samples we should skip.
				if ((sample_mask != null) && (! sample_mask.contains(sample_names[i])))
				{
					continue;
				}

				int n_ref = 0;
				int n_alt = 0;

				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }
				if (g.charAt(1) == record.getReferenceBase()) { n_ref += 1; } else { n_alt += 1; }

				if (n_ref == 2)                    { ref += 1; }
				else if (n_ref == 1 && n_alt == 1) { het += 1; }
				else if (n_alt == 2)               { hom += 1; }

				N += 1;
			}

			double p = (2.0 * ref + het) / (2.0 * (ref + het + hom));
			double q = 1.0 - p;

			//System.out.printf("DBG: p=%f q=%f ref=%d het=%d hom=%d\n", p, q, ref, het, hom);

			double expected_ref = p * p * N; 
			double expected_het = 2.0 * p * q * N;
			double expected_hom = q * q * N;

			double chi_squared = (Math.pow(ref - expected_ref,2)/expected_ref) + (Math.pow(het - expected_het,2)/expected_het) + (Math.pow(hom - expected_hom,2)/expected_hom);

			return chi_squared;
		}

		// This function assumes a 1-degree of freedom chi-squared.
		public static double P_from_Chi(double chi)
		{
			double gamma = 1.772454;
			double a = Math.pow(2,0.5) * gamma;
			double b = Math.pow(chi, 0.5-1.0) * Math.exp((-1.0 * chi)/2.0);
			double ans = (1.0/a) * b;
			return ans;
		}

}

