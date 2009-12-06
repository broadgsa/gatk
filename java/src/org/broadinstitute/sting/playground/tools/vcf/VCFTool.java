package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;

import edu.mit.broad.picard.util.Interval;


import java.io.*;
import java.util.*;


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

			VCFHeader header = reader.getHeader();

			writer = new VCFWriter(header, new File(out_filename));

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

class VCFConcordance extends CommandLineProgram 
{
		@Argument(fullName = "vcf1", shortName = "vcf1", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "vcf2", shortName = "vcf2", doc = "file to open", required = true) public String filename2;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = true) public String output_filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print extremely detailed stats", required = false) public Boolean verbose = false;
		@Argument(fullName = "qual_threshold", shortName = "qual_threshold", doc = "minimum genotype quality to consider", required = false) public long qual_threshold = 1;

		class GenotypeConcordance
		{
			String name;

			protected int[][] counts = {{0,0,0},
							            {0,0,0},
							            {0,0,0}};

			public GenotypeConcordance(String name)
			{
				this.name = name;
			}

			public void add(char ref, String g1, String g2)
			{
				int g1_dosage = 0;
				int g2_dosage = 0;

				if (g1.charAt(0) != ref) { g1_dosage += 1; }
				if (g1.charAt(1) != ref) { g1_dosage += 1; }

				if (g2.charAt(0) != ref) { g2_dosage += 1; }
				if (g2.charAt(1) != ref) { g2_dosage += 1; }

				counts[g1_dosage][g2_dosage] += 1;
			}

			public void add(GenotypeConcordance G)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						counts[i][j] += G.counts[i][j];
					}
				}
			}

			public String toString()
			{
				String s = this.name + "\n";

				int on_diag            = 0;
				int on_diag_not_homref = 0;
				int off_diag           = 0;
				int total              = 0;
				int total_not_homref   = 0;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						s += counts[i][j] + "\t";

						if (i == j) { on_diag += counts[i][j]; }
						if (i == j && i != 0) { on_diag_not_homref += counts[i][j]; }
						if (i != j) { off_diag += counts[i][j]; }
						if (i != 0 || j != 0) { total_not_homref += counts[i][j]; }
						total += counts[i][j];
					}
					s += "\n";
				}

				s += String.format("On-Diagonal                = %.02f\n", 100.0 * (double)on_diag / (double)total);
				s += String.format("On-Diagonal  (not hom-ref) = %.02f\n", 100.0 * (double)on_diag_not_homref / (double)total_not_homref);
				s += String.format("Off-Diagonal               = %.02f\n", 100.0 * (double)off_diag / (double)total_not_homref);
				s += String.format("Total                      = %d\n", total);
				s += String.format("Total (not hom-ref)        = %d\n", total_not_homref);

				s += "\n";
				return s;
			}

			public double errorRate()
			{
				int off_diag           = 0;
				int total_not_homref   = 0;
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						if (i != j) { off_diag += counts[i][j]; }
						if (i != 0 || j != 0) { total_not_homref += counts[i][j]; }
					}
				}
				double error_rate = (double)off_diag / (double)total_not_homref;
				return error_rate;
			}

			public int total()
			{
				int total = 0;
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						total += counts[i][j];
					}
				}
				return total;
			}

			public int totalNonHomRef()
			{
				int total = 0;
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						if (i != 0 || j != 0) { total += counts[i][j]; }
					}
				}
				return total;
			}

		}

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

}

