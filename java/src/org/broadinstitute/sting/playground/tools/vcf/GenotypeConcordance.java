package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;

import edu.mit.broad.picard.util.Interval;


import java.io.*;
import java.util.*;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

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

				//s += String.format("On-Diagonal                = %.02f\n", 100.0 * (double)on_diag / (double)total);
				//s += String.format("On-Diagonal  (not hom-ref) = %.02f\n", 100.0 * (double)on_diag_not_homref / (double)total_not_homref);
				//s += String.format("Off-Diagonal               = %.02f\n", 100.0 * (double)off_diag / (double)total_not_homref);
				s += String.format("Total                      = %d\n", total);
				s += String.format("Total (not hom-ref)        = %d\n", total_not_homref);
				s += String.format("Error Rate                 = %f\n", this.errorRate());

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
