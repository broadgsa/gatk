/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.tools.vcf;

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

			public String toLine()
			{
				int on_diag            = 0;
				int on_diag_not_homref = 0;
				int off_diag           = 0;
				int total              = 0;
				int total_not_homref   = 0;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						if (i == j) { on_diag += counts[i][j]; }
						if (i == j && i != 0) { on_diag_not_homref += counts[i][j]; }
						if (i != j) { off_diag += counts[i][j]; }
						if (i != 0 || j != 0) { total_not_homref += counts[i][j]; }
						total += counts[i][j];
					}
				}

				String s = String.format("SNP %s %d %d %f %f\n", this.name, total, total_not_homref, this.errorRate(), this.hetErrorRate());
				return s;
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

			public double hetErrorRate()
			{
				int true_hets          = 0;
				int correct_hets       = 0;
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						if (j == 1) { true_hets += counts[i][j]; }
					}
				}
				correct_hets = counts[1][1];
				double het_error_rate = 1.0 - ((double)correct_hets / (double)true_hets);
				return het_error_rate;
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
