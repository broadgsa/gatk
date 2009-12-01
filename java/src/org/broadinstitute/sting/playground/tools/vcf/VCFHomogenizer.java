package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;
            
import org.broadinstitute.sting.utils.genotype.vcf.*;
        
import edu.mit.broad.picard.util.Interval;


import java.io.*;
import java.util.*;
import java.util.zip.*;

// Edit a VCF on the fly to be on-spec.

/** 
 * @author jmaguire
 */

class VCFHomogenizer extends InputStream 
{
	private BufferedReader in = null;
	private String currentLine;
	private ByteArrayInputStream stream;

	public VCFHomogenizer(Reader in)
	{
		this.in = (BufferedReader)in;
		currentLine = this.readLine();
		stream = new ByteArrayInputStream(currentLine.getBytes());
	}

	////////////////////////////////////////
	// InputStream methods
	
	public int available()
	{
		return 1;
	}

	public void close() throws java.io.IOException
	{
		this.in.close();
	}

	public void mark(int readlimit)
	{
		System.err.println("not implemented");
		System.exit(-1);
	}

	public void reset()
	{
		System.err.println("not implemented");
		System.exit(-1);
	}

	public boolean markSupported()
	{
		return false; 
	}

	public int read()
	{
		if ((stream == null) || (stream.available() == 0))
		{
			currentLine = this.readLine();
			if (currentLine == null) { return -1; }
			stream = new ByteArrayInputStream(currentLine.getBytes());
		}
		return stream.read();
	}

	// END InputStream methods
	/////////////////////////////////////


	public static VCFHomogenizer create(String filename)
	{
		try
		{
			if (filename.endsWith(".gz"))
			{
				return new VCFHomogenizer(new BufferedReader(
											new InputStreamReader(
												new GZIPInputStream(
													new FileInputStream(filename)))));
			}
			else
			{
				return new VCFHomogenizer(new BufferedReader(
											new InputStreamReader(
												new FileInputStream(filename))));
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(-1);
		}

		return null;
	}

	//my ($chr, $off, $id, $ref, $alt, $qual, $filter, $info, $format, @genotypes) = @tokens;
	private String editLine(String input)
	{
		if (input == null) { return null; }

		//System.out.println("input : " + input);

		/////////
		// Header corrections
		if (input.startsWith("##format=VCFv3.2")) { return "##format=VCRv3.2\n"; }
		if (input.startsWith("#CHROM")) { return input.replaceAll("PROB", "QUAL"); }
		if (input.startsWith("#")) { return input; }

		/////////
		// Line-level corrections
		
		// make "nan" into "NaN"
		input = input.replaceAll("nan", "NaN");
		input = input.replaceAll("DB(\\;|\\s)", "DB=1$1");
		input = input.replaceAll("HM2(\\;|\\s)", "HM2=1$1");
		input = input.replaceAll("HM3(\\;|\\s)", "HM3=1$1");

		String[] tokens = input.split("\\s+");

		/////////
		// Token-level corrections

		// if alt is "N", make it "."
		if (tokens[4].equals("N")) { tokens[4] = "."; }

		String ref = tokens[3];
		String alt = tokens[4];
		String[] alts = alt.split(",");

		for (int i = 9; i < tokens.length; i++)
		{
			if (tokens[i].equals(".")) { tokens[i] = "./.:0"; }

			tokens[i] = tokens[i].replaceAll(ref, "0");
			if (! alt.equals(".")) 
			{ 
				if (alts.length == 1)
				{
					tokens[i] = tokens[i].replaceAll(alt, "1"); 
				}
				else
				{
					for (int j = 0; j < alts.length; j++)
					{
						tokens[i] = tokens[i].replaceAll(alts[j], "1"); 
					}
				}
			}
		}

		/////////
		// Now put it back together and emit.
		String output = tokens[0];
		for (int i = 1; i < tokens.length; i++)
		{
			output = output + "\t" + tokens[i];
		}
		output = output + "\n";

		//System.out.println("output: " + output);

		return output;
	}

	public String readLine()
	{
		try
		{
			String line = in.readLine();
			if (line == null) { return null; }
			else { return editLine(line + "\n"); }
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(-1);
		}

		return null;
	}

	/*
	public int read()
	{
		throw new RuntimeException("VCFHomogenizer.read() not implemented.");
	}

	public byte read(byte[] b)
	{
		throw new RuntimeException("VCFHomogenizer.read(byte[]) not implemented.");
	}

	public int read(byte[] b, int off, int len)
	{
		throw new RuntimeException("VCFHomogenizer.read(byte[], int, int) not implemented.");
	}

	public long skip(long n)
	{
		throw new RuntimeException("VCFHomogenizer.skip(long) not implemented.");
	}

	public boolean markSupported() { return false; }
	*/
}



