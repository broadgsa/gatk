package org.broadinstitute.sting.oneoffprojects.multisamplecaller;

import java.lang.*;
import java.util.*;
import java.io.*;

import org.broadinstitute.sting.utils.*;

public class ConfusionMatrix
{
	private double[][] ILLUMINA;
	private double[][] solid;
	private double[][] LS454;

	public ConfusionMatrix(String file_name)
	{
		//System.out.println("DBG: ConfusionMatrix constructor! (" + file_name + ")");

		ILLUMINA = new double[4][4];
		solid = new double[4][4];
		LS454 = new double[4][4];

		try
		{
			Scanner sc = new Scanner(new File(file_name));
			while (sc.hasNext()) 
			{
				String platform = sc.next();
				char   read     = sc.next().charAt(0);
				char   ref      = sc.next().charAt(0);
				double p        = sc.nextDouble();

				int read_i = BaseUtils.simpleBaseToBaseIndex(read);
				int ref_i  = BaseUtils.simpleBaseToBaseIndex(ref);

				if (platform.equals("ILLUMINA")) { ILLUMINA[read_i][ref_i] = p; }
				if (platform.equals("solid"))    { solid[read_i][ref_i] = p; }
				if (platform.equals("LS454"))    { LS454[read_i][ref_i] = p; }
	
				//System.out.println("DBG: " + key + " " + p);
			}
		} 
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(-1);
		}

	}
	
	double lookup(String platform, char read, char truth)
	{
		int read_i = BaseUtils.simpleBaseToBaseIndex(read);
		int truth_i  = BaseUtils.simpleBaseToBaseIndex(truth);

		double d = 0;

		if (platform.equals("ILLUMINA")) { d = ILLUMINA[read_i][truth_i]; }
		else if (platform.equals("solid"))    { d = solid[read_i][truth_i]; }
		else if (platform.equals("LS454"))    { d = LS454[read_i][truth_i]; }
		else { assert(false); }

		return d;
	}
}
