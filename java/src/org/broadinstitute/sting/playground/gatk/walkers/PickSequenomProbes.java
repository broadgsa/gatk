package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.*;
import java.util.*;


@WalkerName("PickSequenomProbes")
@Requires(value={DataSource.REFERENCE})
@Reference(window=@Window(start=-200,stop=200))
public class PickSequenomProbes extends RefWalker<String, String>
{
    @Argument(required=false, shortName="snp_mask", doc="positions to be masked with N's") public String SNP_MASK = null;

	ArrayList<GenomeLoc> mask = null;
	Object[] mask_array = null;
    public void initialize() 
	{
		System.out.printf("Loading SNP mask...  ");
		if (SNP_MASK != null)
		{
	         mask = new ArrayList<GenomeLoc>();
			 mask.addAll(GenomeLocParser.intervalFileToList(SNP_MASK));

			 Object[] temp_array = mask.toArray();
			 mask_array = new GenomeLoc[temp_array.length];
			 for (int i = 0; i < temp_array.length; i++) { mask_array[i] = (GenomeLoc)temp_array[i]; }
			 Arrays.sort(mask_array);
		}
		System.out.printf("Done.\n");
    }

	private boolean in_mask(GenomeLoc loc)
	{
		return (Arrays.binarySearch(mask_array, loc) >= 0);
	}

    public String map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) 
	{
        String refBase = String.valueOf(ref.getBase());

		System.out.printf("Probing " + ref.getLocus() + " " + ref.getWindow() + "\n");

        Iterator<ReferenceOrderedDatum> rods = rodData.getAllRods().iterator();

		Variation snp = null;
        while (rods.hasNext()) 
		{
            ReferenceOrderedDatum rod = rods.next();
            if (!(rod instanceof Variation))
                continue;

            Variation variant = (Variation) rod;
            if (variant.isSNP())
			{
				snp = variant;
            }
        }

		String contig = context.getLocation().getContig();
		long   offset = context.getLocation().getStart();

		char[] context_bases = ref.getBases();
		long true_offset = offset - 200; 
		for (long i = 0; i < 401; i++)
		{
			GenomeLoc loc = GenomeLocParser.parseGenomeLoc(context.getLocation().getContig() + ":" + true_offset + "-" + true_offset);
			if (in_mask(loc)) { context_bases[(int)i] = 'N'; }
			true_offset += 1;
		}
		char[] leading_bases  = Arrays.copyOfRange(context_bases, 0, 200);
		char[] trailing_bases = Arrays.copyOfRange(context_bases, 201, 401);

		if (snp != null) 
		{
			String snp_string = new String(leading_bases) + new String("[" + refBase + "/" + snp.getAlternativeBaseForSNP() + "]") + new String(trailing_bases);
			String fasta_string = new String(">" + context.getLocation().toString() + "_" + ref.getWindow().toString());
			return fasta_string + "\n" + snp_string + "\n";
		}
		else
		{
			return "";
		}
    }

	public String reduceInit()
	{
		return "";
	}

	public String reduce(String data, String sum)
	{
		out.print(data);
		return "";
	}

    public void onTraversalDone(String sum) 
	{
    }

}
