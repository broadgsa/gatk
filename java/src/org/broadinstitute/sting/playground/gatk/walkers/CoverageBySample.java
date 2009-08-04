
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;

public class CoverageBySample extends LocusWalker<String, String> 
{
    List<String> sample_names = null;

    public boolean requiresReads()     { return true; }    

    private SAMFileHeader header;

    public void initialize() 
    { 
        GenomeAnalysisEngine toolkit = this.getToolkit();
        this.header = toolkit.getSAMFileHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        sample_names    = new ArrayList<String>();
		HashSet<String> unique_sample_names = new HashSet<String>();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
			if (unique_sample_names.contains(sample_name)) { continue; }
			unique_sample_names.add(sample_name);
            sample_names.add(sample_name);
        } 
    }

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
    {
        String line = context.getLocation().getContig() + " " + context.getLocation().getStart() + " " ;

		AlignmentContext[] contexts = filterAlignmentContext(context, sample_names, 0);

		HashMap<String,Integer> counts = countReadsBySample(context);
        for (int i = 0; i < contexts.length; i++)
        {
			List<SAMRecord> reads = contexts[i].getReads();
			List<Integer> offsets = contexts[i].getOffsets();

			out.printf("%s %s ", context.getLocation(), sample_names.get(i));

			int[] forward_counts = new int[4];
			int[] backward_counts = new int[4];

			for (int j = 0; j < reads.size(); j++)
			{
				SAMRecord read = reads.get(j);
				int offset = offsets.get(j);
				boolean backward = read.getReadNegativeStrandFlag();
				char base = Character.toUpperCase((char)(read.getReadBases()[offset]));

				if (BaseUtils.simpleBaseToBaseIndex(base) == -1) { continue; }

				if (backward) { base = Character.toLowerCase(base); }

				if (! backward) { forward_counts[BaseUtils.simpleBaseToBaseIndex(base)]++; }
				else { backward_counts[BaseUtils.simpleBaseToBaseIndex(base)]++; }

				//out.printf("%c", base);
			}
			out.printf("A[%d] C[%d] G[%d] T[%d]   a[%d] c[%d] g[%d] t[%d]",
							forward_counts[0],
							forward_counts[1],
							forward_counts[2],
							forward_counts[3],
							backward_counts[0],
							backward_counts[1],
							backward_counts[2],
							backward_counts[3]);
			out.printf("\n");
        }
        return "";
    }

    private HashMap<String,Integer> countReadsBySample(AlignmentContext context)
    {
		HashMap<String,Integer> counts = new HashMap<String,Integer>();
		for (int i = 0; i < sample_names.size(); i++)
		{
			counts.put(sample_names.get(i), 0);
		}
        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            String RG = (String)(read.getAttribute("RG"));
            String sample = header.getReadGroup(RG).getSample();
			counts.put(sample, counts.get(sample)+1); 
        }
        return counts;
    }

    private AlignmentContext[] filterAlignmentContext(AlignmentContext context, List<String> sample_names, int downsample)
    {
		HashMap<String,Integer> index = new HashMap<String,Integer>();
		for (int i = 0; i < sample_names.size(); i++)
		{
			index.put(sample_names.get(i), i);
		}

		AlignmentContext[] contexts = new AlignmentContext[sample_names.size()];
		ArrayList<SAMRecord>[] reads = new ArrayList[sample_names.size()];
		ArrayList<Integer>[] offsets = new ArrayList[sample_names.size()];

		for (int i = 0; i < sample_names.size(); i++)
		{
			reads[i] = new ArrayList<SAMRecord>();
			offsets[i] = new ArrayList<Integer>();
		}

        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));

            assert(header != null);
            assert(header.getReadGroup(RG) != null);

            String sample = header.getReadGroup(RG).getSample();
			//if (SAMPLE_NAME_REGEX != null) { sample = sample.replaceAll(SAMPLE_NAME_REGEX, "$1"); }
            reads[index.get(sample)].add(read); 
            offsets[index.get(sample)].add(offset); 
        }

        if (downsample != 0)
        {
			for (int j = 0; j < reads.length; j++)
			{
	            List<Integer> perm = new ArrayList<Integer>(); 
	            for (int i = 0; i < reads[j].size(); i++) { perm.add(i); }
	            perm = Utils.RandomSubset(perm, downsample);
	           
	            ArrayList<SAMRecord> downsampled_reads = new ArrayList<SAMRecord>();
	            ArrayList<Integer> downsampled_offsets = new ArrayList<Integer>();
	
	            for (int i = 0; i < perm.size(); i++)
	            {
	                downsampled_reads.add(reads[j].get(perm.get(i)));
	                downsampled_offsets.add(offsets[j].get(perm.get(i)));
	            }
	
	            reads[j] = downsampled_reads;
	            offsets[j] = downsampled_offsets;
				contexts[j] = new AlignmentContext(context.getLocation(), reads[j], offsets[j]);
			}
        }
		else
		{
			for (int j = 0; j < reads.length; j++)
			{
				contexts[j] = new AlignmentContext(context.getLocation(), reads[j], offsets[j]);
			}
		}

        return contexts;
    }

    public void onTraversalDone(String result) 
    {
        return;
    }

    public String reduceInit() 
    { 
		/*
        String header = "contig offset";
        for (int i = 0; i < sample_names.size(); i++)
        {
            header += " " + sample_names.get(i);
        }
        header += "\n";
        out.print(header);
        return header;
		*/
		return "";
    }

    public String reduce(String line, String sum) 
    {
        //out.print(line);
		out.flush();
        return "";
    }

    
}
