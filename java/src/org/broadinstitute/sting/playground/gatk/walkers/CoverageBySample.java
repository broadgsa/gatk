
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

public class CoverageBySample extends LocusWalker<String, String> 
{
    List<String> sample_names = null;

    public boolean requiresReads()     { return true; }    

    private SAMFileHeader header;

    public void initialize() 
    { 
        GenomeAnalysisEngine toolkit = this.getToolkit();
        this.header = toolkit.getEngine().getSAMHeader();
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

    public String map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        String line = context.getLocation().getContig() + " " + context.getLocation().getStart() + " " ;
		HashMap<String,Integer> counts = countReadsBySample(context);
        for (int i = 0; i < sample_names.size(); i++)
        {
            int count = counts.get(sample_names.get(i));
            line += " " + count;
        }
        line += "\n";
        return line;
    }

    private HashMap<String,Integer> countReadsBySample(LocusContext context)
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

    public void onTraversalDone() 
    {
        return;
    }

    public String reduceInit() 
    { 
        String header = "contig offset";
        for (int i = 0; i < sample_names.size(); i++)
        {
            header += " " + sample_names.get(i);
        }
        header += "\n";
        out.print(header);
        return header;
    }

    public String reduce(String line, String sum) 
    {
        out.print(line);
        return "";
    }

    
}
