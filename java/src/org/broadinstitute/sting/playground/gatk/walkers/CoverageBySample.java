
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

    public void initialize() 
    { 
        GenomeAnalysisTK toolkit = this.getToolkit();
        SAMFileHeader header = toolkit.getSamReader().getFileHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        sample_names    = new ArrayList<String>();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
            sample_names.add(sample_name);
        } 
    }

    public String map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        String line = context.getLocation().getContig() + " " + context.getLocation().getStart() + " " ;
        for (int i = 0; i < sample_names.size(); i++)
        {
            int count = countReadsBySample(context, sample_names.get(i));
            line += " " + count;
        }
        line += "\n";
        return line;
    }

    private int countReadsBySample(LocusContext context, String sample_name)
    {
        int count = 0;
        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));
            String sample = read.getHeader().getReadGroup(RG).getSample();
            if (sample == sample_name) 
            { 
                count += 1;
            }
        }
        return count;
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
        return sum + line; 
    }

    
}
