
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;

import java.util.List;

public class ListSampleIds extends LocusWalker<Boolean, Boolean> 
{
    public void initialize() 
    { 
        GenomeAnalysisEngine toolkit = this.getToolkit();
        SAMFileHeader header = toolkit.getSAMFileHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
            out.println(sample_name);
        } 
    }

    public Boolean map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        List<SAMRecord> reads = context.getReads();
        StringBuilder readNames = new StringBuilder();

        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);
            String rg = (String) read.getAttribute("RG");
            SAMFileHeader header = read.getHeader();
            SAMReadGroupRecord readGroup = header.getReadGroup(rg);
            if (readGroup == null) { System.out.printf("."); return false; }
            String sample = readGroup.getSample();
            System.out.printf("FROM_MAP %s\n", sample);
        }

        return true;
    }

    public void onTraversalDone() 
    {
        return;
    }

    public Boolean reduceInit() 
    { 
        return true;
    }

    public Boolean reduce(Boolean mapresult, Boolean sum) 
    {
        out.flush();
        return true;
    }
}
