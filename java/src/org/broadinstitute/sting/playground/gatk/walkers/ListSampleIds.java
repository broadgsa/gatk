
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

public class ListSampleIds extends LocusWalker<Boolean, Boolean> 
{
    public void initialize() 
    { 
        GenomeAnalysisTK toolkit = this.getToolkit();
        SAMFileHeader header = toolkit.getSamReader().getFileHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
            out.println(sample_name);
        } 
    }

    public Boolean map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
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
        System.exit(0);
        return true;
    }
}
