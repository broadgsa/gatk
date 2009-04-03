
package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.List;

public class PoolCallingExperiment extends LocusWalker<AlleleFrequencyEstimate, String> 
{
    List<AlleleFrequencyWalker> deep_callers;
    List<AlleleFrequencyWalker> shallow_callers;
    AlleleFrequencyWalker       pooled_caller;

    @Argument public int DOWNSAMPLE;

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        for (int i = 0; i < context.getReads().size(); i++)
        {
            String read_group = (String)(context.getReads().get(i).getAttribute("RG"));
            String sample     = context.getReads().get(i).getHeader().getReadGroup(read_group).READ_GROUP_SAMPLE_TAG;
            System.out.println("RG: " + read_group + " SAMPLE: " + sample);
        } 
        return null;
    }

    public void onTraversalDone() 
    {
        return;
    }

    public String reduceInit() 
    { 
        return "";
    }

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
    {
        return "";
    }

    
}
