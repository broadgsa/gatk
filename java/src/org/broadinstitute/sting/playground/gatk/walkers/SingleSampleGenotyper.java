package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

import org.broadinstitute.sting.playground.utils.*;

import net.sf.samtools.SAMRecord;

import java.util.List;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, Integer> 
{
    AlleleMetrics metrics;
    
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        return true;    // We are keeping all the reads
    }

    public boolean requiresReads()     { return true; }    

    public void initialize()
    {
        metrics = new AlleleMetrics("metrics.out");
    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context)
    {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        String bases = "";
        String quals = "";

        ref = Character.toUpperCase(ref);

        String rodString = "";
        // Look up dbsnp priors
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) 
        {
            if ( datum != null ) 
            {
                if ( datum instanceof rodDbSNP)
                {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    rodString += dbsnp.toString();
                }
                else 
                {
                    rodString += datum.toSimpleString();
                }
            }
        }
        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        // Accumulate genotype likelihoods
        GenotypeLikelihoods G = new GenotypeLikelihoods(); 
        for ( int i = 0; i < reads.size(); i++ ) 
        {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            bases += read.getReadString().charAt(offset);
            quals += read.getBaseQualityString().charAt(offset);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }
        G.ApplyPrior(ref, Double.NaN);

        System.out.printf("%s %s %s %s\n", context.getLocation(), ref, bases, G.toString(ref), rodString);

        AlleleFrequencyEstimate freq = G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods);

        metrics.nextPosition(freq, tracker);
        metrics.printMetricsAtLocusIntervals(1000);

        return freq;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(AlleleFrequencyEstimate value, Integer sum) 
    {
        return 0;
    }
}
