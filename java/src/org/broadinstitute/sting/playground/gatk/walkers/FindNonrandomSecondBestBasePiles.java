package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pileup;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;

import java.util.List;
import java.util.ArrayList;

public class FindNonrandomSecondBestBasePiles extends LocusWalker<Integer, Integer> {
    @Argument(fullName="verbose",doc="verbose",required=false,defaultValue="false")
    public boolean VERBOSE;
    
    private AlleleFrequencyWalker caller_1b;
    private AlleleFrequencyWalker caller_4b;
    public void initialize() 
    {
        caller_1b = new AlleleFrequencyWalker();
        caller_1b.N = 2;
        caller_1b.DOWNSAMPLE = 0;
        caller_1b.GFF_OUTPUT_FILE = "/dev/null";
        caller_1b.FORCE_1BASE_PROBS = true; 
        caller_1b.initalize();
        caller_1b.reduceInit();

        caller_4b = new AlleleFrequencyWalker();
        caller_4b.N = 2;
        caller_4b.DOWNSAMPLE = 0;
        caller_4b.GFF_OUTPUT_FILE = "/dev/null";
        caller_4b.FORCE_1BASE_PROBS = false; 
        caller_4b.initalize();
        caller_4b.reduceInit();
    }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);

        // First, call the site, because we're mostly interested in hets
        ref = Character.toUpperCase(ref);
        AlleleFrequencyEstimate call_1b = caller_1b.map(tracker, ref, context);
        AlleleFrequencyEstimate call_4b = caller_4b.map(tracker, ref, context);

        // Only output data for sites that are confident disagreements.
        if (Math.abs(call_1b.lodVsRef) < 5.0) { return 0; }
        if (Math.abs(call_4b.lodVsRef) < 5.0) { return 0; }
        if (call_1b.genotype().equals(call_4b.genotype())) { return 0; }

        List<SAMRecord> reads = pileup.getReads();
        List<Integer> offsets = pileup.getOffsets();
        String best_bases = pileup.getBases();
        String second_bases;
        double[] quals = new double[reads.size()];
        double[] second_quals = new double[reads.size()];
        {
            char[] second_bases_array = new char[reads.size()];
	        for ( int i = 0; i < reads.size(); i++ ) 
            {
	            SAMRecord read = reads.get(i);
	            int offset = offsets.get(i);
	            byte qual = (byte)read.getBaseQualities()[offset];
                quals[i] = 1.0 - Math.pow(10,(double)qual/-10.0);

                second_bases_array[i] = BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(((byte[])read.getAttribute("SQ"))[offset]));
                second_quals[i] = QualityUtils.compressedQualityToProb(((byte[])read.getAttribute("SQ"))[offset]);
	        }
            second_bases = new String(second_bases_array);
        }
        
        String rodString = "";
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if ( datum != null && ! (datum instanceof rodDbSNP)) {
                //System.out.printf("rod = %s%n", datum.toSimpleString());
                rodString += datum.toSimpleString();
                //System.out.printf("Rod string %s%n", rodString);
            }
        }
        
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp != null )
            rodString += dbsnp.toMediumString();

        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        char[] bases = {'A', 'C', 'G', 'T'};
        double[][] counts = new double[4][];
        double[] totals = new double[4];
        double[][] fractional_counts = new double[4][];
        double[] fractional_totals = new double[4];
        for (int i = 0; i < 4; i++)
        {
            counts[i] = new double[4];
            fractional_counts[i] = new double[4];
            for (int j = 0; j < best_bases.length(); j++)
            {
                if (best_bases.charAt(j) == bases[i]) 
                {
                    counts[BaseUtils.simpleBaseToBaseIndex(bases[i])][BaseUtils.simpleBaseToBaseIndex(second_bases.charAt(j))] += 1;
                    totals[BaseUtils.simpleBaseToBaseIndex(bases[i])] += 1;

                    fractional_counts[BaseUtils.simpleBaseToBaseIndex(bases[i])][BaseUtils.simpleBaseToBaseIndex(second_bases.charAt(j))] += second_quals[j];
                    fractional_totals[BaseUtils.simpleBaseToBaseIndex(bases[i])] += second_quals[j];
                }
            }
        }

        out.printf("%s\t%c\t%s\n%s\t%c\t%s\n", 
                        pileup.getLocation(), 
                        ref,
                        best_bases,
                        pileup.getLocation(),
                        ref,
                        second_bases);                    
        out.printf("%s\t%s\t%f\t%f\n", 
                        pileup.getLocation(),
                        call_1b.genotype(),
                        call_1b.lodVsRef,
                        call_1b.lodVsNextBest);
        out.printf("%s\t%s\t%f\t%f\n", 
                        pileup.getLocation(),
                        call_4b.genotype(),
                        call_4b.lodVsRef,
                        call_4b.lodVsNextBest);
        for (int i = 0; i < quals.length; i++)
        {
            out.printf("%s\t%c %c %f\t%f\t%f\t%f\t%f\n", 
                        pileup.getLocation(),
                        best_bases.charAt(i), 
                        second_bases.charAt(i), 
                        quals[i],
                        second_quals[i],
                        -10.0 * Math.log10(1.0 - quals[i]),
                        -10.0 * Math.log10(1.0 - second_quals[i]),
                        Math.log10((quals[i])/(second_quals[i])));
        }
        out.printf("\n");
        for (int i = 0; i < 4; i++)
        {
            out.printf("%s\t%c ", pileup.getLocation(), bases[i]);
            for (int j = 0; j < 4; j++)
            {
                out.printf("%.03f ", counts[i][j] / totals[i]);
            }
            out.printf("\n");
        }
        /*
        out.printf("\n");
        for (int i = 0; i < 4; i++)
        {
            out.printf("%s\t%c ", pileup.getLocation(), bases[i]);
            for (int j = 0; j < 4; j++)
            {
                out.printf("%.03f ", fractional_counts[i][j] / fractional_totals[i]);
            }
            out.printf("\n");
        }
        */
        out.printf("\n\n");

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
