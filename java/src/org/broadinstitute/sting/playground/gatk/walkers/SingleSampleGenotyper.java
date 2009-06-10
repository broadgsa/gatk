package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods.IndelCall;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import net.sf.samtools.*;

import java.io.*;
import java.util.*;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, String>  
{
    @Argument(fullName="calls",        shortName="calls",        doc="File to write SNP calls to",  required=true)  public String callsFileName;
    @Argument(fullName="metrics",      shortName="met",          doc="metrics",                     required=false) public String metricsFileName = "/dev/null";
    @Argument(fullName="metInterval",  shortName="mi",           doc="metInterval",                 required=false) public Integer metricsInterval = 50000;
    @Argument(fullName="printMetrics", shortName="printMetrics", doc="printMetrics",                required=false) public Boolean printMetrics = true;
    @Argument(fullName="lodThreshold", shortName="lod",          doc="lodThreshold",                required=false) public Double lodThreshold = 5.0;
    @Argument(fullName="fourBaseMode", shortName="fb",           doc="fourBaseMode",                required=false) public Boolean fourBaseMode = false;
    @Argument(fullName="call_indels",  shortName="call_indels",  doc="Call Indels",                 required=false) public Boolean call_indels = false;
    @Argument(fullName="refPrior", shortName="refPrior", required=false, doc="Prior likelihood of the reference theory") public double REF_PRIOR = 0.999;
    @Argument(fullName="hetVarPrior", shortName="hetVarPrior", required=false, doc="Prior likelihood of a heterozygous variant theory") public double HETVAR_PRIOR = 1e-3;
    @Argument(fullName="homVarPrior", shortName="homVarPrior", required=false, doc="Prior likelihood of the homozygous variant theory") public double HOMVAR_PRIOR = 1e-5;
    @Argument(fullName="geli",         shortName="geli",         doc="If true, output will be in Geli/Picard format", required=false) public boolean GeliOutputFormat = false;
    @Argument(fullName="sample_name_regex", shortName="sample_name_regex", doc="sample_name_regex", required=false) public String SAMPLE_NAME_REGEX = null;

    public AlleleMetrics metrics;
	public PrintStream calls_file;
	public String sample_name;
    
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) { return true; }

    public boolean requiresReads() { return true; }

    public void initialize() 
	{ 
		try 
		{
			sample_name = null; 
			if (metricsFileName != null) { metrics    = new AlleleMetrics(metricsFileName, lodThreshold); }
			if (callsFileName != null)   
			{ 
				calls_file = new PrintStream(callsFileName);
                String header = GeliOutputFormat ? AlleleFrequencyEstimate.geliHeaderString() : AlleleFrequencyEstimate.asTabularStringHeader();
				calls_file.println(header);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace(); 
			System.exit(-1);
		}
	}

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        String rodString = getRodString(tracker);

		ref = Character.toUpperCase(ref);
		if (ref == 'N') { return null; }

		if (context.getReads().size() != 0)
		{
			SAMRecord read = context.getReads().get(0);
	        String RG = (String)(read.getAttribute("RG"));
			SAMReadGroupRecord read_group_record = read.getHeader().getReadGroup(RG);
			if (read_group_record != null)
			{
				String local_sample_name = read.getHeader().getReadGroup(RG).getSample();
				if (SAMPLE_NAME_REGEX != null) { local_sample_name = local_sample_name.replaceAll(SAMPLE_NAME_REGEX, "$1"); }
				if (sample_name == null) { sample_name = local_sample_name; }
				else 
				{ 
					if (! sample_name.equals(local_sample_name)) { System.out.printf("SAMPLE NAME MIXUP: %s vs. %s\n", sample_name, local_sample_name); }
					assert(sample_name.equals(local_sample_name)); 
				}
			}
		}

        AlleleFrequencyEstimate freq = getAlleleFrequency(ref, context, rodString, sample_name);

        if (printMetrics) {
            if (freq != null) { metrics.nextPosition(freq, tracker); }
            metrics.printMetricsAtLocusIntervals(metricsInterval);
        }

        return freq;
    }

    private AlleleFrequencyEstimate getAlleleFrequency(char ref, LocusContext context, String rodString, String sample_name) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();

		if (bases.length() == 0)
		{
	        GenotypeLikelihoods G = new GenotypeLikelihoods();
	        return G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods, sample_name);
		}

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        ref = Character.toUpperCase(ref);

		// Handle indels.
		if (call_indels)
		{
			String[] indels = BasicPileup.indelPileup(reads, offsets);
			IndelCall indel_call = GenotypeLikelihoods.callIndel(indels);
			if (indel_call != null)
			{
				if (! indel_call.type.equals("ref"))
				{ 
					System.out.printf("INDEL %s %s\n", context.getLocation(), indel_call); 
				}
			}
		}
        
		// Handle single-base polymorphisms.
        GenotypeLikelihoods G = new GenotypeLikelihoods();
        for ( int i = 0; i < reads.size(); i++ )  
        {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }
        G.ApplyPrior(ref, this.alt_allele, this.allele_frequency_prior);

        if (fourBaseMode && pileup.getBases().length() < 750) {
            G.applySecondBaseDistributionPrior(pileup.getBases(), pileup.getSecondaryBasePileup());
        }

        return G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods, sample_name);
    }

    private String getRodString(RefMetaDataTracker tracker) {
        String rodString = "";

        for ( ReferenceOrderedDatum datum : tracker.getAllRods() )  {
            if ( datum != null )  {
                if ( datum instanceof rodDbSNP) {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    rodString += dbsnp.toString();
                } else  {
                    rodString += datum.toSimpleString();
                }
            }
        }
        
        if ( rodString != "" ) { rodString = "[ROD: " + rodString + "]"; }
        return rodString;
    }

    double allele_frequency_prior = -1;
	char alt_allele;
    public void setAlleleFrequencyPrior(double freq, char alt)
    {
        this.allele_frequency_prior = freq;
		this.alt_allele = alt;
    }

    // Given result of map function
    private String  confident_ref_interval_contig  = "";
    private long    confident_ref_interval_start   = 0;
    private double  confident_ref_interval_LOD_sum = 0;
    private double  confident_ref_interval_length  = 0;
    private long    last_position_considered       = -1;
    private boolean inside_confident_ref_interval  = false;

    public String reduceInit() 
	{ 
        confident_ref_interval_contig  = "";
        confident_ref_interval_start   = 0;
        confident_ref_interval_LOD_sum = 0;
        confident_ref_interval_length  = 0;
        last_position_considered       = -1;
        inside_confident_ref_interval  = false;
        return "";
    }

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
	{
        if ( alleleFreq.lodVsRef >= lodThreshold ) {
            String line = GeliOutputFormat ? alleleFreq.asGeliString() : alleleFreq.asTabularString();
		    calls_file.println(line);
        }
        return "";
	}

	public void onTraversalDone()
	{
		if (callsFileName != null)   
		{ 
			calls_file.close();
		}
	}


	/*
    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
	{
        // Print RESULT data for confident calls

       long current_offset = alleleFreq.location.getStart(); //Integer.parseInt(tokens[1]);

        if (inside_confident_ref_interval &&
                ((alleleFreq.lodVsRef > -5.0) || (current_offset != last_position_considered + 1)) )
        {
            // No longer hom-ref, so output a ref line.
            double lod = confident_ref_interval_LOD_sum / confident_ref_interval_length;

            calls_file.format("%s\tCALLER\tREFERENCE\t%d\t%d\t%f\t.\t.\tLENGTH %d\n",
                            confident_ref_interval_contig,
                            confident_ref_interval_start,
                            last_position_considered,
                            lod,
                            (int)(confident_ref_interval_length));

            inside_confident_ref_interval = false;
        }
        else if (inside_confident_ref_interval && (alleleFreq.lodVsRef <= -5.0))
        {
            // Still hom-ref so increment the counters.
            confident_ref_interval_LOD_sum += alleleFreq.lodVsRef;
            confident_ref_interval_length  += 1;
        }
        else if ((!inside_confident_ref_interval) && (alleleFreq.lodVsRef > -5.0))
        {
            // do nothing.
        }
        else if ((!inside_confident_ref_interval) && (alleleFreq.lodVsRef <= -5.0))
        {
            // We moved into a hom-ref region so start a new interval.
            confident_ref_interval_contig  = alleleFreq.location.getContig();
            confident_ref_interval_start   = alleleFreq.location.getStart();
            confident_ref_interval_LOD_sum = alleleFreq.lodVsRef;
            confident_ref_interval_length  = 1;
            inside_confident_ref_interval  = true;
        }

        last_position_considered = current_offset;
        
        if (alleleFreq.lodVsRef >= 5) {
            calls_file.print(alleleFreq.asGFFString());

            //String gtype = genotypeTypeString(alleleFreq.qstar, alleleFreq.N);
            //System.out.print("DEBUG " + gtype + " ");
            //if (gtype.contentEquals("het")) {
            //    System.out.println(alleleFreq.ref + "" + alleleFreq.alt);
            //} else if (gtype.contentEquals("hom")) {
            //    System.out.println(alleleFreq.ref + "" + alleleFreq.ref);
            //} else {
            //    System.out.println(alleleFreq.alt + "" + alleleFreq.alt);
            //}
        }

        return "";
    }
	*/
}
