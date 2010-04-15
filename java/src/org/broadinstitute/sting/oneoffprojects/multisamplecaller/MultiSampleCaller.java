
package org.broadinstitute.sting.oneoffprojects.multisamplecaller;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.zip.GZIPOutputStream;

// Beta iterative multi-sample caller
// j.maguire 6-11-2009

public class MultiSampleCaller extends LocusWalker<MultiSampleCaller.MultiSampleCallResult,String>
{
    @Argument(required=false, shortName="fractional_counts", doc="should we use fractional counts?") public boolean FRACTIONAL_COUNTS = false;
    @Argument(required=false, shortName="max_iterations", doc="Maximum number of iterations for EM") public int MAX_ITERATIONS = 10;
    @Argument(fullName="discovery_output", shortName="discovery_output", required=true, doc="file to write SNP discovery output to")       public String DISCOVERY_OUTPUT = null;
    @Argument(fullName="individual_output", shortName="individual_output", required=false, doc="file to write individual SNP calls to") public String INDIVIDUAL_OUTPUT = null;
    @Argument(fullName="stats_output", shortName="stats_output", required=false, doc="file to write stats to") public String STATS_OUTPUT = null;
    @Argument(fullName="sample_name_regex", shortName="sample_name_regex", required=false, doc="sample_name_regex") public String SAMPLE_NAME_REGEX = null;
    @Argument(fullName="call_indels", shortName="call_indels", required=false, doc="call indels?") public boolean CALL_INDELS = false;
    @Argument(fullName="weight_samples", shortName="weight_samples", required=false, doc="rw-weight samples during EM?") public boolean WEIGHT_SAMPLES = false;

    @Argument(fullName="theta", shortName="theta", required=false, doc="rate of sequence divergence") public double THETA = 1e-3;
    @Argument(fullName="allele_frequency_prior", shortName="allele_frequency_prior", required=false, doc="use prior on allele frequencies? (P(f) = theta/(N*f)") public boolean ALLELE_FREQUENCY_PRIOR = false;

    @Argument(fullName="confusion_matrix_file", shortName="confusion_matrix_file", required=false, doc="file containing confusion matrix for all three technologies") public String CONFUSION_MATRIX_FILE = null;

    @Argument(fullName="ALLELE_FREQ_TOLERANCE", shortName="AFT", required=false, doc="") 
    public double ALLELE_FREQ_TOLERANCE = 1e-6;

    @Argument(fullName="append", shortName="append", required=false, doc="if the discovery file already exists, don't re-call sites that are done.") boolean APPEND = false;

    private static final double MIN_LOD_FOR_STRAND = 0.01;

	// Private state.
    protected List<String> sample_names;
    protected SAMFileHeader header;
	protected PrintStream individual_output_file = null;
	protected PrintStream discovery_output_file = null;
	protected PrintStream stats_output_file = null;

    private boolean INCLUDE_STATS = false;
    private boolean INCLUDE_GENOTYPES = false ;

	class MultiSampleCallResult
	{
		GenomeLoc location;
		char ref;
		char alt;
		EM_Result em_result;
		double lod;
		double strand_score;
		double pD;
		double pNull;
		String in_dbsnp;
		int n_ref;
		int n_het;
		int n_hom;
		int EM_N;
		double alt_freq;

		public MultiSampleCallResult(GenomeLoc location, char ref, char alt, EM_Result em_result, double lod, double strand_score, double pD, double pNull, String in_dbsnp, int n_ref, int n_het, int n_hom, int EM_N, double alt_freq)
		{
			this.location = location;
			this.ref = ref;
			this.alt = alt;
			this.em_result = em_result;
			this.lod = lod;
			this.strand_score = strand_score; 
			this.pD = pD;
			this.pNull = pNull;
			this.in_dbsnp = in_dbsnp;
			this.n_ref = n_ref;
			this.n_het = n_het;
			this.n_hom = n_hom;
			this.EM_N = EM_N;
			this.alt_freq = alt_freq;
		}

		public MultiSampleCallResult() { } // this is just so I can do new MultiSampleCallResult().header(). "inner classes cannot have static declarations" :(

		public String header()
		{
			return new String("loc ref alt lod strand_score pD pNull in_dbsnp pA pC pG pT EM_alt_freq EM_N n_ref n_het n_hom");
		}

		public String toString()
		{
			String s = "";
			s = s + String.format("%s %c %c %f %f %f %f %s ", location, ref, alt, lod, strand_score, pD, pNull, in_dbsnp);
			for (int i = 0; i < 4; i++) { s = s + String.format("%f ", em_result.allele_likelihoods[i]); }
			s = s + String.format("%f %d %d %d %d", alt_freq, em_result.EM_N, n_ref, n_het, n_hom);
			return s;
		}
	}

	public static class DepthStats
	{
		public static String Header()
		{
			return "loc ref depth A C G T a c g t mq_min mq_mean mq_median mq_max mq_sd";
		}

		public static String Row(char ref, AlignmentContext context)
		{
			String ans = "";
			List<SAMRecord> reads = context.getReads();
			List<Integer> offsets = context.getOffsets();
			//Pileup pileup = new ReadBackedPileupOld(ref, context);
			
			ans += String.format("%s ", context.getLocation());
			ans += String.format("%c ", ref);
			ans += String.format("%d ", reads.size());
			ans += String.format("%d ", countBase(context, 'A', "+"));
			ans += String.format("%d ", countBase(context, 'C', "+"));
			ans += String.format("%d ", countBase(context, 'G', "+"));
			ans += String.format("%d ", countBase(context, 'T', "+"));
			ans += String.format("%d ", countBase(context, 'A', "-"));
			ans += String.format("%d ", countBase(context, 'C', "-"));
			ans += String.format("%d ", countBase(context, 'G', "-"));
			ans += String.format("%d ", countBase(context, 'T', "-"));

			ans += String.format("%s ", Stats(BasicPileup.mappingQualPileup(reads)));

			return ans;
		}

		static int countBase(AlignmentContext context, char base, String strand)
		{
			int count =  0;
			List<SAMRecord> reads = context.getReads();
			List<Integer> offsets = context.getOffsets();
			for (int i = 0; i < reads.size(); i++)
			{
				if (reads.get(i).getReadString().charAt(offsets.get(i)) == base)
				{
					if (strand.equals("+") && (reads.get(i).getReadNegativeStrandFlag()==false)) { count += 1; }
					else if (strand.equals("-") && (reads.get(i).getReadNegativeStrandFlag()==true)) { count += 1; }
					else if (! (strand.equals("+") || strand.equals("-"))) { count += 1; }
				}
			}
			return count;
		}

		public static String Stats(ArrayList<Byte> X)
		{
			Collections.sort(X);

			long count = 0;
			long sum = 0;
			long min = X.get(0);
			long max = X.get(0);
			long median = X.get(0);
			for (int i = 0; i < X.size(); i++)
			{
				int x = X.get(i);
				if (x < min) { min = x; }
				if (x > max) { max = x; }
				sum += x;
				count += 1;
				if (i == X.size()/2) { median = x; }
			}

			double mean = sum/count;
			for (int i = 0; i < X.size(); i++)
			{
				int x = X.get(i);
				sum += (x-mean)*(x-mean);
				count += 1;
			}
			double sd = Math.sqrt(sum/count);

			return String.format("%d %f %d %d %f", min, mean, median, max, sd);
		}
	}

	GenomeLoc highest_previously_done_loc = null; 
	private boolean in_skip_mask(GenomeLoc loc)
	{
		if (highest_previously_done_loc == null) { return false; }
		if (highest_previously_done_loc.compareTo(loc) < 0) { return false; }
		else { return true; }
	}

    private void maybeInitializeDisoveryOutput() 
	{
        if ( discovery_output_file == null ) 
		{
			File file = new File(DISCOVERY_OUTPUT);
			if ((APPEND == true) && (file.exists()))
			{
				try
				{
					Runtime.getRuntime().exec("cp -v " + DISCOVERY_OUTPUT + " " + DISCOVERY_OUTPUT + ".backup");

					// 1. Read the existing file and record the highest site we've seen.
					Scanner scanner = new Scanner(file);
					while(scanner.hasNext())
					{
						String line = scanner.nextLine();
						String[] tokens = line.split(" +");
						String loc_string = tokens[0];
						if (loc_string.equals("loc")) { continue; }
						highest_previously_done_loc = GenomeLocParser.parseGenomeLoc(loc_string);
					}
	
					// 2. Open the output file for appending.
					discovery_output_file = new PrintStream(new FileOutputStream(DISCOVERY_OUTPUT, true));			
				}
				catch (Exception e) 
				{
	                throw new RuntimeException(e);
	            }
			}
			else
			{
	            try 
				{
	                final String filename = DISCOVERY_OUTPUT;
	                discovery_output_file = new PrintStream(filename);
	                discovery_output_file.println(new MultiSampleCallResult().header());
	            } 
				catch (Exception e) 
				{
	                throw new RuntimeException(e);
	            }
			}
        }
    }

	/////////
	// Walker Interface Functions 
    public void initialize() 
	{ 
		System.out.printf("\n\n\n\n");
		(new ClassicGenotypeLikelihoods()).TEST();

		try
		{
            maybeInitializeDisoveryOutput();

            INCLUDE_GENOTYPES = INDIVIDUAL_OUTPUT != null;
            if ( INCLUDE_GENOTYPES ) {
                individual_output_file = new PrintStream(new GZIPOutputStream(new FileOutputStream(INDIVIDUAL_OUTPUT)));
                individual_output_file.println("loc ref sample_name genotype lodVsNextBest lodVsRef in_dbsnp AA AC AG AT CC CG CT GG GT TT");
            }

            INCLUDE_STATS = STATS_OUTPUT != null;
            if ( INCLUDE_STATS ) {
                stats_output_file = new PrintStream(STATS_OUTPUT);
                stats_output_file.println(DepthStats.Header());
            }
            }
		catch (Exception e)
		{
			e.printStackTrace(); 
			System.exit(-1);
		}

		GenomeAnalysisEngine toolkit = this.getToolkit();
		this.header = toolkit.getSAMFileHeader();
		List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        sample_names    = new ArrayList<String>();

		HashSet<String> unique_sample_names = new HashSet<String>();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
            String platform = (String)(read_groups.get(i).getAttribute(SAMReadGroupRecord.PLATFORM_TAG));

			if (SAMPLE_NAME_REGEX != null) { sample_name = sample_name.replaceAll(SAMPLE_NAME_REGEX, "$1"); }

            //System.out.printf("SAMPLE: %s %s\n", sample_name, platform);

			if (unique_sample_names.contains(sample_name)) { continue; }
			unique_sample_names.add(sample_name);
            sample_names.add(sample_name);

            System.out.printf("UNIQUE_SAMPLE: %s %s\n", sample_name, platform);
        } 

		// Load the confusion matrix if it exists
		if (CONFUSION_MATRIX_FILE != null)
		{
			this.confusion_matrix = new ConfusionMatrix(CONFUSION_MATRIX_FILE);
		}

    }

	
	Date start_time = null;
	int n_sites_processed = 0;

    public MultiSampleCallResult map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
	{
		if (start_time == null) { start_time = new Date(); }

		if (in_skip_mask(context.getLocation()) == true)  { return null; } 

		if ((n_sites_processed % 1000) == 0)
		{
			Date current_time = new Date();
			long elapsed = current_time.getTime() - start_time.getTime();
			out.printf("RUNTIME: %d sites processed in %f seconds; %f seconds per site.\n",
						n_sites_processed,
						(double)elapsed/1000.0,
						((double)elapsed/1000.0)/(double)n_sites_processed);
		}
		n_sites_processed += 1;

		context = filter_each_read(context);

		if (ref.getBase() == 'N') { return null; }
		if (BaseUtils.simpleBaseToBaseIndex(ref.getBase()) == -1) { return null; }
		if (context.getReads().size() <= 0) { return null; }
        if (context.getReads().size() >= 10000) { return null; }     // to deal with big piles -- totally arbitrary threshold

		this.ref = ref.getBase();
		MultiSampleCallResult result = this.MultiSampleCall(tracker, ref.getBase(), context, sample_names);
        if ( INCLUDE_STATS ) stats_output_file.println(DepthStats.Row(ref.getBase(), context));
		return result;
	}

    public void onTraversalDone(String sum) 
	{
		discovery_output_file.flush();
		discovery_output_file.close();

        if ( INCLUDE_STATS ) {
            stats_output_file.flush();
            stats_output_file.close();
        }

		out.println("MultiSampleCaller done.");
		return;
	}

    public String reduceInit() 
	{
		return null;
	}

    public String reduce(MultiSampleCallResult record, String sum) 
	{
		if (record != null)
		{
			discovery_output_file.printf(record.toString() + "\n");
		}
		return null;
	}

	// END Walker Interface Functions
	/////////


	/////////
	// Calling Functions

	char ref;
	protected ConfusionMatrix confusion_matrix;

	ClassicGenotypeLikelihoods reallyMakeGenotypeLikelihood(AlignmentContext context) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

		// Handle single-base polymorphisms.
        ClassicGenotypeLikelihoods G = new ClassicGenotypeLikelihoods();
        for ( int i = 0; i < reads.size(); i++ )  
        {
			//System.out.printf("DBG: %s\n", context.getLocation());

            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
			if (CONFUSION_MATRIX_FILE == null)
			{
            	G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
			}
			else
			{
            	String RG = (String)(read.getAttribute("RG"));

            	assert(header != null);
            	assert(header.getReadGroup(RG) != null);

                String platform = (String)(header.getReadGroup(RG).getAttribute(SAMReadGroupRecord.PLATFORM_TAG));

            	G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset], confusion_matrix, platform);
			}
        }
        
        return G;
	}
	
	HashMap<AlignmentContext, ClassicGenotypeLikelihoods> glCache = new HashMap<AlignmentContext, ClassicGenotypeLikelihoods>();
	int cache_size = 0;

    ClassicGenotypeLikelihoods GenotypeOld(AlignmentContext context, double[] allele_likelihoods, double indel_alt_freq) {
        //ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        //String bases = pileup.getBases();

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        ref = Character.toUpperCase(ref);

        if (reads.size() == 0) { 
            ClassicGenotypeLikelihoods G = new ClassicGenotypeLikelihoods();
            return G;
        }

        // Handle single-base polymorphisms.
        ClassicGenotypeLikelihoods G = new ClassicGenotypeLikelihoods();
        for ( int i = 0; i < reads.size(); i++ )  
        {
                        //System.out.printf("DBG: %s\n", context.getLocation());

            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
                        if (CONFUSION_MATRIX_FILE == null)
                        {
                			G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
                        }
                        else
                        {
			                String RG = (String)(read.getAttribute("RG"));
			
			                assert(header != null);
			                assert(header.getReadGroup(RG) != null);
			
			                String platform = (String)(header.getReadGroup(RG).getAttribute(SAMReadGroupRecord.PLATFORM_TAG));
			
			                G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset], confusion_matrix, platform);
                        }
        }
        G.ApplyPrior(ref, allele_likelihoods);

                // Handle indels
                if (CALL_INDELS)
                {
                        String[] indels = BasicPileup.indelPileup(reads, offsets);
                        IndelLikelihood indel_call = new IndelLikelihood(indels, indel_alt_freq);
                        if (indel_call.getType() != null)
                        {
                                G.addIndelLikelihood(indel_call);
                        }
                        else
                        {
                                G.addIndelLikelihood(null);
                        }
                }

        return G;
    }

    ClassicGenotypeLikelihoods Genotype(AlignmentContext context, double[] allele_likelihoods, double indel_alt_freq) {
        return GenotypeCache(context, allele_likelihoods, indel_alt_freq );
        //return GenotypeOld(context, allele_likelihoods, indel_alt_freq );
    }        


	ClassicGenotypeLikelihoods GenotypeCache(AlignmentContext context, double[] allele_likelihoods, double indel_alt_freq)
	{
        ref = Character.toUpperCase(ref);

		// Handle single-base polymorphisms.
        ClassicGenotypeLikelihoods G = null;
		if ( context.getReads().size() == 0 ) {
	        G = new ClassicGenotypeLikelihoods();
	        return G;
	    } else {
             if ( true && glCache.containsKey(context) ) {
                ClassicGenotypeLikelihoods cached = glCache.get(context); 
                G = (ClassicGenotypeLikelihoods)cached.clone();
            } else {
                G = reallyMakeGenotypeLikelihood(context);
				if (cache_size < 5000)
				{
					//System.out.printf("cache add (%d)\n", cache_size);
                	glCache.put(context, G.clone());
					cache_size += context.getReads().size();
				}
				else
				{
					//System.out.printf("cache skip (%d)\n", cache_size);
				}
            }
            G.ApplyPrior(ref, allele_likelihoods);
        }
        
        return G;
    }

	// This version is a little faster. 
	double[] CountFreqs(ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{
		double[] allele_likelihoods = new double[4];
		
		for (int x = 0; x < genotype_likelihoods.length; x++)
		{
			ClassicGenotypeLikelihoods G = genotype_likelihoods[x];

			if (G.coverage == 0) { continue; }

        	double[] personal_allele_likelihoods = new double[4];
			
			int k = 0;
			for (int i = 0; i < 4; i++)
			{ 
				for (int j = i; j < 4; j++)
				{
					double likelihood = Math.pow(10,G.likelihoods[k]);
					personal_allele_likelihoods[i] += likelihood;
					personal_allele_likelihoods[j] += likelihood;
					k++;
				}
			}
			double sum = 0;
			for (int y = 0; y < 4; y++) { sum += personal_allele_likelihoods[y]; }
			for (int y = 0; y < 4; y++) { personal_allele_likelihoods[y] /= sum; }
			for (int y = 0; y < 4; y++) { allele_likelihoods[y] += personal_allele_likelihoods[y]; }
		}

		double sum = 0;
		for (int i = 0; i < 4; i++) { sum += allele_likelihoods[i]; }
		for (int i = 0; i < 4; i++) { allele_likelihoods[i] /= sum; }

		return allele_likelihoods;
	}


	double CountIndelFreq(ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{ 
		HashMap<String, Double> indel_allele_likelihoods = new HashMap<String, Double>();

		double pRef = 0;
		double pAlt = 0;

		for (int j = 0; j < sample_names.size(); j++)
		{
			double personal_pRef = 0;
			double personal_pAlt = 0;

			IndelLikelihood indel_likelihood = genotype_likelihoods[j].getIndelLikelihood();
			personal_pRef += 2*Math.pow(10, indel_likelihood.pRef()) + Math.pow(10, indel_likelihood.pHet());
			personal_pAlt += 2*Math.pow(10, indel_likelihood.pHom()) + Math.pow(10, indel_likelihood.pHet());

			personal_pRef = personal_pRef / (personal_pAlt + personal_pRef);
			personal_pAlt = personal_pAlt / (personal_pAlt + personal_pRef);

			pRef += personal_pRef;
			pAlt += personal_pAlt;
		}

		pRef = pRef / (pRef + pAlt);
		pAlt = pAlt / (pRef + pAlt);

		return pAlt;
	}

	// Potential precision error here.
	double Compute_pD(ClassicGenotypeLikelihoods[] genotype_likelihoods, double[] sample_weights)
	{
		double pD = 0;
		for (int i = 0; i < sample_names.size(); i++)
		{
			double sum = 0;
			for (int j = 0; j < 10; j++)
			{
				sum += Math.pow(10, genotype_likelihoods[i].likelihoods[j]);
			}
			pD += Math.log10(sample_weights[i] * sum);
		}
		return pD;
	}

	double Compute_pNull(AlignmentContext[] contexts, double[] sample_weights)
	{
		double[] allele_likelihoods = new double[4];
		for (int i = 0; i < 4; i++) { allele_likelihoods[i] = 1e-6/3.0; }
		allele_likelihoods[BaseUtils.simpleBaseToBaseIndex(ref)] = 1.0-1e-6;
		ClassicGenotypeLikelihoods[] G = new ClassicGenotypeLikelihoods[sample_names.size()];
		for (int j = 0; j < sample_names.size(); j++)
		{
			G[j] = Genotype(contexts[j], allele_likelihoods, 1e-6);
		}
		return Compute_pD(G, sample_weights);
	}

	double[] Compute_SampleWeights(ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{
		double[] pD = new double[sample_names.size()];
		double total_pD = 0;
		for (int i = 0; i < sample_names.size(); i++)
		{
			double sum = 0;
			for (int j = 0; j < 10; j++)
			{
				sum += Math.pow(10, genotype_likelihoods[i].likelihoods[j]);
			}
			pD[i] = sum;
			total_pD += pD[i];
		}

		for (int i = 0; i < sample_names.size(); i++)
		{
			pD[i] /= total_pD;
		}

		return pD;
	}

	// Some globals to cache results.
	EM_Result em_result;
	double pD;
	double pNull;
	double lod;
	double LOD(AlignmentContext[] contexts)
	{
		em_result = EM(contexts);
		ClassicGenotypeLikelihoods[] G = em_result.genotype_likelihoods;
		pD = Compute_pD(G, em_result.sample_weights);
		pNull = Compute_pNull(contexts, em_result.sample_weights);

		if (ALLELE_FREQUENCY_PRIOR)
		{
			// Apply p(f).
			double pVar = 0.0;
			for (int i = 1; i < em_result.EM_N; i++) { pVar += THETA/(double)i; }

			double p0 = Math.log10(1 - pVar);
			double pF;

			double MAF = Compute_alt_freq(ref, em_result.allele_likelihoods);

			if (MAF < 1/(2.0*em_result.EM_N)) { pF = p0; }
			else { pF = Math.log10(THETA/(2.0*em_result.EM_N * MAF)); }

			//System.out.printf("DBG %s %c %f %f %f %f (%.20f) %f ", contexts[0].getLocation(), ref, pD, pF, pNull, p0, Compute_alt_freq(ref, em_result.allele_likelihoods), 2.0 * em_result.EM_N);
			//for (int i = 0; i < 4; i++) { System.out.printf("%f ", em_result.allele_likelihoods[i]); }
			//System.out.printf("\n");

			pD = pD + pF;
			pNull = pNull + p0;
		}

		lod = pD - pNull;
		return lod;
	}

	class EM_Result
	{
		String[] sample_names;
		ClassicGenotypeLikelihoods[] genotype_likelihoods;
		double[] allele_likelihoods;
		double[] sample_weights;
		int EM_N;

		public EM_Result(List<String> sample_names, ClassicGenotypeLikelihoods[] genotype_likelihoods, double[] allele_likelihoods, double[] sample_weights)
		{
			this.sample_names = new String[1];
			this.sample_names = sample_names.toArray(this.sample_names);
			this.genotype_likelihoods = genotype_likelihoods;
			this.allele_likelihoods = allele_likelihoods;
			this.sample_weights = sample_weights;

			EM_N = 0;
			for (int i = 0; i < genotype_likelihoods.length; i++) 
			{
				if (genotype_likelihoods[i].coverage > 0) { EM_N += 1; }
			}

		}

	}

    final static double[] sample_weights = new double[1000];
    static {
        for (int i = 0; i < 1000; i++)
        {
            //sample_weights[i] = 1.0/(double)i;
            sample_weights[i] = 1.0;
        }
    }

	EM_Result EM(AlignmentContext[] contexts)
	{
	    final boolean DEBUG_PRINT = false;
		double[] allele_likelihoods = new double[4];

		// These initial conditions should roughly replicate classic SSG. (at least on hets)
		for (int i = 0; i < 4; i++) 
		{ 
			if (i == BaseUtils.simpleBaseToBaseIndex(ref)) { allele_likelihoods[i] = 0.9994999; } //sqrt(0.999) 
			else { allele_likelihoods[i] = 0.0005002502; } // 0.001 / (2 * sqrt(0.999)
		}
		double indel_alt_freq = 1e-4;

		ClassicGenotypeLikelihoods[] G = new ClassicGenotypeLikelihoods[sample_names.size()];
		//ClassicGenotypeLikelihoods[] Weighted_G = new ClassicGenotypeLikelihoods[sample_names.size()];
	    
	    if ( DEBUG_PRINT ) System.out.printf("%n"); 
        
		for (int i = 0; i < MAX_ITERATIONS; i++)
		{
			for (int j = 0; j < sample_names.size(); j++)
			{
				G[j] = Genotype(contexts[j], allele_likelihoods, indel_alt_freq);
				//if (WEIGHT_SAMPLES) { G[j].ApplyWeight(sample_weights[j]); }
			}

            double[] old_allele_likelihoods = allele_likelihoods;
			allele_likelihoods = CountFreqs(G);
		    double alDelta = 0.0;
		    for (int j = 0; j < 4; j++) { alDelta += Math.abs(old_allele_likelihoods[j] - allele_likelihoods[j]); }

			if ( DEBUG_PRINT ) 
			{
			    System.out.printf("%s AL %f %f %f %f => delta=%e < %e == %b%n", 
			        contexts[0].getLocation(), 
			        Math.log10(allele_likelihoods[0]), Math.log10(allele_likelihoods[1]), Math.log10(allele_likelihoods[2]), Math.log10(allele_likelihoods[3]),
			        alDelta, ALLELE_FREQ_TOLERANCE, alDelta < ALLELE_FREQ_TOLERANCE);
			}

            //if ( alDelta < ALLELE_FREQ_TOLERANCE ) {
            //    System.out.printf("Converged after %d iterations%n", i);
            //    break;
            //}

// 			if (CALL_INDELS) 
// 			{
// 				indel_alt_freq = CountIndelFreq(G);
// 			}

// 			if (WEIGHT_SAMPLES)
// 			{
// 				sample_weights = Compute_SampleWeights(G);
// 			}
		}

		return new EM_Result(sample_names, G, allele_likelihoods, sample_weights);
	}

	// Hacky global variables for debugging.
	double StrandScore(AlignmentContext context)
	{
		//AlignmentContext[] contexts = filterAlignmentContext(context, sample_names, 0);

		AlignmentContext fw = filterAlignmentContext(context, "+");
		AlignmentContext bw = filterAlignmentContext(context, "-");
		AlignmentContext[] contexts_fw = filterAlignmentContext(fw, sample_names, 0);
		AlignmentContext[] contexts_bw = filterAlignmentContext(bw, sample_names, 0);

		EM_Result em_fw = EM(contexts_fw);
		EM_Result em_bw = EM(contexts_bw);

		double pNull_fw = Compute_pNull(contexts_fw, em_fw.sample_weights);
		double pNull_bw = Compute_pNull(contexts_bw, em_bw.sample_weights);

		double pD_fw = Compute_pD(em_fw.genotype_likelihoods, em_fw.sample_weights);
		double pD_bw = Compute_pD(em_bw.genotype_likelihoods, em_bw.sample_weights);

		if (ALLELE_FREQUENCY_PRIOR)
		{
			// Apply p(f).
			double pVar = 0.0;
			for (int i = 1; i < em_result.EM_N; i++) { pVar += THETA/(double)i; }

			pD_fw = pD_fw + Math.log10(THETA/(2.0*em_fw.EM_N * Compute_alt_freq(ref, em_fw.allele_likelihoods)));
			pNull_fw = pNull_fw + Math.log10(1 - pVar);

			pD_bw = pD_bw + Math.log10(THETA/(2.0*em_bw.EM_N * Compute_alt_freq(ref, em_bw.allele_likelihoods)));
			pNull_bw = pNull_bw + Math.log10(1 - pVar);
		}

		double EM_alt_freq_fw = Compute_alt_freq(ref, em_fw.allele_likelihoods);
		double EM_alt_freq_bw = Compute_alt_freq(ref, em_bw.allele_likelihoods);

		//double pNull = Compute_pNull(contexts);
		//double lod = LOD(contexts);
	
		double lod_fw = (pD_fw + pNull_bw) - pNull;
		double lod_bw = (pD_bw + pNull_fw) - pNull;
		double strand_score = Math.max(lod_fw - lod, lod_bw - lod);
		return strand_score;
	}

// 	ClassicGenotypeLikelihoods HardyWeinberg(double[] allele_likelihoods)
// 	{
// 		ClassicGenotypeLikelihoods G = new ClassicGenotypeLikelihoods();
// 		int k = 0;
// 		for (int i = 0; i < 4; i++)
// 		{ 
// 			for (int j = i; j < 4; j++)
// 			{
// 				G.likelihoods[k] = allele_likelihoods[i] * allele_likelihoods[j];
// 				k++;
// 			}
// 		}	
// 		return G;
// 	}

	char PickAlt(char ref, double[] allele_likelihoods)
	{
		Integer[] perm = Utils.SortPermutation(allele_likelihoods);
		if (perm[3] != BaseUtils.simpleBaseToBaseIndex(ref)) { return BaseUtils.baseIndexToSimpleBase(perm[3]); }
		else { return BaseUtils.baseIndexToSimpleBase(perm[2]); }
	}

	double Compute_discovery_lod(char ref, ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{
		double pBest = 0;
		double pRef  = 0;
		for (int i = 0; i < genotype_likelihoods.length; i++)
		{
			pBest += genotype_likelihoods[i].BestPosterior();
			pRef  += genotype_likelihoods[i].RefPosterior(ref);
		}
		return pBest - pRef;
	}

	// this one is a bit of a lazy hack.
	double Compute_alt_freq(char ref, double[] allele_likelihoods)
	{
		return allele_likelihoods[BaseUtils.simpleBaseToBaseIndex(PickAlt(ref, allele_likelihoods))];
	}

	int Compute_n_ref(char ref, ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{
		int n = 0;
		for (int i = 0; i < genotype_likelihoods.length; i++)
		{
			if (genotype_likelihoods[i].coverage == 0) { continue; }
			String g = genotype_likelihoods[i].BestGenotype();
			if ((g.charAt(0) == ref) && (g.charAt(1) == ref)) { n += 1; }
		}
		return n;
	}

	int Compute_n_het(char ref, ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{
		int n = 0;
		for (int i = 0; i < genotype_likelihoods.length; i++)
		{
			if (genotype_likelihoods[i].coverage == 0) { continue; }
			String g = genotype_likelihoods[i].BestGenotype();
			if ((g.charAt(0) == ref) && (g.charAt(1) != ref)) { n += 1; }
			if ((g.charAt(0) != ref) && (g.charAt(1) == ref)) { n += 1; }
		}
		return n;
	}

	int Compute_n_hom(char ref, ClassicGenotypeLikelihoods[] genotype_likelihoods)
	{
		int n = 0;
		for (int i = 0; i < genotype_likelihoods.length; i++)
		{
			if (genotype_likelihoods[i].coverage == 0) { continue; }
			String g = genotype_likelihoods[i].BestGenotype();
			if ((g.charAt(0) != ref) && (g.charAt(1) != ref)) { n += 1; }
		}
		return n;
	}

	// This should actually return a GLF Record
	MultiSampleCallResult MultiSampleCall(RefMetaDataTracker tracker, char ref, AlignmentContext context, List<String> sample_names) 
	{
		String in_dbsnp;
        if (tracker.getReferenceMetaData(rodDbSNP.STANDARD_DBSNP_TRACK_NAME).size() > 0) { in_dbsnp = "known"; } else { in_dbsnp = "novel"; }

		AlignmentContext[] contexts = filterAlignmentContext(context, sample_names, 0);
		glCache.clear(); // reset the cache
		cache_size = 0;
		
		double lod = LOD(contexts);		

		int n_ref = Compute_n_ref(ref, em_result.genotype_likelihoods);
		int n_het = Compute_n_het(ref, em_result.genotype_likelihoods);
		int n_hom = Compute_n_hom(ref, em_result.genotype_likelihoods);

		//double strand_score = lod > MIN_LOD_FOR_STRAND ? StrandScore(context) : 0.0;
		double strand_score;
		if (n_het+n_hom > 0) { strand_score = StrandScore(context); }
		else { strand_score = 0; }

		//EM_Result em_result = EM(contexts);
		//ClassicGenotypeLikelihoods population_genotype_likelihoods = HardyWeinberg(em_result.allele_likelihoods);	

		//double pD = Compute_pD(em_result.genotype_likelihoods);
		//double pNull = Compute_pNull(contexts);

		//double discovery_lod = Compute_discovery_lod(ref, em_result.genotype_likelihoods);
		double alt_freq      = Compute_alt_freq(ref, em_result.allele_likelihoods);

		char alt = 'N';
	   	//if (lod > 0.0) { alt = PickAlt(ref, em_result.allele_likelihoods); }
	   	if ((n_het > 0) || (n_hom > 0)) { alt = PickAlt(ref, em_result.allele_likelihoods); }

        if ( INCLUDE_GENOTYPES ) {
        for (int i = 0; i < em_result.genotype_likelihoods.length; i++)
        {
                individual_output_file.printf("%s %c %s ", context.getLocation(), ref, sample_names.get(i));
                individual_output_file.printf("%s %f %f %s ", em_result.genotype_likelihoods[i].BestGenotype(),
                        em_result.genotype_likelihoods[i].LodVsNextBest(),
                        em_result.genotype_likelihoods[i].LodVsRef(ref),
                        in_dbsnp);

            //individual_output.printf("%s ", new ReadBackedPileup(ref, contexts[i]).getBaseCountsString());
            assert(em_result.genotype_likelihoods[i] != null);
            em_result.genotype_likelihoods[i].sort();
            assert(em_result.genotype_likelihoods[i].sorted_likelihoods != null);

            if ( INCLUDE_GENOTYPES ) {
                for (int j = 0; j < em_result.genotype_likelihoods[i].sorted_likelihoods.length; j++)
                {
                    individual_output_file.printf("%f ", em_result.genotype_likelihoods[i].likelihoods[j]);
                }
                individual_output_file.printf("\n");
            }
        }
        }

		return new MultiSampleCallResult(context.getLocation(), ref, alt, em_result, lod, strand_score, pD, pNull, in_dbsnp, n_ref, n_het, n_hom, em_result.EM_N, alt_freq);
	}

	// END Calling Functions
	/////////

	/////////
	// Utility Functions
	
	/// Filter a locus context by forward and backward
	private AlignmentContext filterAlignmentContext(AlignmentContext context, String strand)
	{
		ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
		ArrayList<Integer> offsets = new ArrayList<Integer>();

		for (int i = 0; i < context.getReads().size(); i++)
		{
			SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);

			// Filter for strandedness
			if ((!strand.contains("+")) && (read.getReadNegativeStrandFlag() == false)) { continue; }
			if ((!strand.contains("-")) && (read.getReadNegativeStrandFlag() == true))  { continue; }
			reads.add(read);
			offsets.add(offset);
		}
		return new AlignmentContext(context.getLocation(), reads, offsets);
	}

	// Filter a locus context by sample ID
    protected AlignmentContext[] filterAlignmentContext(AlignmentContext context, List<String> sample_names, int downsample)
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
			if (SAMPLE_NAME_REGEX != null) { sample = sample.replaceAll(SAMPLE_NAME_REGEX, "$1"); }
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

	private AlignmentContext filter_each_read(AlignmentContext L)
	{
		ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
		ArrayList<Integer> offsets = new ArrayList<Integer>();

        for (int i = 0; i < L.getReads().size(); i++)
        {
            SAMRecord read = L.getReads().get(i);
            Integer offset = L.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));

            assert(this.header != null);
            //assert(this.header.getReadGroup(RG) != null);
            if (this.header.getReadGroup(RG) == null) { continue; }

			// skip bogus data
			if (read.getMappingQuality() == 0) { continue; }

            String sample = this.header.getReadGroup(RG).getSample();
			//if (SAMPLE_NAME_REGEX != null) { sample = sample.replaceAll(SAMPLE_NAME_REGEX, "$1"); }
		
            reads.add(read); 
            offsets.add(offset); 
        }

		AlignmentContext ans = new AlignmentContext(L.getLocation(), reads, offsets);

        return ans;
    }

	// END Utility functions
	/////////



    
}
