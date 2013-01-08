package org.broadinstitute.sting.gatk.walkers.coverage;


import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;


import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * print intervals file with all the variant sites that have "most" ( >= 90% by default) of the samples with "good" (>= 10 by default)coverage ("most" and "good" can be set in the command line).
 *
 * <p>
 * CoveredByNSamplesSites is a GATK tool for filter out sites based on their coverage.
 * The sites that pass the filter are printed out to an intervals file.
 *
 * <h2>Input</h2>
 * <p>
 * A variant file and optionally min coverage and sample percentage values.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An intervals file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CoveredByNSamplesSites \
 *   -V input.vcf \
 *   -out output.intervals \
 *   -minCov 15
 * </pre>
 *
 */

@By(DataSource.REFERENCE_ORDERED_DATA)
public class CoveredByNSamplesSites extends RodWalker<Pair<Integer,Integer>, Pair<Integer,Integer>> implements TreeReducible<Pair<Integer,Integer>> {

    @Output(fullName = "OutputIntervals", shortName = "out", doc = "Name of file for output intervals", required = true)
    File intervalsFile;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(fullName = "minCoverage", shortName = "minCov",doc = "only samples that have covarage bigger then minCoverage will be counted",required = false)
    int minCoverage = 10;

    @Argument(fullName = "precentageOfSamples", shortName = "percentage", doc = "only sites where at list percentageOfSamples of the samples have good coverage, will be emited", required = false)
    double percentageOfSamples = 0.9;

    private FileOutputStream outputStream;

    @Override
    public void initialize(){
        if (! intervalsFile.getName().endsWith(".intervals")){
            throw new UserException(String.format("Output interval file %s should be <name>.intervals", intervalsFile));
        }
        try{
            outputStream = new FileOutputStream(intervalsFile);
            if (!intervalsFile.exists()) {
                intervalsFile.createNewFile();
            }
        }
        catch (IOException e){
            System.err.println(String.format("Problems with creating outputStream from %s",intervalsFile));
            e.printStackTrace();
        }
    }

    @Override
    public Pair<Integer,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return new Pair<Integer, Integer>(0,0);

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        if ( VCs.size() == 0 )
            return new Pair<Integer, Integer>(0,0);
        if(VCs.size() != 1)
            throw new RuntimeException("there are more then one vc: "+VCs.size());

        boolean emitSite = false;
        List<GenomeLoc> outputIntervals = new ArrayList<GenomeLoc>();
        for(VariantContext vc : VCs){
            int coveredSamples = 0;
            final GenotypesContext genotypes = vc.getGenotypes();
            final int numOfGenotypes = genotypes.size();
            for(Genotype g : genotypes){
                if(g.getDP() >= minCoverage)
                    coveredSamples++;
            }
            if((double)coveredSamples/numOfGenotypes > percentageOfSamples){
                emitSite = true;
                outputIntervals.add(ref.getLocus());
                //System.out.println(ref.getLocus());
                try{
                    String toPrint = ref.getLocus().toString() + "\n";
                    outputStream.write(toPrint.getBytes());
                }catch (IOException e){
                    e.printStackTrace();
                }
            }

        }
        if (emitSite)
            return new Pair<Integer, Integer>(1,1);
        else
            return new Pair<Integer, Integer>(1,0);
    }

    @Override
    public Pair<Integer,Integer> reduceInit() { return new Pair<Integer, Integer>(0,0); }

    @Override
    public Pair<Integer,Integer> reduce(Pair<Integer,Integer> value, Pair<Integer,Integer> sum) { return new Pair<Integer, Integer>(value.getFirst() + sum.getFirst(),value.getSecond()+sum.getSecond()); }

    @Override
    public Pair<Integer,Integer> treeReduce(Pair<Integer,Integer> lhs, Pair<Integer,Integer> rhs) {
        return new Pair<Integer, Integer>(lhs.getFirst() + rhs.getFirst(),lhs.getSecond() + rhs.getSecond());
    }

    /**
     * Tell the user the number of sites processed and how many passed. Close out the new intervals file.
     *
     * @param result  pair of the number of sites seen and number of sites passed the filter.
     */
    public void onTraversalDone(Pair<Integer,Integer> result) {
        logger.info("Processed " + result.getFirst() + " variant sites and found "+ result.getSecond()+" sites that have "+(percentageOfSamples*100)+"% of the samples with at list "+minCoverage+" coverage.\n");
        logger.info("All these sites were printed to the intervals file: "+intervalsFile.getAbsolutePath());
        try{
            outputStream.close();
        }
        catch (IOException e){
            System.err.println("Couldn't close output steam properly");
            e.printStackTrace();
        }
    }



}
