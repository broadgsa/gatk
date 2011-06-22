package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Read feature association walker -- associates read features between dichotomized, or multi-group cohorts
 * todo -- need a heuristic to stop doing tests where there is certainly no signal
 * todo -- for most features there's a nuisance variable which is the proportion of *paired* reads, perhaps a pair-only setting for read features
 */

@ReadFilters({MaxInsertSizeFilter.class,MappingQualityReadFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckReadFilter.class,NotPrimaryAlignmentReadFilter.class,UnmappedReadFilter.class})
public class RFAWalker extends ReadWalker<SAMRecord,RFWindow> {
    // todo -- this needs to be an argument collection that can get passed around to initialize read features etc
    @ArgumentCollection
    private RFAArgumentCollection collection = new RFAArgumentCollection();

    @Output
    public PrintStream out;

    Map<String,Boolean> caseStatus;

    protected List<ReadFeatureAggregator> aggregators; // no re-instantiation, use a list to ensure ordering

    protected Iterator<GenomeLoc> locusIterator;
    protected GenomeLoc iteratorLoc;
    protected GenomeLoc loc;
    protected String sample;
    private List<String> EMPTY_LIST = new ArrayList<String>(0);

    public void initialize() {
        if ( collection.windowSize % collection.windowJump != 0 ) {
            throw new UserException("Window size is not divisible by window jump.");
        }

        if ( collection.caseFile == null || collection.controlFile == null ) {
            throw new UserException("You must provide both a case file (-case) and a control file (-control) each listing those samples belonging to the cohort");
        }

        caseStatus = new HashMap<String,Boolean>(getToolkit().getSAMFileSamples().size());
        try {
            for ( String sample : new XReadLines(collection.caseFile) ) {
                caseStatus.put(sample,true);
            }
            for ( String sample : new XReadLines(collection.controlFile)) {
                caseStatus.put(sample,false);
            }

            for ( Sample sample : getToolkit().getSAMFileSamples() ) {
                if ( ! caseStatus.containsKey(sample.getId())) {
                    throw new UserException("No case/control status for sample "+sample.getId());
                }
            }

        } catch ( FileNotFoundException e ) {
            throw new UserException("Unable to open a case/control file",e);
        }

        Set<Class<? extends ReadFeatureAggregator>> aggregatorSet = getFeatureAggregators(collection.inputFeatures);
        Set<ReadFeatureAggregator> rfHolder1 = new HashSet<ReadFeatureAggregator>(aggregatorSet.size());
        try {
            for ( Class<? extends ReadFeatureAggregator> featureClass : aggregatorSet ) {
                ReadFeatureAggregator readFeature = featureClass.getConstructor(RFAArgumentCollection.class).newInstance(collection);
                rfHolder1.add(readFeature);
            }
        } catch ( Exception e ) {
            throw new StingException("A read feature instantiation error occurred during initialization",e);
        }

        ReadFeatureAggregator[] rfHolder2 = new ReadFeatureAggregator[rfHolder1.size()];
        int idx = 0;
        for ( ReadFeatureAggregator f : rfHolder1 ) {
            rfHolder2[idx++] = f;
        }
        Arrays.sort(rfHolder2, new Comparator<ReadFeatureAggregator>() {
            @Override
            public int compare(ReadFeatureAggregator a, ReadFeatureAggregator b) {
                return a.getClass().getSimpleName().compareTo(b.getClass().getSimpleName());
            }
        });
        aggregators = Arrays.asList(rfHolder2);

        writeHeader();

        locusIterator = getToolkit().getIntervals().iterator();
        iteratorLoc = locusIterator.hasNext() ? locusIterator.next() : null;

    }

    public RFWindow reduceInit() {
        Set<String> samples = new HashSet<String>(getToolkit().getSamples().size());
        for ( Sample s : getToolkit().getSamples() ) {
            samples.add(s.getId());
        }
        return new RFWindow(aggregators,collection,caseStatus,getToolkit().getGenomeLocParser());
    }

    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if ( ref == null ) { return null; } // unmapped reads have null ref contexts
        //loc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),read.getAlignmentStart());
        GenomeLoc newLoc = ref.getLocus().getStartLocation(); // can be problematic if read aligns prior to start of contig -- should never happen
        if ( newLoc.isPast(iteratorLoc.getStartLocation()) ) {
            loc = newLoc;
        } else {
            loc = iteratorLoc.getStartLocation();
        }
        if ( read == null ) { return null; }
        sample = read.getReadGroup().getSample();
        return read;
    }

    public RFWindow reduce(SAMRecord read, RFWindow prevReduce) {
        if ( iteratorLoc != null && iteratorLoc.isBefore(loc) ) {// test if read is past end of the user interval
            //logger.info(String.format("iteratorLoc: %s    loc: %s",iteratorLoc.toString(),loc.toString()));
            onIntervalDone(prevReduce);
            iteratorLoc = locusIterator.hasNext() ? locusIterator.next() : null;
            if ( loc.startsBefore(iteratorLoc) ) {
                loc = iteratorLoc.getStartLocation();
            }
            reduce(read,prevReduce);
        } else if ( read != null ) {
            // todo -- what happens if first read of an interval is not before or at the start of the interval?\
            List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> completed = prevReduce.inc(read, loc, sample,iteratorLoc);
            // todo -- run tests here; for now just log that a window/multiple windows are complete
            if ( completed.size() > 0 ) {
                // System.out.printf("At %s we have seen %d completed windows%n",loc,completed.size())
                // bed format
                int locShift = 0;
                for ( Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>> samWindow : completed ) {
                    GenomeLoc window = samWindow.first;
                    runWindowTests(samWindow.second, window);
                    locShift += collection.windowJump;
                }
            }
        }
        return prevReduce;
    }


    public RFWindow onIntervalDone(RFWindow rWindow) {
        //logger.info("In onIntervalDone at genome loc "+iteratorLoc.toString()+" with read loc "+loc.toString());
        List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> completed = rWindow.flush(iteratorLoc);
        int locShift = 0;
        for ( Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>> samWindow : completed ) {
            GenomeLoc window = samWindow.first;
            runWindowTests(samWindow.second, window);
            locShift += collection.windowJump;
        }

        return rWindow;
    }

    public Set<Class<? extends ReadFeatureAggregator>> getFeatureAggregators(List<String> requestedFeatures) {
        HashSet<Class<? extends ReadFeatureAggregator>> newFeatureSet = new HashSet<Class<? extends ReadFeatureAggregator>>();
        List<Class<? extends ReadFeatureAggregator>> availableFeatures = new PluginManager<ReadFeatureAggregator>(ReadFeatureAggregator.class).getPlugins();

        if ( collection.inputFeatures == null ) {
            newFeatureSet.addAll(availableFeatures);
            return newFeatureSet;
        }


        Map<String,Class<? extends ReadFeatureAggregator>> classNameToClass = new HashMap<String,Class<? extends ReadFeatureAggregator>>(collection.inputFeatures.size());
        for ( Class<? extends ReadFeatureAggregator> clazz : availableFeatures ) {
            classNameToClass.put(clazz.getSimpleName(),clazz);
        }

        for ( String s : requestedFeatures) {
            if ( classNameToClass.containsKey(s) ) {
                newFeatureSet.add(classNameToClass.get(s));
            } else {
                throw new UserException("The name "+s+" does not correspond to an available read feature class.");
            }
        }

        return newFeatureSet;
    }

    public void runWindowTests(Map<String,List<ReadFeatureAggregator>> window, GenomeLoc loc) {
        // two main tests: fixed-significance shift, and confidence-interval sum
        //System.out.printf("Running tests...%n");
        // todo -- really the aggregators should be iterated over directly (rather than indirectly through the index)
        out.printf("%s\t%d\t%d",loc.getContig(),loc.getStart(),loc.getStop());
        for ( int agIdx = 0; agIdx < aggregators.size(); agIdx ++ ) {
            double fixedDelta = fixedSignificance(window.get("case").get(agIdx),window.get("control").get(agIdx));
            //double weightedDelta = confidenceIntervalSum(window,agIdx);
            //out.printf("\t%.2e\t%.2e",Double.isNaN(fixedDelta) ? 0.0 : fixedDelta,weightedDelta);
            Pair<List<String>,List<String>> caseControlAffected = getAffectedSamples(window,agIdx,fixedDelta);
            List<String> cases = caseControlAffected.getFirst();
            List<String> controls = caseControlAffected.getSecond();
            out.printf("\t%.2e\t%d:%d\t%.2e\t%s;%s",fixedDelta,cases.size(),controls.size(), MathUtils.binomialProbability(cases.size(),cases.size()+controls.size(),0.5),Utils.join(",",cases),Utils.join(",",controls));

        }
        out.printf("%n");
    }

    public Pair<List<String>,List<String>> getAffectedSamples(Map<String,List<ReadFeatureAggregator>> aggregators, int idx, double fixedDelta) {
        if ( fixedDelta == 0.0 ) { return new Pair<List<String>,List<String>>(EMPTY_LIST,EMPTY_LIST); } // todo -- too hacky

        Pair<List<String>,List<String>> ccSampleList = new Pair<List<String>,List<String>>(new ArrayList<String>(), new ArrayList<String>());
        for ( Map.Entry<String,List<ReadFeatureAggregator>> entry : aggregators.entrySet() ) {
            if ( entry.getKey().equals("case") || entry.getKey().equals("control")) { continue; }
            ReadFeatureAggregator aggregator = entry.getValue().get(idx);
            // is this sending a truly significant signal
            double zs = (aggregator.getMean() - collection.EPSILON) * Math.sqrt(aggregator.getnReads())/Math.sqrt(aggregator.getUnbiasedVar());
            if ( zs > collection.sampleZThresh ) {
                if ( caseStatus.get(entry.getKey()) ) {
                    ccSampleList.first.add(entry.getKey());
                } else {
                    ccSampleList.second.add(entry.getKey());
                }
            }
        }

        return ccSampleList;
    }

    public double fixedSignificance(ReadFeatureAggregator caseAg, ReadFeatureAggregator controlAg) {
        if ( caseAg.getnReads() == 0 || controlAg.getnReads() == 0 ) {
            return 0.0;
        }
        double stat_num = caseAg.getMean() - controlAg.getMean();
        double stat_denom = Math.sqrt(caseAg.getUnbiasedVar()/caseAg.getnReads() + controlAg.getUnbiasedVar()/controlAg.getnReads());
        double stat = stat_num/stat_denom;
        //System.out.printf("Mean_dif: %.2e Var: %.2e Stat: %.2f Z: %.2f SS-ZZ: %.2f%n",stat_num,stat_denom,stat, collection.fixedZ, stat*stat-collection.fixedZ*collection.fixedZ);
        if ( stat*stat < collection.fixedZ*collection.fixedZ ) {
            return 0.0;
        } else {
            //System.out.printf("Calculating delta: %.2f%n",(stat < 0) ? stat_denom*(-1*collection.fixedZ-stat) : stat_denom*(stat-collection.fixedZ));
            //return (stat > 0) ? stat_denom*(stat-collection.fixedZ) : stat_denom*(stat+collection.fixedZ);
	    return stat_num;
        }
    }

    public double confidenceIntervalSum(Map<String,List<ReadFeatureAggregator>> window, int offset) {
        // this comment serves as an explicit normality assumption (e.g. that the DF will be large)
        double caseWeightMean = 0.0;
        double caseWeightVar = 0.0;
        double caseSumWeight = 0.0;
        double controlWeightMean = 0.0;
        double controlWeightVar = 0.0;
        double controlSumWeight = 0.0;
        for ( Map.Entry<String,List<ReadFeatureAggregator>> sampleEntry : window.entrySet() ) {
            if ( ! sampleEntry.getKey().equals("case") && ! sampleEntry.getKey().equals("control") ) {
                // first check if the sample is shifted from zero (CI does not include zero)
                // todo -- fixme. This will always be true for insert sizes, clipped reads, mapping quality...should have an avg value

                ReadFeatureAggregator aggregator = sampleEntry.getValue().get(offset);
                if ( aggregator.getnReads() == 0 ) {
                    continue;
                }

                boolean shifted;
                int shiftThresh = 0;/*
                // todo -- this is fucking awful
                if ( aggregator instanceof InsertSize ) {
                    shiftThresh = 150;
                }
                if ( aggregator instanceof ClippedBases ) {
                    shiftThresh = 6;
                }*/
                if ( aggregator.getMean() < shiftThresh ) {
                    // use fixedZ/4 to be a little bit less strenuous -- todo -- make an input?
                    shifted = aggregator.getMean() + collection.fixedZ/4*Math.sqrt(aggregator.getUnbiasedVar())/aggregator.getnReads() < shiftThresh;
                } else {
                    shifted = aggregator.getMean() - collection.fixedZ/4*Math.sqrt(aggregator.getUnbiasedVar())/aggregator.getnReads() > shiftThresh;;
                }
                if ( shifted ) {
                    double twoS2 = 2*aggregator.getUnbiasedVar()*aggregator.getUnbiasedVar();
                    if ( caseStatus.get(sampleEntry.getKey())) {
                        caseWeightMean += (aggregator.getnReads()-1.0)*aggregator.getMean()/(twoS2);
                        caseWeightVar += aggregator.getUnbiasedVar()*Math.pow((aggregator.getnReads()-1.0)/(twoS2),2);
                        caseSumWeight += (aggregator.getnReads()-1.0)/(twoS2);
                    } else {
                        controlWeightMean += (aggregator.getnReads()-1.0)*aggregator.getMean()/(twoS2);
                        controlWeightVar += aggregator.getUnbiasedVar()*Math.pow((aggregator.getnReads()-1.0)/(twoS2),2);
                        controlSumWeight += (aggregator.getnReads()-1.0)/(twoS2);
                    }
                }
            }
        }

        double caseGaussianMean = caseWeightMean/caseSumWeight;
        double controlGaussianMean = controlWeightMean/controlSumWeight;
        double caseGaussianVar = caseWeightVar/(caseSumWeight*caseSumWeight);
        double controlGaussianVar = controlWeightVar/(controlSumWeight*controlSumWeight);
        // todo -- is the z-factor an appropriate statistic?
        //return 1.0 - 3*(caseGaussianVar+controlGaussianVar)/Math.abs(caseGaussianMean-controlGaussianMean);
        if ( caseGaussianMean > controlGaussianMean ) {
            // want to examine the case lower fixedZ*stdev vs the control upper fixedZ*stev
            return (caseGaussianMean-collection.fixedZ/4*Math.sqrt(caseGaussianVar)) - (controlGaussianMean + collection.fixedZ/4*Math.sqrt(controlGaussianVar));
        } else {
            // want to examine the case upper fixedZ*stev vs the control lower fixedZ*stev
            return (controlGaussianMean-collection.fixedZ/4*Math.sqrt(controlGaussianVar)) - ( caseGaussianMean + collection.fixedZ/4*Math.sqrt(caseGaussianVar));
        }
    }

    public void writeHeader() {
        // "%.2e\t%d:%d\t%s,%s\t%.2e"
        StringBuffer buf = new StringBuffer();
        buf.append("description=chr,start,stop");
        for ( ReadFeatureAggregator f : aggregators ) {
            buf.append(",");
            buf.append(f.getClass().getSimpleName());
            buf.append("-d,");
            buf.append(f.getClass().getSimpleName());
            buf.append("-r,");
            buf.append(f.getClass().getSimpleName());
            buf.append("-s,");
            buf.append(f.getClass().getSimpleName());
            buf.append("-p");
        }
        out.printf("track type=bedTable %s%n",buf);
    }
}
