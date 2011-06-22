package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.ClippedBases;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.InsertSize;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.ReadFeatureAggregator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/12/11
 * Time: 1:48 PM
 * To change this template use File | Settings | File Templates.
 */
@ReadFilters({MaxInsertSizeFilter.class,MappingQualityReadFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckReadFilter.class,NotPrimaryAlignmentReadFilter.class,UnmappedReadFilter.class})
@By(DataSource.REFERENCE)
public class RFExtractorWalker extends ReadWalker<SAMRecord,RFWindow> {

    @ArgumentCollection
    public RFAArgumentCollection rfaArgs = new RFAArgumentCollection();

    @Argument(doc="Set the marker threshold to this.",shortName="mmt",fullName="markerModeThreshold",required=false)
    public double markerModeThreshold = 0.05;

    @Argument(doc="Turn on marker mode (1 if sample is significantly greater than the threshold, 0 otherwise)",shortName="mm",fullName="markerMode",required=false)
    public boolean markerMode = false;

    @Argument(doc="Turn on raw count mode: output will be raw aberrant read counts",shortName="c",fullName="count",required=false)
    public boolean countMode = false;


    @Output
    public PrintStream out;

    protected Iterator<GenomeLoc> locusIterator;
    protected GenomeLoc iteratorLoc;
    protected GenomeLoc loc;
    protected String sample;

    public void initialize() {
        if ( rfaArgs.windowSize % rfaArgs.windowJump != 0 ) {
            throw new UserException("Window size is not divisible by window jump.");
        }

        if ( markerMode && markerModeThreshold < 0.0 ) {
            throw new UserException("Cannot have a negative threshold when using marker mode");
        }

        if ( countMode && markerMode ) {
            throw new UserException("Cannot be both in count mode and marker mode");
        }

        locusIterator = getToolkit().getIntervals().iterator();
        iteratorLoc = locusIterator.next();
    }

    public RFWindow reduceInit() {
        Map <String,Boolean> allCase = new HashMap<String,Boolean>(getToolkit().getSamples().size());
        for ( Sample s : getToolkit().getSAMFileSamples() ) {
            allCase.put(s.getId(),true);
            if ( s.getId() == null || s.getId().equals("null") ) {
                throw new StingException("Sample IDs must not be null... " + s.toString() + " " + Boolean.toString(s.hasSAMFileEntry()));
            }
        }

        Set<Class<? extends ReadFeatureAggregator>> aggregatorSet = getFeatureAggregators(rfaArgs.inputFeatures);
        Set<ReadFeatureAggregator> rfHolder1 = new HashSet<ReadFeatureAggregator>(aggregatorSet.size());
        try {
            for ( Class<? extends ReadFeatureAggregator> featureClass : aggregatorSet ) {
                ReadFeatureAggregator readFeature = featureClass.getConstructor(RFAArgumentCollection.class).newInstance(rfaArgs);
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

        List<ReadFeatureAggregator> aggregators = Arrays.asList(rfHolder2);

        out.printf("HEADERchrm:start-stop");
        for ( String s : allCase.keySet() ) {
            for ( ReadFeatureAggregator rfa : aggregators ) {
                out.printf("\t%s.%s",s,rfa.getClass().getSimpleName());
            }
        }
        out.printf("%n");

        return new RFWindow(aggregators,rfaArgs,allCase,getToolkit().getGenomeLocParser());
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
            // todo -- what happens if first read of an interval is not before or at the start of the interval?

            List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> completed = prevReduce.inc(read, loc, sample,iteratorLoc);
            if ( completed.size() > 0 ) {
                // System.out.printf("At %s we have seen %d completed windows%n",loc,completed.size())
                // bed format
                for ( Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>> samWindow : completed ) {
                    GenomeLoc window = samWindow.first;
                    /*if ( prevPrint == null ) {
                        prevPrint = window;
                    } else if ( window.startsBefore(prevPrint) ) {
                        throw new StingException(String.format("Attempting to print at %s after having printed a record at %s",window.toString(),prevPrint.toString()));
                    } else {
                        prevPrint = window;
                    }*/
                    out.printf("%s",window.toString());
                    for ( Map.Entry<String,List<ReadFeatureAggregator>> samEntry : samWindow.second.entrySet() ) {
                        for ( ReadFeatureAggregator aggregator : samEntry.getValue() ) {
                            if ( ! markerMode && ! countMode ) {
                                out.printf("\t%.5e,%d",aggregator.getMean(),aggregator.getnReads());
                            } else if ( markerMode ) {
                                out.printf("\t%d",hasEvent(aggregator,markerModeThreshold,rfaArgs.fixedZ) ? 1 : 0);
                            } else if ( countMode ) {
                                out.printf("\t%d", MathUtils.fastRound(aggregator.getMean()*aggregator.getnReads()));
                            }
                        }
                    }
                    out.printf("%n");
                }
            }
        } else {
            prevReduce.inc(null,loc,null,iteratorLoc);
        }
        return prevReduce;
    }

    public RFWindow onIntervalDone(RFWindow rWindow) {
        //logger.info("In onIntervalDone at genome loc "+iteratorLoc.toString()+" with read loc "+loc.toString());
        List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> completed = rWindow.flush(iteratorLoc);
        for ( Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>> samWindow : completed ) {
            GenomeLoc window = samWindow.first;
            /*if ( prevPrint == null ) {
                        prevPrint = window;
                    } else if ( window.startsBefore(prevPrint) ) {
                        throw new StingException(String.format("Attempting to print at %s after having printed a record at %s",window.toString(),prevPrint.toString()));
                    } else {
                        prevPrint = window;
                    }*/
            out.printf("%s",window.toString());
            for ( Map.Entry<String,List<ReadFeatureAggregator>> samEntry : samWindow.second.entrySet() ) {
                for ( ReadFeatureAggregator aggregator : samEntry.getValue() ) {
                    if ( ! markerMode && ! countMode ) {
                        out.printf("\t%.5e,%d",aggregator.getMean(),aggregator.getnReads());
                    } else if ( markerMode ) {
                        out.printf("\t%d",hasEvent(aggregator,markerModeThreshold,rfaArgs.fixedZ) ? 1 : 0);
                    } else if ( countMode ) {
                        out.printf("\t%d", MathUtils.fastRound(aggregator.getMean()*aggregator.getnReads()));
                    }
                }
            }
            out.printf("%n");
        }

        return rWindow;
    }

    public Set<Class<? extends ReadFeatureAggregator>> getFeatureAggregators(List<String> requestedFeatures) {
        HashSet<Class<? extends ReadFeatureAggregator>> newFeatureSet = new HashSet<Class<? extends ReadFeatureAggregator>>();
        List<Class<? extends ReadFeatureAggregator>> availableFeatures = new PluginManager<ReadFeatureAggregator>(ReadFeatureAggregator.class).getPlugins();

        if ( rfaArgs.inputFeatures == null ) {
            newFeatureSet.addAll(availableFeatures);
            return newFeatureSet;
        }


        Map<String,Class<? extends ReadFeatureAggregator>> classNameToClass = new HashMap<String,Class<? extends ReadFeatureAggregator>>(rfaArgs.inputFeatures.size());
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

    public static boolean hasEvent(ReadFeatureAggregator aggregator, double lowThresh, double sigLevel) {
        return (aggregator.getMean() - lowThresh)*Math.sqrt(aggregator.getnReads())/aggregator.getVar() > sigLevel;
    }
}
