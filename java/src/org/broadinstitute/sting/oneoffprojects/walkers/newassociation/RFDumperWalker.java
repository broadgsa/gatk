package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * A (currently very lame) utility to dump read feature information to a file for validity checking of read feature extractors.
 * Soon to be extended for aggregator and window validation; as well as postspective signal analysis.
 */
public class RFDumperWalker extends ReadWalker<String[],ReadFeatureWindow> {
    @ArgumentCollection
    private RFAArgumentCollection collection = new RFAArgumentCollection();

    List<ReadFeatureAggregator> aggregators;
    GenomeLoc loc;
    boolean paired;
    String sample;
    String name;

    String[] data;

    @Output
    PrintStream out;

    public void initialize() {

        Set<Class<? extends ReadFeatureAggregator>> aggregatorSet = getFeatureAggregators(collection.inputFeatures);
        Set<ReadFeatureAggregator> rfHolder1 = new HashSet<ReadFeatureAggregator>(aggregatorSet.size());
        try {
            for ( Class<? extends ReadFeatureAggregator> featureClass : aggregatorSet ) {
                ReadFeatureAggregator readFeature = featureClass.getConstructor(collection.getClass()).newInstance(collection);
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
        for ( ReadFeatureAggregator ag : aggregators ) {
            logger.info(ag.getClass().getSimpleName());
        }

        data = new String[aggregators.size()];
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

    public ReadFeatureWindow reduceInit() { return null; }

    public String[] map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        // TODO: THIS WILL BREAK IF FEATURE REQUIRES PAIRED READ
        if ( ref == null ) { return null; } // unmapped reads have null ref contexts
        loc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),read.getAlignmentStart());
        sample = read.getReadGroup().getSample();
        name = read.getReadName();
        paired = read.getReadPairedFlag() && ! read.getMateUnmappedFlag();
        int idx = 0;
        for ( ReadFeatureAggregator aggregator : aggregators) {
            data[idx++] = String.format("%s: %s",aggregator.getClass().getSimpleName(),aggregator.parseStr(read));
        }
        return data;
    }

    public ReadFeatureWindow reduce(String[] map, ReadFeatureWindow prevReduce) {
        if ( map == null ) { return null; }
        StringBuffer fStrBuilder = new StringBuffer();
        for ( String f : map ) {
            fStrBuilder.append("\t");
            fStrBuilder.append(f);
        }

        out.printf("%s\t%s\t%s%s%n",loc,name,sample,fStrBuilder);

        return null;
    }
}