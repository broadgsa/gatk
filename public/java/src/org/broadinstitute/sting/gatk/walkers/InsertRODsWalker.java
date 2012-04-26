package org.broadinstitute.sting.gatk.walkers;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 3/30/12
 * Time: 4:47 PM
 * To change this template use File | Settings | File Templates.
 */

import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.db.MongoDB;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.List;

/**
 * Inserts all of the RODs in the input data set. Data is inserted using VariantContext.toMongoDB().
 */
public class InsertRODsWalker extends RodWalker<Integer, Integer> {
    @Input(fullName="input", shortName = "input", doc="The input ROD which should be inserted into the DB.", required=true)
    public RodBinding<Feature> input;

    @Output
    PrintStream out;


    private String RODFileName;

    @Override
    public void initialize() {
        DBCollection mongoAttributes = MongoDB.getAttributesCollection();
        DBCollection mongoSamples = MongoDB.getSamplesCollection();

        RODFileName = input.getSource();
        int lastSep = RODFileName.lastIndexOf(File.separator);
        RODFileName = RODFileName.substring(lastSep + 1);

        // set up indices

        mongoAttributes.ensureIndex("location");
        mongoAttributes.ensureIndex("sourceROD");
        mongoAttributes.ensureIndex("contig");
        mongoAttributes.ensureIndex("start");
        mongoAttributes.ensureIndex("stop");

        mongoSamples.ensureIndex("location");
        mongoSamples.ensureIndex("sample");
        mongoSamples.ensureIndex("sourceROD");
        mongoSamples.ensureIndex("contig");
        mongoSamples.ensureIndex("start");
        mongoSamples.ensureIndex("stop");

        // set up primary keys
        mongoAttributes.ensureIndex(new BasicDBObject("location", 1).append("sourceROD", 1).append("alleles", 1), new BasicDBObject("unique", 1));
        mongoSamples.ensureIndex(new BasicDBObject("location", 1).append("sourceROD", 1).append("alleles", 1).append("sample", 1), new BasicDBObject("unique", 1));
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        DBCollection mongoAttributes = MongoDB.getAttributesCollection();
        DBCollection mongoSamples = MongoDB.getSamplesCollection();

        for ( Feature feature : tracker.getValues(Feature.class, context.getLocation()) ) {
            if ( feature instanceof VariantContext ) {
                VariantContext vc = (VariantContext) feature;

                Pair<BasicDBObject,List<BasicDBObject>> mongoCollections = vc.toMongoDB(RODFileName);
                mongoAttributes.insert(mongoCollections.first);
                for (BasicDBObject sampleForMongo : mongoCollections.second) {
                    mongoSamples.insert(sampleForMongo);
                }
            }
        }

        return 1;
    }

    /**
     * Increment the number of rods processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of rods processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        MongoDB.close();
    }
}