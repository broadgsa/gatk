package org.broadinstitute.sting.gatk.walkers;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 3/30/12
 * Time: 4:47 PM
 * To change this template use File | Settings | File Templates.
 */

import com.mongodb.*;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;

/**
 * Inserts all of the RODs in the input data set. Data is inserted using VariantContext.toMongoDB().
 */
public class InsertRODsWalker extends RodWalker<Integer, Integer> {
    @Input(fullName="input", shortName = "input", doc="The input ROD which should be inserted into the DB.", required=true)
    public RodBinding<Feature> input;

    @Output
    PrintStream out;

    private final static String MONGO_HOST = "gsa4.broadinstitute.org";
    private final static Integer MONGO_PORT = 43054;
    private final static String MONGO_DB_NAME = "bjorn";
    private final static String MONGO_VC_COLLECTION = "vcs";

    protected Mongo mongo;
    protected DBCollection mongoCollection;

    private String RODFileName;

    @Override
    public void initialize()
    {
        try {
            mongo = new Mongo(MONGO_HOST, MONGO_PORT);
            DB mongoDb = mongo.getDB(MONGO_DB_NAME);
            mongoCollection = mongoDb.getCollection(MONGO_VC_COLLECTION);

            RODFileName = input.getSource();
            int lastSep = RODFileName.lastIndexOf(File.separator);
            RODFileName = RODFileName.substring(lastSep + 1);

            // set up indices
            mongoCollection.ensureIndex("location");
            mongoCollection.ensureIndex("sample");
            mongoCollection.ensureIndex("sourceROD");
            mongoCollection.ensureIndex("contig");
            mongoCollection.ensureIndex("start");
            mongoCollection.ensureIndex("stop");

            // set up primary key
            mongoCollection.ensureIndex(new BasicDBObject("location", 1).append("sample", 1).append("sourceROD", 1), new BasicDBObject("unique", 1));

        }
        catch (MongoException e) {}
        catch (java.net.UnknownHostException e) {}
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

        for ( Feature feature : tracker.getValues(Feature.class, context.getLocation()) ) {
            if ( feature instanceof VariantContext ) {
                VariantContext vc = (VariantContext) feature;
                for (BasicDBObject vcForMongo : vc.toMongoDB(RODFileName)) {
                    mongoCollection.insert(vcForMongo);
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
        mongo.close();
    }
}