package org.broadinstitute.sting.utils.db;

import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.Mongo;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.net.UnknownHostException;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 4/26/12
 * Time: 3:01 PM
 * Handles Mongo DB connections
 */
final public class MongoDB {
    private final static String MONGO_HOST = "couchdb.broadinstitute.org";
    private final static Integer MONGO_PORT = 43054;
    private final static String MONGO_DB_NAME = "bjorn";
    private final static String MONGO_ATTRIBUTES_COLLECTION = "attributes";
    private final static String MONGO_SAMPLES_COLLECTION = "samples";

    protected Mongo mongo;
    protected DBCollection mongoAttributes;
    protected DBCollection mongoSamples;

    final private static MongoDB INSTANCE = new MongoDB();

    public static DBCollection getAttributesCollection() {
        return INSTANCE.mongoAttributes;
    }

    public static DBCollection getSamplesCollection() {
        return INSTANCE.mongoSamples;
    }

    public static void close() {
        INSTANCE.mongo.close();
    }

    private MongoDB() {
        try {
            mongo = new Mongo(MONGO_HOST, MONGO_PORT);
            DB mongoDb = mongo.getDB(MONGO_DB_NAME);
            mongoAttributes = mongoDb.getCollection(MONGO_ATTRIBUTES_COLLECTION);
            mongoSamples = mongoDb.getCollection(MONGO_SAMPLES_COLLECTION);
        } catch (UnknownHostException e) {
            throw new StingException(e.getMessage(), e);
        }
    }
}
