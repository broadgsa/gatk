package org.broadinstitute.sting.gatk.datasources.sample;

import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Aug 13, 2010
 * Time: 2:09:43 PM
 */
public class SampleParser {

    private String id;

    private HashMap<String, Object> properties;

    private HashMap<String, String> relationships;

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public HashMap<String, Object> getProperties() {
        return properties;
    }

    public void setProperties(HashMap<String, Object> properties) {
        this.properties = properties;
    }

    public HashMap<String, String> getRelationships() {
        return relationships;
    }

    public void setRelationships(HashMap<String, String> relationships) {
        this.relationships = relationships;
    }

}
