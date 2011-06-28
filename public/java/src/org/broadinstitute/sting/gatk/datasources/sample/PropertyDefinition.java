package org.broadinstitute.sting.gatk.datasources.sample;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Aug 12, 2010
 * Time: 2:09:16 PM
 */
public class PropertyDefinition {

    String property;

    String[] values;

    public String getProperty() {
        return property;
    }

    public void setProperty(String property) {
        this.property = property;
    }

    public String[] getValues() {
        return values;
    }

    public void setValues(String[] values) {
        this.values = values;
    }
}