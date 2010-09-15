package org.broadinstitute.sting.gatk.datasources.sample;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Aug 12, 2010
 * Time: 1:30:44 PM
 */
public class SampleFileParser {

    private SampleAlias[] sampleAliases;

    private String[] allowedProperties;

    private String[] allowedRelationships;

    private PropertyDefinition[] propertyDefinitions;

    private SampleParser[] samples;

    public PropertyDefinition[] getPropertyDefinitions() {
        return propertyDefinitions;
    }

    public void setPropertyDefinitions(PropertyDefinition[] propertyDefinitions) {
        this.propertyDefinitions = propertyDefinitions;
    }

    public SampleFileParser() {

    }

    public String[] getAllowedProperties() {
        return allowedProperties;
    }

    public void setAllowedProperties(String[] allowedProperties) {
        this.allowedProperties = allowedProperties;
    }

    public SampleParser[] getSamples() {
        return samples;
    }

    public void setSamples(SampleParser[] samples) {
        this.samples = samples;
    }

    public String[] getAllowedRelationships() {
        return allowedRelationships;
    }

    public void setAllowedRelationships(String[] allowedRelationships) {
        this.allowedRelationships = allowedRelationships;
    }

    public SampleAlias[] getSampleAliases() {
        return sampleAliases;
    }

    public void setSampleAliases(SampleAlias[] sampleAliases) {
        this.sampleAliases = sampleAliases;
    }

}