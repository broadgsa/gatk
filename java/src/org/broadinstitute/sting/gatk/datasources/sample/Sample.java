package org.broadinstitute.sting.gatk.datasources.sample;


import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Jul 26, 2010
 * Time: 3:31:38 PM
 */
public class Sample implements java.io.Serializable {

    private final String id;

    private boolean hasSampleFileEntry = false; // true if this sample has an entry in a sample file

    private boolean hasSAMFileEntry = false; // true if this sample has an entry in the SAM file

    private HashMap<String, Object> properties = new HashMap<String, Object>();

    private HashMap<String, Sample> relationships = new HashMap<String, Sample>();

    public enum Gender {
        MALE,
        FEMALE,
        UNKNOWN
    }

    public Sample(String id) {
/*        if (id == null) {
            throw new StingException("Error creating sample: sample ID cannot be null");
        }*/
        this.id = id;
    }

    public String getId() {
        return this.id;
    }

    public Map<String, Object> getProperties() {
        return properties;
    }

    public void setProperties(Map<String, Object> properties) {
        this.properties = (HashMap) properties;
    }


    public void setSampleFileEntry(boolean value) {
        this.hasSampleFileEntry = value;
    }

    public boolean hasSAMFileEntry() {
        return this.hasSAMFileEntry;
    }

    public void setSAMFileEntry(boolean value) {
        this.hasSAMFileEntry = value;
    }

    public boolean hasSampleFileEntry() {
        return this.hasSampleFileEntry;
    }

    /**
     * Get one property
     * @param key key of property
     * @return value of property as generic object
     */
    public Object getProperty(String key) {
        return properties.get(key);
    }

    /**
     * Set a property
     * If property already exists, it is overwritten
     * @param key key of property
     * @param value object to be stored in properties array
     */
    public void setProperty(String key, Object value) {

        if (relationships.containsKey(key)) {
            throw new StingException("The same key cannot exist as a property and a relationship");
        }

        if (key.equals("gender") && value.getClass() != Gender.class) {
            throw new StingException("'gender' property must be of type Sample.Gender");
        }

        if (key.equals("population") && value.getClass() != String.class) {
            throw new StingException("'population' property must be of type String");
        }

        properties.put(key, value);
    }

    /**
     * Get one relationship
     * @param key of relationship
     * @return Sample object that this relationship points to
     */
    public Sample getRelationship(String key) {
        return relationships.get(key);
    }

    /**
     * Set one relationship
     * If already set, it is overwritten
     * @param key key of the relationship
     * @param value Sample object this relationship points to
     */
    public void setRelationship(String key, Sample value) {
        if (properties.containsKey(key)) {
            throw new StingException("The same key cannot exist as a property and a relationship");
        }
        relationships.put(key, value);
    }

    /**
     * Get the sample's mother
     * @return sample object with relationship mother, if exists, or null
     */
    public Sample getMother() {
        return getRelationship("mother");
    }

    /**
     * Get the sample's father
     * @return sample object with relationship father, if exists, or null
     */
    public Sample getFather() {
        return getRelationship("father");
    }

    /**
     * Get gender of the sample
     * @return property of key "gender" - must be of type Gender
     */
    public Gender getGender() {
        return (Gender) properties.get("gender");
    }

    public String getPopulation() {
        return (String) properties.get("population");
    }

    public String getFamilyId() {
        return (String) properties.get("familyId");
    }

    /**
     * @return True if sample is male, false if female, unknown, or null
     */
    public boolean isMale() {
        return properties.get("gender") == Gender.MALE;
    }

    /**
     * @return True if sample is female, false if male, unknown or null
     */
    public boolean isFemale() {
        return properties.get("gender") == Gender.MALE;
    }

    /**
     *
     * @param key property key
     * @return true if sample has this property (even if its value is null)
     */
    public boolean hasProperty(String key) {
        return properties.containsKey(key);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Sample sample = (Sample) o;

        if (hasSAMFileEntry != sample.hasSAMFileEntry) return false;
        if (hasSampleFileEntry != sample.hasSampleFileEntry) return false;
        if (id != null ? !id.equals(sample.id) : sample.id != null) return false;
        if (properties != null ? !properties.equals(sample.properties) : sample.properties != null) return false;
        if (relationships != null ? !relationships.equals(sample.relationships) : sample.relationships != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        return id != null ? id.hashCode() : "".hashCode();
    }
}
