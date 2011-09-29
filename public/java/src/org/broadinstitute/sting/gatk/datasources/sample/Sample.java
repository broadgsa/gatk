package org.broadinstitute.sting.gatk.datasources.sample;


import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Jul 26, 2010
 * Time: 3:31:38 PM
 */
public class Sample implements java.io.Serializable {
    private final static String MOTHER = "mother";
    private final static String FATHER = "father";
    private final static String GENDER = "gender";
    private final static String POPULATION = "population";
    private final static String FAMILY = "familyId";
    private final static String AFFECTION = "affection";
    private final static String QUANT_TRAIT = "quantTrait";

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

    public enum Affection {
        /** Status is unknown */
        UNKNOWN,
        /** Suffers from the disease */
        AFFECTED,
        /** Unaffected by the disease */
        UNAFFECTED,
        /** A quantitative trait: value of the trait is stored elsewhere */
        QUANTITATIVE
    }
    public final static double UNSET_QUANTITIATIVE_TRAIT_VALUE = Double.NaN;

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

    @Deprecated
    public void setSampleFileEntry(boolean value) {
        this.hasSampleFileEntry = value;
    }

    @Deprecated
    public boolean hasSAMFileEntry() {
        return this.hasSAMFileEntry;
    }

    @Deprecated
    public void setSAMFileEntry(boolean value) {
        this.hasSAMFileEntry = value;
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

        if (key.equals(GENDER) && value.getClass() != Gender.class) {
            throw new StingException("'gender' property must be of type Sample.Gender");
        }

        if (key.equals(POPULATION) && value.getClass() != String.class) {
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
        return getRelationship(MOTHER);
    }

    /**
     * Get the sample's father
     * @return sample object with relationship father, if exists, or null
     */
    public Sample getFather() {
        return getRelationship(FATHER);
    }

    /**
     * Get gender of the sample
     * @return property of key "gender" - must be of type Gender
     */
    public Gender getGender() {
        return (Gender) properties.get(GENDER);
    }

    public String getPopulation() {
        return (String) properties.get(POPULATION);
    }

    public String getFamilyId() {
        return (String) properties.get(FAMILY);
    }

    /**
     * @return True if sample is male, false if female, unknown, or null
     */
    public boolean isMale() {
        return properties.get(GENDER) == Gender.MALE;
    }

    /**
     * @return True if sample is female, false if male, unknown or null
     */
    public boolean isFemale() {
        return properties.get(GENDER) == Gender.MALE;
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
