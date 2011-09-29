package org.broadinstitute.sting.gatk.samples;


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
    final private String familyID, paternalID, maternalID;
    final private Sample.Gender gender;
    final private double quantitativePhenotype;
    final private Sample.Affection affection;
    final private String population;
    final private String ID;
    final private SampleDataSource dataSource;

    private boolean hasSAMFileEntry = false; // true if this sample has an entry in the SAM file
    private Map<String, Object> properties = new HashMap<String, Object>();

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

    public Sample(final String ID, final SampleDataSource dataSource,
                  final String familyID, final String paternalID, final String maternalID,
                  final Gender gender, final double quantitativePhenotype, final Affection affection,
                  final String population) {
        this.familyID = familyID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.gender = gender;
        this.quantitativePhenotype = quantitativePhenotype;
        this.affection = affection;
        this.population = population;
        this.ID = ID;
        this.dataSource = dataSource;
    }

    public Sample(String id, SampleDataSource dataSource) {
        this(id, dataSource,
                null, null, null,
                Gender.UNKNOWN, UNSET_QUANTITIATIVE_TRAIT_VALUE, Affection.UNKNOWN, null);
    }

    @Deprecated
    public boolean hasSAMFileEntry() {
        return this.hasSAMFileEntry;
    }

    @Deprecated
    public void setSAMFileEntry(boolean value) {
        this.hasSAMFileEntry = value;
    }

    // -------------------------------------------------------------------------------------
    //
    // standard property getters
    //
    // -------------------------------------------------------------------------------------

    public String getID() {
        return ID;
    }


    public String getFamilyID() {
        return familyID;
    }

    public String getPaternalID() {
        return paternalID;
    }

    public String getMaternalID() {
        return maternalID;
    }

    public Affection getAffection() {
        return affection;
    }

    public boolean hasQuantitativeTrait() {
        return affection == Affection.QUANTITATIVE;
    }

    public double getQuantitativePhenotype() {
        return quantitativePhenotype;
    }

    /**
     * Get the sample's mother
     * @return sample object with relationship mother, if exists, or null
     */
    public Sample getMother() {
        return dataSource.getSampleById(maternalID);
    }

    /**
     * Get the sample's father
     * @return sample object with relationship father, if exists, or null
     */
    public Sample getFather() {
        return dataSource.getSampleById(paternalID);
    }

    /**
     * Get gender of the sample
     * @return property of key "gender" - must be of type Gender
     */
    public Gender getGender() {
        return gender;
    }

    public String getPopulation() {
        return population;
    }

    public String getFamilyId() {
        return familyID;
    }

    /**
     * @return True if sample is male, false if female, unknown, or null
     */
    public boolean isMale() {
        return getGender() == Gender.MALE;
    }

    /**
     * @return True if sample is female, false if male, unknown or null
     */
    public boolean isFemale() {
        return getGender() == Gender.MALE;
    }

    // -------------------------------------------------------------------------------------
    //
    // code for working with additional -- none standard -- properites
    //
    // -------------------------------------------------------------------------------------

    public Map<String, Object> getExtraProperties() {
        return Collections.unmodifiableMap(properties);
    }

    /**
     * Get one property
     * @param key key of property
     * @return value of property as generic object
     */
    public Object getExtraPropertyValue(final String key) {
        return properties.get(key);
    }

    /**
     *
     * @param key property key
     * @return true if sample has this property (even if its value is null)
     */
    public boolean hasExtraProperty(String key) {
        return properties.containsKey(key);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Sample sample = (Sample) o;

        if (hasSAMFileEntry != sample.hasSAMFileEntry) return false;
        if (ID != null ? !ID.equals(sample.ID) : sample.ID != null) return false;
        if (properties != null ? !properties.equals(sample.properties) : sample.properties != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return ID != null ? ID.hashCode() : "".hashCode();
    }
}
