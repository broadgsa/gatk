package org.broadinstitute.sting.gatk.datasources.sample;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.yaml.snakeyaml.Loader;
import org.yaml.snakeyaml.TypeDescription;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.Constructor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Jul 26, 2010
 * Time: 3:30:09 PM
 *
 * This class stores and manages sample metadata. This data is encoded in a sample file, which can be included
 * in the GATK by the "--samples" argument. This class reads and parses those files.
 *
 * Although there are a set of public methods for accessing sample data, they aren't used by walkers - they are really
 * only used by GenomeAnalysisEngine. An instance of GenomeAnalysisEngine has one SampleDataSource. When a walker
 * wants to access sample data, it asks GenomeAnalysis to fetch this data from its SampleDataSource.
 *
 */
public class SampleDataSource {

    /**
     * This is where Sample objects are stored. Samples are usually accessed by their ID, which is unique, so
     * this is stored as a HashMap.
     */
    private HashMap<String, Sample> samples = new HashMap<String, Sample>();

    /**
     * Samples can have "aliases", because sometimes the same sample is referenced by different IDs in different
     * datasets. If this is the case, one ID is the "primary ID" and others are "aliases".
     *
     * This maps ID => primary ID for all samples ID strings - both primary IDs and aliases.
     */
    private HashMap<String, String> sampleAliases = new HashMap<String, String>();

    /**
     * While loading sample files, we must be aware of "special" properties and relationships that are always allowed
     */
    public static final String[] specialProperties = new String[] {"familyId", "population", "gender"};
    public static final String[] specialRelationships = new String[] {"mother", "father"};

    /**
     * Constructor takes both a SAM header and sample files because the two must be integrated.
     * @param header SAMFileHeader that has been created for this analysis
     * @param sampleFiles Sample files that were included on the command line
     */
    public SampleDataSource(SAMFileHeader header, List<File> sampleFiles) {

        // create empty sample object for each sample referenced in the SAM header
        for (String sampleName : SampleUtils.getSAMFileSamples(header)) {
            if (!hasSample(sampleName)) {
                Sample newSample = new Sample(sampleName);
                newSample.setSAMFileEntry(true);
                samples.put(sampleName, newSample);
            }
        }

        // add files consecutively
        if (sampleFiles != null) {
            for (File file : sampleFiles) {
                addFile(file);
            }
        }
    }

    /**
     * Hallucinates sample objects for all the samples in the SAM file and stores them
     */
    private void getSamplesFromSAMFile() {
        for (String sampleName : SampleUtils.getSAMFileSamples(GenomeAnalysisEngine.instance.getSAMFileHeader())) {
            if (!hasSample(sampleName)) {
                Sample newSample = new Sample(sampleName);
                newSample.setSAMFileEntry(true);
                samples.put(sampleName, newSample);
            }
        }
    }

    /**
     * Parse one sample file and integrate it with samples that are already there
     * Fail quickly if we find any errors in the file
     */
    private void addFile(File sampleFile) {

        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(sampleFile));
        }
        catch (IOException e) {
            throw new StingException("Could not open sample file " + sampleFile.getAbsolutePath(), e);
        }

        // set up YAML reader - a "Constructor" creates java object from YAML and "Loader" loads the file
        Constructor con = new Constructor(SampleFileParser.class);
        TypeDescription desc = new TypeDescription(SampleFileParser.class);
        desc.putListPropertyType("propertyDefinitions", PropertyDefinition.class);
        desc.putListPropertyType("sampleAliases", SampleAlias.class);
        con.addTypeDescription(desc);
        Loader loader = new Loader(con);
        Yaml yaml = new Yaml(loader);

        // SampleFileParser stores an object representation of a sample file - this is what we'll parse
        SampleFileParser parser;
        try {
            parser = (SampleFileParser) yaml.load(reader);
        }
        catch (Exception e) {    // TODO: should we have more granular exception here?
            throw new StingException("There was a syntactic error with the YAML in sample file " + sampleFile.getAbsolutePath(), e);
        }

        // check to see which validation options were built into the file
        boolean restrictProperties = parser.getAllowedProperties() != null;
        boolean restrictRelationships = parser.getAllowedRelationships() != null;
        boolean restrictPropertyValues = parser.getPropertyDefinitions() != null;

        // propertyValues stores the values that are allowed for a given property
        HashMap<String, HashSet> propertyValues = null;
        if (restrictPropertyValues) {
            propertyValues = new HashMap<String, HashSet>();
            for (PropertyDefinition def : parser.getPropertyDefinitions()) {
                HashSet<String> set = new HashSet<String>();
                for (String value : def.getValues()) {
                    set.add(value);
                }
                propertyValues.put(def.getProperty(), set);
            }
        }

        // make sure the aliases are valid
        validateAliases(parser);

        // loop through each sample in the file - a SampleParser stores an object that will become a Sample
        for (SampleParser sampleParser : parser.getSamples()) {

            // step 1: add the sample if it doesn't already exist
            Sample sample = getSampleById(sampleParser.getId());
            if (sample == null) {
                sample = new Sample(sampleParser.getId());
            }
            addSample(sample);
            sample.setSampleFileEntry(true);

            // step 2: add the properties
            if (sampleParser.getProperties() != null) {
                for (String property : sampleParser.getProperties().keySet()) {

                    // check that property is allowed
                    if (restrictProperties) {
                        if (!isPropertyValid(property, parser.getAllowedProperties())) {
                            throw new StingException(property + " is an invalid property. It is not included in the list " +
                                    "of allowed properties.");
                        }
                    }

                    // next check that the value is allowed
                    if (restrictPropertyValues) {
                        if (!isValueAllowed(property, sampleParser.getProperties().get(property), propertyValues)) {
                            throw new StingException("The value of property '" + property + "' is invalid. " +
                                    "It is not included in the list of allowed values for this property.");
                        }
                    }

                    // next check that there isn't already a conflicting property there
                    if (sample.getProperty(property) != null &&
                            sample.getProperty(property) != sampleParser.getProperties().get(property))
                    {
                        throw new StingException(property + " is a conflicting property!");
                    }

                    // checks are passed - now add the property!
                    saveProperty(sample, property, sampleParser.getProperties().get(property));
                }
            }

            // step 3: add the relationships
            if (sampleParser.getRelationships() != null) {
                for (String relationship : sampleParser.getRelationships().keySet()) {
                    String relativeId = sampleParser.getRelationships().get(relationship);
                    if (relativeId == null) {
                        throw new StingException("The relationship cannot be null");
                    }

                    // first check that it's not invalid
                    if (restrictRelationships) {
                        if (!isRelationshipValid(relationship, parser.getAllowedRelationships())) {
                            throw new StingException(relationship + " is an invalid relationship");
                        }
                    }

                    // next check that there isn't already a conflicting property there
                    if (sample.getRelationship(relationship) != null) {
                        if (sample.getRelationship(relationship).getId() != sampleParser.getProperties().get(relationship)) {
                            throw new StingException(relationship + " is a conflicting relationship!");
                        }
                        // if the relationship is already set - and consistent with what we're reading now - no need to continue
                        else {
                            continue;
                        }
                    }

                    // checks are passed - now save the relationship
                    saveRelationship(sample, relationship, relativeId);
                }
            }
        }

    }

    private boolean isValueAllowed(String key, Object value, HashMap<String, HashSet> valuesList) {

        // if the property values weren't specified for this property, then any value is okay
        if (!valuesList.containsKey(key)) {
            return true;
        }

        // if this property has enumerated values, it must be a string
        else if (value.getClass() != String.class)
            return false;

        // is the value specified or not?
        else if (!valuesList.get(key).contains(value))
            return false;

        return true;
    }

    /**
     * Makes sure that the aliases are valid
     * Checks that 1) no string is used as both a main ID and an alias;
     * 2) no alias is used more than once
     * @param parser
     */
    private void validateAliases(SampleFileParser parser) {

        // no aliases sure validate
        if (parser.getSampleAliases() == null)
            return;

        HashSet<String> mainIds = new HashSet<String>();
        HashSet<String> otherIds = new HashSet<String>();

        for (SampleAlias sampleAlias : parser.getSampleAliases()) {
            mainIds.add(sampleAlias.getMainId());
            for (String otherId : sampleAlias.getOtherIds()) {
                if (mainIds.contains(otherId))
                    throw new StingException(String.format("The aliases in your sample file are invalid - the alias %s cannot " +
                            "be both a main ID and an other ID", otherId));

                if (!otherIds.add(otherId))
                    throw new StingException(String.format("The aliases in your sample file are invalid - %s is listed as an " +
                            "alias more than once.", otherId));
            }
        }
    }

    private boolean isPropertyValid(String property, String[] allowedProperties) {

        // is it a special property that is always allowed?
        for (String allowedProperty : specialProperties) {
            if (property.equals(allowedProperty))
                return true;
        }

        // is it in the allowed properties list?
        for (String allowedProperty : allowedProperties) {
            if (property.equals(allowedProperty))
                return true;
        }

        return false;
    }

    private boolean isRelationshipValid(String relationship, String[] allowedRelationships) {

        // is it a special relationship that is always allowed?
        for (String allowedRelationship : specialRelationships) {
            if (relationship.equals(allowedRelationship))
                return true;
        }

        // is it in the allowed properties list?
        for (String allowedRelationship : allowedRelationships) {
            if (relationship.equals(allowedRelationship))
                return true;
        }

        return false;
    }

    /**
     * Saves a property as the correct type
     * @param key property key
     * @param value property value, as read from YAML parser
     * @return property value to be stored
     */
    private void saveProperty(Sample sample, String key, Object value) {

        // convert gender to the right type, if it was stored as a String
        if (key.equals("gender")) {
            if (((String) value).toLowerCase().equals("male")) {
                value = Sample.Gender.MALE;
            }
            else if (((String) value).toLowerCase().equals("female")) {
                value = Sample.Gender.FEMALE;
            }
            else  if (((String) value).toLowerCase().equals("unknown")) {
                value = Sample.Gender.UNKNOWN;
            }
            else if (value != null) {
                throw new StingException("'gender' property must be male, female, or unknown.");
            }
            value = null;
        }
        sample.setProperty(key, value);
    }

    /**
     * Saves a relationship as the correct type
     * @param key relationship key
     * @param relativeId sample ID string of the relative
     * @return relationship value to be stored
     */
    private void saveRelationship(Sample sample, String key, String relativeId) {

        // get the reference that we'll store as the value
        Sample relative = getSampleById(relativeId);

        // create sample object for the relative, if necessary
        if (relative == null) {
            relative = new Sample(relativeId);
            addSample(relative);
        }
        sample.setRelationship(key, relative);
    }



    /**
     * Filter a sample name in case it is an alias
     * @param sampleId to be filtered
     * @return ID of sample that stores data for this alias
     */
    private String aliasFilter(String sampleId) {
        if (!sampleAliases.containsKey(sampleId))
            return sampleId;
        else
            return sampleAliases.get(sampleId);
    }

    /**
     * Add a sample to the collection
     * @param sample to be added
     */
    private void addSample(Sample sample) {
        samples.put(sample.getId(), sample);
    }

    /**
     * Check if sample with this ID exists
     * Note that this will return true if name passed in is an alias
     * @param id ID of sample to be checked
     * @return true if sample exists; false if not
     */
    public boolean hasSample(String id) {
        return samples.get(aliasFilter(id)) != null;
    }

    /**
     * Get a sample by its ID
     * If an alias is passed in, return the main sample object 
     * @param id
     * @return sample Object with this ID
     */
    public Sample getSampleById(String id) {
        return samples.get(aliasFilter(id));
    }

    /**
     * Get the sample for a given read group
     * Must first look up ID for read group
     * @param readGroup of sample
     * @return sample object with ID from the read group
     */
    public Sample getSampleByReadGroup(SAMReadGroupRecord readGroup) {
        String nameFromReadGroup = readGroup.getSample();
        return getSampleById(nameFromReadGroup);
    }

    /**
     * Get a sample for a given read
     * Must first look up read group, and then sample ID for that read group
     * @param read of sample
     * @return sample object of this read
     */
    public Sample getSampleByRead(SAMRecord read) {
        return getSampleByReadGroup(read.getReadGroup());
    }

    /**
     * Get number of sample objects
     * @return size of samples map
     */
    public int sampleCount() {
        return samples.size();
    }

    /**
     * Return all samples with a given family ID
     * Note that this isn't terribly efficient (linear) - it may be worth adding a new family ID data structure for this
     * @param familyId
     * @return
     */
    public Set<Sample> getFamily(String familyId) {
        HashSet<Sample> familyMembers = new HashSet<Sample>();

        for (Sample sample : samples.values()) {
            if (sample.getFamilyId() != null) {
                if (sample.getFamilyId().equals(familyId))
                    familyMembers.add(sample);
            }
        }
        return familyMembers;
    }

    /**
     * Returns all children of a given sample
     * See note on the efficiency of getFamily() - since this depends on getFamily() it's also not efficient
     * @param sample
     * @return
     */
    public Set<Sample> getChildren(Sample sample) {
        HashSet<Sample> children = new HashSet<Sample>();
        for (Sample familyMember : getFamily(sample.getFamilyId())) {
            if (familyMember.getMother() == sample || familyMember.getFather() == sample) {
                children.add(familyMember);
            }
        }
        return children;
    }


}