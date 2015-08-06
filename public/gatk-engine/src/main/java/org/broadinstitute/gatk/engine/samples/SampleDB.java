/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.samples;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import htsjdk.variant.variantcontext.Genotype;

import java.util.*;

/**
 *
 */
public class SampleDB {
    /**
     * This is where Sample objects are stored. Samples are usually accessed by their ID, which is unique, so
     * this is stored as a HashMap.
     */
    private final HashMap<String, Sample> samples = new HashMap<String, Sample>();

    /**
     * Constructor takes both a SAM header and sample files because the two must be integrated.
     */
    public SampleDB() {

    }

    /**
     * Protected function to add a single sample to the database
     *
     * @param sample to be added
     */
    protected SampleDB addSample(Sample sample) {
        Sample prev = samples.get(sample.getID());
        if ( prev != null )
            sample = Sample.mergeSamples(prev, sample);
        samples.put(sample.getID(), sample);
        return this;
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for getting a sample from the DB
    //
    // --------------------------------------------------------------------------------

    /**
     * Get a sample by its ID
     * If an alias is passed in, return the main sample object 
     * @param id
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(String id) {
        return samples.get(id);
    }

    /**
     *
     * @param read
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(final SAMRecord read) {
        return getSample(read.getReadGroup());
    }

    /**
     *
     * @param rg
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(final SAMReadGroupRecord rg) {
        return getSample(rg.getSample());
    }

    /**
     * @param g Genotype
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(final Genotype g) {
        return getSample(g.getSampleName());
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for accessing samples in the DB
    //
    // --------------------------------------------------------------------------------

    /**
     * Get number of sample objects
     * @return size of samples map
     */
    public int sampleCount() {
        return samples.size();
    }

    public Set<Sample> getSamples() {
        return new LinkedHashSet<>(samples.values());
    }

    public Collection<String> getSampleNames() {
        return Collections.unmodifiableCollection(samples.keySet());
    }


    /**
     * Takes a collection of sample names and returns their corresponding sample objects
     * Note that, since a set is returned, if you pass in a list with duplicates names there will not be any duplicates in the returned set
     * @param sampleNameList Set of sample names
     * @return Corresponding set of samples
     */
    public Set<Sample> getSamples(Collection<String> sampleNameList) {
        HashSet<Sample> samples = new HashSet<Sample>();
        for (String name : sampleNameList) {
            try {
                samples.add(getSample(name));
            }
            catch (Exception e) {
                throw new GATKException("Could not get sample with the following ID: " + name, e);
            }
        }
        return samples;
    }

    // --------------------------------------------------------------------------------
    //
    // Higher level pedigree functions
    //
    // --------------------------------------------------------------------------------

    /**
     * Returns a sorted set of the family IDs in all samples (excluding null ids)
     * @return
     */
    public final Set<String> getFamilyIDs() {
        return getFamilies().keySet();
    }

    /**
     * Returns a map from family ID -> set of family members for all samples with
     * non-null family ids
     *
     * @return
     */
    public final Map<String, Set<Sample>> getFamilies() {
        return getFamilies(null);
    }

    /**
     * Returns a map from family ID -> set of family members for all samples in sampleIds with
     * non-null family ids
     *
     * @param sampleIds - all samples to include. If null is passed then all samples are returned.
     * @return
     */
    public final Map<String, Set<Sample>> getFamilies(Collection<String> sampleIds) {
        final Map<String, Set<Sample>> families = new TreeMap<String, Set<Sample>>();

        for ( final Sample sample : samples.values() ) {
            if(sampleIds == null || sampleIds.contains(sample.getID())){
                final String famID = sample.getFamilyID();
                if ( famID != null ) {
                    if ( ! families.containsKey(famID) )
                        families.put(famID, new TreeSet<Sample>());
                    families.get(famID).add(sample);
                }
            }
        }
        return families;
    }

    /**
     * Returns all the trios present in the sample database. The strictOneChild parameter determines
     * whether multiple children of the same parents resolve to multiple trios, or are excluded
     * @param strictOneChild - exclude pedigrees with >1 child for parental pair
     * @return - all of the mother+father=child triplets, subject to strictOneChild
     */
    public final Set<Trio> getTrios(boolean strictOneChild) {
        Set<Trio> trioSet = new HashSet<Trio>();
        for ( String familyString : getFamilyIDs() ) {
            Set<Sample> family = getFamily(familyString);
            for ( Sample sample : family) {
                if ( sample.getParents().size() == 2 ) {
                    Trio trio = new Trio(sample.getMother(),sample.getFather(),sample);
                    trioSet.add(trio);
                }
            }
        }

        if ( strictOneChild )
            trioSet = removeTriosWithSameParents(trioSet);

        return trioSet;
    }

    /**
     * Returns all the trios present in the db. See getTrios(boolean strictOneChild)
     * @return all the trios present in the samples db.
     */
    public final Set<Trio> getTrios() {
        return getTrios(false);
    }

    /**
     * Subsets a set of trios to only those with nonmatching founders. If two (or more) trio objects have
     * the same mother and father, then both (all) are removed from the returned set.
     * @param trios - a set of Trio objects
     * @return those subset of Trio objects in the input set with nonmatching founders
     */
    private Set<Trio> removeTriosWithSameParents(final Set<Trio> trios) {
        Set<Trio> filteredTrios = new HashSet<Trio>();
        filteredTrios.addAll(trios);
        Set<Trio> triosWithSameParents = new HashSet<Trio>();
        for ( Trio referenceTrio : filteredTrios ) {
            for ( Trio compareTrio : filteredTrios ) {
                if ( referenceTrio != compareTrio &&
                     referenceTrio.getFather().equals(compareTrio.getFather()) &&
                     referenceTrio.getMother().equals(compareTrio.getMother()) ) {
                    triosWithSameParents.add(referenceTrio);
                    triosWithSameParents.add(compareTrio);
                }
            }
        }
        filteredTrios.removeAll(triosWithSameParents);
        return filteredTrios;
    }

    /**
     * Returns the set of all children that have both of their parents.
     * Note that if a family is composed of more than 1 child, each child is
     * returned.
     * @return - all the children that have both of their parents
     * @deprecated - getTrios() replaces this function
     */
    @Deprecated
    public final Set<Sample> getChildrenWithParents(){
        return getChildrenWithParents(false);
    }

    /**
     * Returns the set of all children that have both of their parents.
     * Note that if triosOnly = false, a family is composed of more than 1 child, each child is
     * returned.
     *
     * This method can be used wherever trios are needed
     *
     * @param triosOnly - if set to true, only strict trios are returned
     * @return - all the children that have both of their parents
     * @deprecated - getTrios(boolean strict) replaces this function
     * @bug -- does not work for extracting multiple generations of trios, e.g.
     * ..........Mom1------Dad1
     * ................|
     * ..............Child1--------Mom2
     * .......................|
     * .....................Child2
     */
    @Deprecated
    public final Set<Sample> getChildrenWithParents(boolean triosOnly) {

        Map<String, Set<Sample>> families = getFamilies();
        final Set<Sample> childrenWithParents = new HashSet<Sample>();
        Iterator<Sample> sampleIterator;

        for ( Set<Sample> familyMembers: families.values() ) {
            if(triosOnly && familyMembers.size() != 3)
                continue;

            sampleIterator = familyMembers.iterator();
            Sample sample;
            while(sampleIterator.hasNext()){
                sample = sampleIterator.next();
                if(sample.getParents().size() == 2 && familyMembers.containsAll(sample.getParents()))
                    childrenWithParents.add(sample);
            }

        }
        return childrenWithParents;
    }

    /**
     * Return all samples with a given family ID
     * @param familyId
     * @return
     */
    public Set<Sample> getFamily(String familyId) {
        return getFamilies().get(familyId);
    }

    /**
     * Returns all children of a given sample
     * See note on the efficiency of getFamily() - since this depends on getFamily() it's also not efficient
     * @param sample
     * @return
     */
    public Set<Sample> getChildren(Sample sample) {
        final HashSet<Sample> children = new HashSet<Sample>();
        for ( final Sample familyMember : getFamily(sample.getFamilyID())) {
            if ( familyMember.getMother() == sample || familyMember.getFather() == sample ) {
                children.add(familyMember);
            }
        }
        return children;
    }

    public Set<String> getFounderIds(){
        Set<String> founders = new HashSet<String>();
        for(Sample sample : getSamples()){
            if(sample.getParents().size()<1)
                founders.add(sample.getID());

        }
        return founders;
    }
}
