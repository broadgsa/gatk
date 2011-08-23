/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.NameAwareCodec;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.gatk.refdata.SelfScopingFeatureCodec;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.help.GATKDocUtils;
import org.broadinstitute.sting.utils.help.HelpUtils;

import javax.mail.Header;
import java.io.File;
import java.util.*;


/**
 * Class for managing Tribble Feature readers available to the GATK.  The features
 * are dynamically determined via a PluginManager.  This class provides convenient
 * getter methods for obtaining FeatureDescriptor objects that collect all of the
 * useful information about the Tribble Codec, Feature, and name in one place.
 *
 * @author depristo
 */
public class FeatureManager  {
    public static class FeatureDescriptor implements Comparable<FeatureDescriptor> {
        final String name;
        final FeatureCodec codec;

        public FeatureDescriptor(final String name, final FeatureCodec codec) {
            this.name = name;
            this.codec = codec;
        }

        public String getName() {
            return name;
        }
        public String getSimpleFeatureName() { return getFeatureClass().getSimpleName(); }
        public FeatureCodec getCodec() {
            return codec;
        }
        public Class getCodecClass() { return codec.getClass(); }
        public Class getFeatureClass() { return codec.getFeatureType(); }

        @Override
        public String toString() {
            return String.format("FeatureDescriptor name=%s codec=%s feature=%s",
                    getName(), getCodecClass().getName(), getFeatureClass().getName());
        }

        @Override
        public int compareTo(FeatureDescriptor o) {
            return getName().compareTo(o.getName());
        }
    }

    private final PluginManager<FeatureCodec> pluginManager;
    private final Collection<FeatureDescriptor> featureDescriptors = new TreeSet<FeatureDescriptor>();

    /**
     * Construct a FeatureManager
     */
    public FeatureManager() {
        pluginManager = new PluginManager<FeatureCodec>(FeatureCodec.class, "Codecs", "Codec");

        for (final String rawName: pluginManager.getPluginsByName().keySet()) {
            FeatureCodec codec = pluginManager.createByName(rawName);
            String name = rawName.toUpperCase();
            FeatureDescriptor featureDescriptor = new FeatureDescriptor(name, codec);
            featureDescriptors.add(featureDescriptor);
        }
    }

    /**
     * Return the FeatureDescriptor whose getCodecClass().equals(codecClass).
     *
     * @param codecClass
     * @return A FeatureDescriptor or null if none is found
     */
    @Requires("codecClass != null")
    public FeatureDescriptor getByCodec(Class codecClass) {
        for ( FeatureDescriptor descriptor : featureDescriptors )
            if ( descriptor.getCodecClass().equals(codecClass) )
                return descriptor;
        return null;
    }

    /**
     * Returns a collection of FeatureDescriptors that emit records of type featureClass
     *
     * @param featureClass
     * @return A FeatureDescriptor or null if none is found
     */
    @Requires("featureClass != null")
    public <T extends Feature> Collection<FeatureDescriptor> getByFeature(Class<T> featureClass) {
        Set<FeatureDescriptor> consistentDescriptors = new TreeSet<FeatureDescriptor>();

        if (featureClass == null)
            throw new IllegalArgumentException("trackRecordType value is null, please pass in an actual class object");

        for ( FeatureDescriptor descriptor : featureDescriptors ) {
            if ( featureClass.isAssignableFrom(descriptor.getFeatureClass()))
                consistentDescriptors.add(descriptor);
        }
        return consistentDescriptors;
    }

    /**
     * Return the FeatureDescriptor with getName().equals(name)
     *
     * @param name
     * @return A FeatureDescriptor or null if none is found
     */
    @Requires("name != null")
    public FeatureDescriptor getByName(String name) {
        for ( FeatureDescriptor descriptor : featureDescriptors )
            if ( descriptor.getName().equalsIgnoreCase(name) )
                return descriptor;
        return null;
    }

    /**
     * Returns the FeatureDescriptor that can read the contexts of File file, is one can be determined
     *
     * @param file
     * @return A FeatureDescriptor or null if none is found
     */
    @Requires({"file != null", "file.isFile()", "file.canRead()"})
    public FeatureDescriptor getByFiletype(File file) {
        List<FeatureDescriptor> canParse = new ArrayList<FeatureDescriptor>();
        for ( FeatureDescriptor descriptor : featureDescriptors )
            if ( descriptor.getCodec() instanceof SelfScopingFeatureCodec ) {
                if ( ((SelfScopingFeatureCodec) descriptor.getCodec()).canDecode(file) ) {
                    canParse.add(descriptor);
                }
            }

        if ( canParse.size() == 0 )
            return null;
        else if ( canParse.size() > 1 )
            throw new ReviewedStingException("BUG: multiple feature descriptors can read file " + file + ": " + canParse);
        else
            return canParse.get(0);
    }

    /**
     * Returns the FeatureDescriptor associated with the type described by triplet, or null if none is found
     * @param triplet
     * @return
     */
    @Requires("triplet != null")
    public FeatureDescriptor getByTriplet(RMDTriplet triplet) {
        return getByName(triplet.getType());
    }

    /**
     * @return all of the FeatureDescriptors available to the GATK.  Never null
     */
    @Ensures("result != null")
    public Collection<FeatureDescriptor> getFeatureDescriptors() {
        return Collections.unmodifiableCollection(featureDescriptors);
    }


    /**
     * Returns a list of the available tribble track names (vcf,dbsnp,etc) that we can load
     * @return
     */
    @Ensures("result != null")
    public String userFriendlyListOfAvailableFeatures() {
        return userFriendlyListOfAvailableFeatures(Feature.class);
    }

    /**
     * Returns a list of the available tribble track names (vcf,dbsnp,etc) that we can load
     * restricted to only Codecs producting Features consistent with the requiredFeatureType
     * @return
     */
    @Ensures("result != null")
    public String userFriendlyListOfAvailableFeatures(Class<? extends Feature> requiredFeatureType) {
        final String nameHeader="Name", featureHeader = "FeatureType", docHeader="Documentation";

        int maxNameLen = nameHeader.length(), maxFeatureNameLen = featureHeader.length();
        for ( final FeatureDescriptor descriptor : featureDescriptors ) {
            if ( requiredFeatureType.isAssignableFrom(descriptor.getFeatureClass()) ) {
                maxNameLen = Math.max(maxNameLen, descriptor.getName().length());
                maxFeatureNameLen = Math.max(maxFeatureNameLen, descriptor.getSimpleFeatureName().length());
            }
        }

        StringBuilder docs = new StringBuilder();
        String format = "%" + maxNameLen + "s   %" + maxFeatureNameLen + "s   %s%n";
        docs.append(String.format(format, nameHeader, featureHeader, docHeader));
        for ( final FeatureDescriptor descriptor : featureDescriptors ) {
            if ( requiredFeatureType.isAssignableFrom(descriptor.getFeatureClass()) ) {
                String oneDoc = String.format(format,
                        descriptor.getName(),
                        descriptor.getSimpleFeatureName(),
                        GATKDocUtils.helpLinksToGATKDocs(descriptor.getCodecClass()));
                docs.append(oneDoc);
            }
        }

        return docs.toString();
    }

    /**
     * Create a new FeatureCodec of the type described in descriptor, assigning it the
     * name (if possible) and providing it the genomeLocParser (where necessary)
     *
     * @param descriptor FeatureDescriptor of the Tribble FeatureCodec we want to create
     * @param name the name to assign this codec
     * @return the feature codec itself
     */
    @Requires({"descriptor != null", "name != null", "genomeLocParser != null"})
    @Ensures("result != null")
    public FeatureCodec createCodec(FeatureDescriptor descriptor, String name, GenomeLocParser genomeLocParser) {
        FeatureCodec codex = pluginManager.createByType(descriptor.getCodecClass());
        if ( codex instanceof NameAwareCodec )
            ((NameAwareCodec)codex).setName(name);
        if ( codex instanceof ReferenceDependentFeatureCodec )
            ((ReferenceDependentFeatureCodec)codex).setGenomeLocParser(genomeLocParser);
        return codex;
    }
}
