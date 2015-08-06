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

package org.broadinstitute.gatk.utils.refdata.utils;


import org.broadinstitute.gatk.utils.commandline.Tags;

/**
 * a helper class to manage our triplets of data for the -B command line option (name, type, file)
 * TODO: The presence of four datapoints here suggests that this class' name isn't sufficient to describe its function.  Rename.
 */
public class RMDTriplet {
    public enum RMDStorageType { FILE, STREAM };

    private final String name;
    private final String type;
    private final String file;
    private final RMDStorageType storageType;
    private final Tags tags;

    public RMDTriplet(final String name, final String type, final String file, final RMDStorageType storageType, final Tags tags) {
        this.name = name;
        this.type = type;
        this.file = file;
        this.storageType = storageType;
        this.tags = tags;
    }

    /**
     * Gets the name of this track.  RefMetaDataTrackers can use this identifier to retrieve data of a certain type.
     * @return Name associated with this track.
     */
    public String getName() {
        return name;
    }

    /**
     * Gets the type of this track.  Informs the GATK how to parse this file type.
     * @return Type associated with this track.
     */
    public String getType() {
        return type;
    }

    /**
     * Gets the filename representing this track.  Data is loaded from this file.
     * @return Filename of the RMD.
     */
    public String getFile() {
        return file;
    }

    /**
     * The type of storage being used for this metadata track.  Right now, can be either a
     * file type (can be indexed) or a stream type (can't be indexed).
     * @return Storage type for this RMD 'triplet'.
     */
    public RMDStorageType getStorageType() {
        return storageType;
    }

    /**
     * Gets the key=value tags associated with this track
     * @return Tags associated with this track.
     */
    public Tags getTags() {
        return tags;
    }
}
