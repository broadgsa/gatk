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

package org.broadinstitute.gatk.utils.refdata.tracks;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.MutableIndex;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.utils.SequenceDictionaryUtils;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Utilities for working with Sequence Dictionaries embedded in tribble indices
 *
 * @author Your Name
 * @since Date created
 */
public class IndexDictionaryUtils {
    private final static Logger logger = Logger.getLogger(IndexDictionaryUtils.class);

    // a constant we use for marking sequence dictionary entries in the Tribble index property list
    public static final String SequenceDictionaryPropertyPredicate = "DICT:";

    /**
     * get the sequence dictionary from the track, if available.  If not, make it from the contig list that is always in the index
     * @param index the index file to use
     * @return a SAMSequenceDictionary if available, null if unavailable
     */
    public static SAMSequenceDictionary getSequenceDictionaryFromProperties(Index index) {
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (Map.Entry<String,String> entry : index.getProperties().entrySet()) {
            if (entry.getKey().startsWith(SequenceDictionaryPropertyPredicate))
                dict.addSequence(new SAMSequenceRecord(entry.getKey().substring(SequenceDictionaryPropertyPredicate.length() , entry.getKey().length()),
                        Integer.valueOf(entry.getValue())));
        }
        return dict;
    }

    /**
     * create the sequence dictionary with the contig list; a backup approach
     * @param index the index file to use
     * @param dict the sequence dictionary to add contigs to
     * @return the filled-in sequence dictionary
     */
    static SAMSequenceDictionary createSequenceDictionaryFromContigList(final Index index, final SAMSequenceDictionary dict) {
        final List<String> seqNames = index.getSequenceNames();
        if (seqNames == null) {
            return dict;
        }
        for (final String name : seqNames) {
            SAMSequenceRecord seq = new SAMSequenceRecord(name, 0);
            dict.addSequence(seq);
        }
        return dict;
    }

    /**
     *  Sets the sequence dictionary of the given index.  THE INDEX MUST BE MUTABLE (i.e. not Tabix).
     *
     * @param index the (mutable) index file to use
     * @param dict  the dictionary to use
     */
    public static void setIndexSequenceDictionary(Index index, SAMSequenceDictionary dict) {
        for ( SAMSequenceRecord seq : dict.getSequences() ) {
            final String contig = IndexDictionaryUtils.SequenceDictionaryPropertyPredicate + seq.getSequenceName();
            final String length = String.valueOf(seq.getSequenceLength());
            ((MutableIndex)index).addProperty(contig, length);
        }
    }

    public static void validateTrackSequenceDictionary(final String trackName,
                                                       final SAMSequenceDictionary trackDict,
                                                       final SAMSequenceDictionary referenceDict,
                                                       final ValidationExclusion.TYPE validationExclusionType ) {
        // if the sequence dictionary is empty (as well as null which means it doesn't have a dictionary), skip validation
        if (trackDict == null || trackDict.isEmpty())
            logger.warn("Track " + trackName + " doesn't have a sequence dictionary built in, skipping dictionary validation");
        else {
            Set<String> trackSequences = new TreeSet<String>();
            for (SAMSequenceRecord dictionaryEntry : trackDict.getSequences())
                trackSequences.add(dictionaryEntry.getSequenceName());
            SequenceDictionaryUtils.validateDictionaries(logger, validationExclusionType, trackName, trackDict, "reference", referenceDict, false, null);
        }
    }
}
