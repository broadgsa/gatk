/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class ReadMetaDataTracker
 *         <p/>
 *         a read-based meta data tracker
 */
public class ReadMetaDataTracker {
    private final SAMRecord record;
    private final TreeMap<Long, Set<ReferenceOrderedDatum>> mapping;

    /**
     * create a read meta data tracker, given the read and a queue of RODatum positions
     *
     * @param record the read to create offset from
     * @param mapping  the mapping of reference ordered datum
     */
    public ReadMetaDataTracker(SAMRecord record, TreeMap<Long, Set<ReferenceOrderedDatum>> mapping) {
        this.record = record;
        this.mapping = mapping;
    }

    /**
     * create an alignment of read position to reference ordered datum
     *
     * @param record the SAMRecord
     * @param queue the queue (as a tree set)
     * @param cl the class name, null if not filtered by classname
     * @param name the datum track name, null if not filtered by name
     * 
     * @return a mapping from the position in the read to the reference ordered datum
     */
    private Map<Integer, Set<ReferenceOrderedDatum>> createReadAlignment(SAMRecord record, TreeMap<Long, Set<ReferenceOrderedDatum>> queue, Class cl, String name) {
        Map<Integer, Set<ReferenceOrderedDatum>> ret = new LinkedHashMap<Integer, Set<ReferenceOrderedDatum>>();
        GenomeLoc location = GenomeLocParser.createGenomeLoc(record);
        int length = record.getReadLength();
        for (Long loc : queue.keySet()) {
            //if (location.containsP(loc)) {
                long position = loc - location.getStart();
                if (position >= 0 && position < length) ret.put((int)(position),queue.get(loc));
        }
        return ret;
  
    }

    /**
     * create an alignment of read position to reference ordered datum
     *
     * @return a mapping from the position in the read to the reference ordered datum
     */
    private Map<Long, Set<ReferenceOrderedDatum>> createGenomeLocAlignment(SAMRecord record, TreeMap<Long, Set<ReferenceOrderedDatum>> mapping, Class cl, String name) {
        Map<Long, Set<ReferenceOrderedDatum>> ret = new LinkedHashMap<Long, Set<ReferenceOrderedDatum>>();
        int start = record.getAlignmentStart();
        int stop = record.getAlignmentEnd();
        for (Long location : mapping.keySet()) {
            if (location >= start && location <= stop) ret.put(location,mapping.get(location));
        }
        return ret;
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of read offset to ROD(s)
     */
    public Map<Integer, Set<ReferenceOrderedDatum>> getPositionMapping() {
        return createReadAlignment(record, mapping, null, null);
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of genome loc position to ROD(s)
     */
    public Map<Long, Set<ReferenceOrderedDatum>> getGenomeLocMapping() {
        return createGenomeLocAlignment(record, mapping, null, null);
    }
}
