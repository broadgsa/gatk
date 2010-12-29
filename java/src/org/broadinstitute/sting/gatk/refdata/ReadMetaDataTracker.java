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
import org.broadinstitute.sting.gatk.datasources.providers.RODMetaDataContainer;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
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
    /**
     * The parser, used to create new GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;

    private final SAMRecord record;

    // the buffer of positions and RODs we've stored
    private final TreeMap<Integer, RODMetaDataContainer> mapping;

    /**
     * create a read meta data tracker, given the read and a queue of RODatum positions
     *
     * @param record  the read to create offset from
     * @param mapping the mapping of reference ordered datum
     */
    public ReadMetaDataTracker(GenomeLocParser genomeLocParser, SAMRecord record, TreeMap<Integer, RODMetaDataContainer> mapping) {
        this.genomeLocParser = genomeLocParser;
        this.record = record;
        this.mapping = mapping;
    }

    /**
     * create an alignment of read position to reference ordered datum
     *
     * @param record the SAMRecord
     * @param queue  the queue (as a tree set)
     * @param cl     the class name, null if not filtered by classname
     * @param name   the datum track name, null if not filtered by name
     *
     * @return a mapping from the position in the read to the reference ordered datum
     */
    private Map<Integer, Collection<GATKFeature>> createReadAlignment(SAMRecord record, TreeMap<Integer, RODMetaDataContainer> queue, Class cl, String name) {
        if (name != null && cl != null) throw new IllegalStateException("Both a class and name cannot be specified");
        Map<Integer, Collection<GATKFeature>> ret = new LinkedHashMap<Integer, Collection<GATKFeature>>();
        GenomeLoc location = genomeLocParser.createGenomeLoc(record);
        int length = record.getReadLength();
        for (Integer loc : queue.keySet()) {
            Integer position = loc - location.getStart();
            if (position >= 0 && position < length) {
                Collection<GATKFeature> set;
                if (cl != null)
                    set = queue.get(loc).getSet(cl);
                else
                    set = queue.get(loc).getSet(name);
                if (set != null && set.size() > 0)
                    ret.put(position, set);
            }
        }
        return ret;

    }

    /**
     * create an alignment of read position to reference ordered datum
     *
     * @return a mapping from the position in the read to the reference ordered datum
     */
    private Map<Integer, Collection<GATKFeature>> createGenomeLocAlignment(SAMRecord record, TreeMap<Integer, RODMetaDataContainer> mapping, Class cl, String name) {
        Map<Integer, Collection<GATKFeature>> ret = new LinkedHashMap<Integer, Collection<GATKFeature>>();
        int start = record.getAlignmentStart();
        int stop = record.getAlignmentEnd();
        for (Integer location : mapping.keySet()) {
            if (location >= start && location <= stop)
                if (cl != null)
                    ret.put(location, mapping.get(location).getSet(cl));
                else
                    ret.put(location, mapping.get(location).getSet(name));
        }
        return ret;
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of read offset to ROD(s)
     */
    public Map<Integer, Collection<GATKFeature>> getReadOffsetMapping() {
        return createReadAlignment(record, mapping, null, null);
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of genome loc position to ROD(s)
     */
    public Map<Integer, Collection<GATKFeature>> getContigOffsetMapping() {
        return createGenomeLocAlignment(record, mapping, null, null);
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of read offset to ROD(s)
     */
    public Map<Integer, Collection<GATKFeature>> getReadOffsetMapping(String name) {
        return createReadAlignment(record, mapping, null, name);
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of genome loc position to ROD(s)
     */
    public Map<Integer, Collection<GATKFeature>> getContigOffsetMapping(String name) {
        return createGenomeLocAlignment(record, mapping, null, name);
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of read offset to ROD(s)
     */
    public Map<Integer, Collection<GATKFeature>> getReadOffsetMapping(Class cl) {
        return createReadAlignment(record, mapping, cl, null);
    }

    /**
     * get the position mapping, from read offset to ROD
     *
     * @return a mapping of genome loc position to ROD(s)
     */
    public Map<Integer, Collection<GATKFeature>> getContigOffsetMapping(Class cl) {
        return createGenomeLocAlignment(record, mapping, cl, null);
    }

    /**
     * get the list of all the RODS overlapping this read, without any information about their position
     * @return a Collection (no order guaranteed), of all the RODs covering this read
     */
    public List<GATKFeature> getAllCoveringRods() {
        List<GATKFeature> ret = new ArrayList<GATKFeature>();
        for (Map.Entry<Integer, RODMetaDataContainer> entry : mapping.entrySet())
            ret.addAll(entry.getValue().getSet());
        return ret;
    }
}
