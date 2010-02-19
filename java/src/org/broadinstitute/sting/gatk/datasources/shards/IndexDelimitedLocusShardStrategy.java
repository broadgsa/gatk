package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;

import java.util.*;

import net.sf.samtools.Chunk;
import net.sf.samtools.Bin;
import net.sf.samtools.SAMFileReader2;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */


/**
 * A sharding strategy for loci based on reading of the index.
 */
public class IndexDelimitedLocusShardStrategy implements ShardStrategy {
    /**
     * The data source to use when performing this sharding.
     */
    private final BlockDrivenSAMDataSource blockDrivenDataSource;

    /** our storage of the genomic locations they'd like to shard over */
    private final List<FilePointer> filePointers = new ArrayList<FilePointer>();

    /**
     * An iterator through the available file pointers.
     */
    private final Iterator<FilePointer> filePointerIterator;

    /**
     * construct the shard strategy from a seq dictionary, a shard size, and and genomeLocs
     * @param dataSource Data source from which to load index data.
     * @param locations List of locations for which to load data.
     */
    IndexDelimitedLocusShardStrategy(SAMDataSource dataSource, GenomeLocSortedSet locations) {
        if(!(dataSource instanceof BlockDrivenSAMDataSource))
            throw new StingException("Cannot power an IndexDelimitedLocusShardStrategy with this data source.");

        blockDrivenDataSource = (BlockDrivenSAMDataSource)dataSource;
        final int deepestBinLevel = blockDrivenDataSource.getNumIndexLevels()-1;

        // Create a list of contig name -> genome loc, sorted in INSERTION ORDER.
        LinkedHashMap<String,List<GenomeLoc>> locationToReference = new LinkedHashMap<String,List<GenomeLoc>>();
        for(GenomeLoc location: locations) {
            if(!locationToReference.containsKey(location.getContig()))
                locationToReference.put(location.getContig(),new ArrayList<GenomeLoc>());
            locationToReference.get(location.getContig()).add(location);
        }

        for(String contig: locationToReference.keySet()) {
            // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.
            SortedMap<Bin,List<GenomeLoc>> bins = new TreeMap<Bin,List<GenomeLoc>>();
            for(GenomeLoc location: locationToReference.get(contig)) {
                List<Bin> binsForLocation = blockDrivenDataSource.getOverlappingBins(location);
                for(Bin bin: binsForLocation) {
                    if(blockDrivenDataSource.getLevelForBin(bin) == deepestBinLevel) {
                        final int firstLoc = blockDrivenDataSource.getFirstLocusInBin(bin);
                        final int lastLoc = blockDrivenDataSource.getLastLocusInBin(bin);
                        if(!bins.containsKey(bin))
                            bins.put(bin,new ArrayList<GenomeLoc>());
                        bins.get(bin).add(GenomeLocParser.createGenomeLoc(location.getContig(),
                                                                          Math.max(location.getStart(),firstLoc),
                                                                          Math.min(location.getStop(),lastLoc)));
                    }
                }
            }

            // Add a record of the new bin structure.
            for(SortedMap.Entry<Bin,List<GenomeLoc>> entry: bins.entrySet()) {
                Collections.sort(entry.getValue());                
                filePointers.add(new FilePointer(entry.getKey(),entry.getValue()));
            }
        }

        filePointerIterator = filePointers.iterator();
    }

    /**
     * returns true if there are additional shards
     *
     * @return false if we're done processing shards
     */
    public boolean hasNext() {
        return filePointerIterator.hasNext();
    }

    /**
     * gets the next Shard
     *
     * @return the next shard
     */
    public IndexDelimitedLocusShard next() {
        FilePointer nextFilePointer = filePointerIterator.next();
        Map<SAMFileReader2,List<Chunk>> chunksBounding = blockDrivenDataSource.getFilePointersBounding(nextFilePointer.bin);
        return new IndexDelimitedLocusShard(nextFilePointer.locations,chunksBounding,Shard.ShardType.LOCUS_INTERVAL);
    }

    /** we don't support the remove command */
    public void remove() {
        throw new UnsupportedOperationException("ShardStrategies don't support remove()");
    }

    /**
     * makes the IntervalShard iterable, i.e. usable in a for loop.
     *
     * @return
     */
    public Iterator<Shard> iterator() {
        return this;
    }

    /**
     * Represents a small section of a BAM file, and every associated interval.
     */
    private class FilePointer {
        private final Bin bin;
        private final List<GenomeLoc> locations;

        public FilePointer(Bin bin, List<GenomeLoc> locations) {
            this.bin = bin;
            this.locations = locations;
        }
        
    }

}