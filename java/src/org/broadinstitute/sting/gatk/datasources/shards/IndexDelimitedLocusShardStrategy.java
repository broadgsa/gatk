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

        // Create a list of contig name -> genome loc, sorted in INSERTION ORDER.
        LinkedHashMap<String,List<GenomeLoc>> locationToReference = new LinkedHashMap<String,List<GenomeLoc>>();
        for(GenomeLoc location: locations) {
            if(!locationToReference.containsKey(location.getContig()))
                locationToReference.put(location.getContig(),new ArrayList<GenomeLoc>());
            locationToReference.get(location.getContig()).add(location);
        }

        // TODO: Not sure there's any reason to pre-separate the contigs now that we're using a streaming approach to file pointer allocation.
        for(String contig: locationToReference.keySet()) {
            filePointers.addAll(batchLociIntoBins(locationToReference.get(contig),blockDrivenDataSource.getNumIndexLevels()-1));
        }

        filePointerIterator = filePointers.iterator();
    }

    private List<FilePointer> batchLociIntoBins(final List<GenomeLoc> loci, final int binsDeeperThan) {
        // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.
        List<FilePointer> filePointers = new ArrayList<FilePointer>();
        FilePointer filePointer = null;

        for(GenomeLoc location: loci) {
            int locationStart = (int)location.getStart();
            final int locationStop = (int)location.getStop();

            List<Bin> bins = findBinsAtLeastAsDeepAs(blockDrivenDataSource.getOverlappingBins(location),binsDeeperThan);

            // Recursive stopping condition -- algorithm is at the zero point and no bins have been found.
            if(binsDeeperThan == 0 && bins.size() == 0) {
                filePointers.add(new FilePointer(location));
                continue;
            }

            // No bins found; step up a level and search again.
            if(bins.size() == 0) {
                if(filePointer != null && filePointer.locations.size() > 0) {
                    filePointers.add(filePointer);
                    filePointer = null;
                }                

                filePointers.addAll(batchLociIntoBins(Collections.singletonList(location),binsDeeperThan-1));
                continue;
            }

            // Bins found; try to match bins with locations.
            Collections.sort(bins);
            Iterator<Bin> binIterator = bins.iterator();

            while(locationStop >= locationStart) {
                int binStart = filePointer!=null ? blockDrivenDataSource.getFirstLocusInBin(filePointer.bin) : 0;
                int binStop = filePointer!=null ? blockDrivenDataSource.getLastLocusInBin(filePointer.bin) : 0;

                while(binStop < locationStart && binIterator.hasNext()) {
                    if(filePointer != null && filePointer.locations.size() > 0)
                        filePointers.add(filePointer);

                    filePointer = new FilePointer(binIterator.next());
                    binStart = blockDrivenDataSource.getFirstLocusInBin(filePointer.bin);
                    binStop = blockDrivenDataSource.getLastLocusInBin(filePointer.bin);
                }

                if(locationStart < binStart) {
                    // The region starts before the first bin in the sequence.  Add the region occurring before the sequence.
                    if(filePointer != null && filePointer.locations.size() > 0) {
                        filePointers.add(filePointer);
                        filePointer = null;
                    }

                    final int regionStop = Math.min(locationStop,binStart-1);

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,regionStop);
                    filePointers.addAll(batchLociIntoBins(Collections.singletonList(subset),binsDeeperThan-1));

                    locationStart = regionStop + 1;
                }
                else if(locationStart > binStop) {
                    // The region starts after the last bin in the sequence.  Add the region occurring after the sequence.
                    if(filePointer != null && filePointer.locations.size() > 0) {
                        filePointers.add(filePointer);
                        filePointer = null;
                    }

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,locationStop);
                    filePointers.addAll(batchLociIntoBins(Collections.singletonList(subset),binsDeeperThan-1));

                    locationStart = locationStop + 1;
                }
                else {
                    // The start of the region overlaps the bin.  Add the overlapping subset.
                    final int regionStop = Math.min(locationStop,binStop);
                    filePointer.addLocation(GenomeLocParser.createGenomeLoc(location.getContig(),
                            locationStart,
                            regionStop));
                    locationStart = regionStop + 1;
                }
            }
        }

        if(filePointer != null && filePointer.locations.size() > 0)
            filePointers.add(filePointer);

        return filePointers;
    }

    private List<Bin> findBinsAtLeastAsDeepAs(final List<Bin> bins, final int deepestBinLevel) {
        List<Bin> deepestBins = new ArrayList<Bin>();
        for(Bin bin: bins) {
            if(blockDrivenDataSource.getLevelForBin(bin) >= deepestBinLevel)
                deepestBins.add(bin);
        }
        return deepestBins;
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

        public FilePointer(Bin bin) {
            this.bin = bin;
            this.locations = new ArrayList<GenomeLoc>();
        }

        public FilePointer(GenomeLoc location) {
            bin = null;
            locations = Collections.singletonList(location);
        }

        public void addLocation(GenomeLoc location) {
            locations.add(location);
        }
        
    }

}