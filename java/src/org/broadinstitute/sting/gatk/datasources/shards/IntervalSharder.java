package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import net.sf.samtools.Bin;

/**
 * Shard intervals based on position within the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalSharder {
    protected static List<FilePointer> shardIntervals(final BlockDrivenSAMDataSource dataSource, final List<GenomeLoc> loci, final int binsDeeperThan) {
        // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.
        List<FilePointer> filePointers = new ArrayList<FilePointer>();
        FilePointer filePointer = null;

        for(GenomeLoc location: loci) {
            int locationStart = (int)location.getStart();
            final int locationStop = (int)location.getStop();

            List<Bin> bins = findBinsAtLeastAsDeepAs(dataSource,dataSource.getOverlappingBins(location),binsDeeperThan);

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

                filePointers.addAll(shardIntervals(dataSource,Collections.singletonList(location),binsDeeperThan-1));
                continue;
            }

            // Bins found; try to match bins with locations.
            Collections.sort(bins);
            Iterator<Bin> binIterator = bins.iterator();

            while(locationStop >= locationStart) {
                int binStart = filePointer!=null ? dataSource.getFirstLocusInBin(filePointer.bin) : 0;
                int binStop = filePointer!=null ? dataSource.getLastLocusInBin(filePointer.bin) : 0;

                while(binStop < locationStart && binIterator.hasNext()) {
                    if(filePointer != null && filePointer.locations.size() > 0)
                        filePointers.add(filePointer);

                    filePointer = new FilePointer(binIterator.next());
                    binStart = dataSource.getFirstLocusInBin(filePointer.bin);
                    binStop = dataSource.getLastLocusInBin(filePointer.bin);
                }

                if(locationStart < binStart) {
                    // The region starts before the first bin in the sequence.  Add the region occurring before the sequence.
                    if(filePointer != null && filePointer.locations.size() > 0) {
                        filePointers.add(filePointer);
                        filePointer = null;
                    }

                    final int regionStop = Math.min(locationStop,binStart-1);

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,regionStop);
                    filePointers.addAll(shardIntervals(dataSource,Collections.singletonList(subset),binsDeeperThan-1));

                    locationStart = regionStop + 1;
                }
                else if(locationStart > binStop) {
                    // The region starts after the last bin in the sequence.  Add the region occurring after the sequence.
                    if(filePointer != null && filePointer.locations.size() > 0) {
                        filePointers.add(filePointer);
                        filePointer = null;
                    }

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,locationStop);
                    filePointers.addAll(shardIntervals(dataSource,Collections.singletonList(subset),binsDeeperThan-1));

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

    private static List<Bin> findBinsAtLeastAsDeepAs(final BlockDrivenSAMDataSource dataSource, final List<Bin> bins, final int deepestBinLevel) {
        List<Bin> deepestBins = new ArrayList<Bin>();
        for(Bin bin: bins) {
            if(dataSource.getLevelForBin(bin) >= deepestBinLevel)
                deepestBins.add(bin);
        }
        return deepestBins;
    }
}

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
class FilePointer {
    protected final Bin bin;
    protected final List<GenomeLoc> locations;

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


