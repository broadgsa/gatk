package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

import java.util.*;

import net.sf.samtools.*;
import net.sf.picard.util.PeekableIterator;

/**
 * Shard intervals based on position within the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalSharder {
    protected static List<FilePointer> shardIntervals(final BlockDrivenSAMDataSource dataSource, final List<GenomeLoc> loci) {
        Map<SAMReaderID,List<FilePointer>> filePointersByReader = new HashMap<SAMReaderID,List<FilePointer>>();
        for(SAMReaderID id: dataSource.getReaderIDs()) {
            PreloadedBAMFileIndex index = dataSource.getIndex(id);
            // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.\
            filePointersByReader.put(id,shardIntervalsOverIndex(dataSource,id,index,loci,index.getNumIndexLevels()-1));
            index.close();
        }
        return combineFilePointers(filePointersByReader);
    }

    /**
     * Combine adjacent file pointers into a structure that can be streamed in.
     * @param filePointersByReader File pointers broken down by reader.
     * @return A large structure of file pointers.
     */
    private static List<FilePointer> combineFilePointers(Map<SAMReaderID,List<FilePointer>> filePointersByReader) {
        PeekableIterator<FilePointer> mergingIterator = new PeekableIterator<FilePointer>(new FilePointerMergingIterator(filePointersByReader));

        List<FilePointer> overlappingFilePointers = new ArrayList<FilePointer>();
        List<FilePointer> mergedFilePointers = new ArrayList<FilePointer>();

        while(mergingIterator.hasNext()) {
            GenomeLoc bounds = null;

            // Load up a segment where file pointers overlap
            while(mergingIterator.hasNext() && (overlappingFilePointers.size() == 0 || mergingIterator.peek().getBounds().overlapsP(bounds))) {
                FilePointer filePointer = mergingIterator.next();
                if(bounds != null)
                    bounds = GenomeLocParser.createGenomeLoc(bounds.getContig(),
                            Math.min(bounds.getStart(),filePointer.getBounds().getStart()),
                            Math.max(bounds.getStop(),filePointer.getBounds().getStop()));
                else
                    bounds = filePointer.getBounds();
                overlappingFilePointers.add(filePointer);
            }

            // determine the complete set of unique locations defining this set.
            List<GenomeLoc> overlappingLocations = new ArrayList<GenomeLoc>();
            for(FilePointer filePointer: overlappingFilePointers)
                overlappingLocations.addAll(filePointer.locations);
            Collections.sort(overlappingLocations);
            overlappingLocations = GenomeLocSortedSet.mergeOverlappingLocations(overlappingLocations);

            while(!overlappingLocations.isEmpty()) {
                long overlapStart = overlappingLocations.get(0).getStart();
                long overlapStop = overlappingLocations.get(overlappingLocations.size()-1).getStop();

                for(FilePointer overlappingFilePointer: overlappingFilePointers) {
                    if(overlappingFilePointer.getBounds().getStop() < overlapStart)
                        continue;
                    if(overlappingFilePointer.getBounds().getStart() > overlapStart) overlapStop = Math.min(overlapStop,overlappingFilePointer.getBounds().getStart()-1);
                    if(overlappingFilePointer.getBounds().getStop() < overlapStop) overlapStop = Math.min(overlapStop,overlappingFilePointer.getBounds().getStop());
                }

                // Find the overlapping genome locs.
                List<GenomeLoc> segmentOverlap = new ArrayList<GenomeLoc>();
                for(GenomeLoc overlappingLocation: overlappingLocations) {
                    if(overlappingLocation.getStop() <= overlapStop) {
                        // segment is completely before end of overlap.
                        segmentOverlap.add(overlappingLocation);
                    }
                    else if(overlappingLocation.getStart() <= overlapStop) {
                        // segment is partially before end of overlap.
                        segmentOverlap.add(GenomeLocParser.setStop(overlappingLocation,overlapStop));
                        break;
                    }
                    else {
                        // segment starts after overlap ends.
                        break;
                    }
                }

                // Trim the overlapping genome locs of the overlapping locations list.
                while(!overlappingLocations.isEmpty() && overlappingLocations.get(0).getStart() <= overlapStop) {
                    GenomeLoc location = overlappingLocations.remove(0);
                    if(location.getStop() > overlapStop)
                        overlappingLocations.add(0,GenomeLocParser.setStart(location,overlapStop+1));
                }

                // Merge together all file pointers that overlap with these bounds.
                GenomeLoc overlapBounds = GenomeLocParser.createGenomeLoc(segmentOverlap.get(0).getContigIndex(),overlapStart,overlapStop);
                FilePointer mergedFilePointer = null;
                for(FilePointer overlappingFilePointer: overlappingFilePointers) {
                    if(overlappingFilePointer.getBounds().overlapsP(overlapBounds))
                        mergedFilePointer = overlappingFilePointer.merge(mergedFilePointer,segmentOverlap);
                }

                // Add the resulting file pointer and clear state.
                mergedFilePointers.add(mergedFilePointer);
            }

            // reset
            overlappingFilePointers.clear();
        }

        return mergedFilePointers;
    }

    private static class FilePointerMergingIterator implements Iterator<FilePointer> {
        private PriorityQueue<PeekableIterator<FilePointer>> filePointerQueue;

        public FilePointerMergingIterator(Map<SAMReaderID,List<FilePointer>> filePointers) {
            filePointerQueue = new PriorityQueue<PeekableIterator<FilePointer>>(filePointers.size(),new FilePointerMergingComparator());
            for(List<FilePointer> filePointersByReader: filePointers.values())
                filePointerQueue.add(new PeekableIterator<FilePointer>(filePointersByReader.iterator()));
        }

        public boolean hasNext() {
            return !filePointerQueue.isEmpty();
        }

        public FilePointer next() {
            if(!hasNext()) throw new NoSuchElementException("FilePointerMergingIterator is out of elements");
            PeekableIterator<FilePointer> nextIterator = filePointerQueue.remove();
            FilePointer nextFilePointer = nextIterator.next();
            if(nextIterator.hasNext())
                filePointerQueue.add(nextIterator);
            return nextFilePointer;
        }

        public void remove() { throw new UnsupportedOperationException("Cannot remove from a merging iterator."); }

        private class FilePointerMergingComparator implements Comparator<PeekableIterator<FilePointer>> {
            public int compare(PeekableIterator<FilePointer> lhs, PeekableIterator<FilePointer> rhs) {
                if(!lhs.hasNext() && !rhs.hasNext()) return 0;
                if(!rhs.hasNext()) return -1;
                if(!lhs.hasNext()) return 1;
                return lhs.peek().getBounds().compareTo(rhs.peek().getBounds());
            }
        }
    }

    private static List<FilePointer> shardIntervalsOverIndex(final BlockDrivenSAMDataSource dataSource, final SAMReaderID id, final PreloadedBAMFileIndex index, final List<GenomeLoc> loci, final int binsDeeperThan) {
        // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.
        List<FilePointer> filePointers = new ArrayList<FilePointer>();
        FilePointer lastFilePointer = null;
        Bin lastBin = null;

        for(GenomeLoc location: loci) {
            // If crossing contigs, be sure to reset the filepointer that's been accumulating shard data.
            if(lastFilePointer != null && lastFilePointer.referenceSequence != location.getContigIndex()) {
                filePointers.add(lastFilePointer);
                lastFilePointer = null;
                lastBin = null;
            }

            int locationStart = (int)location.getStart();
            final int locationStop = (int)location.getStop();

            List<Bin> bins = findBinsAtLeastAsDeepAs(index,getOverlappingBins(dataSource,id,index,location),binsDeeperThan);

            // Recursive stopping condition -- algorithm is at the zero point and no bins have been found.
            if(binsDeeperThan == 0 && bins.size() == 0) {
                filePointers.add(new FilePointer(location));
                continue;
            }

            // No bins found; step up a level and search again.
            if(bins.size() == 0) {
                if(lastFilePointer != null && lastFilePointer.locations.size() > 0) {
                    filePointers.add(lastFilePointer);
                    lastFilePointer = null;
                    lastBin = null;
                }

                filePointers.addAll(shardIntervalsOverIndex(dataSource,id,index,Collections.singletonList(location),binsDeeperThan-1));
                continue;
            }

            // Bins found; try to match bins with locations.
            Collections.sort(bins);
            Iterator<Bin> binIterator = bins.iterator();

            while(locationStop >= locationStart) {
                int binStart = lastFilePointer!=null ? index.getFirstLocusInBin(lastBin) : 0;
                int binStop = lastFilePointer!=null ? index.getLastLocusInBin(lastBin) : 0;

                while(binStop < locationStart && binIterator.hasNext()) {
                    if(lastFilePointer != null && lastFilePointer.locations.size() > 0)
                        filePointers.add(lastFilePointer);

                    lastBin = binIterator.next();
                    lastFilePointer = new FilePointer(id,lastBin.referenceSequence,getFilePointersBounding(index,lastBin));
                    binStart = index.getFirstLocusInBin(lastBin);
                    binStop = index.getLastLocusInBin(lastBin);
                }

                if(locationStart < binStart) {
                    // The region starts before the first bin in the sequence.  Add the region occurring before the sequence.
                    if(lastFilePointer != null && lastFilePointer.locations.size() > 0) {
                        filePointers.add(lastFilePointer);
                        lastFilePointer = null;
                        lastBin = null;
                    }

                    final int regionStop = Math.min(locationStop,binStart-1);

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,regionStop);
                    filePointers.addAll(shardIntervalsOverIndex(dataSource,id,index,Collections.singletonList(subset),binsDeeperThan-1));

                    locationStart = regionStop + 1;
                }
                else if(locationStart > binStop) {
                    // The region starts after the last bin in the sequence.  Add the region occurring after the sequence.
                    if(lastFilePointer != null && lastFilePointer.locations.size() > 0) {
                        filePointers.add(lastFilePointer);
                        lastFilePointer = null;
                        lastBin = null;
                    }

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,locationStop);
                    filePointers.addAll(shardIntervalsOverIndex(dataSource,id,index,Collections.singletonList(subset),binsDeeperThan-1));

                    locationStart = locationStop + 1;
                }
                else {
                    // The start of the region overlaps the bin.  Add the overlapping subset.
                    final int regionStop = Math.min(locationStop,binStop);
                    lastFilePointer.addLocation(GenomeLocParser.createGenomeLoc(location.getContig(),
                            locationStart,
                            regionStop));
                    locationStart = regionStop + 1;
                }
            }
        }

        if(lastFilePointer != null && lastFilePointer.locations.size() > 0)
            filePointers.add(lastFilePointer);

        return filePointers;
    }

    private static List<Bin> findBinsAtLeastAsDeepAs(final PreloadedBAMFileIndex index, final List<Bin> bins, final int deepestBinLevel) {
        List<Bin> deepestBins = new ArrayList<Bin>();
        for(Bin bin: bins) {
            if(index.getLevelForBin(bin) >= deepestBinLevel)
                deepestBins.add(bin);
        }
        return deepestBins;
    }

    /**
     * Gets a list of the bins in each BAM file that overlap with the given interval list.
     * @param location Location for which to determine the bin.
     * @return A map of reader back to bin.
     */
    private static List<Bin> getOverlappingBins(final BlockDrivenSAMDataSource dataSource, final SAMReaderID id, final PreloadedBAMFileIndex index, final GenomeLoc location) {
        // All readers will have the same bin structure, so just use the first bin as an example.
        final SAMFileHeader fileHeader = dataSource.getHeader(id);
        int referenceIndex = fileHeader.getSequenceIndex(location.getContig());
        if (referenceIndex != -1) {
            return index.getBinsContaining(referenceIndex,(int)location.getStart(),(int)location.getStop());
        }
        return Collections.emptyList();
    }

    /**
     * Gets the file pointers bounded by this bin, grouped by the reader of origination.
     * @param bin The bin for which to load data.
     * @return A map of the file pointers bounding the bin.
     */
    private static List<Chunk> getFilePointersBounding(final PreloadedBAMFileIndex index, final Bin bin) {
        if(bin != null) {
            List<Chunk> chunks = index.getSearchBins(bin);
            return chunks != null ? chunks : Collections.<Chunk>emptyList();
        }
        else
            return Collections.emptyList();
    }


}

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
class FilePointer {
    protected final Map<SAMReaderID,List<Chunk>> chunks = new HashMap<SAMReaderID,List<Chunk>>();
    protected final int referenceSequence;
    protected final List<GenomeLoc> locations;

    public FilePointer(SAMReaderID id, int referenceSequence, List<Chunk> chunks) {
        this.referenceSequence = referenceSequence;
        this.chunks.put(id,chunks);
        this.locations = new ArrayList<GenomeLoc>();
    }

    public FilePointer(GenomeLoc location) {
        referenceSequence = location.getContigIndex();
        locations = Collections.singletonList(location);
    }

    /**
     * Private constructor for merge operation.
     * @param referenceSequence Sequence to merge.
     * @param locations Merged locations.
     */
    private FilePointer(final int referenceSequence, final List<GenomeLoc> locations) {
        this.referenceSequence = referenceSequence;
        this.locations = locations;
    }

    public FilePointer merge(FilePointer other, List<GenomeLoc> locations) {
        FilePointer merged = new FilePointer(referenceSequence,locations);
        merged.chunks.putAll(this.chunks);
        if(other != null)
            merged.chunks.putAll(other.chunks);
        return merged;
    }

    public void addLocation(GenomeLoc location) {
        locations.add(location);
    }

    public GenomeLoc getBounds() {
        final long boundaryStart = locations.get(0).getStart();
        final long boundaryStop = locations.get(locations.size()-1).getStop();
        return GenomeLocParser.createGenomeLoc(locations.get(0).getContigIndex(),boundaryStart,boundaryStop);    
    }
}


