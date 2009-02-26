package edu.mit.broad.arachne;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.BitSet;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Utility class to read in a set of contig-based genomic intervals in zero-based end inclusive
 * and store them efficiently in memory as a 1-based bit-mask
 */
public class GenomeMask {

    // if memory usage becomes a problem... this could be changed to a SparseBitSet
    // http://java.sun.com/developer/onlineTraining/collections/magercises/BitSet/index.html
    private SortedMap<Integer, BitSet> data = new TreeMap<Integer, BitSet>();


    public GenomeMask(File maskFile) throws IOException {
        BufferedReader baitReader = null;
        try {
            baitReader = new BufferedReader(new FileReader(maskFile));
            String line;
            while ((line = baitReader.readLine()) != null) {
                String[] arr = line.split(" ");
                int contig = Integer.parseInt(arr[0]);

                // covert the coordinates from 0-based, end inclusive to
                // 1-based end inclusive
                int startPos = Integer.parseInt(arr[1]) + 1;
                int endPos = Integer.parseInt(arr[2]) + 1;

                BitSet bits = data.get(contig);
                if (bits == null) {
                    bits = new BitSet(endPos);
                    data.put(contig,bits);
                }

                bits.set(startPos, endPos + 1); // set method is end exclusive
            }
        } finally {
            if (baitReader != null) { baitReader.close(); }
        }
    }

    /**
     * This ctor is useful if initializing a GenomeMask externally.
     */
    public GenomeMask() {
    }

    public boolean get(int contig, int position) {
        BitSet bits = data.get(contig);
        return (bits != null) && bits.get(position);
    }

    public BitSet get(int contig) {
        return data.get(contig);
    }

    /**
     * Get an existing BitSet for the given contig, or create one if not already present.  This is
     * useful when initializing a GenomeMask from an external source.
     * @param contig which BitSet
     * @param numBits if there was not already a BitSet for this contig, one is created and initialized to this size.
     * @return the BitSet for the given contig, creating one if necessary
     */
    public BitSet getOrCreate(int contig, int numBits) {
        BitSet ret = data.get(contig);
        if (ret == null) {
            ret = new BitSet(numBits);
            data.put(contig, ret);
        }
        return ret;
    }

    public int getMaxContig() {
        return data.lastKey();
    }
}
