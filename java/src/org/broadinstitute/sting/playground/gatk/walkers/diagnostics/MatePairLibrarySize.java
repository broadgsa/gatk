package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;

/**
 * For each sequencing library, outputs the distribution of mate pair sizes
 */
public class MatePairLibrarySize extends ReadWalker<Integer, Integer> {
    @Argument(fullName="outdir", shortName="outdir", doc="Directory to output results")
    private File OUT_DIR;

    private HashMap<String, HashMap<Integer, Integer>> matePairSize;

    public void initialize() {
        matePairSize = new HashMap<String, HashMap<Integer, Integer>>();

        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            HashMap<Integer, Integer> mps = new HashMap<Integer, Integer>();

            matePairSize.put(rg.getLibrary(), mps);
        }
    }

    public boolean filter(ReferenceContext ref, SAMRecord read) {
        return (read.getReadPairedFlag() && read.getFirstOfPairFlag());
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        int insert = read.getInferredInsertSize();

        Integer oldcount = matePairSize.get(read.getReadGroup().getLibrary()).get(insert);
        if (oldcount == null) { oldcount = 0; }

        matePairSize.get(read.getReadGroup().getLibrary()).put(insert, oldcount + 1);

        return null;
    }

    public Integer reduceInit() {
        return null;
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        String[] libraries = matePairSize.keySet().toArray(new String[1]);

        for (String library : libraries) {
            try {
                Integer[] sizes = matePairSize.get(library).keySet().toArray(new Integer[1]);

                if (sizes != null && sizes.length > 1) {
                    PrintWriter pw = new PrintWriter(String.format("%s/%s.pairdist", OUT_DIR.getAbsolutePath(), library));
                    Arrays.sort(sizes);

                    pw.printf("%s\t%s%n", "insert", "frequency");

                    for (int insert : sizes) {
                        if (insert >= 0) {
                            pw.printf("%d\t%d%n", insert, matePairSize.get(library).get(insert));
                        }
                    }

                    pw.close();
                }
            } catch (IOException e) {
                throw new StingException("Unable to initialize output files.");
            }
        }
    }
}
