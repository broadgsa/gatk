package org.broadinstitute.sting.oneoffprojects.tools;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;

import java.io.*;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

/**
 * Used to test how long it takes to read through a text files and gzipped files.
 * If the passed-in filename ends with .gz, it will be read using GZIPInputStream.
 * Otherwise, its read using FileReader.
 */
public class CompareBAMAlignments extends CommandLineProgram {

    @Argument(fullName = "input", shortName = "i", doc = "xxx", required = true)
    public List<String> filenames;

    @Argument(fullName = "maxIsize", shortName = "s", doc = "xxx", required=false)
    public int maxISize = -1;

    @Argument(fullName = "incr", shortName = "incr", doc = "xxx", required=false)
    int incr = -1;

    @Override
    protected int execute() {
        try {
            List<Iterator<SAMRecord>> readers = new ArrayList<Iterator<SAMRecord>>();
            for ( String filename : filenames ) {
                final File file = new File(filename);
                SAMFileReader reader = new SAMFileReader(file);
                reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
                readers.add(reader.iterator());
            }

            System.out.println("Reading...");
            int next = incr;
            int counter = 0;

            while( true ) {
                List<SAMRecord> reads = new ArrayList<SAMRecord>();
                for ( Iterator<SAMRecord> reader : readers ) {
                    if ( ! reader.hasNext() ) System.exit(0);
                    reads.add(reader.next());
                }

                // comparing
                SAMRecord read1 = reads.get(0);
                if ( read1.getInferredInsertSize() > maxISize ) {
                    for ( SAMRecord read : reads ) {
                        if(incr > 0 && counter % incr == 0) {
                            next += incr;
                            System.err.println(new Date() + " - counter " + counter);
                            System.err.println("read: " + read.format());
                        }

                        if ( ! read1.getReadName().equals(read.getReadName()) )
                            bad(read1, read, "Names not equal");
                        else {
                            if ( read1.getAlignmentStart() != read.getAlignmentStart() )
                                bad(read1, read, "Alignment starts not equal");
                            if ( ! read1.getCigarString().equals(read.getCigarString()) )
                                bad(read1, read, "Unequal CIGAR strings");
                        }
                    }
                }
                counter++;
            }
        } catch(Exception e) {
            System.err.println("ERROR: " + e);
            e.printStackTrace();
        }

        return 0;
    }

    private void bad(SAMRecord read1, SAMRecord read2, String msg) {
        System.out.printf("%nBAD: %s%n", msg);
        System.out.printf("  read1:  %s %s %s %s%n", read1.getReadName(), read1.getAlignmentStart(), read1.getCigarString(), read1.getInferredInsertSize());
        System.out.printf("  read2:  %s %s %s %s%n", read2.getReadName(), read2.getAlignmentStart(), read2.getCigarString(), read2.getInferredInsertSize());
        // System.exit(1);
    }

    public static void main(String args[])
    {
        try {
            CommandLineProgram.start(new CompareBAMAlignments(), args);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

}
