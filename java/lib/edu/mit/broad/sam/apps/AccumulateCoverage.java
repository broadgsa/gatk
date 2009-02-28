/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.apps;

import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMLocusIterator;
import edu.mit.broad.sam.SAMFileHeader;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.io.FileWriter;
import java.util.List;

public class AccumulateCoverage {

    public static void main(final String[] argv) throws Exception {
        if (argv.length != 1) {
            System.err.println("ERROR: Incorrect number of arguments");
            usage();
            System.exit(1);
        }
        final AccumulateCoverage ac = new AccumulateCoverage(argv[0]);
    }

    private static void usage() {
        System.err.println("USAGE: AccumulateCoverage <SAMFile>");
    }



    public AccumulateCoverage(final String samFile) throws IOException {
        final long startTime = System.currentTimeMillis();
        final Writer writer = new FileWriter("/Users/kcibul/projects/sam/acccov.out");

        final SAMFileReader samReader = new SAMFileReader(new File(samFile));

        // ensure the file is sorted
//TODO: is the SAM reader implementation broken?
        if (samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            System.out.println("SAM Files must be coordinate-sorted, this is " + samReader.getFileHeader().getSortOrder());
            System.exit(1);
        }

        final SAMLocusIterator sli = new SAMLocusIterator(samReader.iterator());

        for (final SAMLocusIterator.LocusInfo li : sli) {

        String chrom = li.getChrom().substring(3);
        if (chrom.equals("M")) { chrom = "0"; }
        if (chrom.equals("X")) { chrom = "23"; }
        if (chrom.equals("Y")) { chrom = "24"; }

        final StringBuilder sb = new StringBuilder();
        sb.append(chrom)
                .append(":")
                .append(li.getPosition()-1)
                .append(" ")
                .append(li.getBases().size())
                .append("\n");

        writer.write(sb.toString());
        //System.out.print(sb);

//        // TODO: zero based or 1 based?
//        System.out.print(li.chrom + "\t" + (li.position-1) + "\t" + li.bases.size() + "\t");
//
//        // TODO: print and capitalize by strand (like pileup)
//        System.out.print(bytesToString(li.bases));
//        System.out.print("\t");
//        System.out.print(phredToFastq(li.qualities));
//        System.out.print("\n");
        }


        writer.flush();
        writer.close();
        final long elapsed = System.currentTimeMillis() - startTime;

        System.out.println("Completed in " + elapsed + "ms");
    }


    static String bytesToString(final List<Byte> data) {
        if (data == null || data.size() == 0) {
            return null;
        }

        final char[] chars = new char[data.size()];
        for (int i = 0; i < data.size(); i++) {
            chars[i] = (char) (data.get(i) & 0xFF);
        }
        return new String(chars);
    }


    static String phredToFastq(final List<Byte> data) {
        final byte[] arrData = new byte[data.size()];
        for(int i=0; i< data.size(); i++) { arrData[i] = data.get(i); }
        return phredToFastq(arrData);
    }

    static String phredToFastq(final byte[] data) {
        if (data == null) {
            return null;
        }
        return phredToFastq(data, 0, data.length);
    }

    static String phredToFastq(final byte[] buffer, final int offset, final int length) {
        final char[] chars = new char[length];
        for (int i = 0; i < length; i++) {
            chars[i] = phredToFastq(buffer[offset+i] & 0xFF);
        }
        return new String(chars);
    }

    static char phredToFastq(final int phredScore) {
        if (phredScore < 0 || phredScore > 63) {
            throw new IllegalArgumentException("Cannot encode phred score: " + phredScore);
        }
        return (char) (33 + phredScore);
    }

}