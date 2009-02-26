/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.arachne;

import java.io.*;

/**
 * Utility to convert fastb to fasta files.
 * More importantly, can be used to extract a subset of the reads.
 */
public class Fastb2Fasta {

    private boolean mVerbose = false;
    private boolean mDebug = false;
    private String mInputPath = null;
    private String mIdListFilePath = null;


    public static void main(String[] args)
        throws Exception {
        new Fastb2Fasta().run(args);
    }

    private void usage() {
        System.out.println("Usage: Fastb2Fasta ... <fastb-file>");
        System.out.println("  -idlist <file-of-read-ids>");
        System.out.println("  -verbose");
        System.out.println("  -debug");
    }

    private boolean parseArguments(String[] args) {

        int argpos = 0;
        int argsleft = 0;

        while (argpos < args.length) {
            argsleft = args.length - argpos;
            String arg = args[argpos];
            if (arg.equals("-idlist") && argsleft > 1) {
                argpos++;
                mIdListFilePath = args[argpos++];
            } else if (arg.equals("-verbose")) {
                argpos++;
                mVerbose = true;
            } else if (arg.equals("-debug")) {
                argpos++;
                mDebug = true;
            } else if (arg.startsWith("-")) {
                usage();
                return false;
            } else {
                break;
            }
        }

        argsleft = args.length - argpos;
        if (argsleft != 1) {
            usage();
            return false;
        }

        mInputPath = args[argpos];
        return true;
    }

    private void run(String[] args)
        throws Exception {

        if (!parseArguments(args)) {
            System.exit(1);
        }

        FastbReader fastbReader = new FastbReader(new File(mInputPath));
        try {
            if (mIdListFilePath != null) {
                LineNumberReader reader = new LineNumberReader(new FileReader(mIdListFilePath));
                while (true) {
                    String line = reader.readLine();
                    if (line == null) {
                        reader.close();
                        break;
                    }
                    Integer id = parseReadId(line);
                    if (id == null) {
                        continue;
                    }
                    if (id < 0 || id >= fastbReader.getSequenceCount()) {
                        System.out.println("ERROR: Illegal sequence id: " + id);
                        System.exit(1);
                    }
                    String sequence = fastbReader.readSequence(id);
                    System.out.println(">" + id);
                    System.out.println(sequence);
                }
            } else {
                int id = 0;
                while (fastbReader.hasNext()) {
                    String sequence = fastbReader.next();
                    System.out.println(">" + id);
                    System.out.println(sequence);
                    id++;
                }
            }
        } finally {
            fastbReader.close();
        }
    }

    private Integer parseReadId(String line) {
        String text = line.trim();
        if (text.length() == 0 || text.charAt(0) == '#') {
            return null;
        }
        String token = text.split("\\s+")[0];
        Integer id = null;
        try {
            id = new Integer(token);
        } catch (NumberFormatException exc) {
            System.out.println("ERROR: Invalid sequence id: " + token);
            System.exit(1);
        }
        return id;
    }
}
