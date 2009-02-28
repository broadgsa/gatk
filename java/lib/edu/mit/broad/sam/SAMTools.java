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
package edu.mit.broad.sam;


import edu.mit.broad.sam.util.CloseableIterator;
import java.io.*;


/**
 * Command line utility for manipulating SAM/BAM files.
 */
public class SAMTools
{
    private String mCommand = null;
    private File mInputFile = null;


    public static void main(final String[] args)
        throws Exception {
        final int status = new SAMTools().run(args);
        if (status != 0) {
            System.exit(status);
        }
    }

    private SAMTools() {
    }

    private void usage() {
        System.out.println();
        System.out.println("SAMTools version 0.1.0");
        System.out.println("Tools for manipulating SAM/BAM files");
        System.out.println();
        System.out.println("Usage: SAMTools <command> <options...>");
        System.out.println();
        System.out.println("Commands:");
        System.out.println("  help");
        System.out.println("  view        <file>");
        System.out.println();
    }

    private boolean parseArguments(final String[] args) {
        if (args.length == 0) {
            usage();
            return true;
        }
        final String command = args[0];
        final int argpos = 1;
        final int argcount = args.length - argpos;
        if (command.equals("help")) {
            usage();
            return true;
        } else if (command.equals("view")) {
            if (argcount != 1) {
                usage();
                return false;
            }
            mInputFile = new File(args[1]);
            if (!mInputFile.exists()) {
                System.out.println("Input file not found: " + mInputFile);
                return false;
            }
        } else {
            System.out.println("Unrecognized command: " + command);
            System.out.println();
            usage();
            return false;
        }
        mCommand = command;
        return true;
    }

    private int run(final String[] args)
        throws Exception {
        if (!parseArguments(args)) {
            return 1;
        }
        if (mCommand == null) {
            return 0;
        }
        if (mCommand.equals("view")) {
            return runView();
        }
        return 1;
    }

    private int runView() {
        final SAMFileReader reader = new SAMFileReader(mInputFile);
        final CloseableIterator<SAMRecord> iterator = reader.iterator();
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            System.out.println(record.format());
        }
        iterator.close();
        return 0;
    }
}
