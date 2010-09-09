package org.broadinstitute.sting.oneoffprojects.tools;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Date;
import java.util.zip.GZIPInputStream;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;

/**
 * Used to test how long it takes to read through a text files and gzipped files.
 * If the passed-in filename ends with .gz, it will be read using GZIPInputStream.
 * Otherwise, its read using FileReader.
 */
public class JavaIOSpeedTest extends CommandLineProgram {

    @Argument(fullName = "input", shortName = "i", doc = "file to read", required = true) public String filename;
    @Argument(fullName = "buffer_size", shortName = "s", doc = "read buffer size in mb", required=false) public int bufferSize = -1;

    @Override
    protected int execute() {
        System.out.println("Filename: " + filename);
        try {
            final File file = new File(filename);
            final Reader reader;
            if(filename.endsWith(".gz")) {
                reader = new InputStreamReader(new GZIPInputStream(new BufferedInputStream(new FileInputStream(file))));
            } else {
                reader = new FileReader(file);
            }

            final BufferedReader br;
            if(bufferSize != -1) {
                br = new BufferedReader(reader, bufferSize * 1000000);
            } else {
                br = new BufferedReader(reader);
            }

            System.out.println("Reading...");
            int incr = 10000000;
            int next = incr;
            int counter = 0;
            while(br.ready()) {
                br.readLine();
                if(++counter == next) {
                    next += incr;
                    System.err.println(new Date() + " - file: " + filename + ",  buffer size: " + bufferSize + "mb,  read " + counter + " lines...");
                }
            }

            System.out.println("Read " + counter + " lines from " + filename);

        } catch(Exception e) {
            System.err.println("ERROR: " + e);
            e.printStackTrace();
        }

        return 0;
    }

    public static void main(String args[])
    {
        try {
            CommandLineProgram.start(new JavaIOSpeedTest(), args);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

}
