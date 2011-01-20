package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Keeps a copy of the processing locks in a file, in addition to tracking in memory via the base class
 */
public class FileBackedGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
    private static Logger logger = Logger.getLogger(FileBackedGenomeLocProcessingTracker.class);
    private static final boolean DEBUG = false;
    private File sharedFile = null;
    private GenomeLocParser parser;
    private RandomAccessFile raFile;
    private long lastReadPosition = 0;

    protected FileBackedGenomeLocProcessingTracker(File sharedFile, RandomAccessFile raFile, GenomeLocParser parser, ClosableReentrantLock lock) {
        super(lock);

        this.sharedFile = sharedFile;
        this.raFile = raFile;
        this.parser = parser;
    }

    protected void close() {
        super.close();
        try {
            raFile.close();
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        }
    }

    @Override
    protected List<ProcessingLoc> readNewLocs() {
        List<ProcessingLoc> newPLocs = new ArrayList<ProcessingLoc>();    // todo -- gratitous object creation
        try {
            //logger.warn(String.format("Reading new locs at: file.length=%d last=%d", raFile.length(), lastReadPosition));
            if ( raFile.length() > lastReadPosition ) {
                raFile.seek(lastReadPosition);

                int counter = 0;
                String line = raFile.readLine();   // Read another line
                while ( line != null ) {
                    String[] parts = line.split(" ");
                    if ( parts.length != 2 ) throw new ReviewedStingException("BUG: bad sharedFile line '" + line + "' at " + raFile.getFilePointer());
                    ProcessingLoc ploc = new ProcessingLoc(parser.parseGenomeLoc(parts[0]), parts[1]);
                    //logger.warn("    Read " + ploc);
                    newPLocs.add(ploc);
                    line = raFile.readLine();
                    counter++;
                }
                lastReadPosition = raFile.getFilePointer();
                if ( DEBUG ) logger.warn(String.format("Read %s locs from file, current pos is %d, # read new locs is %d",
                        counter, lastReadPosition, newPLocs.size()));
            }
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(sharedFile, e);
        } catch (IOException e) {
            throw new ReviewedStingException("Couldn't read sharedFile " + sharedFile, e);
        }

        return newPLocs;
    }

    @Override
    protected void registerNewLoc(ProcessingLoc proc) {
        try {
            String packet = String.format("%s %s%n", proc.getLocation(), proc.getOwner());
            long startPos = raFile.getFilePointer();
            raFile.seek(raFile.length());
            raFile.write(packet.getBytes());
            if ( DEBUG ) logger.warn(String.format("Wrote loc %s to file: %d + %d bytes ending at %d", proc, startPos, packet.length(), raFile.getFilePointer()));
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        }
    }
}


