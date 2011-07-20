package org.broadinstitute.sting.utils.distributedutils;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Keeps a copy of the processing locks in a file
 */
public class FileBackedGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
    private static final Logger logger = Logger.getLogger(FileBackedGenomeLocProcessingTracker.class);
    private static final boolean DEBUG = false;
    private static final String READ_MODE = "r";
    private static final String WRITE_MODE = "rws";

    private final File sharedFile;
    private final GenomeLocParser parser;
    private long lastReadPosition = 0;

    public FileBackedGenomeLocProcessingTracker(File sharedFile, GenomeLocParser parser, ClosableReentrantLock lock, PrintStream status) {
        super(lock, status);

        this.sharedFile = sharedFile;
        this.parser = parser;
    }

    private RandomAccessFile openFile(String mode) {
        try {
            return new RandomAccessFile(sharedFile, mode);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        }
    }

    private void closeFile(RandomAccessFile raFile) {
        try {
            if ( raFile != null ) raFile.close();
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        }
    }

    @Override
    protected List<ProcessingLoc> readNewLocs() {
        List<ProcessingLoc> newPLocs = new ArrayList<ProcessingLoc>();    // todo -- gratitous object creation

        if ( sharedFile.exists() ) {
            RandomAccessFile raFile = null;
            try {
                raFile = openFile(READ_MODE);
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
            } finally {
                closeFile(raFile);
            }
        }

        return newPLocs;
        }

    @Override
    protected void registerNewLocs(Collection<ProcessingLoc> plocs) {
        RandomAccessFile raFile = null;

        try {
            raFile = openFile(WRITE_MODE);
            long startPos = raFile.getFilePointer();
            raFile.seek(raFile.length());
            StringBuffer bytes = new StringBuffer();
            for ( ProcessingLoc ploc : plocs ) {
                String packet = String.format("%s %s%n", ploc.getLocation(), ploc.getOwner());
                bytes.append(packet);
                if ( DEBUG ) logger.warn(String.format("Wrote loc %s to file: %d + %d bytes ending at %d", ploc, startPos, packet.length(), raFile.getFilePointer()));
            }
            raFile.write(bytes.toString().getBytes());
            //raFile.getChannel().force(true);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        } finally {
            closeFile(raFile);
        }
    }
}


