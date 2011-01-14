package org.broadinstitute.sting.utils.threading;

import com.google.common.collect.ArrayListMultimap;
import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.FastaSequenceIndexBuilder;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.file.FSLockWithShared;
import org.broadinstitute.sting.utils.file.FileSystemInabilityToLockException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.*;
import java.nio.channels.ClosedChannelException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 */
public class SharedFileGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
    private static final boolean DEBUG = false;
    private static Logger logger = Logger.getLogger(SharedFileGenomeLocProcessingTracker.class);
    private List<ProcessingLoc> processingLocs = new ArrayList<ProcessingLoc>();
    private File sharedFile = null;
    GenomeLocParser parser;

    // the file lock
    private FileLock lock = null;

    // the file channel we open
    private FileChannel channel = null;

    long lastReadPosition = 0;

    public SharedFileGenomeLocProcessingTracker(File sharedFile, GenomeLocParser parser) {
        try {
            this.sharedFile = sharedFile;
            this.channel = new RandomAccessFile(sharedFile, "rw").getChannel();
            this.parser = parser;
        }
        catch (Exception e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        }
    }

    private void lock() {
        try {
            lock = channel.lock();
        } catch (ClosedChannelException e) {
            throw new ReviewedStingException("Unable to lock file because the file channel is closed. " + sharedFile, e);
        }
//            catch (OverlappingFileLockException e) {
//                logger.debug("Unable to lock file because you already have a lock on this file.");
//                return false;
//            }
        catch (IOException e) {
            throw new ReviewedStingException("Coordination file could not be created because a lock could not be obtained.", e);
        } finally {
            unlock();
        }
    }

    private void unlock() {
        try {
            lock.release();
            //channel.close();
        } catch ( IOException e ) {
            throw new ReviewedStingException("Could not free lock on file " + sharedFile, e);
        }
    }

    private void readLocs() {
        try {
            if ( sharedFile.exists() ) {
                FileInputStream in = new FileInputStream(sharedFile);
                if ( in.getChannel().size() > lastReadPosition ) {
                    in.skip(lastReadPosition);

                    BufferedReader reader = new BufferedReader(new InputStreamReader(in));
                    int counter = 0;
                    String line = reader.readLine();   // Read another line
                    while ( line != null ) {
                        String[] parts = line.split(" ");
                        if ( parts.length != 2 ) throw new ReviewedStingException("BUG: bad sharedFile line " + line);
                        GenomeLoc loc = parser.parseGenomeLoc(parts[0]);
                        String owner = parts[1];
                        processingLocs.add(new ProcessingLoc(loc, owner));
                        line = reader.readLine();
                        counter++;
                    }
                    lastReadPosition = in.getChannel().position();
                    if ( DEBUG ) logger.warn(String.format("Read %s locs from file, current pos is %d, total locs is %d",
                            counter, lastReadPosition, processingLocs.size()));
                }
                in.close();
            }
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(sharedFile, e);
        } catch (IOException e) {
            throw new ReviewedStingException("Couldn't read sharedFile " + sharedFile, e);
        }
    }

    private void writeLoc(ProcessingLoc proc) {
        try {
            PrintStream out = new PrintStream(new FileOutputStream(sharedFile, true));
            out.printf("%s %s%n", proc.getLoc(), proc.getOwner());
            if ( DEBUG ) logger.warn(String.format("Wrote loc %s to file", proc));
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(sharedFile, e);
        }
    }

    public synchronized ProcessingLoc claimOwnership(GenomeLoc loc, String myName) {
        ProcessingLoc owner = null;
        try {
            lock();
            readLocs();
            owner = super.findOwner(loc);
            if ( owner == null ) { // we are unowned
                owner = new ProcessingLoc(loc, myName);
                writeLoc(owner);
            }
        } finally {
            unlock();
        }

        return owner;
    }

    protected synchronized List<ProcessingLoc> getProcessingLocs() {
        try {
            lock();
            readLocs();
        } finally {
            unlock();
        }

        return processingLocs;
    }
}
