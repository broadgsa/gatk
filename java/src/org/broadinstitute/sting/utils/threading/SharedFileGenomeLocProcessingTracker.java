package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.*;
import java.nio.channels.*;
import java.util.ArrayList;
import java.util.List;

/**
 *
 */
//public class SharedFileGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
//    private static final boolean DEBUG = false;
//    private static final boolean REALLY_DEBUG = false;
//
//    private boolean ACTUALLY_USE_FILE_LOCK = true;
//
//    private static Logger logger = Logger.getLogger(SharedFileGenomeLocProcessingTracker.class);
//
//    private Object myLock = new Object();
//    private List<ProcessingLoc> processingLocs;
//    private File sharedFile = null;
//    private GenomeLocParser parser;
//    private FileLock lock = null;
//    private RandomAccessFile raFile;
//    private long lastReadPosition = 0;
//
////    //
////    // TODO -- I CAN'T FOR SOME REASON GET THE FILE LOCK TESTING TO WORK WITH MULTIPLE THREADS IN THE UNIT TEST
////    // TODO -- IT SEEMS THAT SOME LOCKS AREN'T BEING FREED, BUT IT DOESN'T SEEM POSSIBLE GIVEN THE CHECKS
////    // TODO -- IN THE CODE.  I THINK THE LOCK IS SOMEHOW CONTINUING BEYOND THE UNLOCK CALL, OR THAT I NEED
////    // TODO -- TO CLOSE AND REOPEN THE CHANNEL FOR EACH LOCK?
////    //
////
//    public SharedFileGenomeLocProcessingTracker(File sharedFile, GenomeLocParser parser) {
//        this(sharedFile, parser, true);
//    }
//
//    protected SharedFileGenomeLocProcessingTracker(File sharedFile, GenomeLocParser parser, boolean useFileLock) {
//        processingLocs = new ArrayList<ProcessingLoc>();
//        ACTUALLY_USE_FILE_LOCK = false;
//        try {
//            this.sharedFile = sharedFile;
//            this.raFile = new RandomAccessFile(sharedFile, "rws");
//            this.parser = parser;
//        }
//        catch (FileNotFoundException e) {
//            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
//        }
//    }
//
//    public void close() {
//        if ( ACTUALLY_USE_FILE_LOCK ) {
//            try {
//                this.raFile.close();
//            }
//            catch (IOException e) {
//                throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
//            }
//        }
//    }
//
//    private void lock() {
//        if ( ACTUALLY_USE_FILE_LOCK ) {
//
//            // Precondition -- lock is always null while we don't have a lock
//            if ( lock != null )
//                throw new ReviewedStingException("BUG: lock() function called when a lock already is owned!");
//
//            try {
//                lock = raFile.getChannel().lock();
//            } catch (ClosedChannelException e) {
//                throw new ReviewedStingException("Unable to lock file because the file channel is closed. " + sharedFile, e);
//            } catch (FileLockInterruptionException e) {
//                throw new ReviewedStingException("File lock interrupted", e);
//            } catch (NonWritableChannelException e) {
//                throw new ReviewedStingException("File channel not writable", e);
//            } catch (OverlappingFileLockException e) {
//                // this only happens when multiple threads are running, and one is waiting
//                // for the lock above and we come here.
//                throw new ReviewedStingException("BUG: Failed to acquire lock, should never happen.");
//            } catch (IOException e) {
//                throw new ReviewedStingException("Coordination file could not be created because a lock could not be obtained.", e);
//            }
//        }
//    }
//
//    private void unlock(boolean excepting) {
//        if ( ACTUALLY_USE_FILE_LOCK ) {
//
//            // Precondition -- lock is never null while we have a lock
//            if ( lock == null ) {
//                if ( ! excepting )
//                    throw new ReviewedStingException("BUG: call to unlock() when we don't have the lock!");
//            } else {
//                if ( ! lock.isValid() )
//                    throw new ReviewedStingException("BUG: call to unlock() when we don't have a valid lock!");
//                try {
//                    lock.release();
//                    lock = null;
//                    //channel.close();
//                } catch ( IOException e ) {
//                    throw new ReviewedStingException("Could not free lock on file " + sharedFile, e);
//                }
//            }
//        }
//    }
//
//    private List<ProcessingLoc> readLocs() {
//        if ( ACTUALLY_USE_FILE_LOCK ) {
//            // we must have a lock to run this code
//            if ( lock == null || ! lock.isValid() ) throw new ReviewedStingException("File lock must be valid upon entry to readLocs()");
//
//            try {
//                if ( raFile.length() > lastReadPosition ) {
//                    raFile.seek(lastReadPosition);
//
//                    int counter = 0;
//                    String line = raFile.readLine();   // Read another line
//                    while ( line != null ) {
//                        String[] parts = line.split(" ");
//                        if ( parts.length != 2 ) throw new ReviewedStingException("BUG: bad sharedFile line '" + line + "' at " + raFile.getFilePointer());
//                        GenomeLoc loc = parser.parseGenomeLoc(parts[0]);
//                        String owner = parts[1];
//                        processingLocs.add(new ProcessingLoc(loc, owner));
//                        line = raFile.readLine();
//                        counter++;
//                    }
//                    lastReadPosition = raFile.getFilePointer();
//                    if ( DEBUG ) logger.warn(String.format("Read %s locs from file, current pos is %d, total locs is %d",
//                            counter, lastReadPosition, processingLocs.size()));
//                }
//            } catch (FileNotFoundException e) {
//                throw new UserException.CouldNotReadInputFile(sharedFile, e);
//            } catch (IOException e) {
//                throw new ReviewedStingException("Couldn't read sharedFile " + sharedFile, e);
//            }
//        }
//
//        return processingLocs;
//    }
//
//    private void writeLoc(ProcessingLoc proc) {
//        if ( ACTUALLY_USE_FILE_LOCK ) {
//            // we must have a lock to run this code
//            if ( lock == null || ! lock.isValid() )
//                throw new ReviewedStingException("File lock must be valid upon entry to writeLoc()");
//
//            try {
//                String packet = String.format("%s %s%n", proc.getLoc(), proc.getOwner());
//                long startPos = raFile.getFilePointer();
//                raFile.seek(raFile.length());
//                raFile.write(packet.getBytes());
//                if ( DEBUG ) logger.warn(String.format("Wrote loc %s to file: %d + %d bytes ending at %d", proc, startPos, packet.length(), raFile.getFilePointer()));
//            } catch (FileNotFoundException e) {
//                throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
//            } catch (IOException e) {
//                throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
//            }
//        } else {
//            processingLocs.add(proc);
//        }
//    }
//
//    private final void printOwners() {
//        for ( ProcessingLoc proc : processingLocs )
//            System.out.println(proc);
//    }
//
//    public ProcessingLoc claimOwnership(GenomeLoc loc, String myName) {
//        if ( REALLY_DEBUG ) System.out.printf("  claimOwnership %s%n", myName);
//        synchronized (processingLocs) {
//            boolean excepting = true;
//            ProcessingLoc owner = null;
//
//            if ( lock != null ) throw new ReviewedStingException("BUG: into claimOwnership synchronized block while another thread owns the lock");
//
//            if ( REALLY_DEBUG ) System.out.printf("  sync raFile %s %s%n", myName, raFile);
//            try {
//                lock();
//                owner = findOwnerInUnsortedList(loc, readLocs());
//                //owner = super.findOwner(loc);
//                if ( owner == null ) { // we are unowned
//                    owner = new ProcessingLoc(loc, myName);
//                    writeLoc(owner);
//                }
//                excepting = false;
//            } finally {
//                if ( REALLY_DEBUG ) System.out.printf("  claimOwnership unlock %s excepting %s, owner %s%n", myName, excepting, owner);
//                //printOwners();
//                unlock(excepting);
//            }
//
//            if ( lock != null ) throw new ReviewedStingException("BUG: exiting claimOwnership synchronized block without setting lock to null");
//            return owner;
//        }
//    }
//
//    protected List<ProcessingLoc> getProcessingLocs() {
//        synchronized (processingLocs) {
//            boolean excepting = true;
//            if ( lock != null ) throw new ReviewedStingException("BUG: into claimOwnership synchronized block while another thread owns the lock");
//
//            try {
//                lock();
//                readLocs();
//                excepting = false;
//            } finally {
//                unlock(excepting);
//            }
//
//            if ( lock != null ) throw new ReviewedStingException("BUG: exiting getProcessingLocs synchronized block without setting lock to null");
//            return processingLocs;
//        }
//    }
//}
