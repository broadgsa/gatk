package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 25, 2009
 * Time: 8:27:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class PassThroughWriter implements RecordReceiver {
    private SAMFileWriter writer;
    private int reads_written = 0;

    public PassThroughWriter( File f, SAMFileHeader h) {
             writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(h, false, f);
    }

    public PassThroughWriter(String s, SAMFileHeader h) {
        this(new File(s), h);
    }

    public void receive(SAMRecord r) {
        //To change body of implemented methods use File | Settings | File Templates.
        writer.addAlignment(r);
        reads_written++;
    }

    public void close() { writer.close() ; }
    
    /** Returns the number of reads that were so far received by this writer.
     * 
     */
    public int getNumReadsReceived() { return reads_written; }
}
