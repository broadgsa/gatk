/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.io;

import htsjdk.samtools.*;
import htsjdk.samtools.util.ProgressLoggerInterface;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: May 31, 2011
 * Time: 3:52:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class NWaySAMFileWriter implements SAMFileWriter {

    private Map<SAMReaderID,SAMFileWriter> writerMap = null;
    private boolean presorted ;
    GenomeAnalysisEngine toolkit;
    boolean KEEP_ALL_PG_RECORDS = false;

    public NWaySAMFileWriter(GenomeAnalysisEngine toolkit, Map<String,String> in2out, SAMFileHeader.SortOrder order,
                             boolean presorted, boolean indexOnTheFly, boolean generateMD5, SAMProgramRecord pRecord, boolean keep_records) {
        this.presorted = presorted;
        this.toolkit = toolkit;
        this.KEEP_ALL_PG_RECORDS = keep_records;
        writerMap = new HashMap<SAMReaderID,SAMFileWriter>();
        setupByReader(toolkit,in2out,order, presorted, indexOnTheFly, generateMD5, pRecord);
    }

    public NWaySAMFileWriter(GenomeAnalysisEngine toolkit, String ext, SAMFileHeader.SortOrder order,
                              boolean presorted, boolean indexOnTheFly , boolean generateMD5, SAMProgramRecord pRecord, boolean keep_records) {
        this.presorted = presorted;
        this.toolkit = toolkit;
        this.KEEP_ALL_PG_RECORDS = keep_records;
        writerMap = new HashMap<SAMReaderID,SAMFileWriter>();
        setupByReader(toolkit,ext,order, presorted, indexOnTheFly, generateMD5, pRecord);
    }

    public NWaySAMFileWriter(GenomeAnalysisEngine toolkit, Map<String,String> in2out, SAMFileHeader.SortOrder order,
                             boolean presorted, boolean indexOnTheFly, boolean generateMD5) {
        this(toolkit, in2out, order, presorted, indexOnTheFly, generateMD5, null,false);
    }

    public NWaySAMFileWriter(GenomeAnalysisEngine toolkit, String ext, SAMFileHeader.SortOrder order,
                              boolean presorted, boolean indexOnTheFly , boolean generateMD5) {
        this(toolkit, ext, order, presorted, indexOnTheFly, generateMD5, null,false);
    }

    /**
     * Creates a program record for the program, adds it to the list of program records (@PG tags) in the bam file and sets
     * up the writer with the header and presorted status.
     *
     * @param originalHeader      original header
     * @param programRecord       the program record for this program
     */
    public static SAMFileHeader setupWriter(final SAMFileHeader originalHeader, final SAMProgramRecord programRecord) {
        final SAMFileHeader header = originalHeader.clone();
        final List<SAMProgramRecord> oldRecords = header.getProgramRecords();
        final List<SAMProgramRecord> newRecords = new ArrayList<SAMProgramRecord>(oldRecords.size()+1);
        for ( SAMProgramRecord record : oldRecords )
            if ( (programRecord != null && !record.getId().startsWith(programRecord.getId())))
                newRecords.add(record);

        if (programRecord != null) {
            newRecords.add(programRecord);
            header.setProgramRecords(newRecords);
        }
        return header;
    }

    /**
    * Creates a program record for the program, adds it to the list of program records (@PG tags) in the bam file and returns
    * the new header to be added to the BAM writer.
    *
    * @param toolkit             the engine
    * @param walker              the walker object (so we can extract the command line)
    * @param PROGRAM_RECORD_NAME the name for the PG tag
    * @return a pre-filled header for the bam writer
    */
    public static SAMFileHeader setupWriter(final GenomeAnalysisEngine toolkit, final SAMFileHeader originalHeader, final Object walker, final String PROGRAM_RECORD_NAME) {
        final SAMProgramRecord programRecord = createProgramRecord(toolkit, walker, PROGRAM_RECORD_NAME);
        return setupWriter(originalHeader, programRecord);
    }

    /**
     * Creates a program record for the program, adds it to the list of program records (@PG tags) in the bam file and sets
     * up the writer with the header and presorted status.
     *
     * @param writer              BAM file writer
     * @param toolkit             the engine
     * @param preSorted           whether or not the writer can assume reads are going to be added are already sorted
     * @param walker              the walker object (so we can extract the command line)
     * @param PROGRAM_RECORD_NAME the name for the PG tag
     */
    public static void setupWriter(GATKSAMFileWriter writer, GenomeAnalysisEngine toolkit, SAMFileHeader originalHeader, boolean preSorted, Object walker, String PROGRAM_RECORD_NAME) {
        SAMFileHeader header = setupWriter(toolkit, originalHeader, walker, PROGRAM_RECORD_NAME);
        writer.writeHeader(header);
        writer.setPresorted(preSorted);
    }

    /**
     * Creates a program record (@PG) tag
     *
     * @param toolkit             the engine
     * @param walker              the walker object (so we can extract the command line)
     * @param PROGRAM_RECORD_NAME the name for the PG tag
     * @return a program record for the tool
     */
    public static SAMProgramRecord createProgramRecord(GenomeAnalysisEngine toolkit, Object walker, String PROGRAM_RECORD_NAME) {
        final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
        try {
            programRecord.setProgramVersion(CommandLineProgram.getVersionNumber());
        } catch (MissingResourceException e) {
            // couldn't care less if the resource is missing...
        }
        programRecord.setCommandLine(toolkit.createApproximateCommandLineArgumentString(toolkit, walker));
        return programRecord;
    }

    /**
     * Instantiates multiple underlying SAM writes, one per input SAM reader registered with GATK engine (those will be retrieved
     * from <code>toolkit</code>). The <code>in2out</code> map must contain an entry for each input filename and map it
     * onto a unique output file name.
     * @param toolkit
     * @param in2out
     */
    public void setupByReader(GenomeAnalysisEngine toolkit, Map<String,String> in2out, SAMFileHeader.SortOrder order,
                              boolean presorted, boolean indexOnTheFly, boolean generateMD5, SAMProgramRecord pRecord) {
        if ( in2out==null ) throw new GATKException("input-output bam filename map for n-way-out writing is NULL");
        for ( SAMReaderID rid : toolkit.getReadsDataSource().getReaderIDs() ) {

            String fName = toolkit.getReadsDataSource().getSAMFile(rid).getName();

            String outName;
            if ( ! in2out.containsKey(fName) )
                    throw new UserException.BadInput("Input-output bam filename map does not contain an entry for the input file "+fName);
            outName = in2out.get(fName);

            if ( writerMap.containsKey( rid ) )
                throw new GATKException("nWayOut mode: Reader id for input sam file "+fName+" is already registered; "+
                        "map file likely contains multiple entries for this input file");

            addWriter(rid,outName, order, presorted, indexOnTheFly, generateMD5, pRecord);
        }

    }

    /**
     * Instantiates multiple underlying SAM writes, one per input SAM reader registered with GATK engine (those will be retrieved
     * from <code>toolkit</code>). The output file names will be generated automatically by stripping ".sam" or ".bam" off the
     * input file name and adding ext instead (e.g. ".cleaned.bam").
     * onto a unique output file name.
     * @param toolkit
     * @param ext
     */
    public void setupByReader(GenomeAnalysisEngine toolkit, String ext, SAMFileHeader.SortOrder order,
                              boolean presorted, boolean indexOnTheFly, boolean generateMD5, SAMProgramRecord pRecord) {
        for ( SAMReaderID rid : toolkit.getReadsDataSource().getReaderIDs() ) {

            String fName = toolkit.getReadsDataSource().getSAMFile(rid).getName();

            String outName;
            int pos ;
            if ( fName.toUpperCase().endsWith(".BAM") ) pos = fName.toUpperCase().lastIndexOf(".BAM");
            else {
                if ( fName.toUpperCase().endsWith(".SAM") ) pos = fName.toUpperCase().lastIndexOf(".SAM");
                else throw new UserException.BadInput("Input file name "+fName+" does not end with .sam or .bam");
            }
            String prefix = fName.substring(0,pos);
            outName = prefix+ext;

            if ( writerMap.containsKey( rid ) )
                throw new GATKException("nWayOut mode: Reader id for input sam file "+fName+" is already registered");
            addWriter(rid,outName, order, presorted, indexOnTheFly, generateMD5, pRecord);
        }

    }

    private void addWriter(SAMReaderID id , String outName, SAMFileHeader.SortOrder order, boolean presorted,
                           boolean indexOnTheFly, boolean generateMD5, SAMProgramRecord programRecord) {
        File f = new File(outName);
        SAMFileHeader header = setupWriter(toolkit.getSAMFileHeader(id), programRecord);
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCreateIndex(indexOnTheFly);
        factory.setCreateMd5File(generateMD5);
        SAMFileWriter sw = factory.makeSAMOrBAMWriter(header, presorted, f);
        writerMap.put(id,sw);
    }

    public Collection<SAMFileWriter> getWriters() {
        return writerMap.values();
    }

    public void addAlignment(SAMRecord samRecord) {
        final SAMReaderID id = toolkit.getReaderIDForRead(samRecord);
        String rg = samRecord.getStringAttribute("RG");
        if ( rg != null ) {
            String rg_orig = toolkit.getReadsDataSource().getOriginalReadGroupId(rg);
            samRecord.setAttribute("RG",rg_orig);
        }
        addAlignment(samRecord, id);
    }

    public void addAlignment(SAMRecord samRecord, SAMReaderID readerID) {
        writerMap.get(readerID).addAlignment(samRecord);
    }

    public SAMFileHeader getFileHeader() {
        return toolkit.getSAMFileHeader();
    }

    public void close() {
        for ( SAMFileWriter w : writerMap.values() ) w.close();
    }

    @Override
    public void setProgressLogger(final ProgressLoggerInterface logger) {
        for (final SAMFileWriter writer: writerMap.values()) {
            writer.setProgressLogger(logger);
        }
    }
}
