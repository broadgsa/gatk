/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

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

    public NWaySAMFileWriter(GenomeAnalysisEngine toolkit, Map<String,String> in2out, SAMFileHeader.SortOrder order, boolean presorted, boolean indexOnTheFly) {
        this.presorted = presorted;
        this.toolkit = toolkit;
        writerMap = new HashMap<SAMReaderID,SAMFileWriter>();
        setupByReader(toolkit,in2out,order, presorted, indexOnTheFly);
    }

    public NWaySAMFileWriter(GenomeAnalysisEngine toolkit, String ext, SAMFileHeader.SortOrder order, boolean presorted, boolean indexOnTheFly ) {
        this.presorted = presorted;
        this.toolkit = toolkit;
        writerMap = new HashMap<SAMReaderID,SAMFileWriter>();
        setupByReader(toolkit,ext,order, presorted, indexOnTheFly);
    }


    /**
     * Instantiates multiple underlying SAM writes, one per input SAM reader registered with GATK engine (those will be retrieved
     * from <code>toolkit</code>). The <code>in2out</code> map must contain an entry for each input filename and map it
     * onto a unique output file name.
     * @param toolkit
     * @param in2out
     */
    public void setupByReader(GenomeAnalysisEngine toolkit, Map<String,String> in2out, SAMFileHeader.SortOrder order, boolean presorted, boolean indexOnTheFly) {

        if ( in2out==null ) throw new StingException("input-output bam filename map for n-way-out writing is NULL");
        for ( SAMReaderID rid : toolkit.getReadsDataSource().getReaderIDs() ) {

            String fName = toolkit.getReadsDataSource().getSAMFile(rid).getName();

            String outName;
            if ( ! in2out.containsKey(fName) )
                    throw new UserException.BadInput("Input-output bam filename map does not contain an entry for the input file "+fName);
            outName = in2out.get(fName);

            if ( writerMap.containsKey( rid ) )
                throw new StingException("nWayOut mode: Reader id for input sam file "+fName+" is already registered");

            addWriter(rid,outName, order, presorted, indexOnTheFly);
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
    public void setupByReader(GenomeAnalysisEngine toolkit, String ext, SAMFileHeader.SortOrder order, boolean presorted, boolean indexOnTheFly) {
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
                throw new StingException("nWayOut mode: Reader id for input sam file "+fName+" is already registered");

            addWriter(rid,outName, order, presorted, indexOnTheFly);
        }

    }

    private void addWriter(SAMReaderID id , String outName, SAMFileHeader.SortOrder order, boolean presorted, boolean indexOnTheFly) {
        File f = new File(outName);
        SAMFileHeader header = toolkit.getSAMFileHeader(id).clone();
        header.setSortOrder(order);
        SAMFileWriter sw = new SAMFileWriterFactory().setCreateIndex(indexOnTheFly).makeSAMOrBAMWriter(header, presorted, f);
        writerMap.put(id,sw);
    }

    
    public void addAlignment(SAMRecord samRecord) {
        final SAMReaderID id = toolkit.getReaderIDForRead(samRecord);
        writerMap.get(id).addAlignment(samRecord);
    }

    public SAMFileHeader getFileHeader() {
        return toolkit.getSAMFileHeader();  
    }

    public void close() {
        for ( SAMFileWriter w : writerMap.values() ) w.close();
    }
}
