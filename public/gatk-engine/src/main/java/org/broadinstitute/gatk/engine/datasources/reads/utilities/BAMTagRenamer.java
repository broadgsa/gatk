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

package org.broadinstitute.gatk.engine.datasources.reads.utilities;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;

import java.io.File;

/**
 * A simple utility written directly in Picard that will rename tags
 * from one name to another.
 *
 * @author hanna
 * @version 0.1
 */

public class BAMTagRenamer extends CommandLineProgram {
    @Argument(fullName="input",shortName="I",doc="Input file to process",required=true)
    private File input = null;

    @Argument(fullName="output",shortName="O",doc="Output file to create",required=true)
    private File output = null;

    @Argument(fullName="bam_compression",shortName="compress",doc="Compression level to use when writing the BAM file.",required=false)
    private int compressionLevel = 5;

    @Argument(fullName="original_tag_name",shortName="otn",doc="Tag name to be replaced.",required=true)
    private String sourceTagName = null;

    @Argument(fullName="replacement_tag_name",shortName="rtn",doc="Tag name to be used as a replacement.",required=true)
    private String targetTagName = null;

    public int execute() {
        long readsWritten = 0;
        long readsAltered = 0;

        SAMFileReader reader = new SAMFileReader(input);
        SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(),true,output,compressionLevel);

        for(SAMRecord read: reader) {
            Object value = read.getAttribute(sourceTagName);
            if(value != null) {
                read.setAttribute(sourceTagName,null);
                read.setAttribute(targetTagName,value);
                readsAltered++;
            }
            writer.addAlignment(read);
            readsWritten++;
            if(readsWritten % 1000000 == 0)
                System.out.printf("%d reads written.  %d tag names updated from %s to %s.%n",readsWritten,readsAltered,sourceTagName,targetTagName);
        }

        writer.close();
        System.out.printf("%d reads written.  %d tag names updated from %s to %s.%n",readsWritten,readsAltered,sourceTagName,targetTagName);        

        return 0;
    }

    /**
     * Required main method implementation.
     */
    public static void main(String[] argv) {
        BAMTagRenamer instance = new BAMTagRenamer();
        try {
            start(instance, argv);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        System.exit(CommandLineProgram.result);
    }
}
