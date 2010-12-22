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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileReaderArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import net.sf.picard.filter.SamRecordFilter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

/**
 * @author aaron
 */
public abstract class CommandLineExecutable extends CommandLineProgram {
    /**
     * The actual engine which performs the analysis.
     */
    protected GenomeAnalysisEngine engine = new GenomeAnalysisEngine();

    // get the analysis name
    public abstract String getAnalysisName();

    /**
     * Gets the GATK argument bundle.
     * @return A structure consisting of whatever arguments should be used to initialize the GATK engine.
     */
    protected abstract GATKArgumentCollection getArgumentCollection();

    /**
     * A list of all the arguments initially used as sources.
     */
    private final Collection<Object> argumentSources = new ArrayList<Object>();

    /**
     * this is the function that the inheriting class can expect to have called
     * when the command line system has initialized.
     *
     * @return the return code to exit the program with
     */
    protected int execute() throws Exception {
        engine.setParser(parser);
        argumentSources.add(this);

        Walker<?,?> walker = engine.getWalkerByName(getAnalysisName());

        try {
            engine.setArguments(getArgumentCollection());

            // File lists can require a bit of additional expansion.  Set these explicitly by the engine. 
            engine.setSAMFileIDs(unpackBAMFileList(getArgumentCollection()));
            engine.setReferenceMetaDataFiles(unpackRODBindings(getArgumentCollection()));

            engine.setWalker(walker);
            walker.setToolkit(engine);

            Collection<SamRecordFilter> filters = engine.createFilters();
            engine.setFilters(filters);

            // load the arguments into the walker / filters.
            // TODO: The fact that this extra load call exists here when all the parsing happens at the engine
            // TODO: level indicates that we're doing something wrong.  Turn this around so that the GATK can drive
            // TODO: argument processing.
            loadArgumentsIntoObject(walker);
            argumentSources.add(walker);

            for (SamRecordFilter filter: filters) {
                loadArgumentsIntoObject(filter);
                argumentSources.add(filter);
            }

            engine.execute();
            generateGATKRunReport(walker);
        } catch ( Exception e ) {
            generateGATKRunReport(walker, e);
            throw e;
        }

        // always return 0
        return 0;
    }

    /**
     * Generate the GATK run report for this walker using the current GATKEngine, if -et is enabled.
     * This report will be written to either STDOUT or to the run repository, depending on the options
     * for -et.
     *
     * @param e the exception, can be null if no exception occurred
     */
    private void generateGATKRunReport(Walker<?,?> walker, Exception e) {
        if ( getArgumentCollection().phoneHomeType != GATKRunReport.PhoneHomeOption.NO_ET ) {
            GATKRunReport report = new GATKRunReport(walker, e, engine, getArgumentCollection().phoneHomeType );
            if ( getArgumentCollection().phoneHomeType == GATKRunReport.PhoneHomeOption.STDOUT )
                report.postReport(System.out);
            else
                report.postReport();
        }
    }

    /**
     * Convenience method for fully parameterized generateGATKRunReport when an exception has
     * not occurred
     *
     * @param walker
     */
    private void generateGATKRunReport(Walker<?,?> walker) {
        generateGATKRunReport(walker, null);
    }

    /**
     * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
     * @return A collection of type descriptors generating implementation-dependent placeholders.
     */
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        return Arrays.asList( new VCFWriterArgumentTypeDescriptor(engine,System.out,argumentSources),
                              new SAMFileReaderArgumentTypeDescriptor(engine),
                              new SAMFileWriterArgumentTypeDescriptor(engine,System.out),
                              new OutputStreamArgumentTypeDescriptor(engine,System.out) );
    }

    /**
     * GATK can add arguments dynamically based on analysis type.
     *
     * @return true
     */
    @Override
    protected boolean canAddArgumentsDynamically() {
        return true;
    }

    /**
     * GATK provides the walker as an argument source.
     * @return List of walkers to load dynamically.
     */
    @Override
    protected Class[] getArgumentSources() {
        // No walker info?  No plugins.
        if (getAnalysisName() == null) return new Class[] {};

        Collection<Class> argumentSources = new ArrayList<Class>();

        Walker walker = engine.getWalkerByName(getAnalysisName());
        engine.setArguments(getArgumentCollection());
        engine.setWalker(walker);
        walker.setToolkit(engine);
        argumentSources.add(walker.getClass());

        Collection<SamRecordFilter> filters = engine.createFilters();
        for(SamRecordFilter filter: filters)
            argumentSources.add(filter.getClass());

        Class[] argumentSourcesAsArray = new Class[argumentSources.size()];
        return argumentSources.toArray(argumentSourcesAsArray);
    }

    @Override
    protected String getArgumentSourceName( Class argumentSource ) {
        return engine.getWalkerName((Class<Walker>)argumentSource);
    }

    /**
     * Unpack the bam files to be processed, given a list of files.  That list of files can
     * itself contain entries which are lists of other files to be read (note: you cannot have lists of lists of lists)
     *
     * @param argCollection the command-line arguments from which to extract the BAM file list.
     * @return a flattened list of the bam files provided
     */
    private List<SAMReaderID> unpackBAMFileList(GATKArgumentCollection argCollection) {
        List<SAMReaderID> unpackedReads = new ArrayList<SAMReaderID>();
        for( File inputFile: argCollection.samFiles ) {
            if (inputFile.getName().toLowerCase().endsWith(".list") ) {
                try {
                    for(String fileName : new XReadLines(inputFile))
                        unpackedReads.add(new SAMReaderID(new File(fileName),parser.getTags(inputFile)));
                }
                catch( FileNotFoundException ex ) {
                    throw new UserException.CouldNotReadInputFile(inputFile, "Unable to find file while unpacking reads", ex);
                }
            }
            else if(inputFile.getName().toLowerCase().endsWith(".bam")) {
                unpackedReads.add( new SAMReaderID(inputFile,parser.getTags(inputFile)) );
            }
            else if(inputFile.getName().equals("-")) {
                unpackedReads.add(new SAMReaderID(new File("/dev/stdin"),Collections.<String>emptyList()));
            }
            else {
                throw new UserException.CommandLineException(String.format("The GATK reads argument (-I) supports only BAM files with the .bam extension and lists of BAM files " +
                        "with the .list extension, but the file %s has neither extension.  Please ensure that your BAM file or list " +
                        "of BAM files is in the correct format, update the extension, and try again.",inputFile.getName()));
            }
        }
        return unpackedReads;
    }
    /**
     * Convert command-line argument representation of ROD bindings to something more easily understandable by the engine.
     * @param argCollection input arguments to the GATK.
     * @return a list of expanded, bound RODs.
     */
    private Collection<RMDTriplet> unpackRODBindings(GATKArgumentCollection argCollection) {
        Collection<RMDTriplet> rodBindings = new ArrayList<RMDTriplet>();

        for (String binding: argCollection.RODBindings) {
            if(parser.getTags(binding).size() != 2)
                throw new UserException("Invalid syntax for -B (reference-ordered data) input flag.  " +
                        "Please use the following syntax when providing reference-ordered " +
                        "data: -B:<name>,<type> <filename>.");
            // Assume that if tags are present, those tags are name and type.
            // Name is always first, followed by type.
            List<String> parameters = parser.getTags(binding);
            String name = parameters.get(0);
            String type = parameters.get(1);
            rodBindings.add(new RMDTriplet(name,type,binding));
        }

        if (argCollection.DBSNPFile != null) {
            if(argCollection.DBSNPFile.toLowerCase().contains("vcf"))
                throw new UserException("--DBSNP (-D) argument currently does not support VCF.  To use dbSNP in VCF format, please use -B:dbsnp,vcf <filename>.");
            rodBindings.add(new RMDTriplet(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, "dbsnp", argCollection.DBSNPFile));
        }

        return rodBindings;
    }


}
