package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.NoSuchElementException;
import java.rmi.NoSuchObjectException;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 25, 2009
 * Time: 3:47:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class AnalyzePowerWalker extends CoverageAndPowerWalker{
    // runs CoverageAndPowerWalker except compares to Syzygy outputs in a file

    @Argument(fullName = "SyzygyOutputFile", shortName = "sof", doc="Syzygy output file to compare to", required = true)
        String pathToSyzygyFile = null;
    @Argument(fullName = "ColumnOffset", shortName = "co", doc = "Offset of column containing the power in the pf", required = true)
        int colOffset = 0;

    BufferedReader syzyFileReader;
    final String pfFileDelimiter = " ";
    boolean outOfLinesInSyzyFile = false;

    @Override
    public void initialize()
    {
        super.initialize();
        try {
            syzyFileReader = new BufferedReader(new FileReader(pathToSyzygyFile));
            syzyFileReader.readLine();
        } catch (FileNotFoundException e) {
            String newErrMsg = "Syzygy input file " + pathToSyzygyFile + " could be incorrect. File not found.";
            throw new StingException(newErrMsg,e);
        } catch (IOException e) {
            String newErrMsg = "Syzygy input file error: could not read first line of "+pathToSyzygyFile;
            throw new StingException(newErrMsg,e);
        }

    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
        if ( !super.suppress_printing )
        {
            Pair<Double,Byte> powpair = super.boostrapSamplingPowerCalc(context);

            boolean syzyFileIsReady;
            try {
                syzyFileIsReady = syzyFileReader.ready();
            }
            catch(IOException e) {
                syzyFileIsReady = false;
            }

            if(!syzyFileIsReady) {
                throw new StingException("Input file reader was not ready before an attempt to read from it.");
            } else if(!outOfLinesInSyzyFile) {
                double syzyPow = getSyzyPowFromFile();
                out.printf("%s: %d %d %f %f%n", context.getLocation(), context.getReads().size(),powpair.second,powpair.first,syzyPow);
            } else {
                out.printf("%s: %d %d %f%n", context.getLocation(), context.getReads().size(),powpair.second,powpair.first);
            }
        }

        return context.getReads().size();
    }

    public double getSyzyPowFromFile() {
        String thisLine = null;
        try {
            thisLine = syzyFileReader.readLine();
        } catch(IOException e) {
            String newErrMsg = "Ran out of lines in the syzyfile; further output of Syzygy power will be suppressed.";
            outOfLinesInSyzyFile=true;
            logger.warn(newErrMsg + " " + e.toString());
            return -1.1;
        }

        StringTokenizer lineTokenizer = new StringTokenizer(thisLine, pfFileDelimiter);
        try {
            for(int j = 0; j < colOffset; j++) {
                lineTokenizer.nextToken();
            }
            return (Double.valueOf(lineTokenizer.nextToken())/100.0);
        } catch (NoSuchElementException e) {
            String errMsg = "The given column offset for the pool, " + colOffset + " exceeded the number of entries in the file " + pathToSyzygyFile;
            throw new StingException(errMsg);
        }
    }
}
