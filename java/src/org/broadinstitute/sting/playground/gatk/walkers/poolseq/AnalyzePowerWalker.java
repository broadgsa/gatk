package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.playground.utils.PoolUtils;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.NoSuchElementException;
import java.util.List;
import java.rmi.NoSuchObjectException;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 25, 2009
 * Time: 3:47:47 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
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
            System.out.println(syzyFileReader.readLine());
        } catch (FileNotFoundException e) {
            String newErrMsg = "Syzygy input file " + pathToSyzygyFile + " could be incorrect. File not found.";
            throw new StingException(newErrMsg,e);
        } catch (IOException e) {
            String newErrMsg = "Syzygy input file error: could not read first line of "+pathToSyzygyFile;
            throw new StingException(newErrMsg,e);
        }

    }

    @Override
    public Pair<Integer,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {

        Pair<Pair<List<SAMRecord>, List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> splitReads = PoolUtils.splitReadsByReadDirection(context.getReads(),context.getOffsets());
        if ( !super.suppress_printing )
        {
            Pair<double[],byte[]> powPair = super.calculatePower(splitReads,false,context);

            boolean syzyFileIsReady=true;
            try {
                syzyFileIsReady = syzyFileReader.ready();
            }
            catch(IOException e) {
                throw new StingException("Input file reader was not ready before an attempt to read from it", e);
            }

            if(!syzyFileIsReady) {
                throw new StingException("Input file reader was not ready before an attempt to read from it, but there was no IOException");
            } else if(!outOfLinesInSyzyFile) {
                Pair<Double,String> syzyPow = getSyzyPowFromFile();
                out.printf("%s: %d %d %d %d %d %d %f %f %f %f |%s%n", context.getLocation(), splitReads.getFirst().getFirst().size(), splitReads.getFirst().getSecond().size(),
                        context.getReads().size(), powPair.getSecond()[0], powPair.getSecond()[1], powPair.getSecond()[2],
                        powPair.getFirst()[0], powPair.getFirst()[1], powPair.getFirst()[2], syzyPow.getFirst(), syzyPow.getSecond());
            } else {
                out.printf("%s: $d %d %d %d %d %d %d %f %f %f%n", context.getLocation(), splitReads.getFirst().getFirst().size(), splitReads.getFirst().getSecond().size(),
                        context.getReads().size(), powPair.getSecond()[0], powPair.getSecond()[1], powPair.getSecond()[2],
                        powPair.getFirst()[0], powPair.getFirst()[1], powPair.getFirst()[2]);
            }
        }

        return new Pair(splitReads.getFirst().getFirst().size(), splitReads.getFirst().getFirst().size());
    }

    public Pair<Double,String> getSyzyPowFromFile() {
        String thisLine = null;
        try {
            thisLine = syzyFileReader.readLine();
        } catch(IOException e) {
            String newErrMsg = "Ran out of lines in the syzyfile; further output of Syzygy power will be suppressed.";
            outOfLinesInSyzyFile=true;
            logger.warn(newErrMsg + " " + e.toString());
            return new Pair(-1.1, "Printing Stops Here");
        }

        String chromPos = null;
        StringTokenizer lineTokenizer = new StringTokenizer(thisLine, pfFileDelimiter);
        try {
            chromPos = lineTokenizer.nextToken();
            for(int j = 1; j < colOffset; j++) {
                lineTokenizer.nextToken();
            }
            return new Pair((Double.valueOf(lineTokenizer.nextToken())/100.0),chromPos);
        } catch (NoSuchElementException e) {
            String errMsg = "The given column offset for the pool, " + colOffset + " exceeded the number of entries in the file " + pathToSyzygyFile;
            throw new StingException(errMsg);
        }
    }

    public String createHeaderString() {
        return (super.createHeaderString() + "  PowSyz");
    }
}
