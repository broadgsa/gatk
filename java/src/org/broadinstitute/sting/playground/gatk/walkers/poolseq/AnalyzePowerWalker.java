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
import org.broadinstitute.sting.playground.utils.ReadOffsetQuad;

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
    private double[][] syzyPowerTable;
    private String[][] syzySanityCheckTable;

    @Override
    public void initialize()
    {
        super.initialize();
        syzyPowerTable = new double[8][1500000];
        syzySanityCheckTable = new String[8][1500000];
        try {
            syzyFileReader = new BufferedReader(new FileReader(pathToSyzygyFile));
            logger.error("generatePowerTable called");
            generatePowerTable(syzyFileReader);
        } catch (FileNotFoundException e) {
            String newErrMsg = "Syzygy input file " + pathToSyzygyFile + " could be incorrect. File not found.";
            throw new StingException(newErrMsg,e);
        } catch (IOException e) {
            String newErrMsg = "Syzygy input file error: could not read first line of "+pathToSyzygyFile;
            throw new StingException(newErrMsg,e);
        }

    }

    @Override
    public Pair<Integer,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext rawContext)
    {
         AlignmentContext context;
        if (super.getMinQualityScore() <= 0) {
            context = rawContext;
        } else {
            Pair<List<SAMRecord>,List<Integer>> readsFilteredByQuality = filterByQuality(rawContext.getReads(),rawContext.getOffsets(), super.getMinQualityScore());
            context = new AlignmentContext(rawContext.getLocation(),readsFilteredByQuality.getFirst(),readsFilteredByQuality.getSecond());
        }
        ReadOffsetQuad readsByDirection = PoolUtils.splitReadsByReadDirection(context.getReads(),context.getOffsets());
        Pair<Pair<List<SAMRecord>, List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> splitReads = new Pair(new Pair(readsByDirection.getFirstReads(),readsByDirection.getSecondReads()),new Pair(readsByDirection.getFirstOffsets(),readsByDirection.getSecondOffsets()));
        if ( !super.suppress_printing )
        {
            Pair<double[],byte[]> powPair = super.calculatePower(splitReads,false,context);



            int tabIndexChrom = getTableChromIndex(context.getLocation().toString());
            int tabIndexLoc = getTablePosIndex(context.getLocation().toString());
            double syzyPow = 0;
            String syzySanity = "NoLoc";
            if(tabIndexChrom >= 0 && tabIndexLoc >= 0) {
                syzyPow += syzyPowerTable[tabIndexChrom][tabIndexLoc];
                syzySanity = "s " + syzySanityCheckTable[tabIndexChrom][tabIndexLoc];

            }
            out.printf("%s|%s: %d %d %d %d %d %d %f %f %f %f %n", context.getLocation(), syzySanity, splitReads.getFirst().getFirst().size(), splitReads.getFirst().getSecond().size(),
                    context.getReads().size(), powPair.getSecond()[0], powPair.getSecond()[1], powPair.getSecond()[2],
                    powPair.getFirst()[0], powPair.getFirst()[1], powPair.getFirst()[2], syzyPow);

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

        String chromPosTar = null;
        StringTokenizer lineTokenizer = new StringTokenizer(thisLine, pfFileDelimiter);
        try {
            chromPosTar = lineTokenizer.nextToken();
            for(int j = 1; j < colOffset; j++) {
                lineTokenizer.nextToken();
            }
            String chromPos = (new StringTokenizer(chromPosTar, "\t")).nextToken();
            return new Pair((Double.valueOf(lineTokenizer.nextToken())/100.0),chromPos);
        } catch (NoSuchElementException e) {
            String errMsg = "The given column offset for the pool, " + colOffset + " exceeded the number of entries in the file " + pathToSyzygyFile;
            throw new StingException(errMsg);
        }
    }

    public String createHeaderString() {
        return (super.createHeaderString() + "  PowSyz");
    }

    //TODO: rewrite this so it's not a complete hack!!!

    public void generatePowerTable(BufferedReader syzFile) {

        System.out.println("generatePowerTable entered");
        int numreads = 0;

        try{
            while(syzFile.ready()) {
                String line = syzFile.readLine();
                if(line == null || line.equals("")) {
                    break;
                }
                StringTokenizer tok = new StringTokenizer(line, " :\t");
                String chrmNoStrWithChr = tok.nextToken();
                String chrmNoStr = chrmNoStrWithChr.substring(3);
                String locStr = tok.nextToken();
                for(int j = colOffset - 2; j > 0; j --) {
                    tok.nextToken();
                }
                numreads++;
                String syzyPowStr = tok.nextToken();
                if ( chrmNoStr == null || locStr == null || chrmNoStr.equals("") || locStr.equals("")) {
                    break;
                }
                int chrIndex = getTableChromIndex(chrmNoStr,locStr);
                int posIndex = getTablePosIndex(chrmNoStr,locStr);
                syzyPowerTable[chrIndex][posIndex] = Double.valueOf(syzyPowStr)/100.0;
                syzySanityCheckTable[chrIndex][posIndex] = "chrm" + chrmNoStr +":" + locStr;
                if( (numreads % 1000) == 0){
                    System.out.println(numreads + " reads from Syzygy file");
                }
            }
        } catch (IOException e) {
            // do nothing
            System.out.println("IOException caught");
        }
        System.out.println("Table generated.");

   }

    public int getTableChromIndex(String chromNo, String locNo) {
        switch (Integer.valueOf(chromNo)) {
            case 1: return 0;
            case 2: return 1;
            case 3:
                if ( Integer.valueOf(locNo) < 63000000 )
                    return 2;
                else
                    return 3;
            case 7: return 4;
            case 10: return 5;
            case 11: return 6;
            case 12: return 7;
            default:
                System.out.println(chromNo + " " + locNo);
                return -1;
        }
    }

    public int getTableChromIndex(String chromAndPos) {
        StringTokenizer tok = new StringTokenizer(chromAndPos,":");
        return getTableChromIndex(tok.nextToken().substring(3),tok.nextToken());
    }

    public int getTablePosIndex(String chromAndPos) {
        StringTokenizer tok = new StringTokenizer(chromAndPos,":");
        return getTablePosIndex(tok.nextToken().substring(3),tok.nextToken());
    }

    public int getTablePosIndex(String chromNo, String locNo) {
        switch ( getTableChromIndex(chromNo, locNo) ) {
            case 0: return Integer.valueOf(locNo) - 120065274;
            case 1: return Integer.valueOf(locNo) - 43311321;
            case 2: return Integer.valueOf(locNo) - 12157091;
            case 3: return Integer.valueOf(locNo) - 64000000;
            case 4: return Integer.valueOf(locNo) - 27838907;
            case 5: return Integer.valueOf(locNo) - 12003199;
            case 6: return Integer.valueOf(locNo) - 92342389;
            case 7: return Integer.valueOf(locNo) - 69809219;
            default:
                System.out.println(chromNo + " " + locNo);
                return -1;
        }
    }
}
