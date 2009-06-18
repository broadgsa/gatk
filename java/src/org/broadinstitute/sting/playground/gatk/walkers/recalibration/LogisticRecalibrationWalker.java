package org.broadinstitute.sting.playground.gatk.walkers.recalibration;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;

@WalkerName("LogisticRecalibration")
@Requires({DataSource.READS, DataSource.REFERENCE})
public class LogisticRecalibrationWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Argument(shortName="logisticParams", doc="logistic params file", required=true)
    public String logisticParamsFile;

    @Argument(shortName="outputBAM", doc="output BAM file", required=false)
    public String outputBamFile = null;

    @Argument(shortName="useCache", doc="If true, uses high-performance caching of logistic regress results.  Experimental", required=false)
    public boolean useLogisticCache = true;

    Map<Pair<String,String>, LogisticRegressor> regressors = new HashMap<Pair<String,String>, LogisticRegressor>();
    private static Logger logger = Logger.getLogger(LogisticRecalibrationWalker.class);

    // maps from [readGroup] -> [prevBase x base -> [cycle, qual, new qual]]
    HashMap<String, HashMap<String, byte[][]>> cache = new HashMap<String, HashMap<String, byte[][]>>();

    private static byte MAX_Q_SCORE = 64;


    @Argument(shortName="maxReadLen", doc="Maximum allowed read length to allow during recalibration, needed for recalibration table allocation", required=false)
    public static int maxReadLen = 125;

    public void initialize() {
        try {
            List<String> lines = new xReadLines(new File(logisticParamsFile)).readLines();
            ArrayList<Pair<Integer, Integer>> mapping = parseHeader(lines.get(0));
            for ( String line : lines.subList(1,lines.size()) ) {
                // dinuc coeff1 coeff2 ... coeff25
                String[] vals = line.split("\\s+");
                String readGroup = vals[0];
                String dinuc = vals[1];
                LogisticRegressor regressor = new LogisticRegressor(2, 4);
                regressors.put(new Pair<String,String>(readGroup,dinuc), regressor);
                System.out.printf("Vals = %s%n", Utils.join(",", vals));
                for ( int i = 2; i < vals.length; i++ ) {
                    Pair<Integer, Integer> ij = mapping.get(i-2);
                    try {
                        double c = Double.parseDouble(vals[i]);
                        regressor.setCoefficient(ij.first, ij.second, c);
                        System.out.printf("Setting coefficient %s,%s => %s = %1.12f%n", readGroup, dinuc, ij, c);
                    } catch ( NumberFormatException e ) {
                        Utils.scareUser("Badly formed logistic regression header at " + vals[i] + " line: " + line );
                    }                        
                }
            }

            if ( useLogisticCache ) System.out.printf("Building recalibration cache%n");
            
            for ( Map.Entry<Pair<String,String>, LogisticRegressor> e : regressors.entrySet() ) {
                String readGroup = e.getKey().first;
                String dinuc = e.getKey().second;
                LogisticRegressor regressor = e.getValue();
                logger.debug(String.format("Regressor: %s,%s => %s", readGroup, dinuc, regressor));

                if ( useLogisticCache ) {
                    addToLogisticCache(readGroup, dinuc, regressor);
                }
            }
            if ( useLogisticCache ) System.out.printf("done%n");            

            //System.exit(1);
        } catch ( FileNotFoundException e ) {
            Utils.scareUser("Cannot read/find logistic parameters file " + logisticParamsFile);
        }
    }

    // damn I hate parsing lines
    private ArrayList<Pair<Integer, Integer>> parseHeader(String headerLine) {
        ArrayList<Pair<Integer, Integer>> mapping = new ArrayList<Pair<Integer, Integer>>();

        String[] elts = headerLine.split("\\s+");
        if ( ! "rg".equalsIgnoreCase(elts[0]) )
            Utils.scareUser("Badly formatted Logistic regression header, upper left keyword must be rg: " + elts[0] + " line: " + headerLine);        
        if ( ! elts[1].toLowerCase().startsWith("dinuc") ) // checking only start of first field because dinuc will be followed by a version number to be checekde later
            Utils.scareUser("Badly formatted Logistic regression header, second left keyword must be dinuc: " + elts[1] + " line: " + headerLine);

        for ( int k = 2; k < elts.length; k++ ) {
            String paramStr = elts[k];
            String[] ij = paramStr.split(",");
            if ( ij.length != 2 ) {
                Utils.scareUser("Badly formed logistic regression header at " + paramStr + " line: " + headerLine);
            } else {
                try {
                    int i = Integer.parseInt(ij[0]);
                    int j = Integer.parseInt(ij[1]);
                    mapping.add(new Pair<Integer, Integer>(i,j));
                    logger.info(String.format("%d => %d/%d", k, i, j));
                } catch ( NumberFormatException e ) {
                    Utils.scareUser("Badly formed logistic regression header at " + paramStr + " line: " + headerLine );
                }
            }
        }

        return mapping;
    }

    private void addToLogisticCache(final String readGroup, final String dinuc, LogisticRegressor regressor) {
        System.out.printf("%s x %s ", readGroup, dinuc);
        byte[][] dataTable = new byte[maxReadLen][MAX_Q_SCORE];

        for ( int cycle = 1; cycle < maxReadLen; cycle++ ) {
            for ( byte qual = 0; qual < MAX_Q_SCORE; qual++ ) {
                dataTable[cycle][qual] = regressor2newQual(regressor, cycle, qual);
            }
        }

        HashMap<String, byte[][]> lookup1 = cache.containsKey(readGroup) ? cache.get(readGroup) : new HashMap<String, byte[][]>();
        lookup1.put(dinuc, dataTable);
        cache.put(readGroup, lookup1);
    }

    public SAMRecord map(char[] ref, SAMRecord read) {
        if ( useLogisticCache )
            return mapCached(ref, read);
        else
            return mapOriginal(ref, read);
    }

    private byte cache2newQual(final String readGroup, HashMap<String, byte[][]> RGcache, byte prevBase, byte base, LogisticRegressor regressor, int cycle, byte qual) {
        //System.out.printf("Lookup %s %c %c %d %d%n", readGroup, prevBase, base, cycle, qual);
        //String dinuc = String.format("%c%c", (char)prevBase, (char)base);
        byte[] bp = {prevBase, base};
        String dinuc = new String(bp);

        //byte newQualCalc = regressor2newQual(regressor, cycle, qual);
        byte[][] dataTable = RGcache.get(dinuc);

        if ( dataTable == null && prevBase != 'N' && base != 'N' )
            throw new RuntimeException(String.format("Unmapped data table at %s %s", readGroup, dinuc));
        
        byte newQualCached = dataTable != null ? dataTable[cycle][qual] : qual;
        //if ( newQualCached != newQualCalc ) {
        //    throw new RuntimeException(String.format("Inconsistent quals between the cache and calculation for RG=%s: %s %d %d : %d <> %d",
        //            readGroup, dinuc, cycle, qual, newQualCalc, newQualCached));
        //}

        return newQualCached;
    }

    public SAMRecord mapCached(char[] ref, SAMRecord read) {
        if ( read.getReadLength() > maxReadLen ) {
            throw new RuntimeException("Expectedly long read, please increase maxium read len with maxReadLen parameter: " + read.format());
        }

        SAMRecord recalRead = read;
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        byte[] recalQuals = new byte[quals.length];

        // Since we want machine direction reads not corrected positive strand reads, rev comp any negative strand reads
        if (read.getReadNegativeStrandFlag()) {
            bases = BaseUtils.simpleReverseComplement(bases);
            quals = BaseUtils.reverse(quals);
        }

        String readGroup = read.getAttribute("RG").toString();

        HashMap<String, byte[][]> RGcache = cache.get(readGroup);

        int numBases = read.getReadLength();
        recalQuals[0] = quals[0];   // can't change the first -- no dinuc

        for ( int cycle = 1; cycle < numBases; cycle++ ) { // skip first and last base, qual already set because no dinuc
            // Take into account that previous base is the next base in terms of machine chemistry if
            // this is a negative strand
            byte qual = quals[cycle];
            //LogisticRegressor regressor = getLogisticRegressor(readGroup, bases[cycle - 1], bases[cycle]);
            LogisticRegressor regressor = null;            
            byte newQual = cache2newQual(readGroup, RGcache, bases[cycle - 1], bases[cycle], regressor, cycle, qual);
            recalQuals[cycle] = newQual;
        }

        if (read.getReadNegativeStrandFlag())
            recalQuals = BaseUtils.reverse(recalQuals);
        //System.out.printf("OLD: %s%n", read.format());
        read.setBaseQualities(recalQuals);
        //System.out.printf("NEW: %s%n", read.format());
        return recalRead;
    }

    // ----------------------------------------------------------------------------------------------------
    //
    // Old-style, expensive recalibrator
    //
    // ----------------------------------------------------------------------------------------------------
    private LogisticRegressor getLogisticRegressor(final String readGroup, byte prevBase, byte base) {
        String dinuc = String.format("%c%c", (char)prevBase, (char)base);
        return regressors.get(new Pair<String,String>(readGroup,dinuc));
    }

    private byte regressor2newQual(LogisticRegressor regressor, int cycle, byte qual) {
        byte newQual = qual;
        if ( regressor != null ) { // no N or some other unexpected bp in the stream
            double gamma = regressor.regress((double)cycle+1, (double)qual);
            double expGamma = Math.exp(gamma);
            double finalP = expGamma / (1+expGamma);
            newQual = QualityUtils.probToQual(1-finalP);
        }
        return newQual;
    }

    public SAMRecord mapOriginal(char[] ref, SAMRecord read) {
        SAMRecord recalRead = read;
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        byte[] recalQuals = new byte[quals.length];

        // Since we want machine direction reads not corrected positive strand reads, rev comp any negative strand reads
        if (read.getReadNegativeStrandFlag()) {
            bases = BaseUtils.simpleReverseComplement(bases);
            quals = BaseUtils.reverse(quals);
        }

        String readGroup = read.getAttribute("RG").toString();
        int numBases = read.getReadLength();
        recalQuals[0] = quals[0];   // can't change the first -- no dinuc

        for ( int cycle = 1; cycle < numBases; cycle++ ) { // skip first and last base, qual already set because no dinuc
            // Take into account that previous base is the next base in terms of machine chemistry if
            // this is a negative strand
            byte qual = quals[cycle];
            LogisticRegressor regressor = getLogisticRegressor(readGroup, bases[cycle - 1], bases[cycle]);
            byte newQual = regressor2newQual(regressor, cycle, qual);
            recalQuals[cycle] = newQual;
        }

        if (read.getReadNegativeStrandFlag())
            recalQuals = BaseUtils.reverse(quals);
        //System.out.printf("OLD: %s%n", read.format());
        read.setBaseQualities(recalQuals);
        //System.out.printf("NEW: %s%n", read.format());
        return recalRead;
    }

//    public SAMRecord mapOriginalUnmodified(char[] ref, SAMRecord read) {
//        SAMRecord recalRead = read;
//        byte[] bases = read.getReadBases();
//        byte[] quals = read.getBaseQualities();
//        byte[] recalQuals = new byte[quals.length];
//
//        // Since we want machine direction reads not corrected positive strand reads, rev comp any negative strand reads
//        if (read.getReadNegativeStrandFlag()) {
//            bases = BaseUtils.simpleReverseComplement(bases);
//            quals = BaseUtils.reverse(quals);
//        }
//
//        String readGroup = read.getAttribute("RG").toString();
//        int numBases = read.getReadLength();
//        recalQuals[0] = quals[0];   // can't change the first -- no dinuc
//        //recalQuals[numBases-1] = quals[numBases-1]; // can't change last -- no dinuc
//        for ( int cycle = 1; cycle < numBases; cycle++ ) { // skip first and last base, qual already set because no dinuc
//            // Take into account that previous base is the next base in terms of machine chemistry if this is a negative strand
//            //int cycle = i; //read.getReadNegativeStrandFlag() ? numBases - i - 1 : i;
//            String dinuc = String.format("%c%c", bases[cycle - 1], bases[cycle]);
//            byte qual = quals[cycle];
//            LogisticRegressor regressor = regressors.get(new Pair<String,String>(readGroup,dinuc));
//            byte newQual;
//
//            if ( regressor != null ) { // no N or some other unexpected bp in the stream
//                double gamma = regressor.regress((double)cycle+1, (double)qual);
//                double expGamma = Math.exp(gamma);
//                double finalP = expGamma / (1+expGamma);
//                newQual = QualityUtils.probToQual(1-finalP);
//                //newQual = -10 * Math.round(logPOver1minusP)
//                /*double POver1minusP = Math.pow(10, logPOver1minusP);
//                P = POver1minusP / (1 + POver1minusP);*/
//                //newQual = QualityUtils.probToQual(P);
//
//                //newQual = (byte)Math.min(Math.round(-10*logPOver1minusP),63);
//                //System.out.printf("Recal %s %d %d %d%n", dinuc, cycle, qual, newQual);
//            }else{
//                newQual = qual;
//            }
//
//            recalQuals[cycle] = newQual;
//        }
//
//        if (read.getReadNegativeStrandFlag())
//            recalQuals = BaseUtils.reverse(quals);
//        //System.out.printf("OLD: %s%n", read.format());
//        read.setBaseQualities(recalQuals);
//        //System.out.printf("NEW: %s%n", read.format());
//        return recalRead;
//    }


    public void onTraversalDone(SAMFileWriter output) {
        if ( output != null ) {
            output.close();
        }
    }

    public SAMFileWriter reduceInit() {
        if ( outputBamFile != null ) { // ! outputBamFile.equals("") ) {
            SAMFileHeader header = this.getToolkit().getEngine().getSAMHeader();
            return Utils.createSAMFileWriterWithCompression(header, true, outputBamFile, getToolkit().getBAMCompression());
        }
        else {
            return null;
        }
    }

    /**
     * Summarize the error rate data.
     *
     */
    public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
        if ( output != null ) {
            output.addAlignment(read);
        } else {
            out.println(read.format());
        }

        return output;
    }
}