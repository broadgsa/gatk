package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;

@WalkerName("LogisticRecalibration")
public class LogisticRecalibrationWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Argument(shortName="logisticParams", doc="logistic params file", required=true)
    public String logisticParamsFile;

    @Argument(shortName="outputBAM", doc="output BAM file", required=false)
    public String outputBamFile = null;

    Map<Pair<String,String>, LogisticRegressor> regressors = new HashMap<Pair<String,String>, LogisticRegressor>();
    private static Logger logger = Logger.getLogger(LogisticRecalibrationWalker.class);

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
                for ( int i = 2; i <= (vals.length-2); i++ ) {
                    Pair<Integer, Integer> ij = mapping.get(i-1);
                    try {
                        double c = Double.parseDouble(vals[i]);
                        regressor.setCoefficient(ij.first, ij.second, c);
                        System.out.printf("Setting coefficient %s,%s => %s = %f%n", readGroup, dinuc, ij, c);
                    } catch ( NumberFormatException e ) {
                        Utils.scareUser("Badly formed logistic regression header at " + vals[i] + " line: " + line );
                    }                        
                }
            }

            for ( Map.Entry<Pair<String,String>, LogisticRegressor> e : regressors.entrySet() ) {
                String readGroup = e.getKey().first;
                String dinuc = e.getKey().second;
                LogisticRegressor regressor = e.getValue();
                logger.debug(String.format("Regressor: %s,%s => %s", readGroup, dinuc, regressor));
            }

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


    public SAMRecord map(char[] ref, SAMRecord read) {
        System.out.println(String.format("Reading: %s (%s:%d-%d)",read.getReadName(),read.getReferenceName(),read.getAlignmentStart(),read.getAlignmentEnd()));              
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
        //recalQuals[numBases-1] = quals[numBases-1]; // can't change last -- no dinuc
        for ( int cycle = 1; cycle < numBases; cycle++ ) { // skip first and last base, qual already set because no dinuc
            // Take into account that previous base is the next base in terms of machine chemistry if this is a negative strand
            //int cycle = i; //read.getReadNegativeStrandFlag() ? numBases - i - 1 : i;
            String dinuc = String.format("%c%c", bases[cycle - 1], bases[cycle]);
            byte qual = quals[cycle];
            LogisticRegressor regressor = regressors.get(new Pair<String,String>(readGroup,dinuc));
            byte newQual;

            if ( regressor != null ) { // no N or some other unexpected bp in the stream
                double gamma = regressor.regress((double)cycle+1, (double)qual);
                double expGamma = Math.exp(gamma);
                double finalP = expGamma / (1+expGamma);
                newQual = QualityUtils.probToQual(1-finalP);
                //newQual = -10 * Math.round(logPOver1minusP)
                /*double POver1minusP = Math.pow(10, logPOver1minusP);
                P = POver1minusP / (1 + POver1minusP);*/
                //newQual = QualityUtils.probToQual(P);

                //newQual = (byte)Math.min(Math.round(-10*logPOver1minusP),63);
                //System.out.printf("Recal %s %d %d %d%n", dinuc, cycle, qual, newQual);
            }else{
                newQual = qual;
            }

            recalQuals[cycle] = newQual;
        }

        if (read.getReadNegativeStrandFlag())
            recalQuals = BaseUtils.reverse(quals);
        //System.out.printf("OLD: %s%n", read.format());
        read.setBaseQualities(recalQuals);
        //System.out.printf("NEW: %s%n", read.format());
        return recalRead;
    }

    public void onTraversalDone(SAMFileWriter output) {
        if ( output != null ) {
            output.close();
        }
    }

    public SAMFileWriter reduceInit() {
        if ( outputBamFile != null ) { // ! outputBamFile.equals("") ) {
            SAMFileWriterFactory fact = new SAMFileWriterFactory();
            SAMFileHeader header = this.getToolkit().getEngine().getSAMHeader();
            return fact.makeBAMWriter(header, true, new File(outputBamFile));
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