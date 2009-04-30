package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
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
    @Argument(shortName="logisticParams", required=true)
    public String logisticParamsFile;

    @Argument(shortName="outputBAM", required=false, defaultValue="")
    public String outputFile;

    HashMap<String, LogisticRegressor> regressors = new HashMap<String, LogisticRegressor>();
    private static Logger logger = Logger.getLogger(LogisticRecalibrationWalker.class);

    public void initialize() {
        try {
            List<String> lines = new xReadLines(new File(logisticParamsFile)).readLines();
            ArrayList<Pair<Integer, Integer>> mapping = parseHeader(lines.get(0));
            for ( String line : lines.subList(1,lines.size()) ) {
                // dinuc coeff1 coeff2 ... coeff25
                String[] vals = line.split("\\s+");
                String dinuc = vals[0];
                LogisticRegressor regressor = new LogisticRegressor(2, 4);
                regressors.put(dinuc, regressor);
                System.out.printf("Vals = %s%n", Utils.join(",", vals));
                for ( int i = 1; i < vals.length; i++ ) {
                    Pair<Integer, Integer> ij = mapping.get(i-1);
                    try {
                        double c = Double.parseDouble(vals[i]);
                        regressor.setCoefficient(ij.first, ij.second, c);
                        System.out.printf("Setting coefficient %s => %s = %f%n", dinuc, ij, c);
                    } catch ( NumberFormatException e ) {
                        Utils.scareUser("Badly formed logistic regression header at " + vals[i] + " line: " + line );
                    }                        
                }
            }

            for ( Map.Entry<String, LogisticRegressor> e : regressors.entrySet() ) {
                String dinuc = e.getKey();
                LogisticRegressor regressor = e.getValue();
                logger.debug(String.format("Regressor: %s => %s", dinuc, regressor));
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
        if ( ! elts[0].toLowerCase().equals("dinuc") )
            Utils.scareUser("Badly formatted Logistic regression header, upper left keyword must be dinuc: " + elts[0] + " line: " + headerLine);

        for ( int k = 1; k < elts.length; k++ ) {
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


    public SAMRecord map(LocusContext context, SAMRecord read) {
        SAMRecord recalRead = read;
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        byte[] recalQuals = new byte[quals.length];

        int numBases = read.getReadLength();
        recalQuals[0] = quals[0];   // can't change the first -- no dinuc
        recalQuals[numBases-1] = quals[numBases-1]; // can't change last -- no dinuc
        for ( int i = 1; i < numBases-1; i++ ) { // skip first and last base, qual already set because no dinuc
            // Take into account that previous base is the next base in terms of machine chemistry if this is a negative strand
            int cycle = read.getReadNegativeStrandFlag() ? numBases - i - 1 : i;
            String dinuc = String.format("%c%c", bases[i  + (read.getReadNegativeStrandFlag() ? 1 : -1)], bases[i]);
            byte qual = quals[i];
            //System.out.printf("dinuc %c %c%n", bases[i-1], bases[i]);
            LogisticRegressor regressor = regressors.get(dinuc);
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
                //System.out.printf("Recal %s %d %d => %f => %f leads to %d%n", dinuc, cycle, qual, logPOver1minusP, P, newQual);
            }else{
                newQual = qual;
            }

            recalQuals[i] = newQual;
        }

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
        if ( outputFile != null ) { // ! outputFile.equals("") ) {
            SAMFileWriterFactory fact = new SAMFileWriterFactory();
            SAMFileHeader header = this.getToolkit().getSamReader().getFileHeader();
            return fact.makeBAMWriter(header, true, new File(outputFile));
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