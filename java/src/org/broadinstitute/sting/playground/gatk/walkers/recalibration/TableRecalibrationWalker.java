package org.broadinstitute.sting.playground.gatk.walkers.recalibration;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.playground.gatk.walkers.recalibration.RecalData;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;

@WalkerName("TableRecalibration")
public class TableRecalibrationWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Argument(shortName="params", doc="CountCovariates params file", required=true)
    public String paramsFile;

    @Argument(shortName="outputBAM", doc="output BAM file", required=false)
    public String outputBamFile = null;

    private static Logger logger = Logger.getLogger(TableRecalibrationWalker.class);

    private final static boolean DEBUG = false;

    // maps from [readGroup] -> [prevBase x base -> [cycle, qual, new qual]]
    HashMap<String, RecalMapping> cache = new HashMap<String, RecalMapping>();

    @Argument(shortName="serial", doc="", required=false)
    boolean serialRecalibration = false;

    public void initialize() {
        try {
            System.out.printf("Reading data...%n");
            List<RecalData> data = new ArrayList<RecalData>();
            List<String> lines = new xReadLines(new File(paramsFile)).readLines();
            for ( String line : lines ) {
                // rg,dn,logitQ,pos,indicator,count
                // SRR007069,AA,28,1,0,2
                data.add(RecalData.fromCSVString(line));
            }
            initializeCache(data);
        } catch ( FileNotFoundException e ) {
            Utils.scareUser("Cannot read/find parameters file " + paramsFile);
        }
    }

    private void initializeCache(List<RecalData> data) {
        Set<String> readGroups = new HashSet<String>();
        Set<String> dinucs = new HashSet<String>();
        int maxPos = -1;
        int maxQReported = -1;

        logger.info(String.format("No. params       : %d", data.size()));
        // get the maximum data from the file
        for ( RecalData datum : data ) {
            readGroups.add(datum.readGroup);
            dinucs.add(datum.dinuc);
            maxPos = Math.max(maxPos, datum.pos);
            maxQReported = Math.max(maxQReported, datum.qual);
        }
        logger.info(String.format("Read groups      : %d %s", readGroups.size(), readGroups.toString()));
        logger.info(String.format("Dinucs           : %d %s", dinucs.size(), dinucs.toString()));
        logger.info(String.format("Max pos          : %d", maxPos));
        logger.info(String.format("Max Q reported   : %d", maxQReported));

        // initialize the data structure
        HashMap<String, RecalDataManager> managers = new HashMap<String, RecalDataManager>();
        for ( String readGroup : readGroups ) {
            RecalDataManager manager = new RecalDataManager(readGroup,  maxPos, maxQReported, dinucs.size(), true, true);
            managers.put(readGroup, manager);
        }

        // fill in the manager structure
        for ( RecalData datum : data ) {
            managers.get(datum.readGroup).addDatum(datum);
        }

        // fill in the table with mapping objects
        for ( String readGroup : readGroups ) {
            RecalDataManager manager = managers.get(readGroup);
            RecalMapping mapper = null;
            if ( serialRecalibration )
                mapper = new SerialRecalMapping(manager, dinucs, maxPos, maxQReported);
            else
                mapper = new CombinatorialRecalMapping(manager, dinucs, maxPos, maxQReported);
            cache.put(readGroup, mapper);
        }
    }

    public SAMRecord map(char[] ref, SAMRecord read) {
        //if ( read.getReadLength() > maxReadLen ) {
        //    throw new RuntimeException("Expectedly long read, please increase maxium read len with maxReadLen parameter: " + read.format());
        //}

        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        byte[] recalQuals = new byte[quals.length];

        // Since we want machine direction reads not corrected positive strand reads, rev comp any negative strand reads
        if (read.getReadNegativeStrandFlag()) {
            bases = BaseUtils.simpleReverseComplement(bases);
            quals = BaseUtils.reverse(quals);
        }

        String readGroup = read.getAttribute("RG").toString();

        RecalMapping mapper = cache.get(readGroup);

        int numBases = read.getReadLength();
        recalQuals[0] = quals[0];   // can't change the first -- no dinuc

        for ( int cycle = 1; cycle < numBases; cycle++ ) { // skip first and last base, qual already set because no dinuc
            // Take into account that previous base is the next base in terms of machine chemistry if
            // this is a negative strand
            byte qual = quals[cycle];
            byte newQual = mapper.getNewQual(readGroup, bases[cycle - 1], bases[cycle], cycle, qual);
            recalQuals[cycle] = newQual;
            //System.out.printf("Mapping %d => %d%n", qual, newQual);
        }

        if (read.getReadNegativeStrandFlag())
            recalQuals = BaseUtils.reverse(quals);
        //System.out.printf("OLD: %s%n", read.format());
        read.setBaseQualities(recalQuals);
        //System.out.printf("NEW: %s%n", read.format());
        return read;
    }

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

interface RecalMapping {
    public byte getNewQual(final String readGroup, byte prevBase, byte base, int cycle, byte qual);
}

class CombinatorialRecalMapping implements RecalMapping {
    HashMap<String, byte[][]> cache = new HashMap<String, byte[][]>();

    public CombinatorialRecalMapping(RecalDataManager manager, Set<String> dinucs, int maxPos, int maxQReported ) {
        // initialize the data structure
        for ( String dinuc : dinucs ) {
            byte[][] table = new byte[maxPos+1][maxQReported+1];
            cache.put(dinuc, table);
        }

        for ( RecalData datum : manager.getAll() ) {
            //System.out.printf("Adding datum %s%n", datum);
            byte [][] table = cache.get(datum.dinuc);
            if ( table[datum.pos][datum.qual] != 0 )
                throw new RuntimeException(String.format("Duplicate entry discovered: %s", datum));
            //table[datum.pos][datum.qual] = (byte)(1 + datum.empiricalQualByte());
            table[datum.pos][datum.qual] = datum.empiricalQualByte();
        }
    }

    public byte getNewQual(final String readGroup, byte prevBase, byte base, int cycle, byte qual) {
        //System.out.printf("Lookup RG=%s prevBase=%c base=%c cycle=%d qual=%d%n", readGroup, prevBase, base, cycle, qual);
        //String dinuc = String.format("%c%c", (char)prevBase, (char)base);
        byte[] bp = {prevBase, base};
        String dinuc = new String(bp);
        byte[][] dataTable = cache.get(dinuc);

        if ( dataTable == null && prevBase != 'N' && base != 'N' )
            throw new RuntimeException(String.format("Unmapped data table at %s %s", readGroup, dinuc));

        return dataTable != null && cycle < dataTable.length ? dataTable[cycle][qual] : qual;
    }
}

class SerialRecalMapping implements RecalMapping {
    // mapping from dinuc x Q => new Q
    HashMap<String, byte[]> mappingByDinuc;

    // mapping from pos x Q => new Q
    byte[][] mappingByPos;

    public SerialRecalMapping(RecalDataManager manager, Set<String> dinucs, int maxPos, int maxQReported ) {
        mappingByDinuc = new HashMap<String, byte[]>();
        for ( String dinuc : dinucs ) {
            byte[] table = new byte[maxQReported+1];
            mappingByDinuc.put(dinuc, table);
        }
        for ( RecalData datum : manager.getDataByDinuc() ) {
            //System.out.printf("Adding datum %s%n", datum);
            if ( mappingByDinuc.get(datum.dinuc).length <= datum.qual ) {
                throw new RuntimeException(String.format("Unexpectedly massive Q score of %d found, calculated max was %d%n", maxQReported, datum.qual));
            }
            mappingByDinuc.get(datum.dinuc)[datum.qual] = datum.empiricalQualByte();
        }

        // initialize the mapping by position
        mappingByPos = new byte[maxPos+1][maxQReported+1];
        for ( RecalData datum : manager.getDataByPos() ) {
            //System.out.printf("Adding datum %s%n", datum);
            mappingByPos[datum.pos][datum.qual] = datum.empiricalQualByte();
        }
    }

    public byte getNewQual(final String readGroup, byte prevBase, byte base, int cycle, byte qual) {
        //System.out.printf("Lookup RG=%s prevBase=%c base=%c cycle=%d qual=%d%n", readGroup, prevBase, base, cycle, qual);
        //String dinuc = String.format("%c%c", (char)prevBase, (char)base);
        byte[] bp = {prevBase, base};
        String dinuc = new String(bp);

        byte newQualFromDinuc = 0;
        byte newQualFromPos = cycle > 0 && cycle < mappingByPos.length ? mappingByPos[cycle][qual] : qual;
        byte newQual = newQualFromPos;
        if ( prevBase != 'N' && base != 'N' ) {
            // if the qual got mapped too high, assume it's the best we've seen for recalibration purposes
            int newQualIndex = newQual < mappingByDinuc.get(dinuc).length ? newQual : mappingByDinuc.get(dinuc).length - 1;
            newQualFromDinuc = mappingByDinuc.get(dinuc)[newQualIndex];
            newQual = newQualFromDinuc;
        }

        //System.out.printf("Lookup RG=%s prevBase=%c base=%c cycle=%d qual=%d => %d => %d => %d%n",
        //        readGroup, prevBase, base, cycle, qual, newQualFromPos, newQualFromDinuc, newQual);


        return newQual;
    }
}