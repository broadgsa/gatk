/*
 * Copyright (c) 2009 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.File;
import java.io.FileNotFoundException;

@WalkerName("TableRecalibration")
@Requires({DataSource.READS, DataSource.REFERENCE})
public class TableRecalibrationWalker extends ReadWalker<SAMRecord, SAMFileWriter> {
    @Argument(shortName="params", doc="CountCovariates params file", required=true)
    public String paramsFile;

    @Argument(shortName="outputBAM", doc="output BAM file", required=false)
    public String outputBamFile = null;

    @Argument(shortName="rawQempirical", doc="If provide, we will use raw Qempirical scores calculated from the # mismatches and # bases, rather than the more conservative estimate of # mismatches + 1 / # bases + 1", required=false)
    public boolean rawQempirical = false;

    private static Logger logger = Logger.getLogger(TableRecalibrationWalker.class);

    private final static boolean DEBUG = false;

    // maps from [readGroup] -> [prevBase x base -> [cycle, qual, new qual]]
    HashMap<String, RecalMapping> cache = new HashMap<String, RecalMapping>();

    public enum RecalibrationMode {
        COMBINATORIAL,
        JAFFE,
        BY_POS_ONLY,
        BY_DINUC_ONLY,
        SEQUENTIAL,
        ERROR
    }
    
    @Argument(shortName="mode", doc="", required=false)
    public String modeString = RecalibrationMode.SEQUENTIAL.toString();
    public RecalibrationMode mode = RecalibrationMode.ERROR;

    private static Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private static Pattern COLLAPSED_POS_PATTERN = Pattern.compile("^#\\s+collapsed_pos\\s+(\\w+)");
    private static Pattern COLLAPSED_DINUC_PATTERN = Pattern.compile("^#\\s+collapsed_dinuc\\s+(\\w+)");
    private static Pattern HEADER_PATTERN = Pattern.compile("^rg.*");

    public void initialize() {
        //
        // crap hack until Enum arg types are supported
        //
        for ( RecalibrationMode potential : RecalibrationMode.values() ) {
            if ( potential.toString().equals(modeString)) {
                mode = potential;
                break;
            }
        }
        if ( mode == RecalibrationMode.ERROR )
            throw new RuntimeException("Unknown mode requested: " + modeString);


        //
        // Walk over the data file, parsing it into recalibration data
        //
        int lineNumber = 0;
        try {
            System.out.printf("Reading data...%n");
            List<RecalData> data = new ArrayList<RecalData>();
            boolean collapsedPos = false;
            boolean collapsedDinuc = false;

            List<String> lines = new xReadLines(new File(paramsFile)).readLines();
            for ( String line : lines ) {
                lineNumber++;
                if ( HEADER_PATTERN.matcher(line).matches() )
                    continue;
                if ( COMMENT_PATTERN.matcher(line).matches() ) {
                    collapsedPos = parseCommentLine(COLLAPSED_POS_PATTERN, line, collapsedPos);
                    collapsedDinuc = parseCommentLine(COLLAPSED_DINUC_PATTERN, line, collapsedDinuc);
                }
                else {
                    data.add(RecalData.fromCSVString(line));
                }
            }
            initializeCache(data, collapsedPos, collapsedDinuc);
        } catch ( FileNotFoundException e ) {
            Utils.scareUser("Cannot read/find parameters file " + paramsFile);
        } catch ( NumberFormatException e ) {
            throw new RuntimeException("Error parsing recalibration data at line " + lineNumber, e);
        }
    }

    private boolean parseCommentLine(Pattern pat, String line, boolean flag) {
        Matcher m = pat.matcher(line);
        if ( m.matches() ) {
            flag = Boolean.parseBoolean(m.group(1));
        }

        return flag;
    }

    private void initializeCache(List<RecalData> data, boolean collapsedPos, boolean collapsedDinuc ) {
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

        // check the mode
        switch ( mode ) {
            case JAFFE:
                collapsedPos = collapsedDinuc = true;
                break;
            case BY_POS_ONLY:
                if ( collapsedPos )
                    throw new RuntimeException(String.format("Cannot perform position_only recalibration -- data is already partially collapsed by pos=%b and dinuc=%b", collapsedPos, collapsedDinuc));
                collapsedPos = true;
                break;
            case BY_DINUC_ONLY:
                if ( collapsedDinuc )
                    throw new RuntimeException(String.format("Cannot perform dinuc_only recalibration -- data is already partially collapsed by pos=%b and dinuc=%b", collapsedPos, collapsedDinuc));
                collapsedDinuc = true;
                break;
            case COMBINATORIAL:
                if ( collapsedPos || collapsedDinuc )
                    throw new RuntimeException(String.format("Cannot perform combinatorial recalibration -- data is already collapsed by pos=%b and dinuc=%b", collapsedPos, collapsedDinuc));
            case SEQUENTIAL:
                if ( collapsedPos || collapsedDinuc )
                    throw new RuntimeException(String.format("Cannot perform sequential recalibration -- data is already collapsed by pos=%b and dinuc=%b", collapsedPos, collapsedDinuc));
        }

        //
        // initialize the data structures
        //
        HashMap<String, RecalDataManager> managers = new HashMap<String, RecalDataManager>();
        for ( String readGroup : readGroups ) {
            // create a manager for each read group
            RecalDataManager manager = new RecalDataManager(readGroup, ! collapsedPos, ! collapsedDinuc);
            managers.put(readGroup, manager);
        }

        // add the data to their appropriate readGroup managers
        for ( RecalData datum : data ) {
            managers.get(datum.readGroup).addDatum(datum);
        }

        // fill in the table with mapping objects
        for ( String readGroup : readGroups ) {
            RecalDataManager manager = managers.get(readGroup);
            RecalMapping mapper = initializeMapper(manager, rawQempirical, dinucs, maxPos, maxQReported);
            logger.info(String.format("Creating mapper for %s of class %s", readGroup, mapper));            
            cache.put(readGroup, mapper);
        }
    }

    /**
     * Returns a new RecalMapping object appropriate for the requested mode in this.mode for the RecalData
     * in manager.  Uses the dinucs set dinucs, maxPos, and maxQReported to build efficient lookup structures
     * for mapping from old qual and features -> new qual
     *
     * @param manager
     * @param dinucs
     * @param maxPos
     * @param maxQReported
     * @return
     */
    private RecalMapping initializeMapper( RecalDataManager manager, final boolean rawQempirical,
                                           Set<String> dinucs, int maxPos, int maxQReported ) {
        switch ( mode ) {
            case COMBINATORIAL:
            case JAFFE:
                return new CombinatorialRecalMapping(manager, rawQempirical, dinucs, maxPos, maxQReported);
            case SEQUENTIAL:
                return new SerialRecalMapping(manager, rawQempirical, dinucs, maxPos, maxQReported);
            default:
                throw new RuntimeException("Unimplemented recalibration mode: " + mode);
        }
    }

    public SAMRecord map(char[] ref, SAMRecord read) {
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();

        // Since we want machine direction reads not corrected positive strand reads, rev comp any negative strand reads
        if (read.getReadNegativeStrandFlag()) {
            bases = BaseUtils.simpleReverseComplement(bases);
            quals = BaseUtils.reverse(quals);
        }

        byte[] recalQuals = recalibrateBasesAndQuals(read.getAttribute("RG").toString(), bases, quals);

        if (read.getReadNegativeStrandFlag())               // reverse the quals for the neg strand read
            recalQuals = BaseUtils.reverse(recalQuals);
        read.setBaseQualities(recalQuals);
        return read;
    }

    /**
     * Workhorse routine.  Given a read group and an array of bases and quals, returns a new set of recalibrated
     * qualities for each base/qual in bases and quals.  Uses the RecalMapping object associated with readGroup.
     *
     * @param readGroup
     * @param bases
     * @param quals
     * @return
     */
    public byte[] recalibrateBasesAndQuals(final String readGroup, byte[] bases, byte[] quals) {
        byte[] recalQuals = new byte[quals.length];
        RecalMapping mapper = cache.get(readGroup);

        recalQuals[0] = quals[0];               // can't change the first -- no dinuc
        for ( int cycle = 1; cycle < bases.length; cycle++ ) { // skip first and last base, qual already set because no dinuc
            byte qual = quals[cycle];
            byte newQual = mapper.getNewQual(readGroup, bases[cycle - 1], bases[cycle], cycle, qual);
            recalQuals[cycle] = newQual;
            //System.out.printf("Mapping %d => %d%n", qual, newQual);
        }

        return recalQuals;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Standard I/O routines
    //
    // --------------------------------------------------------------------------------------------------------------
    public void onTraversalDone(SAMFileWriter output) {
        if ( output != null ) {
            output.close();
        }
    }

    public SAMFileWriter reduceInit() {
        if ( outputBamFile != null ) {
            SAMFileHeader header = this.getToolkit().getEngine().getSAMHeader();
            return Utils.createSAMFileWriterWithCompression(header, true, outputBamFile, getToolkit().getBAMCompression());
        }
        else {
            return null;
        }
    }

    /**
     * Write out the read
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
    ArrayList<byte[][]> cache;
    RecalDataManager manager;

    public CombinatorialRecalMapping(RecalDataManager manager, final boolean useRawQempirical,
                                     Set<String> dinucs, int maxPos, int maxQReported ) {
        this.manager = manager;

        // initialize the data structure
        cache = new ArrayList<byte[][]>(RecalData.NDINUCS);
        for ( String dinuc : dinucs ) {
            cache.add(new byte[maxPos+1][maxQReported+1]);
        }

        for ( RecalData datum : manager.getAll() ) {
            //System.out.printf("Adding datum %s%n", datum);
            byte [][] table = cache.get(this.manager.getDinucIndex(datum.dinuc));
            int pos = manager.canonicalPos(datum.pos);
            if ( table[pos][datum.qual] != 0 )
                throw new RuntimeException(String.format("Duplicate entry discovered: %s", datum));
            //table[datum.pos][datum.qual] = (byte)(1 + datum.empiricalQualByte());
            table[pos][datum.qual] = datum.empiricalQualByte(useRawQempirical);
            //System.out.printf("Binding %d %d => %d%n", pos, datum.qual, datum.empiricalQualByte(useRawQempirical));
        }
    }

    public byte getNewQual(final String readGroup, byte prevBase, byte base, int cycle, byte qual) {
        //String dinuc = String.format("%c%c", (char)prevBase, (char)base);
        //if ( qual == 2 )
        //    System.out.printf("Qual = 2%n");

        int pos = manager.canonicalPos(cycle);
        int index = this.manager.getDinucIndex(prevBase, base);
        byte[][] dataTable = index == -1 ? null : cache.get(index);

        if ( dataTable == null && prevBase != 'N' && base != 'N' )
            throw new RuntimeException(String.format("Unmapped data table at %s %c%c", readGroup, (char)prevBase, (char)base));

        byte result = dataTable != null && pos < dataTable.length ? dataTable[pos][qual] : qual;

        //if ( result == 2 )
        //    System.out.printf("Lookup RG=%s dinuc=%s cycle=%d pos=%d qual=%d datatable=%s / %d => %d%n",
        //            readGroup, dinuc, cycle, pos, qual, dataTable, dataTable.length, result);

        return result;        
    }
}


/**
 * Implements a serial recalibration of the reads using the combinational table.  First,
 * we perform a positional recalibration, and then a subsequent dinuc correction.
 *
 * Given the full recalibration table, we perform the following preprocessing steps:
 *
 *   - calculate the global quality score shift across all data [DeltaQ]
 *   - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
 *      -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
 *   - The final shift equation is:
 *
 *       Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc)
 */
class SerialRecalMapping implements RecalMapping {
    private double globalDeltaQ = 0.0;
    private double[][] deltaQPosMap, deltaQDinucMap;
    double [] deltaQualMap;
    RecalData [][] qPosSupports, qDinucSupports;

    CombinatorialRecalMapping combiMap;
    RecalDataManager manager;

    String dinucToLookAt = null; // "CC";
    int posToLookAt = 0;
    int qualToLookAt = 25;

    public SerialRecalMapping(RecalDataManager manager, final boolean useRawQempirical,
                              Set<String> dinucs, int maxPos, int maxQReported ) {
        //combiMap = new CombinatorialRecalMapping(manager, dinucs, maxPos, maxQReported );
        this.manager = manager;

        // initialize the data structure
        {
            RecalData datum = new RecalData(0, 0, manager.readGroup, "**").inc(manager.getAll());
            double aggregrateQreported = RecalData.combinedQreported(manager.getAll());
            globalDeltaQ = datum.empiricalQualDouble(useRawQempirical) - aggregrateQreported;
            System.out.printf("Global quality score shift is %.2f - %.2f = %.2f%n", datum.empiricalQualDouble(useRawQempirical), aggregrateQreported, globalDeltaQ);
        }

        for ( RecalData datum : manager.getAll() ) {
            if ( printStateP(datum) )
                System.out.printf("PrintValue: %s%n", datum);
        }

        // Jaffe-style
        deltaQualMap = new double[maxQReported+1];
        for ( RecalData datum : RecalData.sort(manager.combine(true, false, true)) ) {
            deltaQualMap[datum.qual] = datum.empiricalQualDouble(useRawQempirical) - datum.qual - globalDeltaQ;
            System.out.printf("%s => %s%n", datum, deltaQualMap[datum.qual]);            
        }

        // calculate the delta Q pos array
        deltaQPosMap = new double[maxPos+1][maxQReported+1];
        qPosSupports = new RecalData[maxPos+1][maxQReported+1];
        for ( RecalData datumAtPosQual : manager.combineDinucs() ) {
            double offset = globalDeltaQ + deltaQualMap[datumAtPosQual.qual];
            updateCache(qPosSupports, datumAtPosQual, useRawQempirical, deltaQPosMap, datumAtPosQual.pos, datumAtPosQual.qual, offset);
        }

        // calculate the delta Q dinuc array
        deltaQDinucMap = new double[dinucs.size()+1][maxQReported+1];
        qDinucSupports = new RecalData[dinucs.size()+1][maxQReported+1];
        for ( RecalData datumAtDinucQual : manager.combineCycles() ) {
            double offset = globalDeltaQ + deltaQualMap[datumAtDinucQual.qual];
            updateCache(qDinucSupports, datumAtDinucQual, useRawQempirical, deltaQDinucMap, datumAtDinucQual.getDinucIndex(), datumAtDinucQual.qual, offset);
        }

        validateTables(maxPos, maxQReported, dinucs.size(), deltaQPosMap, deltaQDinucMap, useRawQempirical, manager.getAll());
    }


    private void validateTables(int maxPos, int maxQReported, int nDinucs,
                                double[][] deltaQPosMap, double[][] deltaQDinucMap,
                                final boolean useRawQempirical, List<RecalData> data ) {
        for ( int i = 0; i < maxPos; i++ ) {
            for ( int j = 0; j < maxQReported; j++ ) {
                if ( printStateP(i, null, j) )
                    System.out.printf("Mapping: pos=%d qual=%2d delta=%.2f based on %s%n",
                            i, j, deltaQPosMap[i][j], qPosSupports[i][j]);
            }
        }

        for ( int i = 0; i < nDinucs; i++ ) {
            for ( int j = 0; j < maxQReported; j++ ) {
                String dinuc = RecalData.dinucIndex2bases(i);
                if ( printStateP(0, dinuc, j ) )
                    System.out.printf("Mapping: dinuc=%s qual=%2d delta=%.2f based on %s%n",
                            dinuc, j, deltaQDinucMap[i][j], qDinucSupports[i][j]);
            }
        }

        for ( RecalData datum : RecalData.sort(data) ) {
            byte newQual = getNewQual(datum.readGroup, (byte)datum.dinuc.charAt(0), (byte)datum.dinuc.charAt(1), datum.pos, (byte)datum.qual);
            if ( printStateP( datum ) ) {
                System.out.printf("Serial mapping %s => %d => %d vs. %d => delta = %d%n",
                        datum, newQual, datum.empiricalQualByte(useRawQempirical), datum.empiricalQualByte(), newQual - datum.empiricalQualByte(useRawQempirical));
            }
        }
    }

    private void updateCache( RecalData[][] supports, RecalData datum, final boolean useRawQempirical,
                              double[][] table, int i, int j, double meanQ ) {
        if ( table[i][j] != 0 )
            throw new RuntimeException(String.format("Duplicate entry discovered: %s", datum));
        double deltaQ = datum.empiricalQualDouble(useRawQempirical) - datum.qual - meanQ;
        table[i][j] = deltaQ;
        supports[i][j] = datum;
    }

    private boolean printStateP( int cycle, String dinuc, int qual ) {
        return posToLookAt != -1 &&
               ( cycle == 0 || posToLookAt == 0 || cycle == posToLookAt ) &&
               ( dinucToLookAt == null || dinuc == null || manager.getDinucIndex(dinuc) == -1 || dinucToLookAt.equals(dinuc)) &&
               ( qualToLookAt == 0 || qual == qualToLookAt );
    }

    private boolean printStateP( RecalData datum ) {
        return printStateP(datum.pos, datum.dinuc, datum.qual);
    }

    public byte getNewQual(final String readGroup, byte prevBase, byte base, int cycle, byte qual) {
        int pos = manager.canonicalPos(cycle);
        int index = this.manager.getDinucIndex(prevBase, base);

        double deltaQual = deltaQualMap[qual];
        double deltaQPos = pos >= deltaQPosMap.length ? 0.0 : deltaQPosMap[pos][qual];  // == length indices no data for last base
        double deltaQDinuc = index == -1 ? 0.0 : deltaQDinucMap[index][qual];           // -1 indices N in prev or current base
        double newQual = qual + globalDeltaQ + deltaQual + deltaQPos + deltaQDinuc;
        byte newQualByte = QualityUtils.boundQual((int)Math.round(newQual), QualityUtils.MAX_REASONABLE_Q_SCORE);

        //if ( printStateP(pos, RecalData.dinucIndex2bases(index), qual) )
        //    System.out.printf("%s %c%c %d %d => %d + %.2f + %.2f + %.2f + %.2f = %d%n",
        //            readGroup, prevBase, base, cycle, qual,
        //           qual, globalDeltaQ, deltaQual, deltaQPos, deltaQDinuc,
        //            newQualByte);


        //byte combiQualByte = combiMap.getNewQual(readGroup, prevBase, base, cycle, qual);
        //System.out.printf("%s %c%c %d %d => %d + %.2f + %.2f + %.2f = %d vs %d (%d)%n",
        //        readGroup, prevBase, base, cycle, qual,
        //        qual, globalDeltaQ, deltaQPos, deltaQDinuc,
        //        newQualByte, combiQualByte, newQualByte - combiQualByte);

        if ( newQualByte <= 0 && newQualByte >= QualityUtils.MAX_REASONABLE_Q_SCORE )
            throw new RuntimeException(String.format("Illegal base quality score calculated: %s %c%c %d %d => %d + %.2f + %.2f + %.2f = %d",
                    readGroup, prevBase, base, cycle, qual, qual, globalDeltaQ, deltaQPos, deltaQDinuc, newQualByte));

        return newQualByte;
    }
}