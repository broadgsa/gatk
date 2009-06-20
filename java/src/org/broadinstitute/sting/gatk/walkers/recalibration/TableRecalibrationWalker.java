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

    private static Logger logger = Logger.getLogger(TableRecalibrationWalker.class);

    private final static boolean DEBUG = false;

    // maps from [readGroup] -> [prevBase x base -> [cycle, qual, new qual]]
    HashMap<String, RecalMapping> cache = new HashMap<String, RecalMapping>();

    //@Argument(shortName="serial", doc="", required=false)
    //boolean serialRecalibration = false;

    private static Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private static Pattern COLLAPSED_POS_PATTERN = Pattern.compile("^#\\s+collapsed_pos\\s+(\\w+)");
    private static Pattern COLLAPSED_DINUC_PATTERN = Pattern.compile("^#\\s+collapsed_dinuc\\s+(\\w+)");
    private static Pattern HEADER_PATTERN = Pattern.compile("^rg.*");

    public void initialize() {
        try {
            System.out.printf("Reading data...%n");
            List<RecalData> data = new ArrayList<RecalData>();
            boolean collapsedPos = false;
            boolean collapsedDinuc = false;

            List<String> lines = new xReadLines(new File(paramsFile)).readLines();
            for ( String line : lines ) {
                //System.out.printf("Reading line %s%n", line);
                if ( HEADER_PATTERN.matcher(line).matches() )
                    continue;
                if ( COMMENT_PATTERN.matcher(line).matches() ) {
                    collapsedPos = parseCommentLine(COLLAPSED_POS_PATTERN, line, collapsedPos);
                    collapsedDinuc = parseCommentLine(COLLAPSED_DINUC_PATTERN, line, collapsedDinuc);
                    //System.out.printf("Collapsed %b %b%n", collapsedPos, collapsedDinuc);
                }
                else {
                    data.add(RecalData.fromCSVString(line));
                }
            }
            initializeCache(data, collapsedPos, collapsedDinuc);
        } catch ( FileNotFoundException e ) {
            Utils.scareUser("Cannot read/find parameters file " + paramsFile);
        }
    }

    private boolean parseCommentLine(Pattern pat, String line, boolean flag) {
        Matcher m = pat.matcher(line);
        if ( m.matches() ) {
            //System.out.printf("Parsing %s%n", m.group(1));
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

        // initialize the data structure
        HashMap<String, RecalDataManager> managers = new HashMap<String, RecalDataManager>();
        for ( String readGroup : readGroups ) {
            RecalDataManager manager = new RecalDataManager(readGroup, ! collapsedPos, ! collapsedDinuc);
            //RecalDataManager manager = new RecalDataManager(readGroup,  maxPos, maxQReported, dinucs.size(), ! collapsedPos, ! collapsedDinuc);
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
            //if ( serialRecalibration )
            //    mapper = new SerialRecalMapping(manager, dinucs, maxPos, maxQReported);
            //else
            mapper = new CombinatorialRecalMapping(manager, dinucs, maxPos, maxQReported);
            cache.put(readGroup, mapper);
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

        if (read.getReadNegativeStrandFlag())
            recalQuals = BaseUtils.reverse(recalQuals);
        //if ( read.getReadName().equals("30PTAAAXX090126:5:14:132:764#0") )
        //    System.out.printf("OLD: %s%n", read.format());
        read.setBaseQualities(recalQuals);
        //if ( read.getReadName().equals("30PTAAAXX090126:5:14:132:764#0") )
        //    System.out.printf("NEW: %s%n", read.format());
        return read;
    }

    public byte[] recalibrateBasesAndQuals(final String readGroup, byte[] bases, byte[] quals) {
        byte[] recalQuals = new byte[quals.length];
        RecalMapping mapper = cache.get(readGroup);

        recalQuals[0] = quals[0];               // can't change the first -- no dinuc
        for ( int cycle = 1; cycle < bases.length; cycle++ ) { // skip first and last base, qual already set because no dinuc
            byte qual = quals[cycle];
            byte newQual = mapper.getNewQual(readGroup, bases[cycle - 1], bases[cycle], cycle, qual);
            //if ( read.getReadName().equals("30PTAAAXX090126:5:14:132:764#0") )
            //    System.out.printf("Processing cycle=%d qual=%d: neg?=%b => %d at %s%n",
            //            cycle, qual, read.getReadNegativeStrandFlag(), newQual, read.getReadName());
            recalQuals[cycle] = newQual;
            //System.out.printf("Mapping %d => %d%n", qual, newQual);
        }

        return recalQuals;
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
    ArrayList<byte[][]> cache;
    RecalDataManager manager;

    public CombinatorialRecalMapping(RecalDataManager manager, Set<String> dinucs, int maxPos, int maxQReported ) {
        this.manager = manager;

        // initialize the data structure
        cache = new ArrayList<byte[][]>(RecalData.NDINUCS);
        for ( String dinuc : dinucs ) {
            cache.add(new byte[maxPos+1][maxQReported+1]);
        }

        for ( RecalData datum : manager.getAll() ) {
            //System.out.printf("Adding datum %s%n", datum);
            byte [][] table = cache.get(this.manager.getDinucIndex(datum.dinuc));
            if ( table[datum.pos][datum.qual] != 0 )
                throw new RuntimeException(String.format("Duplicate entry discovered: %s", datum));
            //table[datum.pos][datum.qual] = (byte)(1 + datum.empiricalQualByte());
            table[datum.pos][datum.qual] = datum.empiricalQualByte();
            // System.out.printf("Binding %d %d => %d%n", datum.pos, datum.qual, datum.empiricalQualByte());
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

/*class CombinatorialRecalMapping implements RecalMapping {
    HashMap<String, byte[][]> cache = new HashMap<String, byte[][]>();
    RecalDataManager manager;

    public CombinatorialRecalMapping(RecalDataManager manager, Set<String> dinucs, int maxPos, int maxQReported ) {
        this.manager = manager;

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
            // System.out.printf("Binding %d %d => %d%n", datum.pos, datum.qual, datum.empiricalQualByte());
        }
    }

    public byte getNewQual(final String readGroup, byte prevBase, byte base, int cycle, byte qual) {
        //String dinuc = String.format("%c%c", (char)prevBase, (char)base);
        //if ( qual == 2 )
        //    System.out.printf("Qual = 2%n");

        byte[] bp = {prevBase, base};
        String dinuc = manager.canonicalDinuc(new String(bp));
        int pos = manager.canonicalPos(cycle);
        byte[][] dataTable = cache.get(dinuc);

        if ( dataTable == null && prevBase != 'N' && base != 'N' )
            throw new RuntimeException(String.format("Unmapped data table at %s %s", readGroup, dinuc));

        byte result = dataTable != null && pos < dataTable.length ? dataTable[pos][qual] : qual;

        //if ( result == 2 )
        //    System.out.printf("Lookup RG=%s dinuc=%s cycle=%d pos=%d qual=%d datatable=%s / %d => %d%n",
        //            readGroup, dinuc, cycle, pos, qual, dataTable, dataTable.length, result);

        return result;
    }
}*/
/*
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
}*/
