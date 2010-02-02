package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.playground.gatk.walkers.poolseq.PowerBelowFrequencyWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.ConcordanceTruthTable;

import java.io.*;
import java.util.LinkedList;
import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 12, 2009
 * Time: 12:31:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class HapmapPoolAllelicInfoWalker extends LocusWalker<String, PrintWriter> {
    @Argument(fullName="outputFile", shortName="of", doc="File to write to", required=true)
    public String outputFileString = null;
    @Argument(fullName="numIndividualsInPool", shortName="ps",doc="Pool size",required = true)
    public int poolSize = -1;
    @Argument(fullName="sampleNames", shortName="samples", doc="Sample name bindings", required=true)
    public String sampleNameFile = null;
    @Argument(fullName="minCallQuality", shortName="q", doc="Ignore calls with below this quality, defaults to -1")
    public double minCallQ = -1;

    private PrintWriter output;
    private static double EPSILON = Math.pow(10,-4);
    private String[] sampleNames = null;
    private PowerBelowFrequencyWalker powerWalker = null;
    private ConcordanceTruthTable ctt = null;

    public void initialize() {
        sampleNames = generateNameTableFromFile(sampleNameFile);
        powerWalker = new PowerBelowFrequencyWalker();
        powerWalker.initialize();
        powerWalker.setPoolSize(poolSize);
        ctt = new ConcordanceTruthTable(poolSize);
    }

    public PrintWriter reduceInit() {
        try {
            output = new PrintWriter(outputFileString);
        } catch (FileNotFoundException e) {
            throw new StingException("File "+outputFileString+" could not be opened.", e);
        }
        output.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n","Chrom","Pos","Ref","Var","Num_Alleles","Num_Chips","Depth","Power","Support","Called");
        //System.out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n","Chrom","Pos","Ref","Var","Num_Alleles","Depth","Power","Support","Called");
        return output;
    }

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc loc = context.getLocation();
        String chrom = loc.getContig();
        long pos = loc.getStart();
        char refBase = Character.toUpperCase(ref.getBase());
        List<Pair<Genotype, Genotype>> chips = getChips(sampleNames, tracker);
        Pair<Integer,Pair<Integer,Integer>> alleleFreqInfo = ctt.getPooledAlleleFrequency(chips,refBase);
        char alternate;
        if ( alleleFreqInfo.first == ConcordanceTruthTable.VARIANT ) {
            //System.out.println(refBase + " " + alleleFreqInfo.getFirst().getBases());
            alternate = getAlternateBase(chips,refBase);

        } else {
            return null; // early return
        }
        int numVariantAllele = alleleFreqInfo.getSecond().getFirst();
        int numChipsObserved = alleleFreqInfo.getSecond().getSecond();
        int depth = context.size();
        double power = powerWalker.calculatePowerAtFrequency(context,numVariantAllele);
        int called;
        Variation call = (Variation) tracker.lookup("calls",null);
        if ( call == null ) {
            called = 0;
        } else if ( call.isReference() || call.getNegLog10PError() < minCallQ-EPSILON ) {
            called = 0;
        } else {
            called = 1;
        }

        ReadBackedPileup p = context.getPileup();
        int support = p.getBaseCounts()[BaseUtils.simpleBaseToBaseIndex(alternate)];

        // sanity check
        if ( refBase == alternate ) {
            if ( alleleFreqInfo.first == ConcordanceTruthTable.VARIANT ) {
                ;//logger.warn("Called as a variant! Ref: "+ refBase +"Chip data: " + alleleFreqInfo.getFirst().getBases());
            }
        }

        return String.format("%s\t%d\t%c\t%c\t%d\t%d\t%d\t%f\t%d\t%d",chrom,pos,refBase,alternate,numVariantAllele,numChipsObserved,depth,power,support,called);

    }

    public char getAlternateBase(List<Pair<Genotype, Genotype>> chips, char ref) {
        for ( Pair<Genotype, Genotype> chip : chips ) {
            Genotype g = chip.first;
            char[] bases = g.getBases().toCharArray();
            if ( Character.toUpperCase(bases[0]) != ref )
                return bases[0];
            if ( Character.toUpperCase(bases[1]) != ref )
                return bases[1];
        }
        return ref;
    }

    public PrintWriter reduce(String s, PrintWriter p) {
        if ( s == null ) {
            // do nothing
            return p;
        } else {
            //System.out.printf("%s%n",s);
            output.printf("%s%n",s);
            return p;
        }
    }

    public void onTraversalDone(PrintWriter p) {
        output.close();
    }

    private List<Pair<Genotype,Genotype>> getChips(String[] rodNames, RefMetaDataTracker tracker) {
        List<Pair<Genotype, Genotype>> chips = new ArrayList <Pair<Genotype,Genotype>>(rodNames.length);
        for ( String name : rodNames ) {
            RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData(name, null);
            Variation chip = (rods == null ? null : (Variation)rods.getRecords().get(0));
            if ( chip != null ) {
                // chips must be Genotypes
                if ( !(chip instanceof VariantBackedByGenotype) )
                    throw new StingException("Failure: trying to analyze genotypes using non-genotype truth data");
                chips.add(new Pair<Genotype,Genotype>(((VariantBackedByGenotype)chip).getCalledGenotype(),null));
            }
        }

        return chips;
    }
    // private methods for reading in names from a file

    private String[] generateNameTableFromFile(String file) {
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(file));
        } catch( FileNotFoundException e) {
            String errMsg = "Hapmap pool file at "+file+" was not found. Please check filepath.";
            throw new StingException(errMsg, e);
        }

        LinkedList<String> nameList = new LinkedList<String>();

        while(continueReading(reader)) {
            String line = readLine(reader);
            nameList.add(line);
        }

        return nameList.toArray(new String[nameList.size()]);
    }

    private boolean continueReading(BufferedReader reader) {
        boolean continueReading;
        try {
            continueReading = reader.ready();
        } catch(IOException e) {
            continueReading = false;
        }
        return continueReading;
    }

    private String readLine(BufferedReader reader) {
        String line;
        try {
            line = reader.readLine();
        } catch( IOException e) {
            String errMsg = "BufferedReader pointing to "+reader.toString()+" was declared ready but no line could be read from it.";
            throw new StingException(errMsg,e);
        }
        return line;
    }

}
