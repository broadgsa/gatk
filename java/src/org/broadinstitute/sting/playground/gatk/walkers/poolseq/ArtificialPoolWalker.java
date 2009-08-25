/*
 * Copyright (c) 2009 The Broad Institute
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
 * TERMS OF USE: WE ARE NOT LIABLE FOR ANYTHING
 */


package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.SingleSampleGenotyper;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.genotype.GenotypeCall;

import java.util.*;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

// Display the depth of coverage at a given locus and output the power

public class ArtificialPoolWalker extends LocusWalker<List<SAMRecord>[], SAMFileWriter> {
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
    public boolean suppressPrinting = false;
    @Argument(fullName="reductionProportion", shortName="rp", doc="Proportion of total coverage reflected in pool file",required=false)
    public double red_prop = -1;
    @Argument(fullName="auxOutputFile", shortName="af", doc="Filepath for Auxiliary pool information (true genotypes, etc.)", required=true)
    public String auxFilepath = null;
    @Argument(fullName="outputBamFile",shortName="of",doc="Output to this file rather than to standard output")
    SAMFileWriter outputBamFile = null;

    private List<Set<String>> readGroupSets;
    private PrintWriter auxWrite;
    private Pair<String,Double>[] local_genotypes;
    private LinkedList<Integer>[] living_reads;
    private SingleSampleGenotyper ssg;
    private int npeople;
    //@param local_genotypes - holds the genotype (A A/ A C/ etc) for each individual. Updates at each locus.
    //@param auxWrite - the writer to the auxiliary file
    //@param readGroupSets : holds the readgroups (identifiers for individuals from each read)
    //@param ssg - the instance of the single sample genotyper
    //@param living_reads - holds on to the selected reads from each person until the location passes the read endpoint
    //                      Updates at each locus.
    //@param npeople - number of people represented in this pool (size of readGroupSets)


    public void initialize() {
        // initialize the output writer
        try {
            auxWrite = new PrintWriter(new FileOutputStream(auxFilepath));
        }
        catch(FileNotFoundException e) {
            String ermsg = "auxiliary file for Artificial Pool Walker was either a folder or could not be created";
            throw new StingException(ermsg, e);
        }

        // create the file header and write to file

        auxWrite.printf("%s%n",createFileHeader());

        // obtain the readgroup sets

        readGroupSets = getToolkit().getMergedReadGroupsByReaders();
        npeople = readGroupSets.size();

        // if a reduction proportion was unspecified, set it to 1/num_files

        if(red_prop <= 0) {
            red_prop =  1.0/npeople;
        } else {
            // do nothing muhahaha
        }

        // initialize the local genotype array

        local_genotypes = new Pair[npeople];

        // initilize the single sample genotyper

        ssg = new SingleSampleGenotyper();
        ssg.initialize();

        // initialize the living reads array

        living_reads = new LinkedList[npeople];

        for(int iter=0; iter < readGroupSets.size(); iter++) {
            living_reads[iter] = new LinkedList<Integer>();
        }
    }

    public List<SAMRecord>[] map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        updateLiving();
        // each time we move to the next locus, remove from the coverage count those reads that ended
        auxWrite.printf("%s:%s",context.getContig(),context.getPosition());

        return getNewReadsAndGenotypesByGroup(tracker, ref, context);
    }

    public SAMFileWriter reduceInit() {

        return outputBamFile;
    }

    public SAMFileWriter reduce(List<SAMRecord>[] readsByReadGroup, SAMFileWriter outFile) {
        // want to generate a sub sample of the new reads
        int total_new_coverage = 0;

        for(int pers = 0; pers < npeople; pers++) {
            total_new_coverage += readsByReadGroup[pers].size();
        }

        int sought_coverage = (int) Math.ceil(red_prop * (double)total_new_coverage);
        List<SAMRecord>[] randomReadsByGroup = drawReadsRandomlyFromReadsByGroup(readsByReadGroup,sought_coverage);
        printToFileAndAuxFile(randomReadsByGroup,sought_coverage,outFile);
        return outFile;

    }



     // ArtificalPoolWalker helper methods begin below

    private List<SAMRecord>[] drawReadsRandomlyFromReadsByGroup(List<SAMRecord>[] readsByReadGroup, int numToDraw){

        int[] simulatedCoverage = simulateRandomUniformDraws(numToDraw, readsByReadGroup);

        return removeReadsToCoverage(readsByReadGroup,simulatedCoverage);
    }

    private List<SAMRecord>[] removeReadsToCoverage(List<SAMRecord>[] readsByReadGroup, int[] simulatedCoverage) {
        for(int group = 0; group < npeople; group++) {
            int cvgAtGroup = readsByReadGroup[group].size();
            for(int drawDown = simulatedCoverage[group]; drawDown > 0; drawDown --) {
                readsByReadGroup[group].remove((int)Math.floor((double) cvgAtGroup*Math.random()));
                cvgAtGroup--;
            }
        }

        return readsByReadGroup;
    }

    private int[] simulateRandomUniformDraws(int numToDraw, List<SAMRecord>[] readsByReadGroup) {
        int[] drawsByPerson = new int[npeople];

        for(int pers = 0; pers < npeople; pers++) { //array instantiation
            drawsByPerson[pers] = 0;
        }

        for(int draw = 0; draw < numToDraw; draw++) {
            int personToDrawFrom = (int) Math.floor((double) npeople * Math.random());
            drawsByPerson[personToDrawFrom]++;
        }

        return redistributeSimulatedCoverage(drawsByPerson,readsByReadGroup,0, 0, null);
    }

    private int[] redistributeSimulatedCoverage(int[] drawsByPerson, List<SAMRecord>[] readsByReadGroup, int excess, int nOffLimits, boolean[] offLimits) {
        int [] result;
        if(offLimits == null) { // first time calling
            offLimits = new boolean[npeople];
            for(int i = 0; i < npeople; i++) {
                offLimits[i] = false;
            }
            result = redistributeSimulatedCoverage(drawsByPerson,readsByReadGroup,excess,0, offLimits);
        } else if (excess == 0) { // no overt excess, need to check individual persons
            // get excess coverage
            int newExcess = 0;

            for(int pers = 0; pers < npeople; pers++) {
                if(!offLimits[pers] && drawsByPerson[pers] > readsByReadGroup[pers].size()) {
                    newExcess += drawsByPerson[pers] - readsByReadGroup[pers].size();
                    drawsByPerson[pers] = readsByReadGroup[pers].size();
                    offLimits[pers] = true;
                    nOffLimits++;
                } else {
                    // do nothing
                }
            }

            if(newExcess == 0) { // we are done!
                result = drawsByPerson;
            } else { // there is now overt excess
                result = redistributeSimulatedCoverage(drawsByPerson,readsByReadGroup,newExcess,nOffLimits,offLimits);
            }
        } else { // there are overt excess reads to distribute
            int nRemaining = npeople - nOffLimits;
            int[] drawsToAdd = new int[nRemaining];

            for(int j = 0; j < nRemaining; j++) {
                drawsToAdd[j] = 0;
            }

            for(int draw = excess; draw > 0; draw--) {
                drawsToAdd[(int) Math.floor((double)nRemaining * Math.random())]++;
            }

            int notOffLimitsCounter = 0;

            for(int pers = 0; pers < npeople; pers++) {
                if(! offLimits[pers]) {
                    drawsByPerson[pers] += drawsToAdd[notOffLimitsCounter];
                    notOffLimitsCounter++;
                } else {
                    // do nothing
                }
            }
            result = redistributeSimulatedCoverage(drawsByPerson,readsByReadGroup,0,nOffLimits,offLimits);
        }
        return result;
    }

    private int[] createNewPersonIndex(int peopleWithReadsLeft, int[] indexLeft, int indexToRemoveFrom) {
        int[] newPersonIndex = new int[peopleWithReadsLeft-1];
        for(int i = 0; i < peopleWithReadsLeft; i++) {
            if(i < indexToRemoveFrom) {
                newPersonIndex[i] = indexLeft[i];
            } else {
                newPersonIndex[i] = indexLeft[i+1];
            }
        }
        return newPersonIndex;
    }


    private List<SAMRecord>[] getNewReadsAndGenotypesByGroup(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
        // at the same time we will:
        // 1) split reads up by read group
        // 2) get new reads
        // 3) get genotypes for each individual
        // then we will:
        // return the new reads (a subsample), stratified by read group

        List<SAMRecord> allReads = context.getReads();
        int n_reads_total = allReads.size();
        List<Integer> allOffsets = context.getOffsets();
        ListIterator readIterator = allReads.listIterator();
        ListIterator offsetIterator = allOffsets.listIterator();
        LinkedList<SAMRecord>[] reads_by_group = new LinkedList[npeople];
        LinkedList<Integer>[] offsets_by_group = new LinkedList[npeople];
        LinkedList<SAMRecord>[] new_reads_by_group = new LinkedList[npeople];

        for(int j = 0; j < npeople; j++) {
            reads_by_group[j] = new LinkedList<SAMRecord>();
            offsets_by_group[j] = new LinkedList<Integer>();
            new_reads_by_group[j] = new LinkedList<SAMRecord>();
        }

        while(readIterator.hasNext()) {
            int iter_person = 0;
            SAMRecord read = (SAMRecord) readIterator.next();
            int offset = (Integer) offsetIterator.next();
            while(!readGroupSets.get(iter_person).contains((String) read.getAttribute("RG"))) {
                iter_person++;
                // for debug purposes only:
                if(iter_person > npeople) {
                    throw new StingException("Read not found in any read group in ArtificialPoolWalker. Life sucks.");
                }
            }
            // here, iter_person is the person from whom the read came so:
            reads_by_group[iter_person].add(read);
            offsets_by_group[iter_person].add(offset);

            // check to see if the read is new
            if(offset == 0) {
                new_reads_by_group[iter_person].add(read);
                // add this to the new reads set to be randomly selected and printed
            }
        }
        // now we obtain the genotypes
        for(int iter=0; iter<npeople; iter++) {
            AlignmentContext subContext = new AlignmentContext(context.getLocation(), reads_by_group[iter], offsets_by_group[iter]);
            GenotypeCall call = ssg.map(tracker, ref, subContext);
            local_genotypes[iter].set(call.getGenotypes().get(0).getBases(),call.getGenotypes().get(0).getConfidenceScore().getScore());
        }

        return new_reads_by_group;
    }

    private String createFileHeader() {
        String sp = "      ";
        String st1 = "Chrom:Pos" + sp;
        String st2 = "";
        for(int j = 0; j < npeople; j++) {
            st2 = st2 + "Pers " + j + " Gen" + sp; // short for "genotype of person j at this location"
            st2 = st2 + "Pers " + j + " Conf" + sp; // short for "confidence in genotype call of ..."
            st2 = st2 + "Pers " + j + " Cvg" + sp; // short for "coverage of person j at this location"
        }
        String st3 = "TotalCvg";

        return st1+st2+st3;
    }

    private void updateLiving()
    {
       // this function updates the living reads, removing them if the new location goes beyond the end

        int readLengthLeft;
        ListIterator updater;

        for(int j = 0; j < npeople; j++) {
             updater = living_reads[j].listIterator();

            while(updater.hasNext()) {
                readLengthLeft = (Integer) updater.next();
                if(readLengthLeft <= 1) {
                    updater.remove();
                } else {
                    updater.set(readLengthLeft-1);
                }

            } // end while(updater.hasNext())
        } // end for(int j=0; j < npeople; j++)
     }

    private void printToFileAndAuxFile(List<SAMRecord>[] downSampledReadsByGroup, int cvg, SAMFileWriter outFile)
    {
        // NOTE: the chromosome and position were already printed in the auxFile during map.
        for(int pers = 0; pers < npeople; pers++) {
            List<SAMRecord> personReads = downSampledReadsByGroup[pers];
            for(SAMRecord read : personReads) {
                outFile.addAlignment(read);
                living_reads[pers].add(read.getReadLength());
            }
        }
        String auxformstring = "          %s          %f           %d";
        // now print to the aux file
        int totalcvg = 0;
        for(int pers = 0; pers < npeople; pers++) {
            auxWrite.printf(auxformstring, local_genotypes[pers].first, local_genotypes[pers].second,living_reads[pers].size());
            totalcvg += living_reads[pers].size();
        }
        auxWrite.printf("       %d%n", totalcvg);

    }

}
