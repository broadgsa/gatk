package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.playground.utils.ArtificialPoolContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Set;
import java.util.List;
import java.util.ListIterator;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 26, 2009
 * Time: 11:28:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class ArtificialPoolWalker extends LocusWalker<ArtificialPoolContext, ArtificialPoolContext> {

    @Argument(fullName = "AuxOutputFile", shortName = "af", doc = "Auxiliary file for genotyp & coverage output", required = true)
        String auxFilePath = null;
    @Argument(fullName = "OutputBamFile", shortName = "of", doc = "Output to this file rather than standard output", required = false)
        SAMFileWriter outputBamFile = null;

    public void initialize() {
    }

    public ArtificialPoolContext reduceInit() { // try to initialize the file writer
        ArtificialPoolContext apContext = new ArtificialPoolContext();
        apContext.setSingleSampleGenotyper(new UnifiedGenotyper());
        apContext.setReadGroupSets(getToolkit().getMergedReadGroupsByReaders());
        apContext.setAuxWriter(initializeAuxFileWriter(apContext.getTotalNumberOfPeople()));
        apContext.setSAMFileWriter(outputBamFile);
        apContext.initializeUG();

        return apContext;
    }

    public ArtificialPoolContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return new ArtificialPoolContext(tracker,ref,context);
    }

    public ArtificialPoolContext reduce(ArtificialPoolContext mapCon, ArtificialPoolContext redCon){
        /* ArtificialPoolContext updatedContext = ArtificialPoolContext.mapReduceMerge(mapCon,redCon);
        List<SAMRecord>[] newReads = updatedContext.splitReadsByGroup(updatedContext.getNewReads());
        long[] newCvg = updateRunningCoverage(updatedContext.getRunningCoverage(), getCoverageByGroup(newReads));
        updatedContext.setRunningCoverage(newCvg);
        List<SAMRecord>[] sampledReads = ArtificialPoolContext.sampleReads(newReads,runningCoverageToDouble(newCvg));
        printToFiles(sampledReads,updatedContext);*/
        AlignmentContext context = redCon.getAlignmentContext();
        SAMFileWriter writer = redCon.getSAMFileWriter();
        for(SAMRecord read : context.getReads()) {
            writer.addAlignment(read);
        }
        PrintWriter auxWrite = redCon.getWriterToAuxiliaryFile();
        auxWrite.print("This is a test.");
        ArtificialPoolContext updatedContext = redCon;
        return updatedContext;
    }

    // Helper methods follow

    private PrintWriter initializeAuxFileWriter(int nFiles) {
        PrintWriter auxFileWriter;
        try {
            auxFileWriter = new PrintWriter(new FileOutputStream(auxFilePath));
            auxFileWriter.println(createAuxFileHeader(nFiles));
        } catch(FileNotFoundException e) {
            String errmsg = "The filepath you entered "+auxFilePath+" could not be opened. Please double-check that the input is correct.";
            throw new StingException(errmsg, e);
        } catch(IOException e) {
            String errmsg = "The file you entered "+auxFilePath+" could not be written to. Please check your permissions to write to this file.";
            throw new StingException(errmsg,e);
        }

        return auxFileWriter;
    }

    private String createAuxFileHeader(int nFiles) {
        String sp = "      ";
        String st1 = "Chrom:Pos" + sp;
        String st2 = "";
        for(int j = 0; j < nFiles; j++) {
            st2 = st2 + "Pers " + j + " Gen" + sp; // short for "genotype of person j at this location"
            st2 = st2 + "Pers " + j + " Conf" + sp; // short for "confidence in genotype call of ..."
            st2 = st2 + "Pers " + j + " NewCvg" + sp; // short for "coverage of person j at this location"
        }
        String st3 = "TotalCvg";

        return st1+st2+st3;
    }

    private int[] getCoverageByGroup(List<SAMRecord>[] readsByGroup) {
        int[] coverage = new int[readsByGroup.length];
        for(int iterator = 0; iterator < readsByGroup.length; iterator ++) {
            coverage[iterator] = readsByGroup[iterator].size();
        }

        return coverage;
    }

    private long[] updateRunningCoverage(long[] cvgUpToNow, int[] newCvgByGroup) {
        long[] newCvg = new long[cvgUpToNow.length];
        for(int iter = 0; iter < cvgUpToNow.length; iter++) {
            newCvg[iter] = cvgUpToNow[iter] + newCvgByGroup[iter];
        }

        return newCvg;
    }

    private double[] runningCoverageToDouble(long[] cvg) {
        double[] avgProp = new double[cvg.length];
        long sum = 0;
        for(long elem : cvg) {
            sum += elem;
        }
        for(int iter = 0; iter < cvg.length; iter++) {
            avgProp[iter] = cvg[iter]/sum;
        }

        return avgProp;
    }

    private void printToFiles(List<SAMRecord>[] sampledNewReads, ArtificialPoolContext context) {
        SAMFileWriter samWrite = context.getSAMFileWriter();
        String sp = "      ";
        PrintWriter auxWrite = context.getWriterToAuxiliaryFile();
        int readGroupInt = 0;
        for(List<SAMRecord> readGroup : sampledNewReads) {
            for(SAMRecord read : readGroup) {
                samWrite.addAlignment(read);
            }
            auxWrite.print(context.getAlignmentContext().getLocation().toString() + sp);
            auxWrite.print(context.genotypeAndConfidenceToString(readGroupInt,sp));
            readGroupInt++;
        }

    }


}

