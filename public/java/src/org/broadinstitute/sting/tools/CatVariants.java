/*
* Copyright (c) 2012 The Broad Institute
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.sting.tools;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.*;
import java.util.*;


/**
 *
 * Concatenates VCF files of non-overlapped genome intervals, all with the same set of samples
 *
 * <p>
 * The main purpose of this tool is to speed up the gather function when using scatter-gather parallelization.
 * This tool concatenates the scattered output VCF files. It assumes that:
 * - All the input VCFs (or BCFs) contain the same samples in the same order.
 * - The variants in each input file are from non-overlapping (scattered) intervals.
 *
 * When the input files are already sorted based on the intervals start positions, use -assumeSorted.
 *
 * Note: Currently the tool is more efficient when working with VCFs; we will work to make it as efficient for BCFs.
 *
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more variant sets to combine. They should be of non-overlapping genome intervals and with the same samples (in the same order).
 * The input files should be 'name.vcf' or 'name.VCF' or 'name.bcf' or 'name.BCF'.
 * If the files are ordered according to the appearance of intervals in the ref genome, then one can use the -assumeSorted flag.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF. The output file should be 'name.vcf' or 'name.VCF'.
 * <\p>
 *
 * <h3>Important note</h3>
 * <p>This is a command-line utility that bypasses the GATK engine. As a result, the command-line you must use to
 * invoke it is a little different from other GATK tools (see example below), and it does not accept any of the
 * classic "CommandLineGATK" arguments.</p>
 *
 * <h3>Example</h3>
 * <pre>
 * java -cp GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
 *    -R ref.fasta \
 *    -V input1.vcf \
 *    -V input2.vcf \
 *    -out output.vcf \
 *    -assumeSorted
 * </pre>
 *
 * @author Ami Levy Moonshine
 * @since Jan 2012
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP )
public class CatVariants extends CommandLineProgram {
    // setup the logging system, used by some codecs
    private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();

    @Input(fullName = "reference", shortName = "R", doc = "genome reference file <name>.fasta", required = true)
    private File refFile = null;

    /**
     * The VCF or BCF files to merge together
     *
     * CatVariants can take any number of -V arguments on the command line.  Each -V argument
     * will be included in the final merged output VCF. The order of arguments does not matter, but it runs more
     * efficiently if they are sorted based on the intervals and the assumeSorted argument is used.
     *
     */
    @Input(fullName="variant", shortName="V", doc="Input VCF file/s named <name>.vcf or <name>.bcf", required = true)
    private List<File> variant = null;

    @Output(fullName = "outputFile", shortName = "out", doc = "output file name <name>.vcf or <name>.bcf", required = true)
    private File outputFile = null;

    @Argument(fullName = "assumeSorted", shortName = "assumeSorted", doc = "assumeSorted should be true if he input files are already sorted (based on the position of the variants", required = false)
    private Boolean assumeSorted = false;

    /*
     * print usage information
     */
    private static void printUsage() {
        System.err.println("Usage: java -cp dist/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants <reference> <input VCF or BCF files> <outputFile> [sorted (optional)]");
        System.err.println("    The input files can be of type: VCF (ends in .vcf or .VCF)");
        System.err.println("                                    BCF2 (ends in .bcf or .BCF)");
        System.err.println("    Output file must be vcf or bcf file (.vcf or .bcf)");
        System.err.println("    if the input files are already sorted, the last argument can indicate that");
    }

    @Override
    protected int execute() throws Exception {
        //if(help){
        //    printUsage();
        //    return 1;
        //}

        BasicConfigurator.configure();
        logger.setLevel(Level.INFO);

        final ReferenceSequenceFile ref;
        try {
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);
        } catch ( Exception e ) {
            throw new UserException("Couldn't load provided reference sequence file " + refFile, e);
        }

        Comparator<Pair<Integer,File>> positionComparator = new PositionComparator();


        //PriorityQueue<Pair<Integer,FeatureReader<VariantContext>>> queue =
        //        new PriorityQueue<Pair<Integer,FeatureReader<VariantContext>>>(2000, comparator);
        Queue<Pair<Integer,File>> priorityQueue;
        if(assumeSorted)
            priorityQueue = new LinkedList<Pair<Integer,File>>();
        else
            priorityQueue = new PriorityQueue<Pair<Integer,File>>(10000, positionComparator);

        Iterator<File> files = variant.iterator();
        File file;
        while (files.hasNext())   {
            file = files.next();
            if (!(file.getName().endsWith(".vcf") || file.getName().endsWith(".VCF") || file.getName().endsWith(".bcf") || file.getName().endsWith(".BCF"))){
                System.err.println("File " + file.getAbsolutePath() + " should be <name>.vcf or <name>.bcf");
                printUsage();
                return 1;
            }
            if (assumeSorted){
                priorityQueue.add(new Pair<Integer, File>(0,file));
            }
            else{
                if (!file.exists()) {
                    throw new UserException(String.format("File %s doesn't exist",file.getAbsolutePath()));
                }
                FeatureReader<VariantContext> reader;
                boolean useVCF = (file.getName().endsWith(".vcf") || file.getName().endsWith(".VCF"));
                if(useVCF)
                    reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new VCFCodec(), false);
                else
                    reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new BCF2Codec(), false);
                Iterator<VariantContext> it = reader.iterator();
                if(!it.hasNext()){
                    System.err.println(String.format("File %s is empty. This file will be ignored",file.getAbsolutePath()));
                    continue;
                }
                VariantContext vc = it.next();
                int firstPosition = vc.getStart();
                reader.close();
                //queue.add(new Pair<Integer, FeatureReader<VariantContext>>(firstPosition,reader));
                priorityQueue.add(new Pair<Integer, File>(firstPosition,file));
            }

        }

        if (!(outputFile.getName().endsWith(".vcf") || outputFile.getName().endsWith(".VCF"))){
            throw new UserException(String.format("Output file %s should be <name>.vcf", outputFile));
        }

        FileOutputStream outputStream = new FileOutputStream(outputFile);
        EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter outputWriter = VariantContextWriterFactory.create(outputFile, outputStream, ref.getSequenceDictionary(), options);

        boolean firstFile = true;
        int count =0;
        //while(!queue.isEmpty()){
        while(!priorityQueue.isEmpty() ){
            count++;
            //FeatureReader<VariantContext> reader = queue.remove().getSecond();
            file = priorityQueue.remove().getSecond();
            if (!file.exists()) {
                throw new UserException(String.format("File %s doesn't exist",file.getAbsolutePath()));
            }
            FeatureReader<VariantContext> reader;
            boolean useVCF = (file.getName().endsWith(".vcf") || file.getName().endsWith(".VCF"));
            if(useVCF)
                reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new VCFCodec(), false);
            else
                reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new BCF2Codec(), false);

            if(count%10 ==0)
                System.out.print(count);
            else
                System.out.print(".");
            if (firstFile){
                VCFHeader header = (VCFHeader)reader.getHeader();
                outputWriter.writeHeader(header);
                firstFile = false;
            }

            Iterator<VariantContext> it = reader.iterator();

            while (it.hasNext()){
                VariantContext vc = it.next();
                outputWriter.add(vc);
            }

            reader.close();

        }
        System.out.println();

        outputStream.close();
        outputWriter.close();

        return 0;
    }


    public static void main(String[] args){
        try {
            CatVariants instance = new CatVariants();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            printUsage();
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }

    private static class PositionComparator implements Comparator<Pair<Integer,File>> {

        @Override
        public int compare(Pair<Integer,File> p1, Pair<Integer,File> p2) {
            int startPositionP1 = p1.getFirst();
            int startPositionP2 = p2.getFirst();
            if (startPositionP1  == startPositionP2)
                return 0;
            return startPositionP1 < startPositionP2 ? -1 : 1 ;
        }
    }

}
