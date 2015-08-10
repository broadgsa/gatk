/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.index.IndexCreator;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.text.XReadLines;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.bcf2.BCF2Codec;
import org.broadinstitute.gatk.utils.collections.Pair;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.*;
import java.util.*;

/**
 *
 * Concatenate VCF files of non-overlapping genome intervals, all with the same set of samples
 *
 * <p>
 * The main purpose of this tool is to speed up the gather function when using scatter-gather parallelization.
 * This tool concatenates the scattered output VCF files. It assumes that:
 * <ul>
 *     <li>All the input VCFs (or BCFs) contain the same samples in the same order.</li>
 *     <li>The variants in each input file are from non-overlapping (scattered) intervals.</li>
 * </ul>
 * </p>
 * <p>When the input files are already sorted based on the intervals start positions, use -assumeSorted.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more variant sets to combine. They should be of non-overlapping genome intervals and with the same
 * samples (sorted in the same order). If the files are ordered according to the appearance of intervals in the ref
 * genome, then one can use the -assumeSorted flag.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF or BCF. The output file should have the same extension as the input(s).
 * <\p>
 *
 * <h3>Important note</h3>
 * <p>This is a command-line utility that bypasses the GATK engine. As a result, the command-line you must use to
 * invoke it is a little different from other GATK tools (see example below), and it does not accept any of the
 * classic "CommandLineGATK" arguments.</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
 *    -R reference.fasta \
 *    -V input1.vcf \
 *    -V input2.vcf \
 *    -out output.vcf \
 *    -assumeSorted
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>Currently the tool is more efficient when working with VCFs than with BCFs.</p>
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
     * will be included in the final merged output VCF/BCF. The order of arguments does not matter, but it runs more
     * efficiently if they are sorted based on the intervals and the assumeSorted argument is used.
     *
     */
    @Input(fullName="variant", shortName="V", doc="Input VCF file/s", required = true)
    private List<File> variant = null;

    @Output(fullName = "outputFile", shortName = "out", doc = "output file", required = true)
    private File outputFile = null;

    @Argument(fullName = "assumeSorted", shortName = "assumeSorted", doc = "assumeSorted should be true if the input files are already sorted (based on the position of the variants)", required = false)
    private Boolean assumeSorted = false;

    @Argument(fullName = "variant_index_type", doc = "which type of IndexCreator to use for VCF/BCF indices", required = false)
    private GATKVCFIndexType variant_index_type = GATKVCFUtils.DEFAULT_INDEX_TYPE;

    @Argument(fullName = "variant_index_parameter", doc = "the parameter (bin width or features per bin) to pass to the VCF/BCF IndexCreator", required = false)
    private Integer variant_index_parameter = GATKVCFUtils.DEFAULT_INDEX_PARAMETER;

    /*
     * print usage information
     */
    private static void printUsage() {
        System.err.println("Usage: java -cp target/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants --reference <reference> --variant <input VCF or BCF file; can specify --variant multiple times> --outputFile <outputFile> [--assumeSorted]");
        System.err.println("    The output file must be of the same type as all input files.");
        System.err.println("    If the input files are already sorted, then indicate that with --assumeSorted to improve performance.");
    }

    private enum FileType {
        VCF,
        BCF,
        BLOCK_COMPRESSED_VCF,
        INVALID
    }

    private FileType fileExtensionCheck(File inFile, FileType previousFileType) {
        final String inFileName = inFile.toString().toLowerCase();

        if (inFileName.endsWith(".vcf")) {
            if (previousFileType == FileType.VCF || previousFileType == null) {
                return FileType.VCF;
            }
        }

        if (inFileName.endsWith(".bcf")) {
            if (previousFileType == FileType.BCF || previousFileType == null) {
                return FileType.BCF;
            }
        }

        for (String extension : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS) {
            if (inFileName.endsWith(".vcf" + extension)) {
                if (previousFileType == FileType.BLOCK_COMPRESSED_VCF || previousFileType == null) {
                    return FileType.BLOCK_COMPRESSED_VCF;
                }
            }
        }

        System.err.println(String.format("File extension for input file %s is not valid for CatVariants", inFile));
        printUsage();
        return FileType.INVALID;
    }

    private FeatureReader<VariantContext> getFeatureReader(final FileType fileType, final File file) {
        FeatureReader<VariantContext> reader = null;
        switch(fileType) {
            case VCF:
            case BLOCK_COMPRESSED_VCF:
                // getFeatureReader will handle both block-compressed and plain text VCFs
                reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new VCFCodec(), false);
                break;
            case BCF:
                reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), new BCF2Codec(), false);
                break;
        }
        return reader;
    }

    /**
     * Replaces any .list files in rawFileList with the files named in said .list file
     * @param rawFileList the original file list, possibly including .list files
     * @return a new List, with .list files replaced
     */
    private List<File> parseVariantList(final List<File> rawFileList) {
        final List<File> result = new ArrayList<>(rawFileList.size());
        for (final File rawFile : rawFileList) {
            if (rawFile.getName().endsWith(".list")) {
                try {
                    for (final String line : new XReadLines(rawFile, true))
                        result.add(new File(line));
                } catch (IOException e) {
                    throw new UserException.CouldNotReadInputFile(rawFile, e);
                }
            } else {
                result.add(rawFile);
            }
        }
        return result;
    }

    @Override
    protected int execute() throws Exception {
        BasicConfigurator.configure();
        logger.setLevel(Level.INFO);

        final ReferenceSequenceFile ref;
        try {
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);
        } catch ( Exception e ) {
            throw new UserException("Couldn't load provided reference sequence file " + refFile, e);
        }

        variant = parseVariantList(variant);

        Comparator<Pair<VariantContext,File>> positionComparator = new PositionComparator(ref.getSequenceDictionary());

        Queue<Pair<VariantContext,File>> priorityQueue;
        if (assumeSorted)
            priorityQueue = new LinkedList<>();
        else
            priorityQueue = new PriorityQueue<>(10000, positionComparator);

        FileType fileType = null;
        for (File file : variant) {
            // if it returns a valid type, it will be the same for all files
            fileType = fileExtensionCheck(file, fileType);
            if (fileType == FileType.INVALID)
                return 1;

            if (assumeSorted){
                priorityQueue.add(new Pair<VariantContext,File>(null,file));
            }
            else{
                if (!file.exists()) {
                    throw new UserException(String.format("File %s doesn't exist",file.getAbsolutePath()));
                }
                FeatureReader<VariantContext> reader = getFeatureReader(fileType, file);
                Iterator<VariantContext> it = reader.iterator();
                if(!it.hasNext()){
                    System.err.println(String.format("File %s is empty. This file will be ignored",file.getAbsolutePath()));
                    continue;
                }
                VariantContext vc = it.next();
                reader.close();
                priorityQueue.add(new Pair<>(vc,file));
            }
        }

        FileOutputStream outputStream = new FileOutputStream(outputFile);
        EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
        IndexCreator idxCreator = GATKVCFUtils.makeIndexCreator(variant_index_type, variant_index_parameter, outputFile, ref.getSequenceDictionary());
        final VariantContextWriter outputWriter = VariantContextWriterFactory.create(outputFile, outputStream, ref.getSequenceDictionary(), idxCreator, options);

        boolean firstFile = true;
        int count = 0;
        while(!priorityQueue.isEmpty() ){
            count++;
            File file = priorityQueue.remove().getSecond();
            if (!file.exists()) {
                throw new UserException(String.format("File %s doesn't exist",file.getAbsolutePath()));
            }
            FeatureReader<VariantContext> reader = getFeatureReader(fileType, file);

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

    private static class PositionComparator implements Comparator<Pair<VariantContext,File>> {
    	
    	VariantContextComparator comp;
    	
    	public PositionComparator(final SAMSequenceDictionary dict){
    		comp = new VariantContextComparator(dict);
    	}

        @Override
        public int compare(final Pair<VariantContext,File> p1, final Pair<VariantContext,File> p2) {
            final VariantContext startPositionP1 = p1.getFirst();
            final VariantContext startPositionP2 = p2.getFirst();
            return comp.compare(startPositionP1, startPositionP2);
        }
    }
}
