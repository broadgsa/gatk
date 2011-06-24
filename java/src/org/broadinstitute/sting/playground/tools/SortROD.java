/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.tools;

import org.apache.log4j.BasicConfigurator;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.utils.codecs.completegenomics.CGVarCodec;
import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.utils.codecs.soapsnp.SoapSNPCodec;
import org.broad.tribble.gelitext.GeliTextCodec;
import org.broad.tribble.dbsnp.DbSNPCodec;
import org.broad.tribble.bed.BEDCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.playground.gatk.features.maf.MafCodec;

import java.io.*;
import java.util.*;

import net.sf.samtools.util.SortingCollection;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 28, 2011
 * Time: 12:15:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class SortROD {
    // setup the logging system, used by some codecs
    private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();

    /**
     * this class:
     *  1) checks to see that the feature file exists
     *  2) loads an index from disk, if one doesn't exist, it creates it and writes it to disk
     *  3) creates a FeatureSource
     *  4) iterates over the records, emitting a final tally for the number of features seen
     *
     * @param args a single parameter, the file name to load
     */
    public static void main(String[] args) {
        BasicConfigurator.configure();
        logger.setLevel(org.apache.log4j.Level.INFO);
        // check yourself before you wreck yourself - we require one arg, the input file
        if (args.length != 3 )
            printUsage();


        String refarg = args[0];
        if ( ! refarg.endsWith(".fasta")) {
            System.err.println("Reference file name must end with .fasta");
            System.exit(1);
        }

        File refFile = new File(refarg);
        if ( ! refFile.exists() ) {
            System.err.println("Reference file "+refarg+" does not exist");
            System.exit(1);
        }

        String rodType = null;
        String inputArg;
        // our feature file
        int pos = args[1].indexOf(":");
        if ( pos == -1 ) {
            inputArg = args[1];
        } else {
            rodType = args[1].substring(0,pos);
            inputArg = args[1].substring(pos+1);
        }
        File featureFile = new File(inputArg);
        if (!featureFile.exists()) {
            System.err.println("File " + featureFile.getAbsolutePath() + " doesnt' exist");
            printUsage();
        }

        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter(args[2]));
        } catch ( IOException e ) {
            System.err.println("Can not open output file "+args[2]+" for writing");
            System.exit(1);
        }

        // determine the codec
        FeatureCodec featureCodec = getFeatureCodec(featureFile,rodType);
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);

        AsciiLineReader reader = null;
        try {
            reader = new AsciiLineReader(new FileInputStream(featureFile));
        } catch (FileNotFoundException e) {
            System.err.println("File "+featureFile.getAbsolutePath()+" doesn't exist");
            System.exit(1);
        }

        // read the headers
        featureCodec.readHeader(reader);

        GenomeLocParser parser = new GenomeLocParser(ref.getSequenceDictionary());

        SortingCollection<String> sorter = SortingCollection.newInstance(String.class,
                new LineCodec(),
                new FeatureComparator(featureCodec,parser),200000);

        int nLines = 0;
        try {
            String currentLine = reader.readLine();
            while ( currentLine != null ) {
                nLines++;

                // uncomment if null returns should be ignored
                //if ( featureCodec.decodeLoc(currentLine) != null )
                sorter.add(currentLine);

                currentLine = reader.readLine();
            }
            
            for ( String s : sorter ) {
                out.write(s);
                out.write('\n');
            }
            out.close();
        } catch (IOException e) {
            System.err.println("Writing failed to the output file "+args[2]);
            System.exit(1);
        }

        logger.info("Sorting finished. Processed lines: "+nLines);
   //     runWithIndex(featureFile, codec, optimizeIndex);

    }


    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -jar SortROD.jar <reference> [<rodType>:]<inputFile> <outputFile>");
        System.err.println("    Where input can be of type: VCF (ends in .vcf or .VCF)");
        System.err.println("                                Bed (ends in .bed or .bed)");
        System.err.println("                                DbSNP (ends in .snp or .rod)");
        System.err.println("                                MAF (ends in .maf)");
        System.err.println("    If input file has non-standard extension, rodType can be specified");
        System.err.println("    (rodType always takes precedence over file extension, even if the");
        System.err.println("    latter is otherwise recognizable). rodType can be vcf, bed, dbsnp, or maf");
        System.err.println("    Reference is what the input file needs to be sorted against");

        /**
         * you could add others here; also look in the GATK code-base for an example of a dynamic way
         * to load Tribble codecs.
         */
        System.exit(1);
    }


    public static FeatureCodec getFeatureCodec(File featureFile, String rodType) {
        // quickly determine the codec type
        if ( rodType != null ) {
            if (rodType.equals("vcf") ) return new VCFCodec();
            if (rodType.equals("bed") ) return new BEDCodec();
            if (rodType.equals("cgvar") || rodType.equals("CGVar") ) return new CGVarCodec();
            if (rodType.equals("snp") || rodType.equals("dbsnp") ) return new DbSNPCodec();
            if (rodType.equals("geli.calls") || rodType.equals("geli") ) return new GeliTextCodec();
            if (rodType.equals("txt") ) return new SoapSNPCodec();
            if (rodType.equals("maf") ) return new MafCodec();
            throw new StingException("Explicitly specified rod type "+rodType+" is not recognized");
        }
        if ( featureFile.getName().endsWith(".vcf") || featureFile.getName().endsWith(".VCF") )
            return new VCFCodec();
        if (featureFile.getName().endsWith(".bed") || featureFile.getName().endsWith(".BED") )
            return new BEDCodec();
        if ( featureFile.getName().endsWith(".tsv") || featureFile.getName().endsWith(".TSV") )
            return new CGVarCodec();
        if (featureFile.getName().endsWith(".snp") || featureFile.getName().endsWith(".rod") )
            return new DbSNPCodec();
        if (featureFile.getName().endsWith(".geli.calls") || featureFile.getName().endsWith(".geli") )
            return new GeliTextCodec();
        if (featureFile.getName().endsWith(".txt") || featureFile.getName().endsWith(".TXT") )
            return new SoapSNPCodec();
        if (featureFile.getName().endsWith(".maf") || featureFile.getName().endsWith(".MAF") )
            return new MafCodec();
        throw new IllegalArgumentException("Unable to determine correct file type based on the file name, for file -> " + featureFile);
    }

    static class LineCodec implements SortingCollection.Codec<String> {
        OutputStream os;
        InputStream is;

        public void setOutputStream(OutputStream outputStream) {
            os = outputStream;
        }

        public void setInputStream(InputStream inputStream) {
            is = inputStream;
        }

        public void encode(String s) {
            try {
                os.write(s.getBytes());
                os.write('\n');
            } catch (IOException e) {
                throw new StingException("SortingCollection: Write into temporary file failed",e);
            }
        }

        public String decode() {
            List<Byte> l = new ArrayList<Byte>(1024);
            try {
                int c = is.read();
                while ( c != -1 && c != '\n' ) {
                    l.add((byte)c);
                    c = is.read();
                }
            } catch (IOException e) {
                throw new StingException("SortingCollection: Read from temporary file failed",e);
            }
            return new String(toByteArray(l));  //To change body of implemented methods use File | Settings | File Templates.
        }

        public SortingCollection.Codec<String> clone() {
            LineCodec codec = new LineCodec();
            codec.setInputStream(is);
            codec.setOutputStream(os);
            return codec;  //To change body of implemented methods use File | Settings | File Templates.
        }

        private byte [] toByteArray(List<Byte> l) {
            byte[] ret = new byte[l.size()];
            for ( int i = 0 ; i < l.size() ; i++ ) ret[i] = l.get(i);
            return ret;
        }
    }

    static class FeatureComparator implements Comparator<String> {

        FeatureCodec codec ;
        GenomeLocParser parser;

        public FeatureComparator (FeatureCodec codec, GenomeLocParser parser) {
            this.codec = codec;
            this.parser=parser;
        }

        /**
         * Compares its two arguments for order.  Returns a negative integer,
         * zero, or a positive integer as the first argument is less than, equal
         * to, or greater than the second.<p>
         * <p/>
         * In the foregoing description, the notation
         * <tt>sgn(</tt><i>expression</i><tt>)</tt> designates the mathematical
         * <i>signum</i> function, which is defined to return one of <tt>-1</tt>,
         * <tt>0</tt>, or <tt>1</tt> according to whether the value of
         * <i>expression</i> is negative, zero or positive.<p>
         * <p/>
         * The implementor must ensure that <tt>sgn(compare(x, y)) ==
         * -sgn(compare(y, x))</tt> for all <tt>x</tt> and <tt>y</tt>.  (This
         * implies that <tt>compare(x, y)</tt> must throw an exception if and only
         * if <tt>compare(y, x)</tt> throws an exception.)<p>
         * <p/>
         * The implementor must also ensure that the relation is transitive:
         * <tt>((compare(x, y)&gt;0) &amp;&amp; (compare(y, z)&gt;0))</tt> implies
         * <tt>compare(x, z)&gt;0</tt>.<p>
         * <p/>
         * Finally, the implementor must ensure that <tt>compare(x, y)==0</tt>
         * implies that <tt>sgn(compare(x, z))==sgn(compare(y, z))</tt> for all
         * <tt>z</tt>.<p>
         * <p/>
         * It is generally the case, but <i>not</i> strictly required that
         * <tt>(compare(x, y)==0) == (x.equals(y))</tt>.  Generally speaking,
         * any comparator that violates this condition should clearly indicate
         * this fact.  The recommended language is "Note: this comparator
         * imposes orderings that are inconsistent with equals."
         *
         * @param o1 the first object to be compared.
         * @param o2 the second object to be compared.
         * @return a negative integer, zero, or a positive integer as the
         *         first argument is less than, equal to, or greater than the
         *         second.
         * @throws ClassCastException if the arguments' types prevent them from
         *                            being compared by this comparator.
         */
        public int compare(String o1, String o2) {
            Feature f1 = codec.decodeLoc(o1);
            Feature f2 = codec.decodeLoc(o2);
            if ( f1 == null ) {
                if ( f2 == null ) return 0;
                else return -1; // null is less than non-null, this will hopefully push header strings up (but commented out lines will move up too!)
            }
            // f1 is not null
            if ( f2 == null ) return 1;

            GenomeLoc l1 = parser.createGenomeLoc(f1.getChr(),f1.getStart(),f1.getEnd());
            GenomeLoc l2 = parser.createGenomeLoc(f2.getChr(),f2.getStart(),f2.getEnd());
            return l1.compareTo(l2);  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}

