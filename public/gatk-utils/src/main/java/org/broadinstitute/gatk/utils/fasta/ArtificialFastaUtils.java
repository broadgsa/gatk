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

package org.broadinstitute.gatk.utils.fasta;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class ArtificialFastaUtils
 *         <p/>
 *         artificial fasta utility class, for generating fake fastas.
 */
public class ArtificialFastaUtils {
    public enum BASE_PATTERN {
        RANDOM, ALL_A, ALL_T, ALL_C, ALL_G;
    }

    // what bases we support
    public enum BASES {
        A, T, C, G;
    }

    // create an artificial fasta file
    public static void createArtificialFasta(String fileName,
                                             List<String> contigNames,
                                             List<Integer> contigSizes,
                                             BASE_PATTERN pattern) {
        PrintStream s;
        try {
            s = new PrintStream(new FileOutputStream(fileName));
        } catch (FileNotFoundException e) {
            throw new ReviewedGATKException("Filename " + fileName + " passed to the ArtificialFastaUtils generated a FileNotFound exception", e);
        }
        generateFakeFasta(contigNames, contigSizes, pattern, s);
    }

    // create an artificial fasta file
    public static void createArtificialFasta(PrintStream stream,
                                             List<String> contigNames,
                                             List<Integer> contigSizes,
                                             BASE_PATTERN pattern) {

        generateFakeFasta(contigNames, contigSizes, pattern, stream);
    }

    /**
     * create a fake fasta file
     *
     * @param contigNames the pile of contig names
     * @param contigSizes the pile of contig sizes
     * @param pattern     the pattern to use for the base distrobution
     * @param s           the print stream to write to
     */
    private static void generateFakeFasta(List<String> contigNames, List<Integer> contigSizes, BASE_PATTERN pattern, PrintStream s) {
        if (contigNames.size() != contigSizes.size()) {
            throw new ReviewedGATKException("ArtificialContig name and size arrays are not equal sizes");
        }
        for (int x = 0; x < contigNames.size(); x++) {
            ArtificialContig tig = new ArtificialContig(contigNames.get(x), contigSizes.get(x), pattern);
            tig.write(s);
        }
        s.close();
    }

}


/** the fake contig class, a fasta is made up of these */
class ArtificialContig {
    public static final int COLUMN_WIDTH = 80;

    final protected String mName;
    final protected int mSize;
    final protected ArtificialFastaUtils.BASE_PATTERN mPattern;

    public ArtificialContig(String name, int size, ArtificialFastaUtils.BASE_PATTERN pat) {
        this.mName = name;
        this.mSize = size;
        this.mPattern = pat;
    }

    /**
     * write out the contig to a stream
     *
     * @param stream
     */
    public void write(PrintStream stream) {
        stream.println(">" + mName);
        int count = 0;
        while (count < mSize) {
            for (int x = 0; x < COLUMN_WIDTH; x++) {
                stream.print(generateAppropriateBase());
                count++;
                if (count >= mSize) {
                    break;
                }
            }
            stream.println();
        }
    }

    /**
     * generate the appropriate base, given the BASE_PATTERN
     *
     * @return a base, as a string
     */
    public String generateAppropriateBase() {
        switch (mPattern) {
            case RANDOM:
                return (ArtificialFastaUtils.BASES.values()[(int) Math.round(Math.random() * 4)]).toString();
            case ALL_A:
                return "A";
            case ALL_T:
                return "T";
            case ALL_C:
                return "C";
            case ALL_G:
                return "G";
            default:
                throw new ReviewedGATKException("Unknown base pattern");
        }
    }

}