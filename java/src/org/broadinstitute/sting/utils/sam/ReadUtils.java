/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;

import java.util.Map;
import java.util.HashMap;
import java.io.File;

/**
 * A miscellaneous collection of utilities for working with SAM files, headers, etc.
 * Static methods only, please.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadUtils {
    private ReadUtils() { }

    public static SAMFileHeader copySAMFileHeader(SAMFileHeader toCopy) {
        SAMFileHeader copy = new SAMFileHeader();

        copy.setSortOrder(toCopy.getSortOrder());
        copy.setGroupOrder(toCopy.getGroupOrder());
        copy.setProgramRecords(toCopy.getProgramRecords());
        copy.setReadGroups(toCopy.getReadGroups());
        copy.setSequenceDictionary(toCopy.getSequenceDictionary());

        for (Map.Entry<String, Object> e : toCopy.getAttributes())
            copy.setAttribute(e.getKey(), e.getValue());

        return copy;
    }

    public static SAMFileWriter createSAMFileWriterWithCompression(SAMFileHeader header, boolean presorted, String file, int compression) {
        if (file.endsWith(".bam"))
            return new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(file), compression);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, new File(file));
    }

    public static boolean isPlatformRead(SAMRecord read, String name) {
        SAMReadGroupRecord readGroup = read.getReadGroup();
        if (readGroup != null) {
            Object readPlatformAttr = readGroup.getAttribute("PL");
            if (readPlatformAttr != null)
                return readPlatformAttr.toString().toUpperCase().contains(name);
        }
        return false;
    }


    public static boolean is454Read(SAMRecord read) {
        return isPlatformRead(read, "454");
    }

    public static boolean isSOLiDRead(SAMRecord read) {
        return isPlatformRead(read, "SOLID");
    }

    public static boolean isSLXRead(SAMRecord read) {
        return isPlatformRead(read, "ILLUMINA");
    }

    private static final Map<Integer, String> readFlagNames
            = new HashMap<Integer, String>();

    static {
        readFlagNames.put(0x1, "Paired");
        readFlagNames.put(0x2, "Proper");
        readFlagNames.put(0x4, "Unmapped");
        readFlagNames.put(0x8, "MateUnmapped");
        readFlagNames.put(0x10, "Forward");
        //readFlagNames.put(0x20, "MateForward");
        readFlagNames.put(0x4, "FirstOfPair");
        readFlagNames.put(0x8, "SecondOfPair");
        readFlagNames.put(0x100, "NotPrimary");
        readFlagNames.put(0x200, "NON-PF");
        readFlagNames.put(0x400, "Duplicate");
    }

    public static String readFlagsAsString(SAMRecord rec) {
        String flags = "";
        for (int flag : readFlagNames.keySet()) {
            if ((rec.getFlags() & flag) != 0) {
                flags += readFlagNames.get(flag) + " ";
            }
        }
        return flags;
    }    
}
