package org.broadinstitute.sting.utils.codecs.vcf;

import java.util.*;

/**
 * A class for translating between vcf header versions
 */
public class VCFHeaderLineTranslator {
    private static Map<VCFHeaderVersion,VCFLineParser> mapping;

    static {
        mapping = new HashMap<VCFHeaderVersion,VCFLineParser>();
        mapping.put(VCFHeaderVersion.VCF4_0,new VCF4Parser());
        mapping.put(VCFHeaderVersion.VCF4_1,new VCF4Parser());
        mapping.put(VCFHeaderVersion.VCF3_3,new VCF3Parser());
        mapping.put(VCFHeaderVersion.VCF3_2,new VCF3Parser());
    }

    public static Map<String,String> parseLine(VCFHeaderVersion version, String valueLine, List<String> expectedTagOrder) {
        return mapping.get(version).parseLine(valueLine,expectedTagOrder);
    }
}


interface VCFLineParser {
    public Map<String,String> parseLine(String valueLine, List<String> expectedTagOrder);
}


/**
 * a class that handles the to and from disk for VCF 4 lines
 */
class VCF4Parser implements VCFLineParser {
    Set<String> bracketed = new HashSet<String>();

    /**
     * parse a VCF4 line
     * @param valueLine the line
     * @return a mapping of the tags parsed out
     */
    public Map<String, String> parseLine(String valueLine, List<String> expectedTagOrder) {
        // our return map
        Map<String, String> ret = new LinkedHashMap<String, String>();

        // a builder to store up characters as we go
        StringBuilder builder = new StringBuilder();

        // store the key when we're parsing out the values
        String key = "";

        // where are we in the stream of characters?
        int index = 0;

        // are we inside a quotation? we don't special case ',' then
        boolean inQuote = false;

        // a little switch machine to parse out the tags. Regex ended up being really complicated and ugly [yes, but this machine is getting ugly now... MAD]
        for (char c: valueLine.toCharArray()) {
            if ( c == '\"' ) {
                inQuote = ! inQuote;
            } else if ( inQuote ) {
                builder.append(c);
            } else {
                switch (c) {
                    case ('<') : if (index == 0) break; // if we see a open bracket at the beginning, ignore it
                    case ('>') : if (index == valueLine.length()-1) ret.put(key,builder.toString().trim()); break; // if we see a close bracket, and we're at the end, add an entry to our list
                    case ('=') : key = builder.toString().trim(); builder = new StringBuilder(); break; // at an equals, copy the key and reset the builder
                    case (',') : ret.put(key,builder.toString().trim()); builder = new StringBuilder(); break; // drop the current key value to the return map
                    default: builder.append(c); // otherwise simply append to the current string
                }
            }
            
            index++;
        }

        // validate the tags against the expected list
        index = 0;
        if (ret.size() > expectedTagOrder.size()) throw new IllegalArgumentException("Unexpected tag count " + ret.size() + " in string " + expectedTagOrder.size());
        for (String str : ret.keySet()) {
            if (!expectedTagOrder.get(index).equals(str)) throw new IllegalArgumentException("Unexpected tag " + str + " in string " + valueLine);
            index++;
        }
        return ret;
    }
}

class VCF3Parser implements VCFLineParser {

    public Map<String, String> parseLine(String valueLine, List<String> expectedTagOrder) {
        // our return map
        Map<String, String> ret = new LinkedHashMap<String, String>();

        // a builder to store up characters as we go
        StringBuilder builder = new StringBuilder();

        // where are we in the stream of characters?
        int index = 0;
        // where in the expected tag order are we?
        int tagIndex = 0;

        // are we inside a quotation? we don't special case ',' then
        boolean inQuote = false;

        // a little switch machine to parse out the tags. Regex ended up being really complicated and ugly
        for (char c: valueLine.toCharArray()) {
            switch (c) {
                case ('\"') : inQuote = !inQuote; break; // a quote means we ignore ',' in our strings, keep track of it
                case (',') : if (!inQuote) { ret.put(expectedTagOrder.get(tagIndex++),builder.toString()); builder = new StringBuilder(); break; } // drop the current key value to the return map
                default: builder.append(c); // otherwise simply append to the current string
            }
            index++;
        }
        ret.put(expectedTagOrder.get(tagIndex++),builder.toString());
        
        // validate the tags against the expected list
        index = 0;
        if (tagIndex != expectedTagOrder.size()) throw new IllegalArgumentException("Unexpected tag count " + tagIndex + ", we expected " + expectedTagOrder.size());
        for (String str : ret.keySet()){
            if (!expectedTagOrder.get(index).equals(str)) throw new IllegalArgumentException("Unexpected tag " + str + " in string " + valueLine);
            index++;
        }
        return ret;
    }
}