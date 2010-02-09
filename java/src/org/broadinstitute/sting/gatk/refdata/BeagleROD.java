package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.xReadLines;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;

public class BeagleROD extends BasicReferenceOrderedDatum {
    GenomeLoc loc;
    List<String> sampleNames = null;
    Map<String, List<String>> sampleGenotypes = new HashMap<String, List<String>>();

    public BeagleROD(String name) {
        super(name);
    }

    public String toString() { return "BeagleRod"; }

    public String delimiterRegex() {
        return " ";
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public List<String> getSampleNames() {
        return sampleNames;
    }

    public Map<String, List<String>> getGenotypes() {
        return sampleGenotypes;
    }

    public Object initialize(final File source) throws FileNotFoundException {
        String firstLine = new xReadLines(source).next();
        String[] parts = firstLine.split(" ");
        if ( parts[0].equals("I") ) {
            // I id NA12891 NA12891 NA12892 NA12892
            sampleNames = Arrays.asList(parts).subList(2, parts.length);
            return sampleNames;
        } else {
            throw new IllegalStateException("Beagle file " + source + " doesn't have required header line I");
        }
    }

    private static Pattern MARKER_PATTERN = Pattern.compile("c(.+)_p([0-9]+)");

    public static GenomeLoc parseMarkerName(String markerName) {
        Matcher m = MARKER_PATTERN.matcher(markerName);
        if ( m.matches() ) {
            String contig = m.group(1);
            long start = Long.valueOf(m.group(2));
            return GenomeLocParser.createGenomeLoc(contig, start, start);
        } else {
            throw new IllegalArgumentException("Malformatted family structure string: " + markerName + " required format is mom+dad=child");
        }
    }

    public boolean parseLine(final Object header, final String[] parts) throws IOException {
        //System.out.printf("Parsing beagle parts=%s header=%s%n", parts, header);
        List<String> sampleNames = (List<String>)header;

        if ( parts.length == 0 || ! parts[0].equals("M") )
            return false;
        else {
            loc = parseMarkerName(parts[1]);

            for ( int i = 2; i < parts.length; i++ ) {
                String sampleName = sampleNames.get(i-2);
                if ( ! sampleGenotypes.containsKey(sampleName) ) {
                    sampleGenotypes.put(sampleName, new ArrayList<String>());
                }

                sampleGenotypes.get(sampleName).add(parts[i]);
            }

            return true;
        }
    }
}
