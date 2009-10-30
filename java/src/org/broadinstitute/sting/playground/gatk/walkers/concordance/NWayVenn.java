package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map.Entry;

/**
 * Split up N call sets into their various "Venn diagram" sets.
 * Note that to minimize complexity (unlike SimpleVenn), this module does NOT check for discordance
 * (i.e. the intersections contain calls that are present in multiple sets, regardless of whether
 * they agree on the actual variant).
 */
public class NWayVenn implements ConcordanceType {

    // TODO -- change this to use Ryan's generic map object when it's ready
    private HashMap<Integer, PrintWriter> writers = new HashMap<Integer, PrintWriter>();

    private String prefix;
    private int N;

    private PrintWriter union_writer = null;

    public NWayVenn() {}

    public void initialize(String prefix, HashMap<String,String> args) {

        if ( args.get("N") == null )
            throw new StingException("NWayVenn concordance module requires the 'N' argument");
        N = Integer.valueOf(args.get("N"));
        if ( N < 1 )
            throw new StingException("NWayVenn concordance module requires that N be greater than 0");

        this.prefix = prefix;

        // create the list of file names
        HashMap<String, String> files = new HashMap<String, String>();
        files.put("", "");
        for (int i = 1; i <= N; i++) {
            // create 2 new list entries for each of the current file names:
            // one with and one without this set appended to the end
            HashMap<String, String> appendedFiles = new HashMap<String, String>();            
            for ( Entry<String, String> file : files.entrySet() ) {
                appendedFiles.put(file.getKey() + "0", file.getValue());
                appendedFiles.put(file.getKey() + "1", file.getValue() + ".set" + i);
            }
            files = appendedFiles;
        }
        // remove the entry for no hits
        files.remove(Utils.dupString('0', N));

        try {
            for ( Entry<String, String> file : files.entrySet() )
                writers.put(stringToHash(file.getKey()), new PrintWriter(prefix + "." + N + "wayVenn" + file.getValue() + ".calls"));
            union_writer = new PrintWriter(prefix + "." + N + "wayVenn.union.calls");
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }
    }

    public void computeConcordance(RefMetaDataTracker tracker, ReferenceContext ref) {
        Variation[] calls = new Variation[N];
        int goodRods = 0;
        for (int i = 0; i < N; i++) {
            RODRecordList<ReferenceOrderedDatum> callList = tracker.getTrackData("callset" + (i+1), null);
            if ( callList != null ) {
                calls[i] = (Variation)callList.getRecords().get(0);
                goodRods++;
            }
        }

        if ( goodRods == 0 )
            return;

        // union
        printVariant(union_writer, calls);

        // add to the correct set
        StringBuilder hashString = new StringBuilder();
        for (int i = 0; i < N; i++)
            hashString.append(calls[i] != null ? "1" : "0");
        printVariant(writers.get(stringToHash(hashString.toString())), calls);
    }

    private static void printVariant(PrintWriter writer, Variation[] variants) {
        for ( Variation variant : variants ) {
            if ( variant != null ) {
                writer.println(variant.toString());
                return;
            }
        }
    }

    public void cleanup() {
        union_writer.close();
        for ( PrintWriter writer : writers.values() )
            writer.close();
    }

    private int stringToHash(String s) {
        int hash = 0;
        for (int i = 0; i < s.length(); i++) {
            if ( s.charAt(i) == '1' )
                hash += Math.pow(2, i);
        }
        return hash;
    }
}