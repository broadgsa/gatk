package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.StingException;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;

/**
 * Split up two call sets into their various "Venn diagram" sets
 */
public class SimpleVenn implements ConcordanceType {

    private PrintWriter union_writer = null;
    private PrintWriter intersect_writer = null;
    private PrintWriter discord_writer = null;
    private PrintWriter set1_writer = null;
    private PrintWriter set2_writer = null;
    private String prefix;

    public SimpleVenn(String prefix) {
        this.prefix = prefix;
    }

    public void initialize(HashMap<String,String> args) {
        try {
            union_writer = new PrintWriter(prefix + ".venn.union.calls");
            intersect_writer = new PrintWriter(prefix + ".venn.intersection.calls");
            discord_writer = new PrintWriter(prefix + ".venn.discordant.calls");
            set1_writer = new PrintWriter(prefix + ".venn.set1Only.calls");
            set2_writer = new PrintWriter(prefix + ".venn.set2Only.calls");
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }
    }

    public void computeConcordance(RefMetaDataTracker tracker, ReferenceContext ref) {
        AllelicVariant call1 = (AllelicVariant)tracker.lookup("callset1", null);
        AllelicVariant call2 = (AllelicVariant)tracker.lookup("callset2", null);

        if ( call1 == null && call2 == null )
            return;

        // union
        printVariant(union_writer, call1 != null ? call1 : call2);

        // set 1 only
        if ( call2 == null )
            printVariant(set1_writer, call1);

        // set 2 only
        else if ( call1 == null )
            printVariant(set2_writer, call2);

        // intersection (concordant)
        else if ( call1.getAltBasesFWD().equalsIgnoreCase(call2.getAltBasesFWD()) )
            printVariant(intersect_writer, call1);

        // intersection (discordant)
        else
            printVariant(discord_writer, call1);
    }

    private static void printVariant(PrintWriter writer, AllelicVariant variant) {
        writer.println(variant.toString());
    }

    public void cleanup() {
        union_writer.close();
        intersect_writer.close();
        discord_writer.close();
        set1_writer.close();
        set2_writer.close();
    }
}
