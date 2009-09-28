package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Variation;

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
        Variation call1 = (Variation)tracker.lookup("callset1", null);
        Variation call2 = (Variation)tracker.lookup("callset2", null);

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
        else if ( call1.getAlternativeBaseForSNP() == call2.getAlternativeBaseForSNP())
            printVariant(intersect_writer, call1);

        // intersection (discordant)
        else
            printVariant(discord_writer, call1);
    }

    private static void printVariant(PrintWriter writer, Variation variant) {
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
