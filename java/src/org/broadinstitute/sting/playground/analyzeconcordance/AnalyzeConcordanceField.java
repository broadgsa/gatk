package org.broadinstitute.sting.playground.analyzeconcordance;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * Created by IntelliJ IDEA.
 * User: kshakir
 * Date: Feb 11, 2010
 */
public enum AnalyzeConcordanceField {
    N_BASES_COVERED("all_bases", "all,summary,variant_counts", "n bases covered"),
    ALL_DBSNP_RATE("all_dbsnp", "all,summary,db_coverage", "dbsnp_rate"),
    ALL_VARIANT_COUNT("all_variants", "all,summary,variant_counts", "variants"),
    ALL_TITV_RATIO("all_titv", "all,summary,transitions_transversions", "ratio"),
    KNOWN_VARIANT_COUNT("known_variants", "known,summary,variant_counts", "variants"),
    KNOWN_TITV_RATIO("known_titv", "known,summary,transitions_transversions", "ratio"),
    NOVEL_VARIANT_COUNT("novel_variants", "novel,summary,variant_counts", "variants"),
    NOVEL_TITV_RATIO("novel_titv", "novel,summary,transitions_transversions", "ratio");

    private String columnHeader;
    private Pattern pattern;

    private AnalyzeConcordanceField(String columnHeader, String evalHeader, String analysis) {
        this.columnHeader = columnHeader;

        String lineRegex = evalHeader + " {2,}" + analysis + " {2,}([0-9.]+).*";
        this.pattern = Pattern.compile(lineRegex);
    }

    public String getColumnHeader() {
        return this.columnHeader;
    }

    public String parseLine(String line) {
        Matcher matcher = this.pattern.matcher(line);
        if (!matcher.matches())
            return null;
        return matcher.group(1);
    }
}
