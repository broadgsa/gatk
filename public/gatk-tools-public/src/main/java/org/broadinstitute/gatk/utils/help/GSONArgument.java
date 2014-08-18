package org.broadinstitute.gatk.utils.help;

import java.util.List;
import java.util.Map;

/**
 * GSON-friendly version of the argument bindings
 */
public class GSONArgument {

    String summary;
    String name;
    String synonyms;
    String type;
    String required;
    String fulltext;
    String defaultValue;
    String minValue;
    String maxValue;
    String minRecValue;
    String maxRecValue;
    String rodTypes;
    String kind;
    List<Map<String, Object>> options;

    public void populate(   String summary,
                            String name,
                            String synonyms,
                            String type,
                            String required,
                            String fulltext,
                            String defaultValue,
                            String minValue,
                            String maxValue,
                            String minRecValue,
                            String maxRecValue,
                            String rodTypes,
                            String kind,
                            List<Map<String, Object>> options
    ) {
        this.summary = summary;
        this.name = name;
        this.synonyms = synonyms;
        this.type = type;
        this.required = required;
        this.fulltext = fulltext;
        this.defaultValue = defaultValue;
        this.minValue = minValue;
        this.maxValue = maxValue;
        this.minRecValue = minRecValue;
        this.maxRecValue = maxRecValue;
        this.rodTypes = rodTypes;
        this.kind = kind;
        this.options = options;
    }


}
