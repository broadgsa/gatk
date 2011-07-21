/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.help;

import org.broadinstitute.sting.utils.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/21/11
 * Time: 8:51 AM
 *
 * Common documentation information about an GATK capability
 */
public class GATKDoc {
    final DocType type;
    final String name;
    List<String> synonyms;
    String summary, fulltext;
    final Map<String, String> tags;

    public final static String NA_STRING = "None provided";

    public enum DocType {
        WALKER ("Walker"),
        WALKER_ARG ("Walker argument"),
        READ_FILTER ("Read filter"),
        ENGINE_FEATURE ("Engine feature");

        private final String name;
        DocType(String name) {
            this.name = name;
        }

    };

    public GATKDoc(DocType type, String name) {
        this(type, name, new ArrayList<String>(), NA_STRING, NA_STRING, new HashMap<String, String>());
    }

    public GATKDoc(DocType type, String name, List<String> synonyms, String summary, String fulltext, Map<String, String> tags) {
        this.type = type;
        this.name = name;
        this.synonyms = synonyms;
        this.summary = summary;
        this.fulltext = fulltext;
        this.tags = tags;
    }

    public Map<String, String> toDataModel() {
        Map<String, String> model = new HashMap<String, String>();
        model.put("type", type.name);
        model.put("name", name);
        model.put("synonyms", Utils.join(",", synonyms));
        model.put("summary", summary);
        model.put("fulltext", fulltext);
        model.putAll(tags);
        return model;
    }

    public DocType getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    public List<String> getSynonyms() {
        return synonyms;
    }

    public void addSynonym(String synonyms) {
        this.synonyms.add(synonyms);
    }

    public String getSummary() {
        return summary;
    }

    public void setSummary(String summary) {
        this.summary = summary;
    }

    public String getFulltext() {
        return fulltext;
    }

    public void setFulltext(String fulltext) {
        this.fulltext = fulltext;
    }

    public void addTag(String key, String value) {
        this.tags.put(key, value);
    }
}
