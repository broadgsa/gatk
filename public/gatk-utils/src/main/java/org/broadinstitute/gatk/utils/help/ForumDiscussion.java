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

package org.broadinstitute.gatk.utils.help;

import java.util.HashMap;
import java.util.Map;

class ForumDiscussion {

    final private static String POST_TEMPLATE = "<p>A new tool has been released!</p><p>Check out the documentation at <a href='%s'>%s</a>.</p>";

    final int Announce;
    final String Body;
    final String Category;
    final int Closed;
    final String Format;
    final String Name;
    final int Sink;
    final String Tags;
    final String Type;

    public ForumDiscussion(String name, String body, String format, String category,
                           String tagsCSV, String type, int closed, int announce, int sink) {
        this.Name = name;
        this.Body = body;
        this.Format = format;
        this.Category = category;
        this.Tags = tagsCSV;
        this.Type = type;
        this.Closed = closed;
        this.Announce = announce;
        this.Sink = sink;
    }

    public ForumDiscussion(GATKDocWorkUnit tool) {
        this(tool.name,
                String.format(POST_TEMPLATE, GATKDocUtils.URL_ROOT_FOR_RELEASE_GATKDOCS + tool.filename, tool.name),
                "Html", "tool-bulletin", tool.name + "," + tool.group + ",gatkdocs", "Discussion", 0, -1, -1);
    }

    public Map<String, String> getPostData() {
        Map<String, String> output = new HashMap<String, String>();

        output.put("Name", Name);
        output.put("Body", Body);
        output.put("Format", Format);
        output.put("Category", Category);
        if (Tags != null)
            output.put("Tags", Tags);
        if (Type != null)
            output.put("Type", Type);
        if (Closed != -1)
            output.put("Closed", Closed == 1 ? "1" : "0");
        if (Announce != -1)
            output.put("Announce", Announce == 1 ? "1" : "0");
        if (Sink != -1)
            output.put("Sink", Sink == 1 ? "1" : "0");

        return output;
    }
}