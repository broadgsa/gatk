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

import com.sun.javadoc.ClassDoc;
import com.sun.javadoc.RootDoc;
import freemarker.template.Configuration;
import freemarker.template.DefaultObjectWrapper;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.*;
import java.util.*;

/**
 *
 */
public class GATKDoclet extends ResourceBundleExtractorDoclet {
    final protected static Logger logger = Logger.getLogger(GATKDoclet.class);

    public static class DocumentationData implements Comparable<DocumentationData> {
        String name, summary, filename;
        Map<String, Object> forTemplate;
        String group;

        public DocumentationData(String name, String summary, Map<String, Object> forTemplate) {
            this.name = name;
            this.summary = summary;
            this.forTemplate = forTemplate;
        }

        public void setGroup(String group) {
            this.group = group;
        }

        public void setFilename(String filename) {
            this.filename = filename;
        }

        public Map<String, String> toMap() {
            Map<String, String> data = new HashMap<String, String>();
            data.put("name", name);
            data.put("summary", summary);
            data.put("filename", filename);
            data.put("group", group);
            return data;
        }

        public int compareTo(DocumentationData other) {
            return this.name.compareTo(other.name);
        }
    }

    RootDoc rootDoc;

    /**
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws java.io.IOException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws IOException {
        GATKDoclet doclet = new GATKDoclet();
        doclet.processDocs(rootDoc, null);
        return true;
    }

    public static int optionLength(String option) {
        return ResourceBundleExtractorDoclet.optionLength(option);
    }

    @Override
    protected void processDocs(RootDoc rootDoc, PrintStream ignore) {
        // setup the global access to the root
        this.rootDoc = rootDoc;

        try {
            /* ------------------------------------------------------------------- */
            /* You should do this ONLY ONCE in the whole application life-cycle:   */

            Configuration cfg = new Configuration();
            // Specify the data source where the template files come from.
            cfg.setDirectoryForTemplateLoading(new File("settings/helpTemplates/"));
            // Specify how templates will see the data-model. This is an advanced topic...
            cfg.setObjectWrapper(new DefaultObjectWrapper());

            List<DocumentationData> indexData = new ArrayList<DocumentationData>();
            for ( ClassDoc doc : rootDoc.classes() ) {
                System.out.printf("Considering %s%n", doc);
                DocumentedGATKFeatureHandler handler = getHandlerForClassDoc(doc);
                if ( handler != null && handler.shouldBeProcessed(doc) ) {
                    DocumentationData docData = processDocumentationWithHandler(cfg, handler, doc);
                    indexData.add(docData);
                }
            }

            processIndex(cfg, indexData);
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException(e);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }

    private DocumentedGATKFeatureHandler getHandlerForClassDoc(ClassDoc doc) {
        try {
            // todo -- what do I need the ? extends Object to pass the compiler?
            Class<? extends Object> docClass = ResourceBundleExtractorDoclet.getClassForDoc(doc);
            if ( docClass.isAnnotationPresent(DocumentedGATKFeature.class) ) {
                DocumentedGATKFeature feature = docClass.getAnnotation(DocumentedGATKFeature.class);
                DocumentedGATKFeatureHandler handler = feature.handler().newInstance();
                handler.setGroupName(feature.groupName());
                handler.setDoclet(this);
                return handler;
            } else {
                return null; // not annotated so it shouldn't be documented
            }
        } catch ( ClassNotFoundException e) {
            //logger.warn("Couldn't find class for ClassDoc " + doc);
            // we got a classdoc for a class we can't find.  Maybe in a library or something
            return null;
        } catch ( IllegalAccessException e) {
            throw new RuntimeException(e); // the constructor is now private -- this is an error
        } catch ( InstantiationException e) {
            throw new RuntimeException(e); // the constructor is now private -- this is an error
        } catch ( UnsatisfiedLinkError e) {
            return null; // naughty BWA bindings
        }
    }

    private void processIndex(Configuration cfg, List<DocumentationData> indexData) throws IOException {
        /* Get or create a template */
        Template temp = cfg.getTemplate("generic.index.template.html");

        /* Merge data-model with template */
        Writer out = new OutputStreamWriter(new FileOutputStream(new File("testdoc/index.html")));
        try {
            temp.process(groupIndexData(indexData), out);
            out.flush();
        } catch ( TemplateException e ) {
            throw new ReviewedStingException("Failed to create GATK documentation", e);
        }
    }

    private Map<String, Object> groupIndexData(List<DocumentationData> indexData) {
        //
        // root -> data -> { summary -> y, filename -> z }, etc
        //      -> groups -> group1, group2, etc.
        Map<String, Object> root = new HashMap<String, Object>();

        Collections.sort(indexData);

        List<Map<String, String>> data = new ArrayList<Map<String, String>>();
        Set<String> groups = new HashSet<String>();

        for ( DocumentationData indexDatum : indexData ) {
            data.add(indexDatum.toMap());
            groups.add(indexDatum.group);
        }

        root.put("data", data);
        root.put("groups", new ArrayList<String>(groups));

        return root;
    }

    private DocumentationData processDocumentationWithHandler(Configuration cfg,
                                                              DocumentedGATKFeatureHandler handler,
                                                              ClassDoc doc)
            throws IOException {
        System.out.printf("Processing documentation for class %s%n", doc);

        DocumentationData docData = handler.processOne(doc);
        docData.setGroup(handler.getGroupName());

        // Get or create a template
        Template temp = cfg.getTemplate(handler.getTemplateName(doc));

        // Merge data-model with template
        String filename = handler.getDestinationFilename(doc);
        docData.setFilename(filename);
        File outputPath = new File("testdoc/" + filename);
        try {
            Writer out = new OutputStreamWriter(new FileOutputStream(outputPath));
            temp.process(docData.forTemplate, out);
            out.flush();
        } catch ( TemplateException e ) {
            throw new ReviewedStingException("Failed to create GATK documentation", e);
        }

        return docData;
    }
}
