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
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import sun.misc.IOUtils;

import java.io.*;
import java.util.*;

/**
 *
 */
public class GATKDoclet extends ResourceBundleExtractorDoclet {
    final protected static File SETTINGS_DIR = new File("settings/helpTemplates");
    final protected static File DESTINATION_DIR = new File("testdoc");
    final protected static Logger logger = Logger.getLogger(GATKDoclet.class);

    public static class DocWorkUnit implements Comparable<DocWorkUnit> {
        // known at the start
        final String name, filename, group;
        final DocumentedGATKFeatureHandler handler;
        final ClassDoc classDoc;
        final Class clazz;
        final DocumentedGATKFeature annotation;
        final String buildTimestamp, absoluteVersion;

        // set by the handler
        String summary;
        Map<String, Object> forTemplate;

        public DocWorkUnit(String name, String filename, String group,
                           DocumentedGATKFeature annotation, DocumentedGATKFeatureHandler handler,
                           ClassDoc classDoc, Class clazz,
                           String buildTimestamp, String absoluteVersion) {
            this.annotation = annotation;
            this.name = name;
            this.filename = filename;
            this.group = group;
            this.handler = handler;
            this.classDoc = classDoc;
            this.clazz = clazz;
            this.buildTimestamp = buildTimestamp;
            this.absoluteVersion = absoluteVersion;
        }

        public void setHandlerContent(String summary, Map<String, Object> forTemplate) {
            this.summary = summary;
            this.forTemplate = forTemplate;
        }

        public Map<String, String> toMap() {
            Map<String, String> data = new HashMap<String, String>();
            data.put("name", name);
            data.put("summary", summary);
            data.put("filename", filename);
            data.put("group", group);
            return data;
        }

        public int compareTo(DocWorkUnit other) {
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

    public Set<DocWorkUnit> workUnits() {
        TreeSet<DocWorkUnit> m = new TreeSet<DocWorkUnit>();

        for ( ClassDoc doc : rootDoc.classes() ) {
            System.out.printf("Considering %s%n", doc);
            Class clazz = getClassForClassDoc(doc);
            DocumentedGATKFeature feature = getFeatureForClassDoc(doc);
            DocumentedGATKFeatureHandler handler = createHandler(doc, feature);
            if ( handler != null && handler.shouldBeProcessed(doc) ) {
                String filename = handler.getDestinationFilename(doc);
                DocWorkUnit unit = new DocWorkUnit(doc.name(),
                        filename, feature.groupName(),
                        feature, handler, doc, clazz,
                        buildTimestamp, absoluteVersion);
                m.add(unit);
            }
        }

        return m;
    }

    @Override
    protected void processDocs(RootDoc rootDoc, PrintStream ignore) {
        // setup the global access to the root
        this.rootDoc = rootDoc;
        super.loadData(rootDoc, false);

        try {
            // basic setup
            DESTINATION_DIR.mkdirs();
            FileUtils.copyFile(new File(SETTINGS_DIR + "/style.css"), new File(DESTINATION_DIR + "/style.css"));

            /* ------------------------------------------------------------------- */
            /* You should do this ONLY ONCE in the whole application life-cycle:   */

            Configuration cfg = new Configuration();
            // Specify the data source where the template files come from.
            cfg.setDirectoryForTemplateLoading(SETTINGS_DIR);
            // Specify how templates will see the data-model. This is an advanced topic...
            cfg.setObjectWrapper(new DefaultObjectWrapper());

            Set<DocWorkUnit> myWorkUnits = workUnits();
            for ( DocWorkUnit workUnit : myWorkUnits ) {
                processDocWorkUnit(cfg, workUnit, myWorkUnits);
            }

            processIndex(cfg, new ArrayList<DocWorkUnit>(myWorkUnits));
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException(e);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }

    private DocumentedGATKFeatureHandler createHandler(ClassDoc doc, DocumentedGATKFeature feature) {
        try {
            if ( feature != null ) {
                if ( feature.enable() ) {
                    DocumentedGATKFeatureHandler handler = feature.handler().newInstance();
                    handler.setDoclet(this);
                    return handler;
                } else {
                    logger.info("Skipping disabled Documentation for " + doc);
                }
            }
        } catch ( IllegalAccessException e) {
            throw new RuntimeException(e); // the constructor is now private -- this is an error
        } catch ( InstantiationException e) {
            throw new RuntimeException(e); // the constructor is now private -- this is an error
        }

        return null;
    }

    private DocumentedGATKFeature getFeatureForClassDoc(ClassDoc doc) {
        // todo -- what do I need the ? extends Object to pass the compiler?
        Class<? extends Object> docClass = getClassForClassDoc(doc);
        if ( docClass != null && docClass.isAnnotationPresent(DocumentedGATKFeature.class) ) {
            return docClass.getAnnotation(DocumentedGATKFeature.class);
        } else {
            return null; // not annotated so it shouldn't be documented
        }
    }

    private Class<? extends Object> getClassForClassDoc(ClassDoc doc) {
        try {
            // todo -- what do I need the ? extends Object to pass the compiler?
            return (Class<? extends Object>)HelpUtils.getClassForDoc(doc);
        } catch ( ClassNotFoundException e) {
            //logger.warn("Couldn't find class for ClassDoc " + doc);
            // we got a classdoc for a class we can't find.  Maybe in a library or something
            return null;
        } catch ( NoClassDefFoundError e ) {
            return null;
        } catch ( UnsatisfiedLinkError e) {
            return null; // naughty BWA bindings
        }
    }

    public static ClassDoc getClassDocForClass(RootDoc rootDoc, Class clazz) {
        return rootDoc.classNamed(clazz.getName());
    }

    private void processIndex(Configuration cfg, List<DocWorkUnit> indexData) throws IOException {
        /* Get or create a template */
        Template temp = cfg.getTemplate("generic.index.template.html");

        /* Merge data-model with template */
        Writer out = new OutputStreamWriter(new FileOutputStream(new File(DESTINATION_DIR + "/index.html")));
        try {
            temp.process(groupIndexData(indexData), out);
            out.flush();
        } catch ( TemplateException e ) {
            throw new ReviewedStingException("Failed to create GATK documentation", e);
        }
    }

    private Map<String, Object> groupIndexData(List<DocWorkUnit> indexData) {
        //
        // root -> data -> { summary -> y, filename -> z }, etc
        //      -> groups -> group1, group2, etc.
        Map<String, Object> root = new HashMap<String, Object>();

        Collections.sort(indexData);

        Set<DocumentedGATKFeature> docFeatures = new HashSet<DocumentedGATKFeature>();
        List<Map<String, String>> data = new ArrayList<Map<String, String>>();
        for ( DocWorkUnit workUnit : indexData ) {
            data.add(workUnit.toMap());
            docFeatures.add(workUnit.annotation);
        }

        List<Map<String, String>> groups = new ArrayList<Map<String, String>>();
        for ( DocumentedGATKFeature feature : docFeatures ) {
            groups.add(toMap(feature));
        }

        root.put("data", data);
        root.put("groups", groups);
        root.put("timestamp", buildTimestamp);
        root.put("version", absoluteVersion);

        return root;
    }

    private static final Map<String, String> toMap(DocumentedGATKFeature annotation) {
        Map<String, String> root = new HashMap<String, String>();
        root.put("name", annotation.groupName());
        root.put("summary", annotation.summary());
        return root;
    }

    public final static DocWorkUnit findWorkUnitForClass(Class c, Set<DocWorkUnit> all) {
        for ( final DocWorkUnit unit : all )
            if ( unit.clazz.equals(c) )
                return unit;
        return null;
    }

    private void processDocWorkUnit(Configuration cfg, DocWorkUnit unit, Set<DocWorkUnit> all)
            throws IOException {
        System.out.printf("Processing documentation for class %s%n", unit.classDoc);

        unit.handler.processOne(rootDoc, unit, all);

        // Get or create a template
        Template temp = cfg.getTemplate(unit.handler.getTemplateName(unit.classDoc));

        // Merge data-model with template
        File outputPath = new File(DESTINATION_DIR + "/" + unit.filename);
        try {
            Writer out = new OutputStreamWriter(new FileOutputStream(outputPath));
            temp.process(unit.forTemplate, out);
            out.flush();
        } catch ( TemplateException e ) {
            throw new ReviewedStingException("Failed to create GATK documentation", e);
        }
    }
}
