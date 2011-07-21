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

import com.sun.javadoc.*;
import freemarker.template.Configuration;
import freemarker.template.DefaultObjectWrapper;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import scala.reflect.Print;

import java.io.*;
import java.util.*;

/**
 *
 */
public class GATKDoclet extends ResourceBundleExtractorDoclet {
    /**
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws java.io.IOException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws IOException {
        GATKDoclet doclet = new GATKDoclet();
        //PrintStream out = doclet.loadData(rootDoc, false);
        doclet.processDocs(rootDoc, null);
        return true;
    }

    public static int optionLength(String option) {
        return ResourceBundleExtractorDoclet.optionLength(option);
    }

    @Override
    protected void processDocs(RootDoc rootDoc, PrintStream ignore) {
        try {
            /* ------------------------------------------------------------------- */
            /* You should do this ONLY ONCE in the whole application life-cycle:   */

            Configuration cfg = new Configuration();
            // Specify the data source where the template files come from.
            // Here I set a file directory for it:
            cfg.setDirectoryForTemplateLoading(new File("settings/helpTemplates/"));
            // Specify how templates will see the data-model. This is an advanced topic...
            // but just use this:
            cfg.setObjectWrapper(new DefaultObjectWrapper());

            for ( ClassDoc doc : rootDoc.classes() ) {
                if ( ResourceBundleExtractorDoclet.isWalker(doc) ) {
                    System.out.printf("Walker class %s%n", doc);
                    processWalkerDocs(cfg, doc);
                    //return;
                }
//                else
//                    System.out.printf("Excluding non-walker class %s%n", doc);
            }
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException(e);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }

    private void processWalkerDocs(Configuration cfg, ClassDoc doc) throws IOException {
        /* ------------------------------------------------------------------- */
        /* You usually do these for many times in the application life-cycle:  */

        // Create the root hash
        Map root = buildWalkerDataModel(doc);

        /* Get or create a template */
        Template temp = cfg.getTemplate("test.html");

        /* Merge data-model with template */
        Writer out = new OutputStreamWriter(new FileOutputStream(new File("testdoc/" + getClassName(doc).replace(".", "_") + ".html")));
        try {
            temp.process(root, out);
            out.flush();
        } catch ( TemplateException e ) {
            throw new ReviewedStingException("Failed to create GATK documentation", e);
        }
    }


    private Map buildWalkerDataModel(ClassDoc classdoc) {
        Map<String, Object> root = new HashMap<String, Object>();

        root.put("name", classdoc.name());

        // Extract overrides from the doc tags.
        StringBuilder summaryBuilder = new StringBuilder();
        for(Tag tag: classdoc.firstSentenceTags())
             summaryBuilder.append(tag.text());
        root.put("summary", summaryBuilder.toString());
        root.put("description", classdoc.commentText());

        for(Tag tag: classdoc.tags()) {
            root.put(tag.name(), tag.text());
        }

        ParsingEngine parsingEngine = createStandardGATKParsingEngine();
//        for (ArgumentDefinition argumentDefinition : parsingEngine.argumentDefinitions )
//            System.out.println(argumentDefinition);

        Map<String, List<Object>> args = new HashMap<String, List<Object>>();
        root.put("arguments", args);
        args.put("required", new ArrayList<Object>());
        args.put("optional", new ArrayList<Object>());
        args.put("hidden", new ArrayList<Object>());
        args.put("depreciated", new ArrayList<Object>());
        try {
            for ( ArgumentSource argumentSource : parsingEngine.extractArgumentSources(getClassForDoc(classdoc)) ) {
                String kind = "optional";
                if ( argumentSource.isRequired() ) kind = "required";
                else if ( argumentSource.isHidden() ) kind = "hidden";
                else if ( argumentSource.isDeprecated() ) kind = "depreciated";
                args.get(kind).add(argumentDataModel(argumentSource.createArgumentDefinitions().get(0)));
                System.out.printf("Processing %s%n", argumentSource);

//            for(FieldDoc fieldDoc: classdoc.fields()) {
//                //for ( AnnotationDesc desc : fieldDoc.annotations() ) {
//                    System.out.printf("AnnotationDesc %s%n", desc);
//                    if ( implementsInterface(desc.annotationType(), Argument.class, Output.class, Input.class) ) {
//                        (requiredAnnotation(desc) ? requiredArgs : optionalArgs).add(dataModelForArgument(desc));
//                        System.out.printf("Processing %s%n", desc);
//                    } else {
//                        System.out.printf("Skipping %s%n", desc);
//                    }
//                }
//            }
            }
        } catch ( ClassNotFoundException e ) {
            throw new RuntimeException(e);
        }

        System.out.printf("Root is %s%n", root);
        return root;
    }

    protected String withDefault(String val, String def) {
        return val == null ? def : val;
    }

    protected Map<String, Object> argumentDataModel(ArgumentDefinition argumentDefinition) {
        Map<String, Object> root = new HashMap<String, Object>();
        root.put("shortName", withDefault(argumentDefinition.shortName, "None provided"));
        root.put("required", argumentDefinition.required);
        root.put("fullName", withDefault(argumentDefinition.fullName, "None provided"));
        root.put("argumentType", argumentDefinition.argumentType);
        root.put("doc", withDefault(argumentDefinition.doc, "None provided"));
        return root;
    }

    protected ParsingEngine createStandardGATKParsingEngine() {
        CommandLineProgram clp = new CommandLineGATK();
        try {
            CommandLineProgram.start(clp, new String[]{}, true);
            return clp.parser;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    protected Map<String, Object> dataModelForArgument(AnnotationDesc desc) {
        Map<String, Object> root = new HashMap<String, Object>();
        root.put("shortName", "None provided");
        root.put("required", false);
        root.put("fullName", "None provided");
        root.put("doc", "None provided");

        for ( AnnotationDesc.ElementValuePair keyvalue : desc.elementValues() ) {
            root.put(keyvalue.element().name(), keyvalue.value().value());
        }
        return root;
    }

    protected boolean requiredAnnotation(AnnotationDesc desc) {
        for ( AnnotationDesc.ElementValuePair keyvalue : desc.elementValues() ) {
            if ( keyvalue.element().name().equals("required") )
                return keyvalue.value().toString().equals("true");
        }
        return false;
    }
}
