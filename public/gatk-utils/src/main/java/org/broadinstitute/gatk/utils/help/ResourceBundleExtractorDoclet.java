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

import com.sun.javadoc.*;
import org.broadinstitute.gatk.utils.Utils;

import java.io.*;
import java.util.*;

/**
 * Extracts certain types of javadoc (specifically package and class descriptions) and makes them available
 * to applications at runtime.
 *
 * @author mhanna
 * @version 0.1
 */
public class ResourceBundleExtractorDoclet {
    // NOTE: Using log4j during javadoc generation requires
    // a proper Log4J initialization (see CommandLineProgram),
    // or a log4.properties file. This doclet has neither.
    //private static Logger logger = Logger.getLogger(ResourceBundleExtractorDoclet.class);

    /**
     * Taglet for the particular version number.
     */
    public static final String VERSION_TAGLET_NAME = "version";
    public static final String SUMMARY_TAGLET_NAME = "help.summary";
    public static final String DESCRIPTION_TAGLET_NAME = "help.description";

    private final RootDoc rootDoc;
    private final Set<ClassDoc> classDocs;
    private final Set<PackageDoc> packageDocs;
    private final Set<Doc> allDocs;

    protected File outFile = null;
    protected String buildTimestamp = null, absoluteVersion = null;

    /**
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws IOException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws IOException {
        ResourceBundleExtractorDoclet doclet = new ResourceBundleExtractorDoclet(rootDoc);
        doclet.checkUndocumentedClasses();
        if (doclet.isUpToDate()) {
            rootDoc.printNotice("Docs up to date. Not regenerating.");
            return true;
        }
        doclet.processDocs();
        return true;
    }

    /**
     * Validate the given options against options supported by this doclet.
     * @param option Option to validate.
     * @return Number of potential parameters; 0 if not supported.
     */
    @SuppressWarnings("unused") // Used by javadoc system
    public static int optionLength(String option) {
        if(option.equals("-build-timestamp") || option.equals("-out") || option.equals("-absolute-version") ) {
            return 2;
        }
        return 0;
    }

    /**
     * Creates a new resource extractor doclet.
     * @param  rootDoc           the documentation root.
     */
    private ResourceBundleExtractorDoclet(RootDoc rootDoc) {
        this.rootDoc = rootDoc;
        this.classDocs = new TreeSet<>();
        this.packageDocs = new TreeSet<>();
        this.allDocs = new TreeSet<>();
        for (final ClassDoc classDoc: rootDoc.classes()) {
            this.classDocs.add(classDoc);
            // Cache packages as we see them, since there's no direct way to iterate over packages.
            this.packageDocs.add(classDoc.containingPackage());
        }
        this.allDocs.addAll(classDocs);
        this.allDocs.addAll(packageDocs);
        for(final String[] options: rootDoc.options()) {
            if(options[0].equals("-out"))
                this.outFile = new File(options[1]);
            if(options[0].equals("-build-timestamp"))
                this.buildTimestamp = options[1];
            if (options[0].equals("-absolute-version"))
                this.absoluteVersion = options[1];
        }
    }

    private void checkUndocumentedClasses() {
        final Set<String> undocumentedClasses = new TreeSet<>();

        for (final ClassDoc classDoc: classDocs) {
            if(isRequiredJavadocMissing(classDoc) && shouldDocument(classDoc))
                undocumentedClasses.add(classDoc.name());
        }

        if(undocumentedClasses.size() > 0) {
            final String message = String.format("The following are currently undocumented: %s%s%s",
                    Utils.TEXT_BLINK, Utils.join(" ", undocumentedClasses), Utils.TEXT_RESET);
            for (final String line: Utils.warnUserLines(message)) {
                rootDoc.printWarning(line);
            }
        }
    }

    private boolean isUpToDate() {
        if (outFile == null)
            return false;

        final long outFileMillis = outFile.lastModified();

        if (outFileMillis == 0L) {
            return false;
        }

        for (final Doc doc: allDocs) {
            final File docFile = doc.position() == null ? null : doc.position().file();
            if (docFile != null && docFile.lastModified() > outFileMillis) {
                rootDoc.printNotice("At least one item is out of date: " + docFile.getAbsolutePath());
                return false;
            }
        }

        return true;
    }

    protected void processDocs() throws IOException {
        final PrintStream out;
        if (outFile != null) {
            out = new PrintStream(outFile);
        } else {
            out = System.out;
        }
        try {
            // Maintains a collection of resources in memory as they're accumulated.
            final Properties resourceText = new Properties();

            loadExistingResourceFile(resourceText);

            resourceText.setProperty("build.timestamp", buildTimestamp);

            for (final ClassDoc currentClass : classDocs)
                renderHelpText(resourceText, DocletUtils.getClassName(currentClass, false), currentClass);
            for (final PackageDoc currentPackage : packageDocs)
                renderHelpText(resourceText, currentPackage.name(), currentPackage);

            resourceText.store(out, "Strings displayed by the GATK help system");
        } finally {
            if (outFile != null) {
                out.close();
            }
        }
    }

    /**
     * Attempts to load the contents of the resource file named by resourceFileName into
     * our in-memory resource collection resourceText. If the resource file doesn't exist,
     * prints a notice to the user but does not throw an exception back to the calling method,
     * since we'll just create a new resource file from scratch in that case.
     * @throws IOException       if there is an I/O-related error other than FileNotFoundException
     *                           while attempting to read the resource file.
     */
    private void loadExistingResourceFile(final Properties resourceText) throws IOException {
        try {
            try (final BufferedReader resourceFile = new BufferedReader(new FileReader(outFile))) {
                resourceText.load(resourceFile);
            }
        }
        catch ( FileNotFoundException e ) {
            rootDoc.printNotice("Resource file not found -- generating a new one from scratch.");
        }
    }

    /**
     * Determine whether a given class should be documented.
     * @param classDoc the type of the given class.
     * @return True if the class should be documented.  False otherwise.
     */
    protected static boolean shouldDocument(ClassDoc classDoc) {
        if (classDoc.isAbstract()) {
            return false;
        }
        // TODO: Code duplication with GATKDoclet, including DocletUtils.getClassForDoc().
        // TODO: Refactor common methods into DocletUtils, and possibly just use DocumentGATKFeatureObjects.
        final Class<?> docClass;
        try {
            docClass = (Class<?>) DocletUtils.getClassForDoc(classDoc);
        } catch (ClassNotFoundException e) {
            return false;
        } catch (NoClassDefFoundError e) {
            return false;
        } catch (UnsatisfiedLinkError e) {
            return false; // naughty BWA bindings
        }
        if (Throwable.class.isAssignableFrom(docClass)) {
            return false; // UserExceptions
        }
        final DocumentedGATKFeature f = docClass.getAnnotation(DocumentedGATKFeature.class);
        return f != null && f.enable();
    }

    /**
     * Is the javadoc for the given class missing?
     * @param classDoc Class for which to inspect the JavaDoc.
     * @return True if the JavaDoc is missing.  False otherwise.
     */
    private static boolean isRequiredJavadocMissing(ClassDoc classDoc) {
        return classDoc.commentText().length() == 0 || classDoc.commentText().contains("Created by IntelliJ");
    }

    /**
     * Renders all the help text required for a given name.
     * @param resourceText resource text properties
     * @param elementName element name to use as the key
     * @param element Doc element to process.
     */
    private void renderHelpText(final Properties resourceText, final String elementName, final Doc element) {
        StringBuilder summaryBuilder = new StringBuilder();
        for(Tag tag: element.firstSentenceTags())
             summaryBuilder.append(tag.text());
        String summary = summaryBuilder.toString();
        String description = element.commentText();

        // this might seem unnecessary, but the GATK command line program uses this tag to determine the version when running
        if(absoluteVersion != null)
            resourceText.setProperty(String.format("%s.%s",elementName,VERSION_TAGLET_NAME),absoluteVersion);

        // Write out an alternate element summary, if exists.
        resourceText.setProperty(String.format("%s.%s",elementName,SUMMARY_TAGLET_NAME),formatText(summary));

        // Write out an alternate description, if present.
        resourceText.setProperty(String.format("%s.%s",elementName,DESCRIPTION_TAGLET_NAME),formatText(description));
    }

    /**
     * Format text for consumption by the properties file.
     * @param text Text to format.
     * @return Formatted text; string trimmed, newlines removed.
     */
    private static String formatText(String text) {
        Scanner scanner = new Scanner(text);
        StringBuilder output = new StringBuilder();

        while(scanner.hasNextLine()) {
            if(output.length() > 0)
                output.append(' ');
            output.append(scanner.nextLine().trim());
        }

        return output.toString();    
    }
}
