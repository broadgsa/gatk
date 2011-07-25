/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.help;

import com.sun.javadoc.*;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.Utils;

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
    /**
     * Taglet for the particular version number.
     */
    public static final String VERSION_TAGLET_NAME = "version";
    public static final String SUMMARY_TAGLET_NAME = "help.summary";
    public static final String DESCRIPTION_TAGLET_NAME = "help.description";

    /**
     * Maintains a collection of resources in memory as they're accumulated.
     */
    protected final Properties resourceText = new Properties();

    /**
     * Maintains a collection of classes that should really be documented.
     */
    protected final Set<String> undocumentedWalkers = new HashSet<String>();

    protected String buildTimestamp = null, absoluteVersion = null;

    /**
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws IOException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws IOException {
        ResourceBundleExtractorDoclet doclet = new ResourceBundleExtractorDoclet();
        PrintStream out = doclet.loadData(rootDoc, true);
        doclet.processDocs(rootDoc, out);
        return true;
    }

    protected PrintStream loadData(RootDoc rootDoc, boolean overwriteResourcesFile) {
        PrintStream out = System.out;

        for(String[] options: rootDoc.options()) {
            if(options[0].equals("-out")) {
                try {
                    loadExistingResourceFile(options[1], rootDoc);
                    if ( overwriteResourcesFile )
                        out = new PrintStream(options[1]);
                } catch ( FileNotFoundException e ) {
                    throw new RuntimeException(e);
                } catch ( IOException e ) {
                    throw new RuntimeException(e);
                }
            }
            if(options[0].equals("-build-timestamp"))
                buildTimestamp = options[1];
            if (options[0].equals("-absolute-version"))
                absoluteVersion = options[1];
        }

        resourceText.setProperty("build.timestamp",buildTimestamp);
        return out;
    }

    protected void processDocs(RootDoc rootDoc, PrintStream out) {
        // Cache packages as we see them, since there's no direct way to iterate over packages.
        Set<PackageDoc> packages = new HashSet<PackageDoc>();

        for(ClassDoc currentClass: rootDoc.classes()) {
            PackageDoc containingPackage = currentClass.containingPackage();
            packages.add(containingPackage);

            if(isRequiredJavadocMissing(currentClass) && isWalker(currentClass))
                undocumentedWalkers.add(currentClass.name());

            renderHelpText(HelpUtils.getClassName(currentClass),currentClass);
        }

        for(PackageDoc currentPackage: packages)
            renderHelpText(currentPackage.name(),currentPackage);

        try {
            resourceText.store(out,"Strings displayed by the Sting help system");
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException(e);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }

        // ASCII codes for making text blink
        final String blink = "\u001B\u005B\u0035\u006D";
        final String reset = "\u001B\u005B\u006D";

        if(undocumentedWalkers.size() > 0)
            Utils.warnUser(String.format("The following walkers are currently undocumented: %s%s%s", blink, Utils.join(" ",undocumentedWalkers), reset));
    }

    /**
     * Validate the given options against options supported by this doclet.
     * @param option Option to validate.
     * @return Number of potential parameters; 0 if not supported.
     */
    public static int optionLength(String option) {
        if(option.equals("-build-timestamp") || option.equals("-out") || option.equals("-absolute-version") ) {
            return 2;
        }
        return 0;
    }

    /**
     * Attempts to load the contents of the resource file named by resourceFileName into
     * our in-memory resource collection resourceText. If the resource file doesn't exist,
     * prints a notice to the user but does not throw an exception back to the calling method,
     * since we'll just create a new resource file from scratch in that case.
     * @param  resourceFileName  name of the resource file to attempt to load.
     * @param  rootDoc           the documentation root.
     * @throws IOException       if there is an I/O-related error other than FileNotFoundException
     *                           while attempting to read the resource file.
     */
    private void loadExistingResourceFile( String resourceFileName, RootDoc rootDoc ) throws IOException {
        try {
            BufferedReader resourceFile = new BufferedReader(new FileReader(resourceFileName));
            try {
                resourceText.load(resourceFile);
            }
            finally {
                resourceFile.close();
            }
        }
        catch ( FileNotFoundException e ) {
            rootDoc.printNotice("Resource file not found -- generating a new one from scratch.");
        }
    }

    /**
     * Determine whether a given class is a walker.
     * @param classDoc the type of the given class.
     * @return True if the class of the given name is a walker.  False otherwise.
     */
    protected static boolean isWalker(ClassDoc classDoc) {
        return HelpUtils.assignableToClass(classDoc, Walker.class, true);
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
     * @param elementName element name to use as the key
     * @param element Doc element to process.
     */
    private void renderHelpText(String elementName, Doc element) {
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
