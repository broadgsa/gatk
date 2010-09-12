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

import java.util.*;
import java.io.PrintStream;
import java.io.IOException;

import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.gatk.walkers.Walker;

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
    private static final String VERSION_TAGLET_NAME = "version";

    /**
     * Maintains a collection of resources in memory as they're accumulated.
     */
    private static final Properties resourceText = new Properties();

    /**
     * Maintains a collection of classes that should really be documented.
     */
    private static final Set<String> undocumentedWalkers = new HashSet<String>();

    /**
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws IOException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws IOException {
        PrintStream out = System.out;
        String buildTimestamp = null, versionPrefix = null, versionSuffix = null;

        for(String[] options: rootDoc.options()) {
            if(options[0].equals("-out"))
                out = new PrintStream(options[1]);
            if(options[0].equals("-build-timestamp"))
                buildTimestamp = options[1];
            if(options[0].equals("-version-prefix"))
                versionPrefix = options[1];
            if(options[0].equals("-version-suffix"))
                versionSuffix = options[1];
        }

        resourceText.setProperty("build.timestamp",buildTimestamp);

        // Cache packages as we see them, since there's no direct way to iterate over packages.
        Set<PackageDoc> packages = new HashSet<PackageDoc>();

        for(ClassDoc currentClass: rootDoc.classes()) {
            PackageDoc containingPackage = currentClass.containingPackage();
            packages.add(containingPackage);

            if(isRequiredJavadocMissing(currentClass) && isWalker(currentClass))
                undocumentedWalkers.add(currentClass.name());

            renderHelpText(getClassName(currentClass),currentClass,versionPrefix,versionSuffix);
        }

        for(PackageDoc currentPackage: packages)
            renderHelpText(currentPackage.name(),currentPackage,versionPrefix,versionSuffix);

        resourceText.store(out,"Strings displayed by the Sting help system");

        // ASCII codes for making text blink
        final String blink = "\u001B\u005B\u0035\u006D";
        final String reset = "\u001B\u005B\u006D";

        if(undocumentedWalkers.size() > 0)
            Utils.warnUser(String.format("The following walkers are currently undocumented: %s%s%s", blink, Utils.join(" ",undocumentedWalkers), reset));

        return true;
    }

    /**
     * Validate the given options against options supported by this doclet.
     * @param option Option to validate.
     * @return Number of potential parameters; 0 if not supported.
     */
    public static int optionLength(String option) {
        if(option.equals("-build-timestamp") || option.equals("-version-prefix") || option.equals("-version-suffix") || option.equals("-out")) {
            return 2;
        }
        return 0;
    }

    /**
     * Determine whether a given class is a walker.
     * @param classDoc the type of the given class.
     * @return True if the class of the given name is a walker.  False otherwise.
     */
    private static boolean isWalker(ClassDoc classDoc) {
        try {
            Class type = Class.forName(getClassName(classDoc));
            return Walker.class.isAssignableFrom(type) && JVMUtils.isConcrete(type);
        }
        catch(Throwable t) {
            // Ignore errors.
            return false;
        }
    }

    /**
     * Reconstitute the class name from the given class JavaDoc object.
     * @param classDoc the Javadoc model for the given class.
     * @return The (string) class name of the given class.
     */
    private static String getClassName(ClassDoc classDoc) {
        PackageDoc containingPackage = classDoc.containingPackage();
        return containingPackage.name().length() > 0 ?
                String.format("%s.%s",containingPackage.name(),classDoc.name()) :
                String.format("%s",classDoc.name());
    }

    /**
     * Is the javadoc for the given class missing?
     * @param classDoc Class for which to inspect the JavaDoc.
     * @return True if the JavaDoc is missing.  False otherwise.
     */
    private static boolean isRequiredJavadocMissing(ClassDoc classDoc) {
        if(classDoc.containingPackage().name().contains("oneoffprojects"))
            return false;
        return classDoc.commentText().length() == 0 || classDoc.commentText().contains("Created by IntelliJ");
    }

    /**
     * Renders all the help text required for a given name.
     * @param elementName element name to use as the key
     * @param element Doc element to process.
     * @param versionPrefix Text to add to the start of the version string.
     * @param versionSuffix Text to add to the end of the version string.
     */
    private static void renderHelpText(String elementName, Doc element, String versionPrefix, String versionSuffix) {
        // Extract overrides from the doc tags.
        String name = null;
        String version = null;
        StringBuilder summaryBuilder = new StringBuilder();
        for(Tag tag: element.firstSentenceTags())
             summaryBuilder.append(tag.text());
        String summary = summaryBuilder.toString();
        String description = element.commentText();

        for(Tag tag: element.tags()) {
            if(tag.name().equals("@"+DisplayNameTaglet.NAME)) {
                if(name != null)
                    throw new GATKException("Only one display name tag can be used per package / walker.");
                name = tag.text();
            }
            else if(tag.name().equals("@"+VERSION_TAGLET_NAME))
                version = String.format("%s%s%s", (versionPrefix != null) ? versionPrefix : "",
                                                  tag.text(),
                                                  (versionSuffix != null) ? versionSuffix : "");
            else if(tag.name().equals("@"+SummaryTaglet.NAME))
                summary = tag.text();
            else if(tag.name().equals("@"+DescriptionTaglet.NAME))
                description = tag.text();
        }

        // Write out an alternate element name, if exists.
        if(name != null)
            resourceText.setProperty(String.format("%s.%s",elementName,DisplayNameTaglet.NAME),name);

        if(version != null)
            resourceText.setProperty(String.format("%s.%s",elementName,VERSION_TAGLET_NAME),version);

        // Write out an alternate element summary, if exists.
        resourceText.setProperty(String.format("%s.%s",elementName,SummaryTaglet.NAME),formatText(summary));

        // Write out an alternate description, if present.
        resourceText.setProperty(String.format("%s.%s",elementName,DescriptionTaglet.NAME),formatText(description));
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
