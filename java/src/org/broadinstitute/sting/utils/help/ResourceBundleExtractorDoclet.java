package org.broadinstitute.sting.utils.help;

import com.sun.javadoc.*;

import java.util.HashSet;
import java.util.Set;
import java.util.Scanner;
import java.io.PrintStream;
import java.io.FileNotFoundException;

import org.broadinstitute.sting.utils.StingException;

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
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws FileNotFoundException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws FileNotFoundException {
        PrintStream out = System.out;
        String versionPrefix = null, versionSuffix = null;

        for(String[] options: rootDoc.options()) {
            if(options[0].equals("-out"))
                out = new PrintStream(options[1]);
            if(options[0].equals("-version-prefix"))
                versionPrefix = options[1];
            if(options[0].equals("-version-suffix"))
                versionSuffix = options[1];
        }

        // Cache packages as we see them, since there's no direct way to iterate over packages.
        Set<PackageDoc> packages = new HashSet<PackageDoc>();

        for(ClassDoc currentClass: rootDoc.classes()) {
            PackageDoc containingPackage = currentClass.containingPackage();
            packages.add(containingPackage);
            String className = containingPackage.name().length() > 0 ?
                    String.format("%s.%s",containingPackage.name(),currentClass.name()) :
                    String.format("%s",currentClass.name());

            renderHelpText(className,currentClass,out,versionPrefix,versionSuffix);
        }

        for(PackageDoc currentPackage: packages)
            renderHelpText(currentPackage.name(),currentPackage,out,versionPrefix,versionSuffix);

        return true;
    }

    /**
     * Validate the given options against options supported by this doclet.
     * @param option Option to validate.
     * @return Number of potential parameters; 0 if not supported.
     */
    public static int optionLength(String option) {
        if(option.equals("-out") || option.equals("-version-prefix") || option.equals("-version-suffix")) {
            return 2;
        }
        return 0;
    }

    /**
     * Renders all the help text required for a given name.
     * @param elementName element name to use as the key
     * @param element Doc element to process.
     * @param out Output stream to which to write.
     * @param versionSuffix Text to add to the end of the version string.
     */
    private static void renderHelpText(String elementName, Doc element, PrintStream out, String versionPrefix, String versionSuffix) {
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
                    throw new StingException("Only one display name tag can be used per package / walker.");
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
            out.printf("%s.%s=%s%n",elementName,DisplayNameTaglet.NAME,name);

        if(version != null)
            out.printf("%s.%s=%s%n",elementName,VERSION_TAGLET_NAME,version);

        // Write out an alternate element summary, if exists.
        out.printf("%s.%s=%s%n",elementName,SummaryTaglet.NAME,formatText(summary));

        // Write out an alternate description, if present.
        out.printf("%s.%s=%s%n",elementName,DescriptionTaglet.NAME,formatText(description));
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
