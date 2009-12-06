package org.broadinstitute.sting.utils.doc;

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
public class HelpExtractorDoclet {
    /**
     * Extracts the contents of certain types of javadoc and adds them to an XML file.
     * @param rootDoc The documentation root.
     * @return Whether the JavaDoc run succeeded.
     * @throws FileNotFoundException if output can't be written.
     */
    public static boolean start(RootDoc rootDoc) throws FileNotFoundException {
        PrintStream out = System.out;

        for(String[] options: rootDoc.options()) {
            if(options[0].equals("-out"))
                out = new PrintStream(options[1]);
        }

        // Cache packages as we see them, since there's no direct way to iterate over packages.
        Set<PackageDoc> packages = new HashSet<PackageDoc>();

        for(ClassDoc currentClass: rootDoc.classes()) {
            PackageDoc containingPackage = currentClass.containingPackage();
            packages.add(containingPackage);
            String className = containingPackage.name().length() > 0 ?
                    String.format("%s.%s",containingPackage.name(),currentClass.name()) :
                    String.format("%s",currentClass.name());

            renderHelpText(className,currentClass,out);
        }

        for(PackageDoc currentPackage: packages)
            renderHelpText(currentPackage.name(),currentPackage,out);

        return true;
    }

    /**
     * Validate the given options against options supported by this doclet.
     * @param option Option to validate.
     * @return Number of potential parameters; 0 if not supported.
     */
    public static int optionLength(String option) {
        if(option.equals("-out")) {
            return 2;
        }
        return 0;
    }

    /**
     * Renders all the help text required for a given name.
     * @param elementName element name to use as the key
     * @param element Doc element to process.
     * @param out Output stream to which to write.
     */
    private static void renderHelpText(String elementName, Doc element, PrintStream out) {
        // Provide a new display name if provided by the user.
        Tag[] tags = element.tags("@"+DisplayNameTaglet.NAME);
        if(tags.length > 1)
            throw new StingException("Cannot provide multiple display names for a single tag.");
        if(tags.length == 1)
            out.printf("%s.%s=%s%n",elementName,DisplayNameTaglet.NAME,tags[0].text());

        String commentText = formatText(element.commentText());
        if(commentText.length() > 0)
            out.printf("%s=%s%n",elementName,commentText);
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
