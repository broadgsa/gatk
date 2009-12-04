package org.broadinstitute.sting.utils.doc;

import com.sun.javadoc.RootDoc;
import com.sun.javadoc.PackageDoc;
import com.sun.javadoc.ClassDoc;

import java.util.HashSet;
import java.util.Set;
import java.util.Scanner;
import java.io.PrintStream;
import java.io.FileNotFoundException;

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
            String commentText = formatText(currentClass.commentText());

            if(commentText.length() > 0)
                out.printf("%s=%s%n",className,commentText);
        }

        for(PackageDoc currentPackage: packages) {
            String commentText = formatText(currentPackage.commentText());
            if(commentText.length() > 0)
                out.printf("%s=%s%n",currentPackage.name(),commentText);
        }

        return true;
    }

    /**
     * Validate the given options against 
     * @param option
     * @return
     */
    public static int optionLength(String option) {
        if(option.equals("-out")) {
            return 2;
        }
        return 0;
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
