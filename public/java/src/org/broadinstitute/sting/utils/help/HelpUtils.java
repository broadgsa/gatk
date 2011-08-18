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

import com.sun.javadoc.FieldDoc;
import com.sun.javadoc.PackageDoc;
import com.sun.javadoc.ProgramElementDoc;
import org.broadinstitute.sting.utils.classloader.JVMUtils;

import java.lang.reflect.Field;

public class HelpUtils {
    public final static String URL_ROOT_FOR_RELEASE_GATKDOCS = "http://www.broadinstitute.org/gsa/gatkdocs/release/";
    public final static String URL_ROOT_FOR_STABLE_GATKDOCS = "http://iwww.broadinstitute.org/gsa/gatkdocs/stable/";
    public final static String URL_ROOT_FOR_UNSTABLE_GATKDOCS = "http://iwww.broadinstitute.org/gsa/gatkdocs/unstable/";

    protected static boolean implementsInterface(ProgramElementDoc classDoc, Class... interfaceClasses) {
        for (Class interfaceClass : interfaceClasses)
            if (assignableToClass(classDoc, interfaceClass, false))
                return true;
        return false;
    }

    protected static boolean assignableToClass(ProgramElementDoc classDoc, Class lhsClass, boolean requireConcrete) {
        try {
            Class type = getClassForDoc(classDoc);
            return lhsClass.isAssignableFrom(type) && (!requireConcrete || JVMUtils.isConcrete(type));
        } catch (Throwable t) {
            // Ignore errors.
            return false;
        }
    }

    protected static Class getClassForDoc(ProgramElementDoc doc) throws ClassNotFoundException {
        return Class.forName(getClassName(doc));
    }

    protected static Field getFieldForFieldDoc(FieldDoc fieldDoc) {
        try {
            Class clazz = getClassForDoc(fieldDoc.containingClass());
            return JVMUtils.findField(clazz, fieldDoc.name());
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Reconstitute the class name from the given class JavaDoc object.
     *
     * @param doc the Javadoc model for the given class.
     * @return The (string) class name of the given class.
     */
    protected static String getClassName(ProgramElementDoc doc) {
        PackageDoc containingPackage = doc.containingPackage();
        return containingPackage.name().length() > 0 ?
                String.format("%s.%s", containingPackage.name(), doc.name()) :
                String.format("%s", doc.name());
    }

    public static String htmlFilenameForClass(Class c) {
        return c.getName().replace(".", "_") + ".html";
    }

    public static String helpLinksToGATKDocs(Class c) {
        String classPath = htmlFilenameForClass(c);
        StringBuilder b = new StringBuilder();
        b.append(URL_ROOT_FOR_RELEASE_GATKDOCS).append(classPath);
        //b.append("stable   version: ").append(URL_ROOT_FOR_STABLE_GATKDOCS).append(classPath).append("\n");
        //b.append("unstable version: ").append(URL_ROOT_FOR_UNSTABLE_GATKDOCS).append(classPath).append("\n");
        return b.toString();
    }
}