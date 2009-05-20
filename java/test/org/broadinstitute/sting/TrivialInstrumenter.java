package org.broadinstitute.sting;

import javassist.*;
import org.junit.Ignore;

import java.io.IOException;
import java.lang.instrument.ClassFileTransformer;
import java.lang.instrument.IllegalClassFormatException;
import java.lang.instrument.Instrumentation;

/**
 *
 * User: aaron
 * Date: May 5, 2009
 * Time: 4:27:04 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date May 5, 2009
 * <p/>
 * Class TrivalInstrumenter
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */

/** A trivial example program that basically just says hello! */
@Ignore
public class TrivialInstrumenter implements ClassFileTransformer {
    public static void premain(String options, Instrumentation ins) {
        if (options != null) {
            System.out.printf("  I've been called with options: \"%s\"\n", options);
        } else
            ins.addTransformer(new TrivialInstrumenter());
    }

    public byte[] transform(ClassLoader loader,
                            String className,
                            Class cBR, java.security.ProtectionDomain pD,
                            byte[] classfileBuffer)
            throws IllegalClassFormatException {
        int size = classfileBuffer.length;

        if (className.contains("broadinstitute") &&
                className.endsWith("Test") &&
                !(className.endsWith("BaseTest"))) {
            ClassPool pool = ClassPool.getDefault();
            CtClass cl = null;

            try {
                cl = pool.makeClass(new java.io.ByteArrayInputStream(classfileBuffer));
                if (cl.isInterface() == false) {
                    for (CtBehavior meth : cl.getDeclaredMethods()) {

                        if (meth.isEmpty() == false) {
                            Object anns[] = meth.getAvailableAnnotations();
                            boolean weAreAJunitTest = false;
                            for (Object obj : anns) {
                                if (obj instanceof org.junit.Test) {
                                    weAreAJunitTest = true;
                                }
                            }
                            if (weAreAJunitTest) {
                                addAnnouncement(meth, cl);
                            }
                        }
                    }
                    classfileBuffer = cl.toBytecode();
                    return classfileBuffer;
                }

                // baseTearDown
            } catch (NotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (CannotCompileException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } finally {
                if (cl != null) {
                    cl.detach();
                }
            }

        }
        return null;
    }

    private static void addTiming(CtClass clas, String mname)
            throws NotFoundException, CannotCompileException {

        //  get the method information (throws exception if method with
        //  given name is not declared directly by this class, returns
        //  arbitrary choice if more than one with the given name)
        CtMethod mold = clas.getDeclaredMethod(mname);

        //  rename old method to synthetic name, then duplicate the
        //  method with original name for use as interceptor
        String nname = mname + "$impl";
        mold.setName(nname);
        CtMethod mnew = CtNewMethod.copy(mold, mname, clas, null);

        //  start the body text generation by saving the start time
        //  to a local variable, then call the timed method; the
        //  actual code generated needs to depend on whether the
        //  timed method returns a value
        String type = mold.getReturnType().getName();
        StringBuffer body = new StringBuffer();
        body.append("{\nlong start = System.currentTimeMillis();\n");
        if (!"void".equals(type)) {
            body.append(type + " result = ");
        }
        body.append(nname + "($$);\n");

        //  finish body text generation with call to print the timing
        //  information, and return saved value (if not void)
        body.append("System.out.println(\"Call to method " + mname +
                " took \" +\n (System.currentTimeMillis()-start) + " +
                "\" ms.\");\n");
        if (!"void".equals(type)) {
            body.append("return result;\n");
        }
        body.append("}");

        //  replace the body of the interceptor method with generated
        //  code block and add it to class
        mnew.setBody(body.toString());
        clas.addMethod(mnew);

        //  print the generated code block just to show what was done
        //System.out.println("Interceptor method body:");
        //System.out.println(body.toString());
    }

    private void addAnnouncement(CtBehavior method, CtClass cl)
            throws NotFoundException, CannotCompileException {
        String name = method.getName();
        method.insertAfter("logger.warn(\"\");");
    }
}


