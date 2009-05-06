package org.broadinstitute.sting;

import org.apache.bcel.Constants;
import org.apache.bcel.Repository;
import org.apache.bcel.classfile.JavaClass;
import org.apache.bcel.generic.*;
import org.junit.Ignore;

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
        if (className.contains("broadinstitute") && !(className.endsWith("BaseTest"))) {
            JavaClass jclas = null;
            try {
                jclas = Repository.lookupClass(className);
                if (jclas == null) {
                    return null;
                }
            } catch (Exception e) {
                return null;
            }
            if (!(jclas.getSuperClass().getClassName().contains("BaseTest"))) {
                return null;
            }
            ClassGen cgen = new ClassGen(jclas);
            ConstantPoolGen pgen = cgen.getConstantPool();
            InstructionFactory fact = new InstructionFactory(cgen, pgen);
            createFields(cgen, pgen);
            createBeforeMethod(cgen, pgen, fact);
            createAfterMethod(cgen, pgen, fact);


            return cgen.getJavaClass().getBytes();
        }
        return null; // No transformation required
    }

    /**
     * create any fields we need here
     *
     * @param cgen our classgen for the class to add to
     * @param pgen our constant pool generator, check out the JVM spec about this
     */
    private void createFields(ClassGen cgen, ConstantPoolGen pgen) {
        FieldGen field;

        field = new FieldGen(Constants.ACC_PRIVATE | Constants.ACC_STATIC, Type.LONG, "startTime", pgen);
        cgen.addField(field.getField());
    }

    /**
     * create the before method
     *
     * @param cgen our classgen for the class to add to
     * @param pgen our constant pool generator
     * @param fact the instruction factory we're using
     */
    private void createBeforeMethod(ClassGen cgen, ConstantPoolGen pgen, InstructionFactory fact) {
        InstructionList il = new InstructionList();
        MethodGen method = new MethodGen(Constants.ACC_PUBLIC | Constants.ACC_FINAL, Type.VOID, Type.NO_ARGS, new String[]{}, "baseSetup", cgen.getClassName(), il, pgen);

        InstructionHandle ih_0 = il.append(fact.createInvoke("java.lang.System", "currentTimeMillis", Type.LONG, Type.NO_ARGS, Constants.INVOKESTATIC));
        il.append(fact.createFieldAccess(cgen.getClassName(), "startTime", Type.LONG, Constants.PUTSTATIC));
        InstructionHandle ih_6 = il.append(fact.createReturn(Type.VOID));

        method.setMaxStack();
        method.setMaxLocals();
        cgen.addMethod(method.getMethod());
        il.dispose();
    }

    /**
     * create the after method
     *
     * @param cgen our classgen for the class to add to
     * @param pgen our constant pool generator
     * @param fact the instruction factory we're using
     */
    private void createAfterMethod(ClassGen cgen, ConstantPoolGen pgen, InstructionFactory fact) {
        InstructionList il = new InstructionList();
        MethodGen method = new MethodGen(Constants.ACC_PUBLIC, Type.VOID, Type.NO_ARGS, new String[]{}, "baseTearDown", cgen.getClassName(), il, pgen);

        InstructionHandle ih_0 = il.append(fact.createInvoke("java.lang.System", "currentTimeMillis", Type.LONG, Type.NO_ARGS, Constants.INVOKESTATIC));
        il.append(fact.createStore(Type.LONG, 1));
        InstructionHandle ih_4 = il.append(fact.createFieldAccess(cgen.getClassName(), "logger", new ObjectType("org.apache.log4j.Logger"), Constants.GETSTATIC));
        il.append(new PUSH(pgen, cgen.getClassName() + " runtime: %dms"));
        il.append(new PUSH(pgen, 1));
        il.append(fact.createNewArray(Type.OBJECT, (short) 1));
        il.append(InstructionConstants.DUP);
        il.append(new PUSH(pgen, 0));
        il.append(fact.createLoad(Type.LONG, 1));
        il.append(fact.createFieldAccess(cgen.getClassName(), "startTime", Type.LONG, Constants.GETSTATIC));
        il.append(InstructionConstants.LSUB);
        il.append(fact.createInvoke("java.lang.Long", "valueOf", new ObjectType("java.lang.Long"), new Type[]{Type.LONG}, Constants.INVOKESTATIC));
        il.append(InstructionConstants.AASTORE);
        il.append(fact.createInvoke("java.lang.String", "format", Type.STRING, new Type[]{Type.STRING, new ArrayType(Type.OBJECT, 1)}, Constants.INVOKESTATIC));
        il.append(fact.createInvoke("org.apache.log4j.Logger", "warn", Type.VOID, new Type[]{Type.OBJECT}, Constants.INVOKEVIRTUAL));
        InstructionHandle ih_30 = il.append(fact.createReturn(Type.VOID));
        method.setMaxStack();
        method.setMaxLocals();
        cgen.addMethod(method.getMethod());
        il.dispose();
    }


}
