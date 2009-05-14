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
            System.err.println("looking at " + className);
            ClassGen cgen = new ClassGen(jclas);
            ConstantPoolGen pgen = cgen.getConstantPool();
            InstructionFactory fact = new InstructionFactory(cgen, pgen);
            createFields(cgen, pgen);
            /*for (Method m : cgen.getMethods()) {
                System.err.println("looking at " + m.getName());
                addStringOutputToMethod(jclas, cgen, pgen, m, fact);
            }*/
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
        field = new FieldGen(Constants.ACC_PRIVATE, Type.STRING, "currentTestName", pgen);
        cgen.addField(field.getField());
    }

    /*
    private void addStringOutputToMethod(JavaClass classname, ClassGen cgen, ConstantPoolGen pgen, Method meth, InstructionFactory fact) {

        if(true) {return;}
        if (meth.getName().contains("<")) {
            System.err.println("Nope -> " + meth.getName());
            return;
        }
        //if (meth.isPublic()) {
        boolean outputInstead = true;
        MethodGen g = new MethodGen(meth, cgen.getClassName(), pgen);
        InstructionList il = g.getInstructionList();
        //if (outputInstead) {
        BufferedWriter outputStream = null;
        BufferedWriter outputStream2 = null;
        //}
        Instruction returnInstruction = null;
        InstructionHandle[] iHandles = il.getInstructionHandles();
        for (int f = 0; f < iHandles.length; f++) {
            if (iHandles[f].getInstruction() instanceof ReturnInstruction) {
                returnInstruction = iHandles[f].getInstruction();
                //System.out.println("found the invoke virtual");
                break;
            }
        }
        if (outputInstead) {
            try {
                outputStream =
                        new BufferedWriter(new FileWriter("one.txt"));
                outputStream2 =
                        new BufferedWriter(new FileWriter("two.txt"));

                outputStream.write(meth.getName() + " of " + meth.getClass());

                for (Instruction i : il.getInstructions()) {

                    outputStream.write(i.getName() + " <code> " + i.getOpcode() + " <toString> " + i.toString() + "\n");

                }
                outputStream.close();
            }
            catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        //InstructionHandle handle = il.getEnd();
        //il.s
        //il.insert(getFieldInstruction, fact.createLoad(Type.OBJECT, 0));
        //il.insert(getFieldInstruction, fact.createNew("java.lang.String"));
        //il.insert(getFieldInstruction, InstructionConstants.DUP);
        //il.insert(getFieldInstruction, new PUSH(pgen, meth.getName()));
        //il.insert(getFieldInstruction, fact.createInvoke("java.lang.String", "<init>", Type.VOID, new Type[]{Type.STRING}, Constants.INVOKESPECIAL));
        //il.insert(getFieldInstruction, fact.createFieldAccess(cgen.getClassName(), "currentTestName", Type.STRING, Constants.PUTFIELD));
        //il.insert(getFieldInstruction,fact.createPrintln("Hello World"));
        /*il.insert(returnInstruction, new ALOAD(0));
        il.insert(returnInstruction, fact.createNew("java.lang.String"));
        il.insert(returnInstruction, InstructionConstants.DUP);
        il.insert(returnInstruction, new PUSH(pgen, meth.getName()));
        il.insert(returnInstruction, fact.createInvoke("java.lang.String", "<init>", Type.VOID, new Type[]{Type.STRING}, Constants.INVOKESPECIAL));
        il.insert(returnInstruction, fact.createFieldAccess(classname.replace("/","."), "currentTestName", Type.STRING, Constants.PUTFIELD));*/
        /*il.setPositions();
        g.setMaxStack();
        g.setMaxLocals();
        g.removeLineNumbers();
        //org.apache.bcel.classfile.LocalVariableTypeTable table;
        InstructionList inst = g.getInstructionList();
        if (outputInstead) {
            try {
                outputStream2.write(meth.getName() + " of " + meth.getClass() + " classname: " + classname.getClassName() + "\n");

                for (Instruction i : inst.getInstructions()) {

                    outputStream2.write(i.getName() + " <code> " + i.getOpcode() + " <toString> " + i.toString() + "\n");

                }
                outputStream2.close();
            }
            catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        //.
        cgen.replaceMethod(meth, g.getMethod());
        il.dispose();*/
        //}
        /*if (meth.isPublic()) {
           InstructionList il = new InstructionList();
           MethodGen method = new MethodGen(Constants.ACC_PUBLIC, Type.VOID, Type.NO_ARGS, new String[]{}, meth.getName(), cgen.getClassName(), il, pgen);

           InstructionHandle ih_0 = il.append(fact.createLoad(Type.OBJECT, 0));
           il.append(fact.createNew("java.lang.String"));
           il.append(InstructionConstants.DUP);
           il.append(new PUSH(pgen, "grrr"));
           il.append(fact.createInvoke("java.lang.String", "<init>", Type.VOID, new Type[]{Type.STRING}, Constants.INVOKESPECIAL));
           il.append(fact.createFieldAccess(cgen.getClassName(), "currentTestName", Type.STRING, Constants.PUTFIELD));
           InstructionHandle ih_13 = il.append(fact.createReturn(Type.VOID));
           method.setMaxStack();
           method.setMaxLocals();
           cgen.removeMethod(meth);
           cgen.addMethod(method.getMethod());
           il.dispose();
       }

    }*/


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
        il.append(fact.createNew("java.lang.StringBuilder"));
        il.append(InstructionConstants.DUP);
        il.append(fact.createInvoke("java.lang.StringBuilder", "<init>", Type.VOID, Type.NO_ARGS, Constants.INVOKESPECIAL));
        il.append(new PUSH(pgen, "Test Name: "));
        il.append(fact.createInvoke("java.lang.StringBuilder", "append", new ObjectType("java.lang.StringBuilder"), new Type[]{Type.STRING}, Constants.INVOKEVIRTUAL));
        il.append(fact.createLoad(Type.OBJECT, 0));
        il.append(fact.createFieldAccess(cgen.getClassName(), "currentTestName", Type.STRING, Constants.GETFIELD));
        il.append(fact.createInvoke("java.lang.StringBuilder", "append", new ObjectType("java.lang.StringBuilder"), new Type[]{Type.STRING}, Constants.INVOKEVIRTUAL));
        il.append(new PUSH(pgen, " runtime: %dms"));
        il.append(fact.createInvoke("java.lang.StringBuilder", "append", new ObjectType("java.lang.StringBuilder"), new Type[]{Type.STRING}, Constants.INVOKEVIRTUAL));
        il.append(fact.createInvoke("java.lang.StringBuilder", "toString", Type.STRING, Type.NO_ARGS, Constants.INVOKEVIRTUAL));
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
        InstructionHandle ih_55 = il.append(fact.createReturn(Type.VOID));
        method.setMaxStack();
        method.setMaxLocals();
        cgen.addMethod(method.getMethod());
        il.dispose();
    }


}
