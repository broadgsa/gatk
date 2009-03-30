package org.broadinstitute.sting.utils.cmdLine;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 24, 2009
 * Time: 11:11:36 AM
 * To change this template use File | Settings | File Templates.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.FIELD})
public @interface Argument {
    String fullName() default "";
    String shortName() default "";
    String doc() default "";
    boolean required() default true;
    String exclusive() default "";
    String defaultValue() default "";
}
