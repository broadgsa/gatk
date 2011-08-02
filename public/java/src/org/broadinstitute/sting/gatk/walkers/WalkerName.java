package org.broadinstitute.sting.gatk.walkers;

import java.lang.annotation.*;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 26, 2009
 * Time: 3:00:16 PM
 * To change this template use File | Settings | File Templates.
 */
@Documented
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface WalkerName {
    public String value() default "";
}
