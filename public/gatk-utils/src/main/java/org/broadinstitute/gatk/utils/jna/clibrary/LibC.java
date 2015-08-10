/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.jna.clibrary;

import com.sun.jna.LastErrorException;
import com.sun.jna.Native;
import com.sun.jna.NativeLong;
import com.sun.jna.Structure;
import com.sun.jna.ptr.NativeLongByReference;

/**
 * Sparse port of the Standard C Library libc -lc.
 */
@SuppressWarnings("unused")
public class LibC {

    static {
        Native.register("c");
    }

    /** Operation not permitted */
    public static final int EPERM = 1;

    /** No such file or directory */
    public static final int ENOENT = 2;

    /** No such process */
    public static final int ESRCH = 3;

    /** Interrupted system call */
    public static final int EINTR = 4;

    /** I/O error */
    public static final int EIO = 5;

    /** No such device or address */
    public static final int ENXIO = 6;

    /** Argument list too long */
    public static final int E2BIG = 7;

    /** Exec format error */
    public static final int ENOEXEC = 8;

    /** Bad file number */
    public static final int EBADF = 9;

    /** No child processes */
    public static final int ECHILD = 10;

    /** Try again */
    public static final int EAGAIN = 11;

    /** Out of memory */
    public static final int ENOMEM = 12;

    /** Permission denied */
    public static final int EACCES = 13;

    /** Bad address */
    public static final int EFAULT = 14;

    /** Block device required */
    public static final int ENOTBLK = 15;

    /** Device or resource busy */
    public static final int EBUSY = 16;

    /** File exists */
    public static final int EEXIST = 17;

    /** Cross-device link */
    public static final int EXDEV = 18;

    /** No such device */
    public static final int ENODEV = 19;

    /** Not a directory */
    public static final int ENOTDIR = 20;

    /** Is a directory */
    public static final int EISDIR = 21;

    /** Invalid argument */
    public static final int EINVAL = 22;

    /** File table overflow */
    public static final int ENFILE = 23;

    /** Too many open files */
    public static final int EMFILE = 24;

    /** Not a typewriter */
    public static final int ENOTTY = 25;

    /** Text file busy */
    public static final int ETXTBSY = 26;

    /** File too large */
    public static final int EFBIG = 27;

    /** No space left on device */
    public static final int ENOSPC = 28;

    /** Illegal seek */
    public static final int ESPIPE = 29;

    /** Read-only file system */
    public static final int EROFS = 30;

    /** Too many links */
    public static final int EMLINK = 31;

    /** Broken pipe */
    public static final int EPIPE = 32;

    /** Math argument out of domain of func */
    public static final int EDOM = 33;

    /** Math result not representable */
    public static final int ERANGE = 34;

    /**
     * Inserts or resets the environment variable name in the current environment list.  If the variable name does not exist
     * in the list, it is inserted with the given value.  If the variable does exist, the argument overwrite is tested; if overwrite is zero, the
     * variable is not reset, otherwise it is reset to the given value.
     * @param name the environment variable name
     * @param value the given value
     * @param overwrite if overwrite is zero, the variable is not reset, otherwise it is reset to the given value
     * @return the value 0 if successful; otherwise the value -1 is returned and the global variable errno is set to indicate the error.
     * @throws LastErrorException [ENOMEM] The function failed because it was unable to allocate memory for the environment.
     */
    public static native int setenv(String name, String value, int overwrite) throws LastErrorException;

    /**
     * Obtains the current value of the environment variable, name.
     * @param name the environment variable name
     * @return the value of the environment variable as a NUL-terminated string.  If the variable name is not in the current environment, NULL is returned.
     */
    public static native String getenv(String name);

    /**
     * The unsetenv() function deletes all instances of the variable name pointed to by name from the list.  Note that only the variable name
     * (e.g., "NAME") should be given; "NAME=value" will not work.
     * @param name the environment variable name
     * @return the value 0 if successful; otherwise the value -1 is returned and the global variable errno is set to indicate the error.
     * @throws LastErrorException The function failed.
     */
    public static native int unsetenv(String name) throws LastErrorException;

    public static class timeval extends Structure {
        public static class ByReference extends timeval implements Structure.ByReference {
        }

        public static class ByValue extends timeval implements Structure.ByValue {
        }

        public NativeLong tv_sec;
        public NativeLong tv_usec;
    }

    /**
     * The time() function returns the value of time in seconds since 0 hours, 0 minutes, 0 seconds, January 1, 1970, Coordinated Universal Time, without including leap seconds.  If an error occurs, time() returns the value (time_t)-1.
     * The return value is also stored in *tloc, provided that t is non-null.
     * @param t the value of time in seconds,  provided that t is non-null.
     * @return the value of time in seconds
     */
    public static native NativeLong time(NativeLongByReference t);

    /**
     * Returns the difference between two calendar times, (time1 - time0), expressed in seconds.
     * @param time1 Time 1
     * @param time0 Time 0
     * @return the difference between two calendar times, (time1 - time0), expressed in seconds.
     */
    public static native double difftime(NativeLong time1, NativeLong time0);
}
