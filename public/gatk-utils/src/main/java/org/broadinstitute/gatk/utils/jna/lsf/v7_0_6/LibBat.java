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

package org.broadinstitute.gatk.utils.jna.lsf.v7_0_6;

import com.sun.jna.*;
import com.sun.jna.ptr.*;
import org.broadinstitute.gatk.utils.jna.clibrary.JNAUtils;
import org.broadinstitute.gatk.utils.jna.clibrary.LibC;

/*
  NOTE: This library uses Pointer for some Struct.ByReference members going
  against the JNA recommendations at http://jna.java.net/#structure_use
  Instead stuct arrays are Pointers and each structure contains a
  constructor that can accept the Pointer iff the size of the array is
  known to be greater than zero.

  This was especially problematic in jobInfoEnt->items->resName. When
  jobInfo->reserveCnt was zero jobInfoItems->items was not necessarily null.

  LSF will often reuse memory for structure arrays but will set the
  array size / count (reserveCnt above) to zero when the array should
  not be accessed. When LSF has reused memory and points to a non-null
  structure pointer (items) the inner structure may contain further
  garbage pointers (especially items->resName).

  When JNA sees a non-null Structure.ByReference it will autoRead() the
  member. When autoRead() eventually gets to the items->resName trying
  to run strlen on the bad memory address causes a SIGSEGV.

  By using a Pointer instead of the Structure.ByReference JNA will not
  automatically autoRead(), and the API user will have to pass the
  pointer to the Structure on their own.
*/

/**
 * JNA wrappers for LSF's lsbatch.h and -lbat
 *
 * $Id: lsbatch.h,v 2.1043 2009/08/06 16:50:49 bxia Exp $
 * -----------------------------------------------------------------
 *
 *  Lsbatch Distributed Batch Utility --
 *
 *  Header file for all lsbatch components: applications, lsblib,
 *                                          mbatchd and sbatchd
 *
 * ------------------------------------------------------------------
 */
@SuppressWarnings("unused")
public class LibBat {

    static {
        // via Platform LSF Configuration Reference, by default quiet the BSUB output.
        if ("Y".equals(System.getProperty("BSUB_QUIET", "Y")))
            LibC.setenv("BSUB_QUIET", "Y", 1);
        String lsfLibDir = System.getenv("LSF_LIBDIR");
        if (lsfLibDir != null) {
            NativeLibrary.addSearchPath("lsf", lsfLibDir);
            NativeLibrary.addSearchPath("bat", lsfLibDir);
        }
        /*
        LSF 7.0.6 on the mac is missing the unsatisfied exported symbol for environ which was removed on MacOS X 10.5+.
        nm $LSF_LIBDIR/liblsf.dylib | grep environ
        See "man environ" for more info, along with http://lists.apple.com/archives/java-dev/2007/Dec/msg00096.html
        For now, we export environ ourselves using libenvironhack.dylib available in c/libenvironhack.
        */
        if (Platform.isMac())
            NativeLibrary.getInstance("environhack");
        NativeLibrary liblsf = NativeLibrary.getInstance("lsf");
        Native.register("bat");
        // HACK: Running into a weird error:
        //   java.lang.UnsatisfiedLinkError: Unable to load library 'bat': <$LSF_LIBDIR>/libbat.so: undefined symbol: xdr_resourceInfoReq
        // This function is very clearly unsatisfied by running 'nm $LSF_LIBDIR/libbat.so | grep xdr_resourceInfoReq' but is
        // found in liblsf.so when running 'nm $LSF_LIBDIR/liblsf.so | grep xdr_resourceInfoReq'. For now holding on to a reference
        // to the LSF lib just in case this is a problem with the NativeLibrary's internal WeakReferences and the library being unloaded?
        liblsf.getFunction("xdr_resourceInfoReq").getName();
    }

    // Via support@platform.com:
    //    For equivalent api of bsub -a "xxx aaa qqq", option -a is not in struct submit, we
    //    have to use setOption_ to set it. setOption_ can be used in user program by including
    //    cmd.h or opensource.h of LSF opensource. You can refer to cmd.sub.c in opensource.
    //
    //    Here is a demonstration on the api for bsub -a
    //    =========================================================================
    //    /*define external setOption_ function*/
    //    extern int setOption_(int argc, char **argv, char *template,
    //    struct submit *req, int mask, int mask2, char **errMsg);
    //
    //    int setEsub(char *esub, struct submit *req) {
    //    int x;
    //    char *template, *arg[3];
    //    /*set esub with the following strings and set array length*/
    //    arg[0] = "blah";
    //    arg[1] = "-a";
    //    arg[2] = test;
    //    /* -a "test", You can add additional esubs in here.  Just make sure they're space delimited.  ie. "test mpich lammpi" */
    //    x=3;
    //    /*set template*/
    //    template = "a:"
    //    /*run setOption_()*/
    //    if (setOption_(x, arg, template, req, ~0, ~0, ~0, NULL) == -1) {
    //    return(-1);
    //    }
    //    else {
    //    return(0);
    //    }
    //    }
    //    =========================================================================

    /**
     * Used for setting esub and other options not in struct submit.
     * Via support@platform.com
     *
     * @param argc number of args
     * @param argv arguments including a first argument that will not be used
     * @param template a colon delimited list of arguments in getopt format
     * @param jobSubReq the lsf submit
     * @param mask unknown
     * @param mask2 unknown
     * @param mask3 unknown
     * @param errMsg unknown
     * @return -1 if the option setting failed
     */
    public static native int setOption_(int argc, Pointer argv, String template, submit jobSubReq, int mask, int mask2, int mask3, Pointer errMsg);

    /** Max job name length as defined by 'man bsub'. */
    public static final int MAX_JOB_NAME_LEN = 4094;

/* if only everyone had <paths.h> */
    public static final String _PATH_NULL = "/dev/null";

    //public static int SKIP_SPACES (int word)  { while (word[0] == ' ' )  word++; }

    //public static void FREEUP_ARRAY(int num, Pointer vector) { FREE_STRING_VECTOR_ENTRIES(num, vector);  FREEUP(vector); }


/* event log version:
*  each new major release requires to add a new line
 */
    public static final float LSB_EVENT_VERSION3_0 = 3.0f;
    public static final float LSB_EVENT_VERSION3_1 = 3.1f;
    public static final float LSB_EVENT_VERSION3_2 = 3.2f;
    public static final float LSB_EVENT_VERSION4_0 = 4.0f;
    public static final float LSB_EVENT_VERSION4_1 = 4.1f;
    public static final float LSB_EVENT_VERSION4_2 = 4.2f;
    public static final float LSB_EVENT_VERSION5_0 = 5.0f;
    public static final float LSB_EVENT_VERSION5_1 = 5.1f;
    public static final float LSB_EVENT_VERSION6_0 = 6.0f;
    public static final float LSB_EVENT_VERSION6_1 = 6.1f;
    public static final float LSB_EVENT_VERSION6_2 = 6.2f;
    public static final float LSB_EVENT_VERSION7_0 = 7.0f;
    public static final float LSB_EVENT_VERSION7_0_1 = 7.01f;
    public static final float LSB_EVENT_VERSION7_0_2 = 7.02f;
    public static final float LSB_EVENT_VERSION7_0_3 = 7.03f;
    public static final float LSB_EVENT_VERSION7_0_4 = 7.04f;
    public static final float LSB_EVENT_VERSION7_0_5 = 7.05f;
    public static final float LSB_EVENT_VERSION7_0_6 = 7.06f;

/* current event version number of the mbatchd */
    public static final String THIS_VERSION = "7.06";

    public static final int MAX_VERSION_LEN = 12;

/* num of users per host partition */
    public static final int MAX_HPART_USERS = 100;

/* max byte limit, OS independent */
    public static final int MAX_CHARLEN = 20;

/* the max length of name */
    public static final int MAX_LSB_NAME_LEN = 60;

/*the max length of user group*/
    public static final int MAX_LSB_UG_NAME_LEN = 512;

/*Maximum levels that a user group hierachy can have*/
    public static final int MAX_LSB_UG_HIERDEPTH = 64;

/* the max length of command */
    public static final int MAX_CMD_DESC_LEN = 512;

/* for the local cluster */
    public static final int MAX_CALENDARS = 256;

/* max num of user equivalent entries */
    public static final int MAX_USER_EQUIVALENT = 128;

/* max num of user mapping entries */
    public static final int MAX_USER_MAPPING = 128;

/* max external msg's description length */
    public static final int MAXDESCLEN = 20 * 512;

/* num of user or host group */
    public static final int MAX_GROUPS = 1024;

/*
*  RFC #725
 */

/* max len. of a filename */
    public static final int MAXFULLFILENAMELEN = 4096;
    public static final int MAXFULLPATHNAMELEN = 2 * MAXFULLFILENAMELEN;
    public static final int FILENAMEPADDING = 128;

    public static final String DEFAULT_MSG_DESC = "no description";

    public static final int MSGSIZE = 4096;

/* RFC #725
*  extend the MSG size to 4*max filename len
 */
    public static final int MAXFULLMSGSIZE = 4 * MAXFULLFILENAMELEN;

/* host status (hStatus) bits */
    /**
     *  \addtogroup host_status host_status
     *  The status of the host. It is the bitwise inclusive OR of some of the following:
     */

    /**
     * < Ready to accept and run jobs
     */
    public static final int HOST_STAT_OK = 0x0;

/* Load is not good enough */
    public static final int HOST_STAT_BUSY = 0x01;
    /**
     * < The host load is greater than a scheduling threshold. In this status, no new job will be scheduled to run on this host.
     */

/* Run windows are closed */
    public static final int HOST_STAT_WIND = 0x02;
    /**
     * < The host dispatch window is closed. In this status, no new job will be accepted.
     */

/* Disabled by admin */
    public static final int HOST_STAT_DISABLED = 0x04;
    /**
     * < The host has been disabled by the LSF administrator and will not accept jobs. In this status, no new job will be scheduled to  run on this host.
     */

/* Lim locked by admin */
    public static final int HOST_STAT_LOCKED = 0x08;
    /**< The host is locked by a exclusive task. In this status, no new job will be scheduled to run on this host.*/

    /**
     * < Great than job limit
     */
    public static final int HOST_STAT_FULL = 0x10;
    /**< The host has reached its job limit. In this status, no new job will be scheduled to run on this host.*/

    /**
     * < The sbatchd on this host is unreachable.
     */
    public static final int HOST_STAT_UNREACH = 0x20;

    /**
     * < The LIM and sbatchd on this host are unavailable.
     */
    public static final int HOST_STAT_UNAVAIL = 0x40;

    /**
     * < The host does not have an LSF license.
     */
    public static final int HOST_STAT_UNLICENSED = 0x80;

    /**
     * < The host is running an sbatchd but not a LIM.
     */
    public static final int HOST_STAT_NO_LIM = 0x100;

    /**
     * < Running exclusive job
     */
    public static final int HOST_STAT_EXCLUSIVE = 0x200;

    /**
     * < Lim locked by master LIM
     */
    public static final int HOST_STAT_LOCKED_MASTER = 0x400;

    /**
     * < Close a remote lease host. This flag is  used together with HOST_STAT_DISABLED.
     */
    public static final int HOST_STAT_REMOTE_DISABLED = 0x800;

    /**
     * < Close a remote lease host due to the  lease is renewing or terminating.
     */
    public static final int HOST_STAT_LEASE_INACTIVE = 0x1000;

/* if LSF_HPC_EXTENTIONS="LSB_HCLOSE_BY_RES" is set in lsf.conf
*  host will be closed if RES is unavailable.
 */

    /**
     * < Host is disabled by RES
     */
    public static final int HOST_STAT_DISABLED_RES = 0x4000;

/* Kite#29531 a bit set in hData->hStatus
*  to show whether the host is closed by
*  admin or closed because RMS is not available.
 */

    /**
     * < Host is disabled by RMS
     */
    public static final int HOST_STAT_DISABLED_RMS = 0x8000;

/* lsf70 project scheduling, a removed host from mbatchd move into
*  a new status HOST_STAT_LOCKED_EGO
 */

    /**
     * < The host is disabled by EGO
     */
    public static final int HOST_STAT_LOCKED_EGO = 0x10000;

    /**
     * < If none of the above hold, hStatus is set to HOST_STAT_OK to indicate that the host is ready to accept and run jobs.
     */
    public static final int HOST_CLOSED_BY_ADMIN = 0x20000;

    /**
     * < Running cu exclusive job
     */
    public static final int HOST_STAT_CU_EXCLUSIVE = 0x40000;

/* host is ok */

    public static boolean LSB_HOST_OK(int status) {
        return (status == HOST_STAT_OK);
    }

/* host is busy */

    public static boolean LSB_HOST_BUSY(int status) {
        return ((status & HOST_STAT_BUSY) != 0);
    }

/* host is closed */

    public static boolean LSB_HOST_CLOSED(int status) {
        return ((status & (HOST_STAT_WIND | HOST_STAT_DISABLED | HOST_STAT_LOCKED | HOST_STAT_LOCKED_MASTER | HOST_STAT_FULL | HOST_STAT_CU_EXCLUSIVE | HOST_STAT_EXCLUSIVE | HOST_STAT_LEASE_INACTIVE | HOST_STAT_NO_LIM)) != 0);
    }

/* host is full */

    public static boolean LSB_HOST_FULL(int status) {
        return ((status & HOST_STAT_FULL) != 0);
    }

/* host is unlicensed */

    public static boolean LSB_HOST_UNLICENSED(int status) {
        return ((status & HOST_STAT_UNLICENSED) != 0);
    }

/* host is unreach */

    public static boolean LSB_HOST_UNREACH(int status) {
        return ((status & HOST_STAT_UNREACH) != 0);
    }

/* host is unavail */

    public static boolean LSB_HOST_UNAVAIL(int status) {
        return ((status & HOST_STAT_UNAVAIL) != 0);
    }


    /* host busy reason bits */
    /**
     *  \addtogroup host_load_BusyReason host_load_BusyReason
     *  If hStatus is HOST_STAT_BUSY, these indicate the host loadSched or loadStop
     *  busy reason. If none of the thresholds have been exceeded, the value is
     *  HOST_BUSY_NOT. Otherwise the value is the bitwise inclusive OR of some of the
     *  following:
     */

    /**
     * < Host not busy
     */
    public static final int HOST_BUSY_NOT = 0x000;

    /**
     * < The 15 second average CPU run queue length is too high.
     */
    public static final int HOST_BUSY_R15S = 0x001;

    /**
     * < The 1 minute average CPU run queue length is too high.
     */
    public static final int HOST_BUSY_R1M = 0x002;

    /**
     * < The 15 minute average CPU run queue length is too high.
     */
    public static final int HOST_BUSY_R15M = 0x004;

    /**
     * < The CPU utilization is too high.
     */
    public static final int HOST_BUSY_UT = 0x008;

    /**
     * < The paging rate is too high.
     */
    public static final int HOST_BUSY_PG = 0x010;

    /**
     * < The I/O rate is too high.
     */
    public static final int HOST_BUSY_IO = 0x020;

    /**
     * < There are too many login sessions.
     */
    public static final int HOST_BUSY_LS = 0x040;

    /**
     * < Host has not been idle long enough.
     */
    public static final int HOST_BUSY_IT = 0x080;

    /**
     * < There is not enough free space in the file  system containing /tmp.
     */
    public static final int HOST_BUSY_TMP = 0x100;

    /**
     * < There is not enough free swap space.
     */
    public static final int HOST_BUSY_SWP = 0x200;

    /**
     * < There is not enough free memory.
     */
    public static final int HOST_BUSY_MEM = 0x400;

/* host is busy */

    public static boolean LSB_ISBUSYON(int[] status, int index) {
        return (((status[(index) / LibLsf.INTEGER_BITS]) & (1 << (index) % LibLsf.INTEGER_BITS)) != 0);
    }


/* queue status (qStatus) bits */
    /**
     *  \addtogroup queue_status queue_status
     *  queue status (qStatus) bits
     */

    /**
     * < The queue is open to accept newly submitted jobs.
     */
    public static final int QUEUE_STAT_OPEN = 0x01;

    /**
     * < The queue is actively dispatching jobs. The queue can be inactivated and  reactivated by the LSF administrator using  \ref lsb_queuecontrol. The queue will also be inactivated when its run or dispatch window  is closed. In this case it cannot be reactivated manually; it will be reactivated by the LSF system when its run and dispatch windows reopen.
     */
    public static final int QUEUE_STAT_ACTIVE = 0x02;

    /**
     * < The queue run and dispatch windows are open. The initial state of a queue at LSF boot time is open and either active or inactive, depending on its run and dispatch windows.
     */
    public static final int QUEUE_STAT_RUN = 0x04;

    /**
     * < Remote queue rejecting jobs.
     */
    public static final int QUEUE_STAT_NOPERM = 0x08;

    /**
     * < Remote queue status is disconnected.
     */
    public static final int QUEUE_STAT_DISC = 0x10;

    /**
     * < Queue run windows are closed.
     */
    public static final int QUEUE_STAT_RUNWIN_CLOSE = 0x20;

/* queue attribute (QAttrib) bits */
    /**
     *  \addtogroup queue_attribute queue_attribute
     *  queue attribute (QAttrib) bits.
     */

    /**
     * < This queue accepts jobs which request exclusive execution.
     */
    public static final int Q_ATTRIB_EXCLUSIVE = 0x01;

    /**
     * < This queue is a default LSF queue.
     */
    public static final int Q_ATTRIB_DEFAULT = 0x02;

    /**
     * < This queue uses the FAIRSHARE scheduling policy. The user shares  are given in userShares.
     */
    public static final int Q_ATTRIB_FAIRSHARE = 0x04;

    /**
     * < This queue uses the PREEMPTIVE scheduling policy.
     */
    public static final int Q_ATTRIB_PREEMPTIVE = 0x08;

    /**
     * < This is an NQS forward queue. The target NQS queues are given in nqsQueues. For NQS forward queues, the hostList, procJobLimit, windows, mig and windowsD fields are meaningless.
     */
    public static final int Q_ATTRIB_NQS = 0x10;

    /**
     * < This queue can receive jobs from other clusters
     */
    public static final int Q_ATTRIB_RECEIVE = 0x20;

    /**
     * < This queue uses a preemptable scheduling policy.
     */
    public static final int Q_ATTRIB_PREEMPTABLE = 0x40;

    /**
     * < This queue uses a backfilling policy.
     */
    public static final int Q_ATTRIB_BACKFILL = 0x80;

    /**
     * < This queue uses a host preference policy.
     */
    public static final int Q_ATTRIB_HOST_PREFER = 0x100;

    /**
     * < This queue can't preempt any other another queue.
     */
    public static final int Q_ATTRIB_NONPREEMPTIVE = 0x200;

    /**
     * < This queue can't be preempted from any queue.
     */
    public static final int Q_ATTRIB_NONPREEMPTABLE = 0x400;

    /**
     * < This queue does not accept batch interactive jobs.
     */
    public static final int Q_ATTRIB_NO_INTERACTIVE = 0x800;

    /**
     * < This queue only accepts batch interactive jobs.
     */
    public static final int Q_ATTRIB_ONLY_INTERACTIVE = 0x1000;

    /**
     * < No host type related resource name specified in resource requirement.
     */
    public static final int Q_ATTRIB_NO_HOST_TYPE = 0x2000;

    /**
     * < This queue disables deadline constrained resource scheduling.
     */
    public static final int Q_ATTRIB_IGNORE_DEADLINE = 0x4000;

    /**
     * < Jobs may run as chkpntable.
     */
    public static final int Q_ATTRIB_CHKPNT = 0x8000;

    /**
     * < Jobs may run as rerunnable.
     */
    public static final int Q_ATTRIB_RERUNNABLE = 0x10000;

    /**
     * < Excluding remote jobs when local jobs are present in the queue.
     */
    public static final int Q_ATTRIB_EXCL_RMTJOB = 0x20000;

    /**
     * < Turn on a multicluster fast scheduling policy.
     */
    public static final int Q_ATTRIB_MC_FAST_SCHEDULE = 0x40000;

    /**
     * < Push interactive jobs in front of other jobs in queue.
     */
    public static final int Q_ATTRIB_ENQUE_INTERACTIVE_AHEAD = 0x80000;

/* Only one of the following four flags could be TRUE. By default, the queue
*  is a local queue only(none of them is set.)
*      0x100000 - 0xf00000 is used for MC attribute
 */


    /**
     * < Flags used by MultiCluster.
     */
    public static final int Q_MC_FLAG = 0xf00000;

    /**
     * < Lease and local.
     */
    public static final int Q_ATTRIB_LEASE_LOCAL = 0x100000;

    /**
     * < Lease only; no local.
     */
    public static final int Q_ATTRIB_LEASE_ONLY = 0x200000;

    /**
     * < Remote batch and local.
     */
    public static final int Q_ATTRIB_RMT_BATCH_LOCAL = 0x300000;

    /**
     * < Remote batch only.
     */
    public static final int Q_ATTRIB_RMT_BATCH_ONLY = 0x400000;


    /**
     * < Memory reservation.
     */
    public static final int Q_ATTRIB_RESOURCE_RESERVE = 0x1000000;

    /**
     * < Cross-queue fairshare.
     */
    public static final int Q_ATTRIB_FS_DISPATCH_ORDER_QUEUE = 0x2000000;

    /**
     * < Batch queue/partition
     */
    public static final int Q_ATTRIB_BATCH = 0x4000000;

    /**
     * < Online partition
     */
    public static final int Q_ATTRIB_ONLINE = 0x8000000;

    /**
     * < Interruptible backfill
     */
    public static final int Q_ATTRIB_INTERRUPTIBLE_BACKFILL = 0x10000000;

    /**
     * < Absolute Priority scheduling (APS) value.
     */
    public static final int Q_ATTRIB_APS = 0x20000000;

    /**
     * < No queue with RESOURCE_RESERVE or SLOT_RESERVE has higher priority than this queue.
     */
    public static final int Q_ATTRIB_NO_HIGHER_RESERVE = 0x40000000;

    /**
     * < No host valid
     */
    public static final int Q_ATTRIB_NO_HOST_VALID = 0x80000000;


/* macros to check queue near real time attributes */

    public static int IS_ONLINE_QUEUE(queueInfoEnt Q) {
        return (Q.qAttrib & Q_ATTRIB_ONLINE);
    }

    public static int IS_BATCH_QUEUE(queueInfoEnt Q) {
        return (Q.qAttrib & Q_ATTRIB_BATCH);
    }

/* macros to check queue remote attributes */

    public static boolean IS_LEASE_LOCAL_QUEUE(queueInfoEnt Q) {
        return ((Q.qAttrib & Q_MC_FLAG) == Q_ATTRIB_LEASE_LOCAL);
    }

    public static boolean IS_LEASE_ONLY_QUEUE(queueInfoEnt Q) {
        return ((Q.qAttrib & Q_MC_FLAG) == Q_ATTRIB_LEASE_ONLY);
    }

    public static boolean IS_RMT_BATCH_LOCAL_QUEUE(queueInfoEnt Q) {
        return ((Q.qAttrib & Q_MC_FLAG) == Q_ATTRIB_RMT_BATCH_LOCAL);
    }

    public static boolean IS_RMT_BATCH_ONLY_QUEUE(queueInfoEnt Q) {
        return ((Q.qAttrib & Q_MC_FLAG) == Q_ATTRIB_RMT_BATCH_ONLY);
    }

    public static boolean IS_LEASE_QUEUE(queueInfoEnt Q) {
        return (IS_LEASE_LOCAL_QUEUE(Q) || IS_LEASE_ONLY_QUEUE(Q));
    }

    public static boolean IS_RMT_BATCH_QUEUE(queueInfoEnt Q) {
        return (IS_RMT_BATCH_LOCAL_QUEUE(Q) || IS_RMT_BATCH_ONLY_QUEUE(Q));
    }

    public static boolean IS_MC_QUEUE(queueInfoEnt Q) {
        return (IS_LEASE_QUEUE(Q) || IS_RMT_BATCH_QUEUE(Q));
    }

    public static int SET_LEASE_LOCAL_QUEUE(queueInfoEnt Q) {
        return (Q.qAttrib |= Q_ATTRIB_LEASE_LOCAL);
    }

    public static int SET_LEASE_ONLY_QUEUE(queueInfoEnt Q) {
        return (Q.qAttrib |= Q_ATTRIB_LEASE_ONLY);
    }

    public static int SET_RMT_BATCH_LOCAL_QUEUE(queueInfoEnt Q) {
        return (Q.qAttrib |= Q_ATTRIB_RMT_BATCH_LOCAL);
    }

    public static int SET_RMT_BATCH_ONLY_QUEUE(queueInfoEnt Q) {
        return (Q.qAttrib |= Q_ATTRIB_RMT_BATCH_ONLY);
    }

    public static int CLR_MC_QUEUE_FLAG(queueInfoEnt Q) {
        return (Q.qAttrib &= ~Q_MC_FLAG);
    }


/* the bits 0x10000000 to 0x80000000 is reserved for internal use (daemons.h) */

/* exit code for mbatchd */
    public static final int MASTER_NULL = 200;
    public static final int MASTER_RESIGN = 201;
    public static final int MASTER_RECONFIG = 202;
    public static final int MASTER_FATAL = 203;
    public static final int MASTER_MEM = 204;
    public static final int MASTER_CONF = 205;
    public static final int MASTER_EVENT = 206;
    public static final int MASTER_DISABLE = 207;

/* sub type of mbatchd die */
    public static final int MBD_USER_CMD = 1;
    public static final int MBD_NON_USER_CMD = 2;

    /**
     *  \addtogroup job_states job_states
     *  define job states
     */

    /**
     * < State null
     */
    public static final int JOB_STAT_NULL = 0x00;

    /**
     * < The job is pending, i.e., it  has not been dispatched yet.
     */
    public static final int JOB_STAT_PEND = 0x01;

    /**
     * < The pending job was suspended by its owner or the LSF system administrator.
     */
    public static final int JOB_STAT_PSUSP = 0x02;

    /**
     * < The job is running.
     */
    public static final int JOB_STAT_RUN = 0x04;

    /**
     * < The running job was suspended  by the system because an execution  host was overloaded or the queue run  window closed. (see \ref lsb_queueinfo,  \ref lsb_hostinfo, and lsb.queues.)
     */
    public static final int JOB_STAT_SSUSP = 0x08;

    /**
     * < The running job was suspended by its owner or the LSF systemadministrator.
     */
    public static final int JOB_STAT_USUSP = 0x10;

    /**
     * < The job has terminated with a non-zero status - it may have been aborted due  to an error in its execution, or  killed by its owner or by the  LSF system administrator.
     */
    public static final int JOB_STAT_EXIT = 0x20;

    /**
     * < The job has terminated with status 0.
     */
    public static final int JOB_STAT_DONE = 0x40;

    /**
     * < Post job process done successfully
     */
    public static final int JOB_STAT_PDONE = (0x80);

    /**
     * < Post job process has error
     */
    public static final int JOB_STAT_PERR = (0x100);

    /**
     * < Chunk job waiting its turn to exec
     */
    public static final int JOB_STAT_WAIT = (0x200);

    /**
     * < The slave batch daemon (sbatchd) on  the host on which the job is processed  has lost contact with the master batch  daemon (mbatchd).
     */
    public static final int JOB_STAT_UNKWN = 0x10000;

    /**
     *  \addtogroup event_types event_types
     *  define statements used by \ref lsb_geteventrec. Events logged in lsb.events file
     */

    /**
     * < Submit a new job
     */
    public static final int EVENT_JOB_NEW = 1;

    /**
     * < mbatchd is trying to start a job
     */
    public static final int EVENT_JOB_START = 2;

    /**
     * < Job's status change event
     */
    public static final int EVENT_JOB_STATUS = 3;

    /**
     * < Job switched to another queue
     */
    public static final int EVENT_JOB_SWITCH = 4;

    /**
     * < Move a pending job's position within a queue
     */
    public static final int EVENT_JOB_MOVE = 5;

    /**
     * < Queue's status changed by Platform LSF  administrator (bhc operation)
     */
    public static final int EVENT_QUEUE_CTRL = 6;

    /**
     * < Host status changed by Platform LSF  administrator (bhc operation)
     */
    public static final int EVENT_HOST_CTRL = 7;

    /**
     * < Log parameters before mbatchd died
     */
    public static final int EVENT_MBD_DIE = 8;

    /**
     * < Action that was not taken because the  mbatchd was unable to contact the sbatchd on the job's execution host
     */
    public static final int EVENT_MBD_UNFULFILL = 9;

    /**
     * < Job finished (Logged in lsb.acct)
     */
    public static final int EVENT_JOB_FINISH = 10;

    /**
     * < The complete list of load indices, including external load indices
     */
    public static final int EVENT_LOAD_INDEX = 11;

    /**
     * < Job checkpointed.
     */
    public static final int EVENT_CHKPNT = 12;

    /**
     * < Job migrated
     */
    public static final int EVENT_MIG = 13;

    /**
     * < The pre-execution command started
     */
    public static final int EVENT_PRE_EXEC_START = 14;

    /**
     * < New mbatchd start event
     */
    public static final int EVENT_MBD_START = 15;

    /**
     * < The job has been routed to NQS
     */
    public static final int EVENT_JOB_ROUTE = 16;

    /**
     * < Job modification request
     */
    public static final int EVENT_JOB_MODIFY = 17;

    /**
     * < Signal/delete a job
     */
    public static final int EVENT_JOB_SIGNAL = 18;

    /**
     * < Add new calendar to the system
     */
    public static final int EVENT_CAL_NEW = 19;

    /**
     * < Calendar modified
     */
    public static final int EVENT_CAL_MODIFY = 20;

    /**
     * < Delete a calendar in the system
     */
    public static final int EVENT_CAL_DELETE = 21;

    /**
     * < Job forwarded to another cluster
     */
    public static final int EVENT_JOB_FORWARD = 22;

    /**
     * < Job from a remote cluster dispatched
     */
    public static final int EVENT_JOB_ACCEPT = 23;

    /**
     * < Job status successfully sent to  submission cluster
     */
    public static final int EVENT_STATUS_ACK = 24;

    /**
     * < Job started successfully on the  execution host
     */
    public static final int EVENT_JOB_EXECUTE = 25;

    /**
     * < Send a message to a job
     */
    public static final int EVENT_JOB_MSG = 26;

    /**
     * < The message has been delivered
     */
    public static final int EVENT_JOB_MSG_ACK = 27;

    /**
     * < Job is requeued
     */
    public static final int EVENT_JOB_REQUEUE = 28;

    /**
     * < Submission mbatchd logs this after sending  an occupy request to execution mbatchd
     */
    public static final int EVENT_JOB_OCCUPY_REQ = 29;

    /**
     * < Submission mbatchd logs this event after  all execution mbatchds have vacated the occupied hosts for the job
     */
    public static final int EVENT_JOB_VACATED = 30;

    /**
     * < A signal action on a job has been  initiated or finished
     */
    public static final int EVENT_JOB_SIGACT = 32;

    /**
     * < sbatchd's new job status
     */
    public static final int EVENT_SBD_JOB_STATUS = 34;

    /**
     * < sbatchd accepts job start
     */
    public static final int EVENT_JOB_START_ACCEPT = 35;

    /**
     * < Undelete a calendar in the system
     */
    public static final int EVENT_CAL_UNDELETE = 36;

    /**
     * < Job is cleaned out of the core
     */
    public static final int EVENT_JOB_CLEAN = 37;

    /**
     * < Job exception was detected
     */
    public static final int EVENT_JOB_EXCEPTION = 38;

    /**
     * < Adding a new job group
     */
    public static final int EVENT_JGRP_ADD = 39;

    /**
     * < Modifying a job group
     */
    public static final int EVENT_JGRP_MOD = 40;

    /**
     * < Controlling a job group
     */
    public static final int EVENT_JGRP_CTRL = 41;

    /**
     * < Forcing a job to start on specified  hosts (brun operation)
     */
    public static final int EVENT_JOB_FORCE = 42;

    /**
     * < Switching the event file lsb.events
     */
    public static final int EVENT_LOG_SWITCH = 43;

    /**
     * < Job modification request
     */
    public static final int EVENT_JOB_MODIFY2 = 44;

    /**
     * < Log job group status
     */
    public static final int EVENT_JGRP_STATUS = 45;

    /**
     * < Job attributes have been set
     */
    public static final int EVENT_JOB_ATTR_SET = 46;

    /**
     * < Send an external message to a job
     */
    public static final int EVENT_JOB_EXT_MSG = 47;

    /**
     * < Update data status of a message for a job
     */
    public static final int EVENT_JOB_ATTA_DATA = 48;

    /**
     * < Insert one job to a chunk
     */
    public static final int EVENT_JOB_CHUNK = 49;

    /**
     * < Save unreported sbatchd status
     */
    public static final int EVENT_SBD_UNREPORTED_STATUS = 50;

    /**
     * < Reservation finished
     */
    public static final int EVENT_ADRSV_FINISH = 51;

    /**
     * < Dynamic host group control
     */
    public static final int EVENT_HGHOST_CTRL = 52;

    /**
     * < Saved current CPU allocation on service partition
     */
    public static final int EVENT_CPUPROFILE_STATUS = 53;

    /**
     * < Write out data logging file
     */
    public static final int EVENT_DATA_LOGGING = 54;

    /**
     * < Write job rusage in lsb.stream
     */
    public static final int EVENT_JOB_RUN_RUSAGE = 55;

    /**
     * < Stream closed and new stream opened.
     */
    public static final int EVENT_END_OF_STREAM = 56;

    /**
     * < SLA goal is reavaluated
     */
    public static final int EVENT_SLA_RECOMPUTE = 57;

    /**
     * < Write performance metrics to lsb.stream
     */
    public static final int EVENT_METRIC_LOG = 58;

    /**
     * < Write task finish log to ssched.acct
     */
    public static final int EVENT_TASK_FINISH = 59;

    /**
     * < Resize allocation is made
     */
    public static final int EVENT_JOB_RESIZE_NOTIFY_START = 60;

    /**
     * < Resize notification action initialized
     */
    public static final int EVENT_JOB_RESIZE_NOTIFY_ACCEPT = 61;

    /**
     * < Resize notification action completed
     */
    public static final int EVENT_JOB_RESIZE_NOTIFY_DONE = 62;

    /**
     * < Job resize release request is received
     */
    public static final int EVENT_JOB_RESIZE_RELEASE = 63;

    /**
     * < Job resize cancel request is received
     */
    public static final int EVENT_JOB_RESIZE_CANCEL = 64;

    /**
     * < Job resize event for lsb.acct
     */
    public static final int EVENT_JOB_RESIZE = 65;

    /**
     * < Saved array element's resource consumption  for LSF simulator
     */
    public static final int EVENT_JOB_ARRAY_ELEMENT = 66;

    /**
     * < Saved LSF simulator status
     */
    public static final int EVENT_MBD_SIM_STATUS = 67;

/* event kind
 */

    /**
     * < it is a job related event
     */
    public static final int EVENT_JOB_RELATED = 1;

    /**
     * < it is a non job related event
     */
    public static final int EVENT_NON_JOB_RELATED = 0;

    /*
   *  EXCLUSIVE PENDING REASONS
   *  a job must stay pending as long as ONE of the exclusive reasons exists
    */

/* Job Related Reasons (001 - 300)
 */
    /**
     * \addtogroup pending_reasons pending_reasons
     * \brief          Each entry in the table contains one of the following pending reasons
     */

    /**
     * < Virtual code; not a reason
     */
    public static final int PEND_JOB_REASON = 0;

    /**
     * < A new job is waiting to be scheduled
     */
    public static final int PEND_JOB_NEW = 1;

    /**
     * < The job is held until its specified start time
     */
    public static final int PEND_JOB_START_TIME = 2;

    /**
     * < The job is waiting for its dependency condition(s) to be satisfied
     */
    public static final int PEND_JOB_DEPEND = 3;

    /**
     * < The dependency condition is invalid or never satisfied
     */
    public static final int PEND_JOB_DEP_INVALID = 4;

    /**
     * < The migrating job is waiting to be rescheduled
     */
    public static final int PEND_JOB_MIG = 5;

    /**
     * < The job's pre-exec command exited with non-zero status
     */
    public static final int PEND_JOB_PRE_EXEC = 6;

    /**
     * < Unable to access jobfile
     */
    public static final int PEND_JOB_NO_FILE = 7;

    /**
     * < Unable to set job's environment variables
     */
    public static final int PEND_JOB_ENV = 8;

    /**
     * < Unable to determine the job's home or working directories
     */
    public static final int PEND_JOB_PATHS = 9;

    /**
     * < Unable to open the job's input and output files
     */
    public static final int PEND_JOB_OPEN_FILES = 10;

    /**
     * < Job execution initialization failed
     */
    public static final int PEND_JOB_EXEC_INIT = 11;

    /**
     * < Unable to copy restarting job's checkpoint files
     */
    public static final int PEND_JOB_RESTART_FILE = 12;

    /**
     * < Scheduling of the job is delayed
     */
    public static final int PEND_JOB_DELAY_SCHED = 13;

    /**
     * < Waiting for the re-scheduling of the job after switching queues
     */
    public static final int PEND_JOB_SWITCH = 14;

    /**
     * < An event is rejected by eeventd due to a syntax error
     */
    public static final int PEND_JOB_DEP_REJECT = 15;

    /**
     * < A JobScheduler feature is not enabled
     */
    public static final int PEND_JOB_JS_DISABLED = 16;

    /**
     * < Failed to get user password
     */
    public static final int PEND_JOB_NO_PASSWD = 17;

    /**
     * < The job is pending due to logon failure
     */
    public static final int PEND_JOB_LOGON_FAIL = 18;

    /**
     * < The job is waiting to be re-scheduled after its parameters have been changed
     */
    public static final int PEND_JOB_MODIFY = 19;

    /**
     * < The job time event is invalid
     */
    public static final int PEND_JOB_TIME_INVALID = 20;

    /**
     * < The job time event has expired
     */
    public static final int PEND_TIME_EXPIRED = 21;

    /**
     * < The job has been requeued
     */
    public static final int PEND_JOB_REQUEUED = 23;

    /**
     * < Waiting for the next time event
     */
    public static final int PEND_WAIT_NEXT = 24;

    /**
     * < The parent group is held
     */
    public static final int PEND_JGRP_HOLD = 25;

    /**
     * < The parent group is inactive
     */
    public static final int PEND_JGRP_INACT = 26;

    /**
     * < The group is waiting for scheduling
     */
    public static final int PEND_JGRP_WAIT = 27;

    /**
     * < The remote cluster(s) are unreachable
     */
    public static final int PEND_JOB_RCLUS_UNREACH = 28;

    /**
     * < SNDJOBS_TO queue rejected by remote  clusters
     */
    public static final int PEND_JOB_QUE_REJECT = 29;

    /**
     * < Waiting for new remote scheduling  session
     */
    public static final int PEND_JOB_RSCHED_START = 30;

    /**
     * < Waiting for allocation reply from remote clusters
     */
    public static final int PEND_JOB_RSCHED_ALLOC = 31;

    /**
     * < The job is forwarded to a remote cluster
     */
    public static final int PEND_JOB_FORWARDED = 32;

    /**
     * < The job running remotely is in a zombie state
     */
    public static final int PEND_JOB_RMT_ZOMBIE = 33;

    /**
     * < Job's enforced user group share account not selected
     */
    public static final int PEND_JOB_ENFUGRP = 34;

    /**
     * < The system is unable to schedule the job
     */
    public static final int PEND_SYS_UNABLE = 35;

    /**
     * < The parent group has just been released
     */
    public static final int PEND_JGRP_RELEASE = 36;

    /**
     * < The job has run since group active
     */
    public static final int PEND_HAS_RUN = 37;

    /**
     * < The job has reached its running element limit
     */
    public static final int PEND_JOB_ARRAY_JLIMIT = 38;

    /**
     * < Checkpoint directory is invalid
     */
    public static final int PEND_CHKPNT_DIR = 39;

    /**
     * < The first job in the chunk failed  (all other jobs in the chunk are set to PEND)
     */
    public static final int PEND_CHUNK_FAIL = 40;

    /**
     * < Optimum number of running jobs for SLA has been reached
     */
    public static final int PEND_JOB_SLA_MET = 41;

    /**
     * < Specified application profile does not exist
     */
    public static final int PEND_JOB_APP_NOEXIST = 42;

    /**
     * < Job no longer satisfies application  PROCLIMIT configuration
     */
    public static final int PEND_APP_PROCLIMIT = 43;

    /**
     * < No hosts for the job from EGO
     */
    public static final int PEND_EGO_NO_HOSTS = 44;

    /**
     * < The specified job group has reached its job limit
     */
    public static final int PEND_JGRP_JLIMIT = 45;

    /**
     * < Job pre-exec retry limit
     */
    public static final int PEND_PREEXEC_LIMIT = 46;

    /**
     * < Job re-queue limit
     */
    public static final int PEND_REQUEUE_LIMIT = 47;

    /**
     * < Job has bad res req
     */
    public static final int PEND_BAD_RESREQ = 48;

    /**
     * < Job's reservation is inactive
     */
    public static final int PEND_RSV_INACTIVE = 49;

    /**
     * < Job was in PSUSP with bad res req, after successful bmod  waiting for the user to bresume
     */
    public static final int PEND_WAITING_RESUME = 50;

    /**
     * < Job slot request cannot satisfy compound  resource requirement
     */
    public static final int PEND_SLOT_COMPOUND = 51;

/*
*  Queue and System Related Reasons (301 - 599)
 */

    /**
     * < The queue is inactivated by the administrator
     */
    public static final int PEND_QUE_INACT = 301;

    /**
     * < The queue is inactivated by its time windows
     */
    public static final int PEND_QUE_WINDOW = 302;

    /**
     * < The queue has reached its job slot limit
     */
    public static final int PEND_QUE_JOB_LIMIT = 303;

    /**
     * < The user has reached the per-user job slot limit of the queue
     */
    public static final int PEND_QUE_USR_JLIMIT = 304;

    /**
     * < Not enough per-user job slots of the queue for the parallel job
     */
    public static final int PEND_QUE_USR_PJLIMIT = 305;

    /**
     * < The queue's pre-exec command exited with non-zero status
     */
    public static final int PEND_QUE_PRE_FAIL = 306;

    /**
     * < The job was not accepted by the NQS host,  Attempt again later
     */
    public static final int PEND_NQS_RETRY = 307;

    /**
     * < Unable to send the job to an NQS host
     */
    public static final int PEND_NQS_REASONS = 308;

    /**
     * < Unable to contact NQS host
     */
    public static final int PEND_NQS_FUN_OFF = 309;

    /**
     * < The system is not ready for scheduling after reconfiguration
     */
    public static final int PEND_SYS_NOT_READY = 310;

    /**
     * < The requeued job is waiting for rescheduling
     */
    public static final int PEND_SBD_JOB_REQUEUE = 311;

    /**
     * < Not enough hosts to meet the job's spanning requirement
     */
    public static final int PEND_JOB_SPREAD_TASK = 312;

    /**
     * < Not enough hosts to meet the queue's spanning requirement
     */
    public static final int PEND_QUE_SPREAD_TASK = 313;

    /**
     * < The queue has not enough job slots for the parallel job
     */
    public static final int PEND_QUE_PJOB_LIMIT = 314;

    /**
     * < The job will not finish before queue's run window is closed
     */
    public static final int PEND_QUE_WINDOW_WILL_CLOSE = 315;

    /**
     * < Job no longer satisfies queue  PROCLIMIT configuration
     */
    public static final int PEND_QUE_PROCLIMIT = 316;

    /**
     * < Job requeued due to plug-in failure
     */
    public static final int PEND_SBD_PLUGIN = 317;

    /**
     * < Waiting for lease signing
     */
    public static final int PEND_WAIT_SIGN_LEASE = 318;

/* waitint for scheduling for SLOT_SHARE*/
    public static final int PEND_WAIT_SLOT_SHARE = 319;

/*
*  User Related Reasons (601 - 800)
 */

    /**
     * < The job slot limit is reached
     */
    public static final int PEND_USER_JOB_LIMIT = 601;

    /**
     * < A user group has reached its job slot limit
     */
    public static final int PEND_UGRP_JOB_LIMIT = 602;

    /**
     * < The job slot limit for the parallel job is reached
     */
    public static final int PEND_USER_PJOB_LIMIT = 603;

    /**
     * < A user group has reached its job slot limit for the parallel job
     */
    public static final int PEND_UGRP_PJOB_LIMIT = 604;

    /**
     * < Waiting for scheduling after resumed by user
     */
    public static final int PEND_USER_RESUME = 605;

    /**
     * < The job was suspended by the user while pending
     */
    public static final int PEND_USER_STOP = 607;

    /**
     * < Unable to determine user account for execution
     */
    public static final int PEND_NO_MAPPING = 608;

    /**
     * < The user has no permission to run the job on remote host/cluster
     */
    public static final int PEND_RMT_PERMISSION = 609;

    /**
     * < The job was suspended by LSF admin or root while pending
     */
    public static final int PEND_ADMIN_STOP = 610;

    /**
     * < The requested label is not valid
     */
    public static final int PEND_MLS_INVALID = 611;

    /**
     * < The requested label is above user allowed range
     */
    public static final int PEND_MLS_CLEARANCE = 612;

    /**
     * < The requested label rejected by /etc/rhost.conf
     */
    public static final int PEND_MLS_RHOST = 613;

    /**
     * < The requested label does not dominate current label
     */
    public static final int PEND_MLS_DOMINATE = 614;

    /**
     * < The requested label problem
     */
    public static final int PEND_MLS_FATAL = 615;

    /**
     * < LSF internally bstoped a pending job
     */
    public static final int PEND_INTERNAL_STOP = 616;

/*
*  NON-EXCLUSIVE PENDING REASONS
*  A job may still start even though non-exclusive reasons exist.
 */

/*
*  Host(sbatchd)-Job Related Reasons (1001 - 1300)
 */

    /**
     * < The job's resource requirements not satisfied
     */
    public static final int PEND_HOST_RES_REQ = 1001;

    /**
     * < The job's requirement for exclusive execution not satisfied
     */
    public static final int PEND_HOST_NONEXCLUSIVE = 1002;

    /**
     * < Higher or equal priority jobs already suspended by system
     */
    public static final int PEND_HOST_JOB_SSUSP = 1003;

    /**
     * < The job failed to compete with other jobs on host partition
     */
    public static final int PEND_HOST_PART_PRIO = 1004;

    /**
     * < Unable to get the PID of the restarting job
     */
    public static final int PEND_SBD_GETPID = 1005;

    /**
     * < Unable to lock the host for exclusively executing the job
     */
    public static final int PEND_SBD_LOCK = 1006;

    /**
     * < Cleaning up zombie job
     */
    public static final int PEND_SBD_ZOMBIE = 1007;

    /**
     * < Can't run jobs submitted by root.  The job is rejected by the sbatchd
     */
    public static final int PEND_SBD_ROOT = 1008;

    /**
     * < Job can't finish on the host before queue's run window is closed
     */
    public static final int PEND_HOST_WIN_WILL_CLOSE = 1009;

    /**
     * < Job can't finish on the host before job's termination deadline
     */
    public static final int PEND_HOST_MISS_DEADLINE = 1010;

    /**
     * < The specified first execution host is  not eligible for this job at this time
     */
    public static final int PEND_FIRST_HOST_INELIGIBLE = 1011;

    /**
     * < Exclusive job reserves slots on host
     */
    public static final int PEND_HOST_EXCLUSIVE_RESERVE = 1012;

    /**
     * < Resized shadow job  or non-first resReq of a compound resReq job try to reuse the first execution host
     */
    public static final int PEND_FIRST_HOST_REUSE = 1013;
/*
*  Host Related Reasons (1301 - 1600)
 */

    /**
     * < The host is closed by the LSF administrator
     */
    public static final int PEND_HOST_DISABLED = 1301;

    /**
     * < The host is locked by the LSF administrator
     */
    public static final int PEND_HOST_LOCKED = 1302;

    /**
     * < Not enough job slots for the parallel job
     */
    public static final int PEND_HOST_LESS_SLOTS = 1303;

    /**
     * < Dispatch windows are closed
     */
    public static final int PEND_HOST_WINDOW = 1304;

    /**
     * < The job slot limit reached
     */
    public static final int PEND_HOST_JOB_LIMIT = 1305;

    /**
     * < The queue's per-CPU job slot limit is reached
     */
    public static final int PEND_QUE_PROC_JLIMIT = 1306;

    /**
     * < The queue's per-host job slot limit is reached
     */
    public static final int PEND_QUE_HOST_JLIMIT = 1307;

    /**
     * < The user's per-CPU job slot limit is reached
     */
    public static final int PEND_USER_PROC_JLIMIT = 1308;

    /**
     * < The host's per-user job slot limit is reached
     */
    public static final int PEND_HOST_USR_JLIMIT = 1309;

    /**
     * < Not a member of the queue
     */
    public static final int PEND_HOST_QUE_MEMB = 1310;

    /**
     * < Not a user-specified host
     */
    public static final int PEND_HOST_USR_SPEC = 1311;

    /**
     * < The user has no access to the host partition
     */
    public static final int PEND_HOST_PART_USER = 1312;

    /**
     * < There is no such user account
     */
    public static final int PEND_HOST_NO_USER = 1313;

    /**
     * < Just started a job recently
     */
    public static final int PEND_HOST_ACCPT_ONE = 1314;

    /**
     * < Load info unavailable
     */
    public static final int PEND_LOAD_UNAVAIL = 1315;

    /**
     * < The LIM is unreachable by the sbatchd
     */
    public static final int PEND_HOST_NO_LIM = 1316;

    /**
     * < The host does not have a valid LSF software license
     */
    public static final int PEND_HOST_UNLICENSED = 1317;

    /**
     * < The queue's resource requirements are not satisfied
     */
    public static final int PEND_HOST_QUE_RESREQ = 1318;

    /**
     * < The submission host type is not the same
     */
    public static final int PEND_HOST_SCHED_TYPE = 1319;

    /**
     * < There are not enough processors to meet the job's spanning requirement.  The job level locality is unsatisfied.
     */
    public static final int PEND_JOB_NO_SPAN = 1320;

    /**
     * < There are not enough processors to meet the queue's spanning requirement.  The queue level locality is unsatisfied.
     */
    public static final int PEND_QUE_NO_SPAN = 1321;

    /**
     * < An exclusive job is running
     */
    public static final int PEND_HOST_EXCLUSIVE = 1322;

    /**
     * < Job Scheduler is disabled on the host.  It is not licensed to accept repetitive jobs.
     */
    public static final int PEND_HOST_JS_DISABLED = 1323;

    /**
     * < The user group's per-CPU job slot limit is reached
     */
    public static final int PEND_UGRP_PROC_JLIMIT = 1324;

    /**
     * < Incorrect host, group or cluster name
     */
    public static final int PEND_BAD_HOST = 1325;

    /**
     * < Host is not used by the queue
     */
    public static final int PEND_QUEUE_HOST = 1326;

    /**
     * < Host is locked by master LIM
     */
    public static final int PEND_HOST_LOCKED_MASTER = 1327;

    /**
     * < Not enough reserved job slots at this time for specified reservation ID
     */
    public static final int PEND_HOST_LESS_RSVSLOTS = 1328;

    /**
     * < Not enough slots or resources for whole duration of the job
     */
    public static final int PEND_HOST_LESS_DURATION = 1329;

    /**
     * < Specified reservation has expired or has been deleted
     */
    public static final int PEND_HOST_NO_RSVID = 1330;

    /**
     * < The host is closed due to lease is inactive
     */
    public static final int PEND_HOST_LEASE_INACTIVE = 1331;

    /**
     * < Not enough job slot(s) while advance reservation is active
     */
    public static final int PEND_HOST_ADRSV_ACTIVE = 1332;

    /**
     * < This queue is not configured to send jobs to the cluster specified in the advance
     */
    public static final int PEND_QUE_RSVID_NOMATCH = 1333;

    /**
     * < Individual host based reasons
     */
    public static final int PEND_HOST_GENERAL = 1334;

    /**
     * < Host does not belong to the specified  advance reservation
     */
    public static final int PEND_HOST_RSV = 1335;

    /**
     * < Host does not belong to a compute unit  of the required type
     */
    public static final int PEND_HOST_NOT_CU = 1336;

    /**
     * < A compute unit containing the host is  used exclusively
     */
    public static final int PEND_HOST_CU_EXCL = 1337;

    /**
     * < CU-level excl. job cannot start since CU  is occupied
     */
    public static final int PEND_HOST_CU_OCCUPIED = 1338;

    /**
     * < Insufficiently many usable slots on the  host's compute unit
     */
    public static final int PEND_HOST_USABLE_CU = 1339;

    /**
     * < No first execution compute unit satisfies CU 'usablepercu' requirement.
     */
    public static final int PEND_JOB_FIRST_CU = 1340;

    /**
     * < A CU containing the host is reserved  exclusively
     */
    public static final int PEND_HOST_CU_EXCL_RSV = 1341;

    /**
     * < Maxcus cannot be satisfied
     */
    public static final int PEND_JOB_CU_MAXCUS = 1342;

    /**
     * < Balance cannot be satisfied
     */
    public static final int PEND_JOB_CU_BALANCE = 1343;

    /**
     * < Cu not supported on toplib integration hosts
     */
    public static final int PEND_CU_TOPLIB_HOST = 1344;

/*
*  sbatchd Related Reasons (1601 - 1900)
 */

    /**
     * < Cannot reach sbatchd
     */
    public static final int PEND_SBD_UNREACH = 1601;

    /**
     * < Number of jobs exceed quota
     */
    public static final int PEND_SBD_JOB_QUOTA = 1602;

    /**
     * < The job failed in talking to the server to start the job
     */
    public static final int PEND_JOB_START_FAIL = 1603;

    /**
     * < Failed in receiving the reply from the server when starting the job
     */
    public static final int PEND_JOB_START_UNKNWN = 1604;

    /**
     * < Unable to allocate memory to run job.  There is no memory on the sbatchd.
     */
    public static final int PEND_SBD_NO_MEM = 1605;

    /**
     * < Unable to fork process to run the job.  There are no more processes on the sbatchd.
     */
    public static final int PEND_SBD_NO_PROCESS = 1606;

    /**
     * < Unable to communicate with the job process
     */
    public static final int PEND_SBD_SOCKETPAIR = 1607;

    /**
     * < The slave batch server failed to accept the job
     */
    public static final int PEND_SBD_JOB_ACCEPT = 1608;

    /**
     * < Lease job remote dispatch failed
     */
    public static final int PEND_LEASE_JOB_REMOTE_DISPATCH = 1609;

    /**
     * < Failed to restart job from last checkpoint
     */
    public static final int PEND_JOB_RESTART_FAIL = 1610;
/*
*  Load Related Reasons (2001 - 2300)
 */

    /**
     * < The load threshold is reached
     */
    public static final int PEND_HOST_LOAD = 2001;

/*
*  Queue Resource Reservation Related Reasons (2301 - 2600)
 */

    /**
     * < The queue's requirements for resource  reservation are not satisfied.
     */
    public static final int PEND_HOST_QUE_RUSAGE = 2301;

/*
*  Jobs Resource Reservation Related Reasons (2601 - 2900)
 */

    /**
     * < The job's requirements for resource  reservation are not satisfied.
     */
    public static final int PEND_HOST_JOB_RUSAGE = 2601;

/*
*  Remote Forwarding Related Reasons (2901 - 3200)
 */

    /**
     * < Remote job not recongized by remote cluster, waiting for rescheduling
     */
    public static final int PEND_RMT_JOB_FORGOTTEN = 2901;

    /**
     * < Remote import limit reached, waiting  for rescheduling
     */
    public static final int PEND_RMT_IMPT_JOBBKLG = 2902;

    /**
     * < Remote schedule time reached,  waiting for rescheduling
     */
    public static final int PEND_RMT_MAX_RSCHED_TIME = 2903;

    /**
     * < Remote pre-exec retry limit reached, waiting for rescheduling
     */
    public static final int PEND_RMT_MAX_PREEXEC_RETRY = 2904;

    /**
     * < Remote queue is closed
     */
    public static final int PEND_RMT_QUEUE_CLOSED = 2905;

    /**
     * < Remote queue is inactive
     */
    public static final int PEND_RMT_QUEUE_INACTIVE = 2906;

    /**
     * < Remote queue is congested
     */
    public static final int PEND_RMT_QUEUE_CONGESTED = 2907;

    /**
     * < Remote queue is disconnected
     */
    public static final int PEND_RMT_QUEUE_DISCONNECT = 2908;

    /**
     * < Remote queue is not configured to accept jobs from this cluster
     */
    public static final int PEND_RMT_QUEUE_NOPERMISSION = 2909;

    /**
     * < Job's termination time exceeds the job creation time on remote cluster
     */
    public static final int PEND_RMT_BAD_TIME = 2910;

    /**
     * < Permission denied on the execution cluster
     */
    public static final int PEND_RMT_PERMISSIONS = 2911;

    /**
     * < Job's required on number of processors cannot be satisfied on the remote cluster
     */
    public static final int PEND_RMT_PROC_NUM = 2912;

    /**
     * < User is not defined in the fairshare policy of the remote queue
     */
    public static final int PEND_RMT_QUEUE_USE = 2913;

    /**
     * < Remote queue is a non-interactive queue
     */
    public static final int PEND_RMT_NO_INTERACTIVE = 2914;

    /**
     * < Remote queue is an interactive-only queue
     */
    public static final int PEND_RMT_ONLY_INTERACTIVE = 2915;

    /**
     * < Job's required maximum number of  processors is less then the minimum number
     */
    public static final int PEND_RMT_PROC_LESS = 2916;

    /**
     * < Job's required resource limit exceeds that of the remote queue
     */
    public static final int PEND_RMT_OVER_LIMIT = 2917;

    /**
     * < Job's resource requirements do not match with those of the remote queue
     */
    public static final int PEND_RMT_BAD_RESREQ = 2918;

    /**
     * < Job failed to be created on the remote cluster
     */
    public static final int PEND_RMT_CREATE_JOB = 2919;

    /**
     * < Job is requeued for rerun on the execution cluster
     */
    public static final int PEND_RMT_RERUN = 2920;

    /**
     * < Job is requeued on the execution cluster due to exit value
     */
    public static final int PEND_RMT_EXIT_REQUEUE = 2921;

    /**
     * < Job was killed and requeued on the execution cluster
     */
    public static final int PEND_RMT_REQUEUE = 2922;

    /**
     * < Job was forwarded to remote cluster
     */
    public static final int PEND_RMT_JOB_FORWARDING = 2923;

    /**
     * < Remote import queue defined for the job in lsb.queues is either not ready or not valid
     */
    public static final int PEND_RMT_QUEUE_INVALID = 2924;

    /**
     * < Remote queue is a non-exclusive queue
     */
    public static final int PEND_RMT_QUEUE_NO_EXCLUSIVE = 2925;

    /**
     * < Job was rejected; submitter does not belong to the specified User Group in the remote cluster or the user group does not exist in the remote cluster
     */
    public static final int PEND_RMT_UGROUP_MEMBER = 2926;

    /**
     * < Remote queue is rerunnable: can not accept interactive jobs
     */
    public static final int PEND_RMT_INTERACTIVE_RERUN = 2927;

    /**
     * < Remote cluster failed in talking to server to start the job
     */
    public static final int PEND_RMT_JOB_START_FAIL = 2928;

    /**
     * < Job was rejected; submitter does not belong to the specified User Group in the remote cluster or the user group does not exist in the remote cluster
     */
    public static final int PEND_RMT_FORWARD_FAIL_UGROUP_MEMBER = 2930;

    /**
     * < Specified remote reservation has expired or has been deleted
     */
    public static final int PEND_RMT_HOST_NO_RSVID = 2931;

    /**
     * < Application profile could not be found in the remote cluster.
     */
    public static final int PEND_RMT_APP_NULL = 2932;

    /**
     * < Job's required RUNLIMIT exceeds  RUNTIME*  JOB_RUNLIMIT_RATIO of the remote cluster.
     */
    public static final int PEND_RMT_BAD_RUNLIMIT = 2933;

    /**
     * < Job's required RUNTIME exceeds the hard runtime limit in the remote queue.
     */
    public static final int PEND_RMT_OVER_QUEUE_LIMIT = 2934;

    /**
     * < Job will be pend when no slots available among remote queues.
     */
    public static final int PEND_RMT_WHEN_NO_SLOTS = 2935;
/* SUSPENDING REASONS */

/*
*  General Resource Limits Related Reasons ( 3201 - 4800)
 */

    /**
     * < Resource limit defined on user  or user group has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_USER = 3201;

    /**
     * < Resource (%s) limit defined on queue has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_QUEUE = 3501;

    /**
     * < Resource limit defined on project has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_PROJECT = 3801;

    /**
     * < Resource (%s) limit defined cluster wide has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_CLUSTER = 4101;

    /**
     * < Resource (%s) limit defined on host and/or host group has  been reached.
     */
    public static final int PEND_GENERAL_LIMIT_HOST = 4401;

    /**
     * < JOBS limit defined for the user or user group has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_JOBS_USER = 4701;

    /**
     * < JOBS limit defined for the queue has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_JOBS_QUEUE = 4702;

    /**
     * < JOBS limit defined for the project has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_JOBS_PROJECT = 4703;

    /**
     * < JOBS limit defined cluster-wide has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_JOBS_CLUSTER = 4704;

    /**
     * < JOBS limit defined on host or host group has been reached.
     */
    public static final int PEND_GENERAL_LIMIT_JOBS_HOST = 4705;

/* LSF2 Presto RLA-related reasons    (4900 - 4989) */

    /**
     * < RMS scheduler plugin  internal error.
     */
    public static final int PEND_RMS_PLUGIN_INTERNAL = 4900;

    /**
     * < RLA communication failure.
     */
    public static final int PEND_RMS_PLUGIN_RLA_COMM = 4901;

    /**
     * < RMS is not available.
     */
    public static final int PEND_RMS_NOT_AVAILABLE = 4902;

    /**
     * < Cannot satisfy the topology  requirement.
     */
    public static final int PEND_RMS_FAIL_TOPOLOGY = 4903;

    /**
     * < Cannot allocate an RMS resource.
     */
    public static final int PEND_RMS_FAIL_ALLOC = 4904;

    /**
     * < RMS job with special topology requirements cannot be preemptive or backfill job.
     */
    public static final int PEND_RMS_SPECIAL_NO_PREEMPT_BACKFILL = 4905;

    /**
     * < RMS job with special topology requirements cannot reserve slots.
     */
    public static final int PEND_RMS_SPECIAL_NO_RESERVE = 4906;

    /**
     * < RLA internal error.
     */
    public static final int PEND_RMS_RLA_INTERNAL = 4907;

    /**
     * < Not enough slots for job.  Job with RMS topology requirements cannot reserve slots, be preemptive, or be a backfill job.
     */
    public static final int PEND_RMS_NO_SLOTS_SPECIAL = 4908;

    /**
     * < User account does not exist on the execution host.
     */
    public static final int PEND_RMS_RLA_NO_SUCH_USER = 4909;

    /**
     * < Unknown host and/or partition unavailable.
     */
    public static final int PEND_RMS_RLA_NO_SUCH_HOST = 4910;

    /**
     * < Cannot schedule chunk jobs to RMS hosts.
     */
    public static final int PEND_RMS_CHUNKJOB = 4911;

    /**
     * < RLA protocol mismatch.
     */
    public static final int PEND_RLA_PROTOMISMATCH = 4912;

    /**
     * < Contradictory topology requirements specified.
     */
    public static final int PEND_RMS_BAD_TOPOLOGY = 4913;

    /**
     * < Not enough slots to satisfy manditory contiguous requirement.
     */
    public static final int PEND_RMS_RESREQ_MCONT = 4914;

    /**
     * < Not enough slots to satisfy RMS ptile requirement.
     */
    public static final int PEND_RMS_RESREQ_PTILE = 4915;

    /**
     * < Not enough slots to satisfy RMS nodes requirement.
     */
    public static final int PEND_RMS_RESREQ_NODES = 4916;

    /**
     * < Cannot satisfy RMS base node requirement.
     */
    public static final int PEND_RMS_RESREQ_BASE = 4917;

    /**
     * < Cannot satisfy RMS rails requirement.
     */
    public static final int PEND_RMS_RESREQ_RAILS = 4918;

    /**
     * < Cannot satisfy RMS railmask requirement.
     */
    public static final int PEND_RMS_RESREQ_RAILMASK = 4919;


/*
*  Maui Integration Related Reasons ( 5000 - 5100)
 */

    /**
     * < Unable to communicate with external Maui scheduler.
     */
    public static final int PEND_MAUI_UNREACH = 5000;

    /**
     * < Job is pending at external Maui scheduler.
     */
    public static final int PEND_MAUI_FORWARD = 5001;

    /**
     * < External Maui scheduler sets detail reason.
     */
    public static final int PEND_MAUI_REASON = 5030;

/*
*  SGI CPUSET Integration Related Reasons ( 5200 - 5299)
 */

    /**
     * < CPUSET attach failed.  Job requeued
     */
    public static final int PEND_CPUSET_ATTACH = 5200;

    /**
     * < Not a cpuset host
     */
    public static final int PEND_CPUSET_NOT_CPUSETHOST = 5201;

    /**
     * < Topd initialization failed
     */
    public static final int PEND_CPUSET_TOPD_INIT = 5202;

    /**
     * < Topd communication timeout
     */
    public static final int PEND_CPUSET_TOPD_TIME_OUT = 5203;

    /**
     * < Cannot satisfy the cpuset  allocation requirement
     */
    public static final int PEND_CPUSET_TOPD_FAIL_ALLOC = 5204;

    /**
     * < Bad cpuset allocation request
     */
    public static final int PEND_CPUSET_TOPD_BAD_REQUEST = 5205;

    /**
     * < Topd internal error
     */
    public static final int PEND_CPUSET_TOPD_INTERNAL = 5206;

    /**
     * < Cpuset system API failure
     */
    public static final int PEND_CPUSET_TOPD_SYSAPI_ERR = 5207;

    /**
     * < Specified static cpuset does  not exist on the host
     */
    public static final int PEND_CPUSET_TOPD_NOSUCH_NAME = 5208;

    /**
     * < Cpuset is already allocated   for this job
     */
    public static final int PEND_CPUSET_TOPD_JOB_EXIST = 5209;

    /**
     * < Topd malloc failure
     */
    public static final int PEND_CPUSET_TOPD_NO_MEMORY = 5210;

    /**
     * < User account does not exist   on the cpuset host
     */
    public static final int PEND_CPUSET_TOPD_INVALID_USER = 5211;

    /**
     * < User does not have permission   to run job within cpuset
     */
    public static final int PEND_CPUSET_TOPD_PERM_DENY = 5212;

    /**
     * < Topd is not available
     */
    public static final int PEND_CPUSET_TOPD_UNREACH = 5213;

    /**
     * < Topd communication failure
     */
    public static final int PEND_CPUSET_TOPD_COMM_ERR = 5214;


    /**
     * < CPUSET scheduler plugin internal error
     */
    public static final int PEND_CPUSET_PLUGIN_INTERNAL = 5215;

    /**
     * < Cannot schedule chunk jobs to cpuset hosts
     */
    public static final int PEND_CPUSET_CHUNKJOB = 5216;

    /**
     * < Can't satisfy CPU_LIST   requirement
     */
    public static final int PEND_CPUSET_CPULIST = 5217;

    /**
     * < Cannot satisfy CPUSET MAX_RADIUS requirement
     */
    public static final int PEND_CPUSET_MAXRADIUS = 5218;

/* Bproc integration related reasons (5300 - 5320)
 */

    /**
     * < Node allocation failed
     */
    public static final int PEND_NODE_ALLOC_FAIL = 5300;

/* Eagle pending reasons  (5400 - 5449) */

    /**
     * < RMS resource is not available
     */
    public static final int PEND_RMSRID_UNAVAIL = 5400;


    /**
     * < Not enough free cpus to satisfy job requirements
     */
    public static final int PEND_NO_FREE_CPUS = 5450;

    /**
     * < Topology unknown or recently changed
     */
    public static final int PEND_TOPOLOGY_UNKNOWN = 5451;

    /**
     * < Contradictory topology requirement specified
     */
    public static final int PEND_BAD_TOPOLOGY = 5452;

    /**
     * < RLA communications failure
     */
    public static final int PEND_RLA_COMM = 5453;

    /**
     * < User account does not exist on execution host
     */
    public static final int PEND_RLA_NO_SUCH_USER = 5454;

    /**
     * < RLA internal error
     */
    public static final int PEND_RLA_INTERNAL = 5455;

    /**
     * < Unknown host and/or partition unavailable
     */
    public static final int PEND_RLA_NO_SUCH_HOST = 5456;

    /**
     * < Too few slots for specified topology requirement
     */
    public static final int PEND_RESREQ_TOOFEWSLOTS = 5457;

/* PSET pending reasons (5500 - 5549) */

    /**
     * < PSET scheduler plugin internal error
     */
    public static final int PEND_PSET_PLUGIN_INTERNAL = 5500;

    /**
     * < Cannot satisfy PSET ptile requirement
     */
    public static final int PEND_PSET_RESREQ_PTILE = 5501;

    /**
     * < Cannot satisfy PSET cells requirement
     */
    public static final int PEND_PSET_RESREQ_CELLS = 5502;

    /**
     * < Cannot schedule chunk jobs to PSET hosts
     */
    public static final int PEND_PSET_CHUNKJOB = 5503;

    /**
     * < Host does not support processor set functionality
     */
    public static final int PEND_PSET_NOTSUPPORT = 5504;

    /**
     * < PSET bind failed. Job requeued
     */
    public static final int PEND_PSET_BIND_FAIL = 5505;

    /**
     * < Cannot satisfy PSET CELL_LIST  requirement
     */
    public static final int PEND_PSET_RESREQ_CELLLIST = 5506;


/* SLURM pending reasons (5550 - 5599) */

    /**
     * < SLURM scheduler plugin internal error
     */
    public static final int PEND_SLURM_PLUGIN_INTERNAL = 5550;

    /**
     * < Not enough resource to satisfy  SLURM nodes requirment
     */
    public static final int PEND_SLURM_RESREQ_NODES = 5551;

    /**
     * < Not enough resource to satisfy  SLURM node attributes requirment.
     */
    public static final int PEND_SLURM_RESREQ_NODE_ATTR = 5552;

    /**
     * < Not enough resource to satisfy SLURM exclude requirment.
     */
    public static final int PEND_SLURM_RESREQ_EXCLUDE = 5553;

    /**
     * < Not enough resource to satisfy SLURM nodelist requirment.
     */
    public static final int PEND_SLURM_RESREQ_NODELIST = 5554;

    /**
     * < Not enough resource to satisfy SLURM contiguous requirment.
     */
    public static final int PEND_SLURM_RESREQ_CONTIGUOUS = 5555;

    /**
     * < SLURM allocation is not available. Job requeued.
     */
    public static final int PEND_SLURM_ALLOC_UNAVAIL = 5556;

    /**
     * < Invalid grammar in SLURM constraints option, job will never run.
     */
    public static final int PEND_SLURM_RESREQ_BAD_CONSTRAINT = 5557;

/* Cray X1 pending reasons (5600 - 5649) */

    /**
     * < Not enough SSPs for job
     */
    public static final int PEND_CRAYX1_SSP = 5600;

    /**
     * < Not enough MSPs for job
     */
    public static final int PEND_CRAYX1_MSP = 5601;

    /**
     * < Unable to pass limit information to psched.
     */
    public static final int PEND_CRAYX1_PASS_LIMIT = 5602;

/* Cray XT3 pending reasons (5650 - 5699) */

    /**
     * < Unable to create or assign a  partition by CPA
     */
    public static final int PEND_CRAYXT3_ASSIGN_FAIL = 5650;

/* BlueGene pending reasons (5700 - 5749) */

    /**
     * < BG/L: Scheduler plug-in internal error.
     */
    public static final int PEND_BLUEGENE_PLUGIN_INTERNAL = 5700;

    /**
     * < BG/L: Allocation is not available. Job requeued.
     */
    public static final int PEND_BLUEGENE_ALLOC_UNAVAIL = 5701;

    /**
     * < BG/L: No free base partitions available for a full block allocation.
     */
    public static final int PEND_BLUEGENE_NOFREEMIDPLANES = 5702;

    /**
     * < BG/L: No free quarters available for a small block allocation.
     */
    public static final int PEND_BLUEGENE_NOFREEQUARTERS = 5703;

    /**
     * < BG/L: No free node cards available for a small block allocation.
     */
    public static final int PEND_BLUEGENE_NOFREENODECARDS = 5704;

/* resize enhancement releated pending reasons */

    /**
     * < First execution host unavailable
     */
    public static final int PEND_RESIZE_FIRSTHOSTUNAVAIL = 5705;

    /**
     * < Master is not in the RUN state
     */
    public static final int PEND_RESIZE_MASTERSUSP = 5706;

    /**
     * < Host is not same as for master
     */
    public static final int PEND_RESIZE_MASTER_SAME = 5707;

    /**
     * < Host already used by master
     */
    public static final int PEND_RESIZE_SPAN_PTILE = 5708;

    /**
     * < The job can only use first host
     */
    public static final int PEND_RESIZE_SPAN_HOSTS = 5709;

    /**
     * < The job cannot get slots on remote hosts
     */
    public static final int PEND_RESIZE_LEASE_HOST = 5710;

/* compound resreq related pending reasons (5800 - ??) */

    /**
     * < The job cannot get slots on  pre-7Update5 remote hosts
     */
    public static final int PEND_COMPOUND_RESREQ_OLD_LEASE_HOST = 5800;

    /**
     * < Hosts using LSF HPC system  integrations do not support compound resource requirements.
     */
    public static final int PEND_COMPOUND_RESREQ_TOPLIB_HOST = 5801;
/* multi-phase resreq related pending reasons (5900 - ??) */

    /**
     * < The job cannot get slots on  pre-7Update6 remote hosts
     */
    public static final int PEND_MULTIPHASE_RESREQ_OLD_LEASE_HOST = 5900;

/* EGO-Enabled SLA pending reasons (5750 - 5799) */

    /**
     * < Host does not have enough slots for this SLA job.
     */
    public static final int PEND_PS_PLUGIN_INTERNAL = 5750;

    /**
     * < EGO SLA: Failed to synchronize resource with MBD.
     */
    public static final int PEND_PS_MBD_SYNC = 5751;


/* PLATFORM reserves pending reason number from 1 - 20000.
*  External plugin is suggested to use platform's reserved pending reason
*  number. However, they can use pending reason number between 20001 - 25000
*  as customer specific pending reasons. In this case, bjobs -p will only show
*  the reason number without detailed message
 */

    /**
     * < Customized pending reason number between min and max.
     */
    public static final int PEND_CUSTOMER_MIN = 20001;

    /**
     * < Customized pending reason number between min and max.
     */
    public static final int PEND_CUSTOMER_MAX = 25000;


    /**
     * < The maximum number of reasons
     */
    public static final int PEND_MAX_REASONS = 25001;

    /**
     * \addtogroup suspending_reasons  suspending_reasons
     * suspending_reasons is part of pending_reasons
     */
/* SUSPENDING REASONS */

/* User related reasons */

    /**
     * < Virtual code. Not a reason
     */
    public static final int SUSP_USER_REASON = 0x00000000;

    /**
     * < The job is waiting to be re-scheduled after being resumed by the user.
     */
    public static final int SUSP_USER_RESUME = 0x00000001;

    /**
     * < The user suspended the job.
     */
    public static final int SUSP_USER_STOP = 0x00000002;

/* Queue and system related reasons */

    /**
     * < Virtual code. Not a reason
     */
    public static final int SUSP_QUEUE_REASON = 0x00000004;

    /**
     * < The run window of the queue is closed.
     */
    public static final int SUSP_QUEUE_WINDOW = 0x00000008;

    /**
     * < Suspended after preemption. The system needs to re-allocate CPU utilization by job priority.
     */
    public static final int SUSP_RESCHED_PREEMPT = 0x00000010;

    /**
     * < The LSF administrator has locked the execution host.
     */
    public static final int SUSP_HOST_LOCK = 0x00000020;

    /**
     * < A load index exceeds its threshold. The subreasons field indicates which indices.
     */
    public static final int SUSP_LOAD_REASON = 0x00000040;

    /**
     * < The job was preempted by mbatchd because of a higher priorty job.
     */
    public static final int SUSP_MBD_PREEMPT = 0x00000080;

    /**
     * < Preempted by sbatchd. The job limit of the host/user has been reached.
     */
    public static final int SUSP_SBD_PREEMPT = 0x00000100;

    /**
     * < The suspend conditions of the queue,  as specified by the STOP_COND parameter in lsb.queues, are true.
     */
    public static final int SUSP_QUE_STOP_COND = 0x00000200;

    /**
     * < The resume conditions of the queue, as specified by the RESUME_COND parameter in lsb.queues, are false.
     */
    public static final int SUSP_QUE_RESUME_COND = 0x00000400;

    /**
     * < The job was suspended due to the paging rate and the host is not idle yet.
     */
    public static final int SUSP_PG_IT = 0x00000800;

    /**
     * < Resets the previous reason.
     */
    public static final int SUSP_REASON_RESET = 0x00001000;

    /**
     * < Load information on the execution hosts is unavailable.
     */
    public static final int SUSP_LOAD_UNAVAIL = 0x00002000;

    /**
     * < The job was suspened by root or the LSF administrator.
     */
    public static final int SUSP_ADMIN_STOP = 0x00004000;

    /**
     * < The job is terminated due to resource limit.
     */
    public static final int SUSP_RES_RESERVE = 0x00008000;

    /**
     * < The job is locked by the mbatchd.
     */
    public static final int SUSP_MBD_LOCK = 0x00010000;

    /**
     * < The job's requirements for resource  reservation are not satisfied.
     */
    public static final int SUSP_RES_LIMIT = 0x00020000;

    /**
     * < The job is suspended while the sbatchd is restarting.
     */
    public static final int SUSP_SBD_STARTUP = 0x00040000;

    /**
     * < The execution host is locked by the master LIM.
     */
    public static final int SUSP_HOST_LOCK_MASTER = 0x00080000;

    /**
     * < An advance reservation using the  host is active
     */
    public static final int SUSP_HOST_RSVACTIVE = 0x00100000;

    /**
     * < There is a detailed reason in the subreason field
     */
    public static final int SUSP_DETAILED_SUBREASON = 0x00200000;
    /* GLB suspending reason */

    /**
     * < The job is preempted by glb
     */
    public static final int SUSP_GLB_LICENSE_PREEMPT = 0x00400000;

    /* Cray X1 suspend reasons */

    /**
     * < Job not placed by Cray X1  psched
     */
    public static final int SUSP_CRAYX1_POSTED = 0x00800000;

    /**
     * < Job suspended when its advance  reservation expired
     */
    public static final int SUSP_ADVRSV_EXPIRED = 0x01000000;

    /**
     * \addtogroup suspending_subreasons  suspending_subreasons
     * suspending_subreasons has the following options:
     */

    /**
     * < Sub reason of SUSP_RES_LIMIT: RUNLIMIT is reached.
     */
    public static final int SUB_REASON_RUNLIMIT = 0x00000001;

    /**
     * < Sub reason of SUSP_RES_LIMIT: DEADLINE is reached.
     */
    public static final int SUB_REASON_DEADLINE = 0x00000002;

    /**
     * < Sub reason of SUSP_RES_LIMIT: PROCESSLIMIT is reached.
     */
    public static final int SUB_REASON_PROCESSLIMIT = 0x00000004;

    /**
     * < Sub reason of SUSP_RES_LIMIT: CPULIMIT is reached.
     */
    public static final int SUB_REASON_CPULIMIT = 0x00000008;

    /**
     * < Sub reason of SUSP_RES_LIMIT: MEMLIMIT is reached.
     */
    public static final int SUB_REASON_MEMLIMIT = 0x00000010;

    /**
     * < Sub reason of SUSP_RES_LIMIT: THREADLIMIT is reached.
     */
    public static final int SUB_REASON_THREADLIMIT = 0x00000020;

    /**
     * < Sub reason of SUSP_RES_LIMIT: SWAPLIMIT is reached.
     */
    public static final int SUB_REASON_SWAPLIMIT = 0x00000040;

    /**
     * < Account ID does not match those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_ACCOUNTID = 0x00000001;

    /**
     * < Attribute does not match  those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_ATTRIBUTE = 0x00000002;

    /**
     * < Blocked by one or more gates
     */
    public static final int SUB_REASON_CRAYX1_BLOCKED = 0x00000004;

    /**
     * < Application is in the process of being restarted  and it is under the control  of CPR
     */
    public static final int SUB_REASON_CRAYX1_RESTART = 0x00000008;

    /**
     * < Depth does not match those  allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_DEPTH = 0x00000010;

    /**
     * < GID does not match those  allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_GID = 0x00000020;

    /**
     * < No GASID is available
     */
    public static final int SUB_REASON_CRAYX1_GASID = 0x00000040;

    /**
     * < Hard label does not match  those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_HARDLABEL = 0x00000080;

    /**
     * < Limit exceeded in regions   or domains
     */
    public static final int SUB_REASON_CRAYX1_LIMIT = 0x00000100;

    /**
     * < Memory size does not match  those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_MEMORY = 0x00000200;

    /**
     * < Soft label does not match   those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_SOFTLABEL = 0x00000400;

    /**
     * < Size gate (width times  depth larger than gate  allows)
     */
    public static final int SUB_REASON_CRAYX1_SIZE = 0x00000800;

    /**
     * < Time limit does not match those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_TIME = 0x00001000;

    /**
     * < UID does not match those  allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_UID = 0x00002000;

    /**
     * < Width does not match those allowed by the gate
     */
    public static final int SUB_REASON_CRAYX1_WIDTH = 0x00004000;
/*
*  EXITING REASONS: currently only to indicate exited due to
*  1) rerunnable job being restart from last chkpnt;
*  2) being killed while execution host is unavailable
 */

    /** Job finished normally */
    public static final int EXIT_NORMAL = 0x00000000;

    /** Rerunnable job to be restarted */
    public static final int EXIT_RESTART = 0x00000001;

    /** Job killed while host unavailable */
    public static final int EXIT_ZOMBIE = 0x00000002;

    /** Job is finished and put into pend list */
    public static final int FINISH_PEND = 0x00000004;

    /** The job is killed while the execution host is unreach */
    public static final int EXIT_KILL_ZOMBIE = 0x00000008;

    /** The job in ZOMBIE is removed */
    public static final int EXIT_ZOMBIE_JOB = 0x00000010;

    /** Rerun a job without creating a ZOMBIE job */
    public static final int EXIT_RERUN = 0x00000020;

    /** Remote job has no mapping user name here */
    public static final int EXIT_NO_MAPPING = 0x00000040;

    /** Remote job has no permission running here */
    public static final int EXIT_REMOTE_PERMISSION = 0x00000080;

    /** Remote job cannot run locally because of environment problem */
    public static final int EXIT_INIT_ENVIRON = 0x00000100;

    /** Remote job failed in pre_exec command */
    public static final int EXIT_PRE_EXEC = 0x00000200;

    /** The job is killed and will be later requeued */
    public static final int EXIT_REQUEUE = 0x00000400;

    /** Job could not be killed but was removed from system */
    public static final int EXIT_REMOVE = 0x00000800;

    /** Requeue by exit value */
    public static final int EXIT_VALUE_REQUEUE = 0x00001000;

    /** Cancel request received from remote cluster. */
    public static final int EXIT_CANCEL = 0x00002000;

    /** MED killed job on web server */
    public static final int EXIT_MED_KILLED = 0x00004000;

    /** Remote lease job exit on execution, side, return to pend on submission */
    public static final int EXIT_REMOTE_LEASE_JOB = 0x00008000;

    /** Exit when cwd does not exist*/
    public static final int EXIT_CWD_NOTEXIST = 0x00010000;


    /** Mode indicating running in batch, js, or batch-js mode */
    public static final int LSB_MODE_BATCH = 0x1;
    public static final int LSB_MODE_JS = 0x2;
    public static final int LSB_MODE_BATCH_RD = 0x4;

    public static final int RLIMIT_CPU = 0;
    public static final int RLIMIT_FSIZE = 1;
    public static final int RLIMIT_DATA = 2;
    public static final int RLIMIT_STACK = 3;
    public static final int RLIMIT_CORE = 4;
    public static final int RLIMIT_RSS = 5;
    public static final int RLIM_INFINITY = 0x7fffffff;

/*
*  Error codes for lsblib calls
*  Each error code has its corresponding error message defined in lsb.err.c
*  The code number is just the position number of its message.
*  Adding a new code here must add its message there in the corresponding
*  position.  Changing any code number here must change the position there.
 */
/* Error codes related to job */

    /** No error at all */
    public static final int LSBE_NO_ERROR = 0;

    /** No matching job found */
    public static final int LSBE_NO_JOB = 1;

    /** Job not started yet */
    public static final int LSBE_NOT_STARTED = 2;

    /** Job already started */
    public static final int LSBE_JOB_STARTED = 3;

    /** Job already finished */
    public static final int LSBE_JOB_FINISH = 4;

    /** Ask sbatchd to stop the wrong job */
    public static final int LSBE_STOP_JOB = 5;

    /** Depend_cond syntax error */
    public static final int LSBE_DEPEND_SYNTAX = 6;

    /** Queue doesn't accept EXCLUSIVE job */
    public static final int LSBE_EXCLUSIVE = 7;

    /** Root is not allowed to submit jobs */
    public static final int LSBE_ROOT = 8;

    /** Job is already being migrated */
    public static final int LSBE_MIGRATION = 9;

    /** Job is not chkpntable */
    public static final int LSBE_J_UNCHKPNTABLE = 10;

    /** Job has no output so far */
    public static final int LSBE_NO_OUTPUT = 11;

    /** No jobId can be used now */
    public static final int LSBE_NO_JOBID = 12;

    /** Queue only accepts bsub -I job */
    public static final int LSBE_ONLY_INTERACTIVE = 13;

    /** Queue doesn't accept bsub -I job */
    public static final int LSBE_NO_INTERACTIVE = 14;

/** Error codes related to user, queue and host */

    /** No user defined in lsb.users file */
    public static final int LSBE_NO_USER = 15;

    /** Bad user name */
    public static final int LSBE_BAD_USER = 16;

    /** User permission denied */
    public static final int LSBE_PERMISSION = 17;

    /** No such queue in the system */
    public static final int LSBE_BAD_QUEUE = 18;

    /** Queue name should be given */
    public static final int LSBE_QUEUE_NAME = 19;

    /** Queue has been closed */
    public static final int LSBE_QUEUE_CLOSED = 20;

    /** Queue windows are closed */
    public static final int LSBE_QUEUE_WINDOW = 21;

    /** User cannot use the queue */
    public static final int LSBE_QUEUE_USE = 22;

    /** Bad host name or host group name" */
    public static final int LSBE_BAD_HOST = 23;

    /** Too many processors requested */
    public static final int LSBE_PROC_NUM = 24;

    /** No host partition in the system */
    public static final int LSBE_NO_HPART = 25;

    /** Bad host partition name */
    public static final int LSBE_BAD_HPART = 26;

    /** No group defined in the system */
    public static final int LSBE_NO_GROUP = 27;

    /** Bad host/user group name */
    public static final int LSBE_BAD_GROUP = 28;

    /** Host is not used by the queue */
    public static final int LSBE_QUEUE_HOST = 29;

    /** User reach UJOB_LIMIT of the queue */
    public static final int LSBE_UJOB_LIMIT = 30;

    /** No host available for migration */
    public static final int LSBE_NO_HOST = 31;


    /** chklog is corrupted */
    public static final int LSBE_BAD_CHKLOG = 32;

    /** User reach PJOB_LIMIT of the queue */
    public static final int LSBE_PJOB_LIMIT = 33;

    /** request from non LSF host rejected*/
    public static final int LSBE_NOLSF_HOST = 34;

/** Error codes related to input arguments of lsblib call */

    /** Bad argument for lsblib call */
    public static final int LSBE_BAD_ARG = 35;

    /** Bad time spec for lsblib call */
    public static final int LSBE_BAD_TIME = 36;

    /** Start time is later than end time */
    public static final int LSBE_START_TIME = 37;

    /** Bad CPU limit specification */
    public static final int LSBE_BAD_LIMIT = 38;

    /** Over hard limit of queue */
    public static final int LSBE_OVER_LIMIT = 39;

    /** Empty job (command) */
    public static final int LSBE_BAD_CMD = 40;

    /** Bad signal value; not supported */
    public static final int LSBE_BAD_SIGNAL = 41;

    /** Bad job name */
    public static final int LSBE_BAD_JOB = 42;

    /** Queue reach QJOB_LIMIT of the queue */
    public static final int LSBE_QJOB_LIMIT = 43;

    /** Expired job terminate time*/
    public static final int LSBE_BAD_TERM = 44;
/** 44 is reserved for future use */

/** Error codes related to lsb.events file */

    /** Unknown event in event log file */
    public static final int LSBE_UNKNOWN_EVENT = 45;

    /** bad event format in event log file */
    public static final int LSBE_EVENT_FORMAT = 46;

    /** End of file */
    public static final int LSBE_EOF = 47;
/** 48-49 are reserved for future use */

/** Error codes related to system failure */

    /** mbatchd internal error */
    public static final int LSBE_MBATCHD = 50;

    /** sbatchd internal error */
    public static final int LSBE_SBATCHD = 51;

    /** lsbatch lib internal error */
    public static final int LSBE_LSBLIB = 52;

    /** LSLIB call fails */
    public static final int LSBE_LSLIB = 53;

    /** System call fails */
    public static final int LSBE_SYS_CALL = 54;

    /** Cannot alloc memory */
    public static final int LSBE_NO_MEM = 55;

    /** Lsbatch service not registered */
    public static final int LSBE_SERVICE = 56;

    /** LSB_SHAREDIR not defined */
    public static final int LSBE_NO_ENV = 57;

    /** chkpnt system call fail */
    public static final int LSBE_CHKPNT_CALL = 58;

    /** mbatchd cannot fork */
    public static final int LSBE_NO_FORK = 59;

/** Error codes related to communication between mbatchd/lsblib/sbatchd */

    /** LSBATCH protocol error */
    public static final int LSBE_PROTOCOL = 60;

    /** XDR en/decode error */
    public static final int LSBE_XDR = 61;

    /** No appropriate port can be bound */
    public static final int LSBE_PORT = 62;

    /** Timeout in contacting mbatchd */
    public static final int LSBE_TIME_OUT = 63;

    /** Timeout on connect() call */
    public static final int LSBE_CONN_TIMEOUT = 64;

    /** Connection refused by server */
    public static final int LSBE_CONN_REFUSED = 65;

    /** server connection already exists */
    public static final int LSBE_CONN_EXIST = 66;

    /** server is not connected */
    public static final int LSBE_CONN_NONEXIST = 67;

    /** sbatchd cannot be reached */
    public static final int LSBE_SBD_UNREACH = 68;

    // Search for any ; \s+ /** and fix the comments
    /** Operation cannot be performed right now, op. will be retried. */
    public static final int LSBE_OP_RETRY = 69;

    /** user has no enough job slots */
    public static final int LSBE_USER_JLIMIT = 70;
/** 71 is reserved for future use */

/** Error codes related to NQS */

    /** Bad specification for a NQS job */
    public static final int LSBE_NQS_BAD_PAR = 72;


    /** Client host has no license */
    public static final int LSBE_NO_LICENSE = 73;

/** Error codes related to calendar */

    /** Bad calendar name */
    public static final int LSBE_BAD_CALENDAR = 74;

    /** No calendar found */
    public static final int LSBE_NOMATCH_CALENDAR = 75;

    /** No calendar in system */
    public static final int LSBE_NO_CALENDAR = 76;

    /** Bad calendar time events */
    public static final int LSBE_BAD_TIMEEVENT = 77;

    /** Calendar exist already */
    public static final int LSBE_CAL_EXIST = 78;

    /** Calendar function is not enabled*/
    public static final int LSBE_CAL_DISABLED = 79;

/** Error codes related to modify job's parameters */

    /** the job's params cannot be changed */
    public static final int LSBE_JOB_MODIFY = 80;
    /** the changed once parameters are not used */
    public static final int LSBE_JOB_MODIFY_ONCE = 81;


    /** the job is not a repetitive job */
    public static final int LSBE_J_UNREPETITIVE = 82;

    /** bad cluster name */
    public static final int LSBE_BAD_CLUSTER = 83;

/** Error codes related jobs driven by calendar */

    /** Job can not be killed in pending */
    public static final int LSBE_PEND_CAL_JOB = 84;
    /** This Running turn is being terminated */
    public static final int LSBE_RUN_CAL_JOB = 85;


    /** Modified parameters are being used */
    public static final int LSBE_JOB_MODIFY_USED = 86;

    /** Can not get user's token */
    public static final int LSBE_AFS_TOKENS = 87;

/** Error codes related to event */

    /** Bad event name */
    public static final int LSBE_BAD_EVENT = 88;

    /** No event found */
    public static final int LSBE_NOMATCH_EVENT = 89;

    /** No event in system */
    public static final int LSBE_NO_EVENT = 90;

/** Error codes related to user, queue and host */

    /** User reach HJOB_LIMIT of the queue */
    public static final int LSBE_HJOB_LIMIT = 91;

/** Error codes related to bmsg */

    /** Message delivered */
    public static final int LSBE_MSG_DELIVERED = 92;
    /** MBD could not find the message that SBD mentions about */
    public static final int LSBE_NO_JOBMSG = 93;

    /** x */
    public static final int LSBE_MSG_RETRY = 94;

/** Error codes related to resource requirement */

    /** Bad resource requirement */
    public static final int LSBE_BAD_RESREQ = 95;


    /** No enough hosts */
    public static final int LSBE_NO_ENOUGH_HOST = 96;

/** Error codes related to configuration lsblib call */

    /** Fatal error in reading conf files */
    public static final int LSBE_CONF_FATAL = 97;

    /** Warning error in reading conf files */
    public static final int LSBE_CONF_WARNING = 98;


    /** CONF used calendar cannot be modified */
    public static final int LSBE_CAL_MODIFY = 99;

    /** Job created calendar cannot be modified */
    public static final int LSBE_JOB_CAL_MODIFY = 100;
    /** FAIRSHARE queue or HPART defined */
    public static final int LSBE_HP_FAIRSHARE_DEF = 101;

    /** No resource specified */
    public static final int LSBE_NO_RESOURCE = 102;

    /** Bad resource name */
    public static final int LSBE_BAD_RESOURCE = 103;
    /** Calendar not allowed for interactive job */
    public static final int LSBE_INTERACTIVE_CAL = 104;
    /** Interactive job cannot be rerunnable */
    public static final int LSBE_INTERACTIVE_RERUN = 105;

    /** PTY and infile specified */
    public static final int LSBE_PTY_INFILE = 106;

    /** JobScheduler is disabled */
    public static final int LSBE_JS_DISABLED = 107;

    /** Submission host and its host type can not be found any more */
    public static final int LSBE_BAD_SUBMISSION_HOST = 108;
    /** Lock the job so that it cann't be resume by sbatchd */
    public static final int LSBE_LOCK_JOB = 109;

    /** user not in the user group */
    public static final int LSBE_UGROUP_MEMBER = 110;
    /** Operation not supported for a Multicluster job */
    public static final int LSBE_UNSUPPORTED_MC = 111;
    /** Operation permission denied for a Multicluster job */
    public static final int LSBE_PERMISSION_MC = 112;

    /** System Calendar exist already */
    public static final int LSBE_SYSCAL_EXIST = 113;

    /** exceed q's resource reservation */
    public static final int LSBE_OVER_RUSAGE = 114;

    /** bad host spec of run/cpu limits */
    public static final int LSBE_BAD_HOST_SPEC = 115;

    /** calendar syntax error */
    public static final int LSBE_SYNTAX_CALENDAR = 116;

    /** delete a used calendar */
    public static final int LSBE_CAL_USED = 117;

    /** cyclic calednar dependence */
    public static final int LSBE_CAL_CYC = 118;

    /** bad user group name */
    public static final int LSBE_BAD_UGROUP = 119;

    /** esub aborted request */
    public static final int LSBE_ESUB_ABORT = 120;

    /** Bad exception handler syntax */
    public static final int LSBE_EXCEPT_SYNTAX = 121;
    /** Bad exception condition specification */
    public static final int LSBE_EXCEPT_COND = 122;
    /** Bad or invalid action specification */
    public static final int LSBE_EXCEPT_ACTION = 123;

    /** job dependence, not deleted immed */
    public static final int LSBE_JOB_DEP = 124;
/** error codes for job group */

    /** the job group exists */
    public static final int LSBE_JGRP_EXIST = 125;

    /** the job group doesn't exist */
    public static final int LSBE_JGRP_NULL = 126;

    /** the group contains jobs */
    public static final int LSBE_JGRP_HASJOB = 127;

    /** the unknown group control signal */
    public static final int LSBE_JGRP_CTRL_UNKWN = 128;

    /** Bad Job Group name */
    public static final int LSBE_JGRP_BAD = 129;

    /** Job Array */
    public static final int LSBE_JOB_ARRAY = 130;

    /** Suspended job not supported */
    public static final int LSBE_JOB_SUSP = 131;

    /** Forwarded job not suported */
    public static final int LSBE_JOB_FORW = 132;

    /** parent group is held */
    public static final int LSBE_JGRP_HOLD = 133;

    /** bad index */
    public static final int LSBE_BAD_IDX = 134;

    /** index too big */
    public static final int LSBE_BIG_IDX = 135;

    /** job array not exist*/
    public static final int LSBE_ARRAY_NULL = 136;

    /** Void calendar */
    public static final int LSBE_CAL_VOID = 137;

    /** the job exists */
    public static final int LSBE_JOB_EXIST = 138;

    /** Job Element fail */
    public static final int LSBE_JOB_ELEMENT = 139;

    /** Bad jobId */
    public static final int LSBE_BAD_JOBID = 140;

    /** cannot change job name */
    public static final int LSBE_MOD_JOB_NAME = 141;

/** error codes for frame job */

    /** Bad frame expression */
    public static final int LSBE_BAD_FRAME = 142;

    /** Frame index too long */
    public static final int LSBE_FRAME_BIG_IDX = 143;

    /** Frame index syntax error */
    public static final int LSBE_FRAME_BAD_IDX = 144;


    /** child process died */
    public static final int LSBE_PREMATURE = 145;

/** error code for user not in project group */

    /** Invoker is not in project group */
    public static final int LSBE_BAD_PROJECT_GROUP = 146;

/** error code for user group / host group */

    /** No host group defined in the system */
    public static final int LSBE_NO_HOST_GROUP = 147;

    /** No user group defined in the system */
    public static final int LSBE_NO_USER_GROUP = 148;

    /** Bad jobid index file format */
    public static final int LSBE_INDEX_FORMAT = 149;

/** error codes for IO_SPOOL facility */

    /** source file does not exist */
    public static final int LSBE_SP_SRC_NOT_SEEN = 150;

    /** Number of failed spool hosts reached max */
    public static final int LSBE_SP_FAILED_HOSTS_LIM = 151;

    /** spool copy failed for this host*/
    public static final int LSBE_SP_COPY_FAILED = 152;

    /** fork failed */
    public static final int LSBE_SP_FORK_FAILED = 153;

    /** status of child is not available */
    public static final int LSBE_SP_CHILD_DIES = 154;

    /** child terminated with failure */
    public static final int LSBE_SP_CHILD_FAILED = 155;

    /** Unable to find a host for spooling */
    public static final int LSBE_SP_FIND_HOST_FAILED = 156;

    /** Cannot get $JOB_SPOOLDIR for this host */
    public static final int LSBE_SP_SPOOLDIR_FAILED = 157;

    /** Cannot delete spool file for this host */
    public static final int LSBE_SP_DELETE_FAILED = 158;


    /** Bad user priority */
    public static final int LSBE_BAD_USER_PRIORITY = 159;

    /** Job priority control undefined */
    public static final int LSBE_NO_JOB_PRIORITY = 160;

    /** Job has been killed & requeued */
    public static final int LSBE_JOB_REQUEUED = 161;

    /** Remote job cannot kill-requeued */
    public static final int LSBE_JOB_REQUEUE_REMOTE = 162;

    /** Cannot submit job array to a NQS queue */
    public static final int LSBE_NQS_NO_ARRJOB = 163;

/** error codes for EXT_JOB_STATUS */

    /** No message available */
    public static final int LSBE_BAD_EXT_MSGID = 164;

    /** Not a regular file */
    public static final int LSBE_NO_IFREG = 165;

    /** MBD fail to create files in the directory*/
    public static final int LSBE_BAD_ATTA_DIR = 166;

    /** Fail to transfer data */
    public static final int LSBE_COPY_DATA = 167;

    /** exceed the limit on data transferring of a msg*/
    public static final int LSBE_JOB_ATTA_LIMIT = 168;

    /** cannot resize a chunk job, cannot bswitch a run/wait job */
    public static final int LSBE_CHUNK_JOB = 169;

/** Error code used in communications with dlogd */


    /** dlogd is already connected */
    public static final int LSBE_DLOGD_ISCONN = 170;

/** Error code for LANL3_1ST_HOST */

    /** Multiple first execution host */
    public static final int LSBE_MULTI_FIRST_HOST = 171;

    /** Host group as first execution host */
    public static final int LSBE_HG_FIRST_HOST = 172;

    /** Host partition as first execution host */
    public static final int LSBE_HP_FIRST_HOST = 173;

    /** "others" as first execution host */
    public static final int LSBE_OTHERS_FIRST_HOST = 174;

/** error code for multi-cluster: remote only queue */

    /** cannot specify exec host */
    public static final int LSBE_MC_HOST = 175;

    /** cannot specify repetitive job */
    public static final int LSBE_MC_REPETITIVE = 176;

    /** cannot be a chkpnt job */
    public static final int LSBE_MC_CHKPNT = 177;

    /** cannot specify exception */
    public static final int LSBE_MC_EXCEPTION = 178;

    /** cannot specify time event */
    public static final int LSBE_MC_TIMEEVENT = 179;

    /** Too few processors requested */
    public static final int LSBE_PROC_LESS = 180;
    /** bmod pending options and running options together towards running job */
    public static final int LSBE_MOD_MIX_OPTS = 181;

    /** cannot bmod remote running job */
    public static final int LSBE_MOD_REMOTE = 182;
    /** cannot bmod cpulimit without LSB_JOB_CPULIMIT defined */
    public static final int LSBE_MOD_CPULIMIT = 183;
    /** cannot bmod memlimit without LSB_JOB_MEMLIMIT defined */
    public static final int LSBE_MOD_MEMLIMIT = 184;

    /** cannot bmod err file name */
    public static final int LSBE_MOD_ERRFILE = 185;

    /** host is locked by master LIM*/
    public static final int LSBE_LOCKED_MASTER = 186;
    /** warning time period is invalid */
    public static final int LSBE_WARNING_INVALID_TIME_PERIOD = 187;
    /** either warning time period or warning action is not specified */
    public static final int LSBE_WARNING_MISSING = 188;
    /** The job arrays involved in  one to one dependency do not  have the same size. */
    public static final int LSBE_DEP_ARRAY_SIZE = 189;

    /** Not enough processors to be reserved (lsb_addreservation()) */
    public static final int LSBE_FEWER_PROCS = 190;

    /** Bad reservation ID */
    public static final int LSBE_BAD_RSVID = 191;

    /** No more reservation IDs can be used now */
    public static final int LSBE_NO_RSVID = 192;

    /** No hosts are exported */
    public static final int LSBE_NO_EXPORT_HOST = 193;

    /** Trying to control remote hosts*/
    public static final int LSBE_REMOTE_HOST_CONTROL = 194;

/*Can't open a remote host closed by the remote cluster admin */
    public static final int LSBE_REMOTE_CLOSED = 195;

    /** User suspended job */
    public static final int LSBE_USER_SUSPENDED = 196;

    /** Admin suspended job */
    public static final int LSBE_ADMIN_SUSPENDED = 197;

    /** Not a local host name in  bhost -e command */
    public static final int LSBE_NOT_LOCAL_HOST = 198;

    /** The host's lease is not active. */
    public static final int LSBE_LEASE_INACTIVE = 199;

    /** The advance reserved host is not on queue. */
    public static final int LSBE_QUEUE_ADRSV = 200;

    /** The specified host(s) is not exported. */
    public static final int LSBE_HOST_NOT_EXPORTED = 201;

    /** The user specified host is not inn advance reservation */
    public static final int LSBE_HOST_ADRSV = 202;

    /** The remote cluster is not connected */
    public static final int LSBE_MC_CONN_NONEXIST = 203;

    /** The general resource limit broken */
    public static final int LSBE_RL_BREAK = 204;

/** ---- The following RMS errors are obsoleted in Eagle */

    /** cannot submit a job with special topology requirement to a preemptive queue*/
    public static final int LSBE_LSF2TP_PREEMPT = 205;

    /** cannot submit a job with special topology requirement to a queue with slot reservation*/
    public static final int LSBE_LSF2TP_RESERVE = 206;
    /** cannot submit a job with special topology requirement to a queue with backill */
    public static final int LSBE_LSF2TP_BACKFILL = 207;
    /** ---- The above RMS errors are obsoleted in Eagle */

    /** none existed policy name */
    public static final int LSBE_RSV_POLICY_NAME_BAD = 208;

    /** All normal user has no privilege */
    public static final int LSBE_RSV_POLICY_PERMISSION_DENIED = 209;

    /** user has no privilege */
    public static final int LSBE_RSV_POLICY_USER = 210;

    /** user has no privilege to create reservation on host */
    public static final int LSBE_RSV_POLICY_HOST = 211;

    /** time window is not allowed by policy */
    public static final int LSBE_RSV_POLICY_TIMEWINDOW = 212;

    /** the feature is disabled */
    public static final int LSBE_RSV_POLICY_DISABLED = 213;
    /** the general limit related errors */

    /** There are no general limit defined */
    public static final int LSBE_LIM_NO_GENERAL_LIMIT = 214;

    /** There are no resource usage */
    public static final int LSBE_LIM_NO_RSRC_USAGE = 215;

    /** Convert data error */
    public static final int LSBE_LIM_CONVERT_ERROR = 216;

    /** There are no qualified host found in cluster*/
    public static final int LSBE_RSV_NO_HOST = 217;

    /** Cannot modify job group on element of job array */
    public static final int LSBE_MOD_JGRP_ARRAY = 218;

    /** Cannot combine modify job group or service class option with others */
    public static final int LSBE_MOD_MIX = 219;

    /** the service class doesn't exist */
    public static final int LSBE_SLA_NULL = 220;

    /** Modify job group for job in service class is not supported*/
    public static final int LSBE_MOD_JGRP_SLA = 221;

    /** User or user group is not a member of the specified service class */
    public static final int LSBE_SLA_MEMBER = 222;

    /** There is no exceptional host found */
    public static final int LSBE_NO_EXCEPTIONAL_HOST = 223;

    /** warning action (signal) is invalid */
    public static final int LSBE_WARNING_INVALID_ACTION = 224;


    /** Extsched option syntax error */
    public static final int LSBE_EXTSCHED_SYNTAX = 225;

    /** SLA doesn't work with remote only queues */
    public static final int LSBE_SLA_RMT_ONLY_QUEUE = 226;

    /** Cannot modify service class on element of job array */
    public static final int LSBE_MOD_SLA_ARRAY = 227;

    /** Modify service class for job in job group is not supported*/
    public static final int LSBE_MOD_SLA_JGRP = 228;

    /** Max. Pending job error */
    public static final int LSBE_MAX_PEND = 229;

    /** System concurrent query exceeded */
    public static final int LSBE_CONCURRENT = 230;

    /** Requested feature not enabled */
    public static final int LSBE_FEATURE_NULL = 231;


    /** Host is already member of group */
    public static final int LSBE_DYNGRP_MEMBER = 232;

    /** Host is not a dynamic host */
    public static final int LSBE_BAD_DYN_HOST = 233;

    /** Host was not added with badmin hghostadd */
    public static final int LSBE_NO_GRP_MEMBER = 234;

    /** Cannot create job info file */
    public static final int LSBE_JOB_INFO_FILE = 235;

    /** Cannot modify rusage to a new || (or) expression after the job is dispatched */
    public static final int LSBE_MOD_OR_RUSAGE = 236;

    /** Bad host group name */
    public static final int LSBE_BAD_GROUP_NAME = 237;

    /** Bad host name */
    public static final int LSBE_BAD_HOST_NAME = 238;

    /** Bsub is not permitted on DT cluster */
    public static final int LSBE_DT_BSUB = 239;


    /** The parent symphony job/group was  gone when submitting jobs*/
    public static final int LSBE_PARENT_SYM_JOB = 240;

    /** The partition has no cpu alllocated */
    public static final int LSBE_PARTITION_NO_CPU = 241;

    /** batch partition does not accept online jobs: obsolete */
    public static final int LSBE_PARTITION_BATCH = 242;

    /** online partition does not accept batch jobs */
    public static final int LSBE_PARTITION_ONLINE = 243;

    /** no batch licenses */
    public static final int LSBE_NOLICENSE_BATCH = 244;

    /** no online licenses */
    public static final int LSBE_NOLICENSE_ONLINE = 245;

    /** signal is not supported for service job */
    public static final int LSBE_SIGNAL_SRVJOB = 246;

    /** the begin time is not later than current time. */
    public static final int LSBE_BEGIN_TIME_INVALID = 247;

    /** the end time is not later than current time. */
    public static final int LSBE_END_TIME_INVALID = 248;

    /** Bad regular expression */
    public static final int LSBE_BAD_REG_EXPR = 249;


    /** Host group has regular expression */
    public static final int LSBE_GRP_REG_EXPR = 250;

    /** Host group have no member */
    public static final int LSBE_GRP_HAVE_NO_MEMB = 251;

    /** the application doesn't exist */
    public static final int LSBE_APP_NULL = 252;

    /** job's proclimit rejected by App */
    public static final int LSBE_PROC_JOB_APP = 253;

    /** app's proclimit rejected by Queue */
    public static final int LSBE_PROC_APP_QUE = 254;

    /** application name is too long */
    public static final int LSBE_BAD_APPNAME = 255;

    /** Over hard limit of queue */
    public static final int LSBE_APP_OVER_LIMIT = 256;

    /** Cannot remove default application */
    public static final int LSBE_REMOVE_DEF_APP = 257;

    /** Host is disabled by EGO */
    public static final int LSBE_EGO_DISABLED = 258;

    /** Host is a remote host. Remote hosts cannot be added to a local host group. */
    public static final int LSBE_REMOTE_HOST = 259;

    /** SLA is exclusive, only accept exclusive job. */
    public static final int LSBE_SLA_EXCLUSIVE = 260;

    /** SLA is non-exclusive, only accept non-exclusive job */
    public static final int LSBE_SLA_NONEXCLUSIVE = 261;

    /** The feature has already been started */
    public static final int LSBE_PERFMON_STARTED = 262;

    /** The Featurn has already been turn down */
    public static final int LSBE_PERFMON_STOPED = 263;

    /** Current sampling period is already set to %%s,seconds. Ignored*/
    public static final int LSBE_PERFMON_PERIOD_SET = 264;

    /** Default spool dir is disabled */
    public static final int LSBE_DEFAULT_SPOOL_DIR_DISABLED = 265;

    /** job belongs to an APS queue and cannot be moved */
    public static final int LSBE_APS_QUEUE_JOB = 266;

    /** job is not in an absolute priority enabled queue */
    public static final int LSBE_BAD_APS_JOB = 267;

    /** Wrong aps admin value */
    public static final int LSBE_BAD_APS_VAL = 268;

    /** Trying to delete a non-existent APS string */
    public static final int LSBE_APS_STRING_UNDEF = 269;

    /** A job cannot be assigned an SLA and an APS queue with factor FS */
    public static final int LSBE_SLA_JOB_APS_QUEUE = 270;

    /** bmod -aps | -apsn option cannot be mixed with other option */
    public static final int LSBE_MOD_MIX_APS = 271;

    /** specified ADMIN factor/system APS value out of range */
    public static final int LSBE_APS_RANGE = 272;

    /** specified ADMIN factor/system APS value is zero */
    public static final int LSBE_APS_ZERO = 273;


    /** res port is unknown */
    public static final int LSBE_DJOB_RES_PORT_UNKNOWN = 274;

    /** timeout on res communication */
    public static final int LSBE_DJOB_RES_TIMEOUT = 275;

    /** I/O error on remote stream */
    public static final int LSBE_DJOB_RES_IOERR = 276;

    /** res internal failure */
    public static final int LSBE_DJOB_RES_INTERNAL_FAILURE = 277;


    /** can not run outside LSF */
    public static final int LSBE_DJOB_CAN_NOT_RUN = 278;

    /** distributed job's validation failed due to incorrect job ID or index */
    public static final int LSBE_DJOB_VALIDATION_BAD_JOBID = 279;

    /** distributed job's validation failed due to incorrect host selection */
    public static final int LSBE_DJOB_VALIDATION_BAD_HOST = 280;

    /** distributed job's validation failed due to incorrect user */
    public static final int LSBE_DJOB_VALIDATION_BAD_USER = 281;

    /** failed while executing tasks */
    public static final int LSBE_DJOB_EXECUTE_TASK = 282;

    /** failed while waiting for tasks to finish*/
    public static final int LSBE_DJOB_WAIT_TASK = 283;


    /** HPC License not exist */
    public static final int LSBE_APS_HPC = 284;

    /** Integrity check of bsub command failed */
    public static final int LSBE_DIGEST_CHECK_BSUB = 285;

    /** Distributed Application Framework disabled */
    public static final int LSBE_DJOB_DISABLED = 286;

/** Error codes related to runtime estimation and cwd */

    /** Bad runtime specification */
    public static final int LSBE_BAD_RUNTIME = 287;

    /** RUNLIMIT: Cannot exceed RUNTIME*JOB_RUNLIMIT_RATIO */
    public static final int LSBE_BAD_RUNLIMIT = 288;

    /** RUNTIME: Cannot exceed the hard runtime limit in the queue */
    public static final int LSBE_OVER_QUEUE_LIMIT = 289;

    /** RUNLIMIT: Is not set by command line */
    public static final int LSBE_SET_BY_RATIO = 290;

    /** current working directory name too long */
    public static final int LSBE_BAD_CWD = 291;


    /** Job group limit is greater than its parent group */
    public static final int LSBE_JGRP_LIMIT_GRTR_THAN_PARENT = 292;

    /** Job group limit is less than its children groups */
    public static final int LSBE_JGRP_LIMIT_LESS_THAN_CHILDREN = 293;

    /** Job Array end index should be specified explicitly */
    public static final int LSBE_NO_ARRAY_END_INDEX = 294;

    /** cannot bmod runtime without LSB_MOD_ALL_JOBS=y defined */
    public static final int LSBE_MOD_RUNTIME = 295;

    /** EP3 */
    public static final int LSBE_BAD_SUCCESS_EXIT_VALUES = 296;
    public static final int LSBE_DUP_SUCCESS_EXIT_VALUES = 297;
    public static final int LSBE_NO_SUCCESS_EXIT_VALUES = 298;

    public static final int LSBE_JOB_REQUEUE_BADARG = 299;
    public static final int LSBE_JOB_REQUEUE_DUPLICATED = 300;

    /** "all" with number */
    public static final int LSBE_JOB_REQUEUE_INVALID_DIGIT = 301;

    /** ~digit without "all" */
    public static final int LSBE_JOB_REQUEUE_INVALID_TILDE = 302;
    public static final int LSBE_JOB_REQUEUE_NOVALID = 303;


    /** No matching job group found */
    public static final int LSBE_NO_JGRP = 304;
    public static final int LSBE_NOT_CONSUMABLE = 305;

/** AR pre/post */

    /** Cannot parse an Advance Reservation -exec string */
    public static final int LSBE_RSV_BAD_EXEC = 306;

    /** Unknown AR event type */
    public static final int LSBE_RSV_EVENTTYPE = 307;

    /** pre/post cannot have postive offset */
    public static final int LSBE_RSV_SHIFT = 308;

    /** pre-AR command cannot have offset < 0 in user-created AR */
    public static final int LSBE_RSV_USHIFT = 309;

    /** only one pre- and one post- cmd permitted per AR */
    public static final int LSBE_RSV_NUMEVENTS = 310;

/*Error codes related to AR Modification*/

    /** ID does not correspond to a known AR. */
    public static final int LSBE_ADRSV_ID_VALID = 311;

    /** disable non-recurrent AR. */
    public static final int LSBE_ADRSV_DISABLE_NONRECUR = 312;

    /** modification is rejected because AR is activated. */
    public static final int LSBE_ADRSV_MOD_ACTINSTANCE = 313;

    /** modification is rejected because host slots is not available. */
    public static final int LSBE_ADRSV_HOST_NOTAVAIL = 314;

    /** the  time of the AR cannot be modified since resource is not available. */
    public static final int LSBE_ADRSV_TIME_MOD_FAIL = 315;

    /** resource requirement (-R) must be followed a slot requirment (-n) */
    public static final int LSBE_ADRSV_R_AND_N = 316;

/*modification is rejected because trying to empty the AR. */
    public static final int LSBE_ADRSV_EMPTY = 317;

/*modification is rejected because switching AR type. */
    public static final int LSBE_ADRSV_SWITCHTYPE = 318;

/*modification is rejected because specifying -n for system AR. */
    public static final int LSBE_ADRSV_SYS_N = 319;

    /** disable string is not valid. */
    public static final int LSBE_ADRSV_DISABLE = 320;

    /** Unique AR ID required */
    public static final int LSBE_ADRSV_ID_UNIQUE = 321;

    /** Bad reservation name */
    public static final int LSBE_BAD_RSVNAME = 322;

    /** Cannot change the start time of an active reservation. */
    public static final int LSBE_ADVRSV_ACTIVESTART = 323;

    /** AR ID is refernced by a job */
    public static final int LSBE_ADRSV_ID_USED = 324;

    /** the disable period has already been disabled */
    public static final int LSBE_ADRSV_PREVDISABLED = 325;

    /** an active period of a recurring reservation cannot be disabled */
    public static final int LSBE_ADRSV_DISABLECURR = 326;

    /** modification is rejected because specified hosts or host groups do not belong to the reservation */
    public static final int LSBE_ADRSV_NOT_RSV_HOST = 327;

/*new parser */

/*checking resreq return ok */
    public static final int LSBE_RESREQ_OK = 328;

/*checking resreq return error */
    public static final int LSBE_RESREQ_ERR = 329;


    /** modification is rejected because reservation has running jobs on the specified hosts or host groups */
    public static final int LSBE_ADRSV_HOST_USED = 330;


    /** The checkpoint directory is too long */
    public static final int LSBE_BAD_CHKPNTDIR = 331;

    /** trying to modify in a remote cluster */
    public static final int LSBE_ADRSV_MOD_REMOTE = 332;
    public static final int LSBE_JOB_REQUEUE_BADEXCLUDE = 333;

    /** trying to disable for a date in the past */
    public static final int LSBE_ADRSV_DISABLE_DATE = 334;

    /** cannot mix the -Un option with others for started jobs */
    public static final int LSBE_ADRSV_DETACH_MIX = 335;

    /** cannot detach a started job when the reservation is active */
    public static final int LSBE_ADRSV_DETACH_ACTIVE = 336;

    /** invalid time expression: must specify day for both start and end time */
    public static final int LSBE_MISSING_START_END_TIME = 337;

    /** Queue level limitation */
    public static final int LSBE_JOB_RUSAGE_EXCEED_LIMIT = 338;

    /** Queue level limitation */
    public static final int LSBE_APP_RUSAGE_EXCEED_LIMIT = 339;

    /** Hosts and host groups specified by -m are not used by the queue */
    public static final int LSBE_CANDIDATE_HOST_EMPTY = 340;

    /** An int must follow an open bracket */
    public static final int LSBE_HS_BAD_AFTER_BRACKT = 341;

    /** An end index must follow a dash */
    public static final int LSBE_HS_NO_END_INDEX = 342;

    /** Integers must come before and after the comma */
    public static final int LSBE_HS_BAD_COMMA = 343;

    /** Incorrect condensed host specification */
    public static final int LSBE_HS_BAD_FORMAT = 344;

    /** The start index must be less than end index */
    public static final int LSBE_HS_BAD_ORDER = 345;

    /** The end index must be less than 10 digits */
    public static final int LSBE_HS_BAD_MANY_DIGITS = 346;

    /** Number of digits in the start index must be less than that of end index */
    public static final int LSBE_HS_BAD_NUM_DIGITS = 347;

    /** The end index cannot start with zero (0) */
    public static final int LSBE_HS_BAD_END_INDEX = 348;

    /** Index must be an integer or a range */
    public static final int LSBE_HS_BAD_INDEX = 349;

/** host group admin*/

    /** When a Host Group Admin (badmin hclose or hopen) closes or opens a host,  the usage of the -C "message" option must be compulsory, as is the logging  of the name of the person performing the action. */
    public static final int LSBE_COMMENTS = 350;


    /** First hosts specified by -m are not used by the queue */
    public static final int LSBE_FIRST_HOSTS_NOT_IN_QUEUE = 351;


    /** The job is not started */
    public static final int LSBE_JOB_NOTSTART = 352;

    /** Accumulated runtime of the job is not available */
    public static final int LSBE_RUNTIME_INVAL = 353;

    /** SSH feature can only be used for interactive job */
    public static final int LSBE_SSH_NOT_INTERACTIVE = 354;

    /** Run time specification is less than the accumulated run time */
    public static final int LSBE_LESS_RUNTIME = 355;

    /** Resize job notification command */
    public static final int LSBE_RESIZE_NOTIFY_CMD_LEN = 356;

    /** Job is not resizable */
    public static final int LSBE_JOB_RESIZABLE = 357;

    /** Bad bresize release host spec */
    public static final int LSBE_RESIZE_RELEASE_HOSTSPEC = 358;

    /** no resize notify matches in mbatchd*/
    public static final int LSBE_NO_RESIZE_NOTIFY = 359;

    /** Can't release first exec host */
    public static final int LSBE_RESIZE_RELEASE_FRISTHOST = 360;

    /** resize event in progress */
    public static final int LSBE_RESIZE_EVENT_INPROGRESS = 361;

    /** too few or too many slots */
    public static final int LSBE_RESIZE_BAD_SLOTS = 362;

    /** No active resize request */
    public static final int LSBE_RESIZE_NO_ACTIVE_REQUEST = 363;

    /** specified host not part of the job's allocation*/
    public static final int LSBE_HOST_NOT_IN_ALLOC = 364;

    /** nothing released */
    public static final int LSBE_RESIZE_RELEASE_NOOP = 365;

    /** Can't resize a brun job */
    public static final int LSBE_RESIZE_URGENT_JOB = 366;
    public static final int LSBE_RESIZE_EGO_SLA_COEXIST = 367;

    /** hpc jobs can't be resized */
    public static final int LSBE_HOST_NOT_SUPPORT_RESIZE = 368;

    /** Application doesn't allow resizable */
    public static final int LSBE_APP_RESIZABLE = 369;

    /** can't operate on lost & found hosts*/
    public static final int LSBE_RESIZE_LOST_AND_FOUND = 370;

    /** can't resize while the first host is lost & found*/
    public static final int LSBE_RESIZE_FIRSTHOST_LOST_AND_FOUND = 371;

    /** bad host name (for resize) */
    public static final int LSBE_RESIZE_BAD_HOST = 372;

    /** proper app is required by an auto-resizable job */
    public static final int LSBE_AUTORESIZE_APP = 373;

    /** cannot resize job because there is a pedning resize request */
    public static final int LSBE_RESIZE_PENDING_REQUEST = 374;

    /** number of hosts specified by -m exceeding configuration */
    public static final int LSBE_ASKED_HOSTS_NUMBER = 375;

    /** All hosts reserved by advanced reservation are invalid in intersected hosts */
    public static final int LSBE_AR_HOST_EMPTY = 376;

    /** First hosts specified by -m are not used by advanced reservation */
    public static final int LSBE_AR_FIRST_HOST_EMPTY = 377;

    /** Internal jobbroker error */
    public static final int LSBE_JB = 378;

    /** Internal jobbroker database library error */
    public static final int LSBE_JB_DBLIB = 379;

    /** Jobbroker cannot reach database */
    public static final int LSBE_JB_DB_UNREACH = 380;

    /** Jobbroker cannot reach mbatchd */
    public static final int LSBE_JB_MBD_UNREACH = 381;

    /** BES server returned an error */
    public static final int LSBE_JB_BES = 382;

    /** Unsupported BES operation */
    public static final int LSBE_JB_BES_UNSUPPORTED_OP = 383;

    /** invalid LS project name*/
    public static final int LSBE_LS_PROJECT_NAME = 384;

    /** the end time is not later than start  time. */
    public static final int LSBE_END_TIME_INVALID_COMPARE_START = 385;

    /** one host cannot be defined in more than one host partition.*/
    public static final int LSBE_HP_REDUNDANT_HOST = 386;

    /** The application level compound resreq causes slots requirements conflict */
    public static final int LSBE_COMPOUND_APP_SLOTS = 387;

    /** The queue level compound resreq causes slots requirements conflict */
    public static final int LSBE_COMPOUND_QUEUE_SLOTS = 388;

    /** Resizable job cannot work with compound resreq */
    public static final int LSBE_COMPOUND_RESIZE = 389;
/** compute unit support */

    /** Compute units cannot have overlapping hosts */
    public static final int LSBE_CU_OVERLAPPING_HOST = 390;

    /** The compute unit cannot contain other compute units */
    public static final int LSBE_CU_BAD_HOST = 391;

    /** The compute unit cannot contain host or host group as a member */
    public static final int LSBE_CU_HOST_NOT_ALLOWED = 392;

    /** Only lowest level compute units are allowed to add hosts as a member */
    public static final int LSBE_CU_NOT_LOWEST_LEVEL = 393;

    /** You cannot modify a compute unit resource requirement when a job is already running */
    public static final int LSBE_CU_MOD_RESREQ = 394;

    /** A compute unit resource requirement cannot be specified for auto resizable jobs */
    public static final int LSBE_CU_AUTORESIZE = 395;

    /** No COMPUTE_UNIT_TYPES are specified in lsb.params */
    public static final int LSBE_NO_COMPUTE_UNIT_TYPES = 396;

    /** No compute unit defined in the system */
    public static final int LSBE_NO_COMPUTE_UNIT = 397;

    /** No such compute unit defined in the system */
    public static final int LSBE_BAD_COMPUTE_UNIT = 398;

    /** The queue is not configured to accept exclusive compute unit jobs */
    public static final int LSBE_CU_EXCLUSIVE = 399;

    /** The queue is not configured to accept higher level of exclusive compute unit jobs */
    public static final int LSBE_CU_EXCLUSIVE_LEVEL = 400;

    /** Job cannot be switched due to the exclusive compute unit reqirement */
    public static final int LSBE_CU_SWITCH = 401;

    /** Job level compound resreq causes slots requirements conflict */
    public static final int LSBE_COMPOUND_JOB_SLOTS = 402;

    /** "||" used in rusage[] of queue resource requirement. It's conflict with job level compound resource requirement */
    public static final int LSBE_COMPOUND_QUEUE_RUSAGE_OR = 403;

    /** balance and usablecuslots cannot both be used in a compute unit resource requirement */
    public static final int LSBE_CU_BALANCE_USABLECUSLOTS = 404;

    /** TS jobs cannot use compound resource requirement (application level) */
    public static final int LSBE_COMPOUND_TSJOB_APP = 405;

    /** TS jobs cannot use compound resource requirement (queue level) */
    public static final int LSBE_COMPOUND_TSJOB_QUEUE = 406;
    /** Job dependency conditions using a job name or job name wild-card exceed limitation set by MAX_JOB_NAME_DEP in lsb.params */
    public static final int LSBE_EXCEED_MAX_JOB_NAME_DEP = 407;

    /** "is waiting for the remote cluster to synchronize." */
    public static final int LSBE_WAIT_FOR_MC_SYNC = 408;

    /** Job cannot exceed queue level RESRSV_LIMIT limitation */
    public static final int LSBE_RUSAGE_EXCEED_RESRSV_LIMIT = 409;

    /** job description too long */
    public static final int LSBE_JOB_DESCRIPTION_LEN = 410;

    /** Cannot use simulation options */
    public static final int LSBE_NOT_IN_SIMMODE = 411;

    /** Value of runtime simulation is incorrect */
    public static final int LSBE_SIM_OPT_RUNTIME = 412;

    /** Value of cputime simulation is incorrect */
    public static final int LSBE_SIM_OPT_CPUTIME = 413;

    /** Incorrect maxmem simulation opt */
    public static final int LSBE_SIM_OPT_MAXMEM = 414;

    /** Incorrect job exitstatus simulation opt */
    public static final int LSBE_SIM_OPT_EXITSTATUS = 415;

    /** Incorrect job simulation option syntax */
    public static final int LSBE_SIM_OPT_SYNTAX = 416;

    /** Number of the above error codes */
    public static final int LSBE_NUM_ERR = 417;

    /**
     * *****************************************************
     */

/* op codes for hand shake protocol between client/server */
    public static final int PREPARE_FOR_OP = 1024;
    public static final int READY_FOR_OP = 1023;

/*
*  Data structures for lsblib interface
 */


    /**
     * \addtogroup lsb_submit_options lsb_submit_options
     * define statements used by lsb_submit.
     */

/* lsb_submit() options */
    /**
     * < Flag to indicate jobName parameter has data. Equivalent to bsub -J command line option existence.
     */
    public static final int SUB_JOB_NAME = 0x01;
    /**
     * < Flag to indicate queue parameter has data. Equivalent to bsub -q command line option existence.
     */
    public static final int SUB_QUEUE = 0x02;
    /**
     * < Flat to indicate numAskedHosts parameter has data. Equivalent to bsub -m command line option existence.
     */
    public static final int SUB_HOST = 0x04;
    /**
     * < Flag to indicate inFile parameter has data. Equivalent to bsub -i command line option existence.
     */
    public static final int SUB_IN_FILE = 0x08;
    /**
     * < Flag to indicate outFile parameter has data. Equivalent to bsub -o command line option existence.
     */
    public static final int SUB_OUT_FILE = 0x10;
    /**
     * < Flag to indicate errFile parameter has data. Equivalent to bsub -e command line option existence.
     */
    public static final int SUB_ERR_FILE = 0x20;
    /**
     * < Flag to indicate execution of a job on a host by itself requested. Equivalent to bsub -x command line option existence.
     */
    public static final int SUB_EXCLUSIVE = 0x40;
    /**
     * < Flag to indicate whether to send mail to the user when the job finishes. Equivalent to bsub -N command line option existence.
     */
    public static final int SUB_NOTIFY_END = 0x80;
    /**
     * < Flag to indicate whether to send mail to the user when the job is dispatched. Equivalent to bsub -B command line option existence.
     */
    public static final int SUB_NOTIFY_BEGIN = 0x100;
    /**
     * < Flag to indicate userGroup name parameter has data. Equivalent to bsub -G command line option existence.
     */
    public static final int SUB_USER_GROUP = 0x200;
    /**
     * < Flag to indicatechkpntPeriod parameter has data . Equivalent to bsub -k command line option existence.
     */
    public static final int SUB_CHKPNT_PERIOD = 0x400;
    /**
     * < Flag to indicate chkpntDir parameter has data. Equivalent to bsub -k command line option existence.
     */
    public static final int SUB_CHKPNT_DIR = 0x800;
    /**
     * < Indicates the job is checkpointable. Equivalent to bsub -k command line option.
     */
    public static final int SUB_CHKPNTABLE = SUB_CHKPNT_DIR;
    /**
     * < Flag to indicate whether to force the job to restart even if non-restartable conditions exist. These conditions are operating system specific. Equivalent to brestart() -f command line option existence.
     */
    public static final int SUB_RESTART_FORCE = 0x1000;
    /**
     * < Flag to indicate restart of a
     * checkpointed job. Only jobs that have been successfully checkpointed
     * can be restarted. Jobs are re-submitted and assigned a new job ID.
     * By default, jobs are restarted with the same output file, file
     * transfer specifications, job name, window signal value, checkpoint
     * directory and period, and rerun options as the original job. To
     * restart a job on another host, both hosts must be binary compatible,
     * run the same OS version, have access to the executable, have access
     * to all open files (LSF must locate them with an absolute path name),
     * and have access to the checkpoint directory. Equivalent to bsub -k
     * command line option existence.
     */
    public static final int SUB_RESTART = 0x2000;
    /**
     * < Indicates the job is re-runnable.
     * If the execution host of the job is considered down, the batch
     * system will re-queue this job in the same job queue, and re-run
     * it from the beginning when a suitable host is found. Everything
     * will be as if it were submitted as a new job, and a new job ID will
     * be assigned. The user who submitted the failed job will receive a
     * mail notice of the job failure, requeueing of the job, and the
     * new job ID.
     * <p/>
     * For a job that was checkpointed before the execution host went down,
     * the job will be restarted from the last checkpoint. Equivalent to
     * bsub -r command line option existence.
     */
    public static final int SUB_RERUNNABLE = 0x4000;
    /**
     * < Flag to indicate sigValue parameter
     * has data. Sends a signal as the queue window closes.
     */
    public static final int SUB_WINDOW_SIG = 0x8000;
    /**
     * < Flag to indicate hostSpec parameter
     * has data.
     */
    public static final int SUB_HOST_SPEC = 0x10000;
    /**
     * < Flag to indicate dependCond parameter
     * has data. Equivalent to bsub -w command line option existence.
     */
    public static final int SUB_DEPEND_COND = 0x20000;
    /**
     * < Flag to indicate resReq parameter
     * has data. Equivalent to bsub -R command line option existence.
     */
    public static final int SUB_RES_REQ = 0x40000;
    /**
     * < Flag to indicate nxf parameter and structure xf have data.
     */
    public static final int SUB_OTHER_FILES = 0x80000;
    /**
     * < Flag to indicate preExecCmd
     * parameter has data. Equivalent to bsub -E command line option
     * existence.
     */
    public static final int SUB_PRE_EXEC = 0x100000;
    /**
     * < Equivalent to bsub -L command line option existence.
     */
    public static final int SUB_LOGIN_SHELL = 0x200000;
    /**
     * < Flag to indicate mailUser parameter has data.
     */
    public static final int SUB_MAIL_USER = 0x400000;
    /**
     * < Flag to indicate newCommand parameter has data. Equivalent to bmod bsub_options existence.
     */
    public static final int SUB_MODIFY = 0x800000;
    /**
     * < Flag to indicate modify option once.
     */
    public static final int SUB_MODIFY_ONCE = 0x1000000;
    /**
     * < Flag to indicate ProjectName
     * parameter has data . Equivalent to bsub -P command line option
     * existence.
     */
    public static final int SUB_PROJECT_NAME = 0x2000000;
    /**
     * < Indicates that the job is submitted
     * as a batch interactive job. When this flag is given, \ref lsb_submit
     * does not return unless an error occurs during the submission process.
     * When the job is started, the user can interact with the job's
     * standard input and output via the terminal. See the -I option
     * in bsub for the description of a batch interactive job. Unless
     * the SUB_PTY flag is specified, the job will run without a
     * pseudo-terminal. Equivalent to bsub -I command line option.
     */
    public static final int SUB_INTERACTIVE = 0x4000000;
    /**
     * < Requests pseudo-terminal support
     * for a job submitted with the SUB_INTERACTIVE flag. This flag is
     * ignored if SUB_INTERACTIVE is not specified. A pseudo-terminal
     * is required to run some applications (such as: vi). Equivalent to
     * bsub -Ip command line option.
     */
    public static final int SUB_PTY = 0x8000000;
    /**< Requests pseudo-terminal shell
     *  mode support for a job submitted with the SUB_INTERACTIVE and
     *  SUB_PTY flags. This flag is ignored if SUB_INTERACTIVE and SUB_PTY
     *  are not specified. This flag should be specified for submitting
     *  interactive shells, or applications which redefine the ctrl-C and
     *  ctrl-Z keys (such as: jove). Equivalent to bsub -Is
     *  command line option. */
    public static final int SUB_PTY_SHELL = 0x10000000;

    /**
     * < Exception handler for job.
     */
    public static final int SUB_EXCEPT = 0x20000000;

    /**
     * < Specifies time_event.
     */
    public static final int SUB_TIME_EVENT = 0x40000000;
/* the last bit 0x80000000 is reserved for internal use */

    /**
     * \addtogroup lsb_submit_options2 lsb_submit_options2
     * define statements used by \ref lsb_submit.
     */

    /**< Hold the job after it is submitted. The job will be in PSUSP status. Equivalent to bsub -H command line option. */
    public static final int SUB2_HOLD = 0x01;

    /**
     * < New cmd for bmod.
     */
    public static final int SUB2_MODIFY_CMD = 0x02;

    /**//* Removed access to SUB2_BSUB_BLOCK since it exits the process (including the JVM) with the exit code of the submitted job. -kshakir December 14, 2010
     * < Submit a job in a synchronous
     * mode so that submission does not return until the job terminates.
     * Note once this flag is set, the \ref lsb_submit will never return if
     * the job is accepted by LSF. Programs that wishes to know the status
     * of the submission needs to fork, with the child process invoking the
     * API call in the blocking mode and the parent process wait on the
     * child process (see wait() for details.
     */
    //public static final int SUB2_BSUB_BLOCK = 0x04;

    /**
     * < Submit from NT.
     */
    public static final int SUB2_HOST_NT = 0x08;

    /**
     * < Submit fom UNIX.
     */
    public static final int SUB2_HOST_UX = 0x10;

    /**
     * < Submit to a chkpntable queue.
     */
    public static final int SUB2_QUEUE_CHKPNT = 0x20;

    /**
     * < Submit to a rerunnable queue.
     */
    public static final int SUB2_QUEUE_RERUNNABLE = 0x40;

    /**
     * < Spool job command.
     */
    public static final int SUB2_IN_FILE_SPOOL = 0x80;

    /**
     * < Inputs the specified file with spooling
     */
    public static final int SUB2_JOB_CMD_SPOOL = 0x100;

    /**
     * < Submits job with priority.
     */
    public static final int SUB2_JOB_PRIORITY = 0x200;

    /**
     * < Job submitted without -n, use queue's default proclimit
     */
    public static final int SUB2_USE_DEF_PROCLIMIT = 0x400;

    /**
     * < bmod -c/-M/-W/-o/-e
     */
    public static final int SUB2_MODIFY_RUN_JOB = 0x800;

    /**
     * < bmod options only to pending jobs
     */
    public static final int SUB2_MODIFY_PEND_JOB = 0x1000;

    /**
     * < Job action warning time. Equivalent to bsub or bmod -wt.
     */
    public static final int SUB2_WARNING_TIME_PERIOD = 0x2000;

    /**
     * < Job action to be taken before a job control action occurs. Equivalent to bsub or bmod -wa.
     */
    public static final int SUB2_WARNING_ACTION = 0x4000;

    /**
     * < Use an advance reservation created with the brsvadd command. Equivalent to bsub -U.
     */
    public static final int SUB2_USE_RSV = 0x8000;

    /**
     * < Windows Terminal Services job
     */
    public static final int SUB2_TSJOB = 0x10000;

/* SUB2_LSF2TP is obsolete in Eagle. We keep it here for backward
*  compatibility */

    /**
     * < Parameter is deprecated
     */
    public static final int SUB2_LSF2TP = 0x20000;

    /**
     * < Submit into a job group
     */
    public static final int SUB2_JOB_GROUP = 0x40000;

    /**
     * < Submit into a service class
     */
    public static final int SUB2_SLA = 0x80000;

    /**
     * < Submit with -extsched options
     */
    public static final int SUB2_EXTSCHED = 0x100000;

    /**
     * < License Scheduler project
     */
    public static final int SUB2_LICENSE_PROJECT = 0x200000;

    /**
     * < Overwrite the standard output of the job. Equivalent to bsub -oo.
     */
    public static final int SUB2_OVERWRITE_OUT_FILE = 0x400000;

    /**
     * < Overwrites the standard error output of the job. Equivalent to bsub -eo.
     */
    public static final int SUB2_OVERWRITE_ERR_FILE = 0x800000;

/* Following are for symphony submission definition.
*  Note that SYM_GRP is an LSF job, which represents a symphony group.
 */

    /**
     * < (symphony) session job
     */
    public static final int SUB2_SSM_JOB = 0x1000000;

    /**
     * < (symphony) symphony job
     */
    public static final int SUB2_SYM_JOB = 0x2000000;

    /**
     * < (symphony) service(LSF) job
     */
    public static final int SUB2_SRV_JOB = 0x4000000;

    /**
     * < (symphony) "group" job
     */
    public static final int SUB2_SYM_GRP = 0x8000000;

    /**
     * < (symphony) symphony job has child symphony job
     */
    public static final int SUB2_SYM_JOB_PARENT = 0x10000000;

    /**
     * < (symphony) symphony job has real time feature
     */
    public static final int SUB2_SYM_JOB_REALTIME = 0x20000000;

    /**
     * < (symphony) symphony job has dummy feature to hold all persistent service jobs.
     */
    public static final int SUB2_SYM_JOB_PERSIST_SRV = 0x40000000;

    /**
     * < Persistent session job
     */
    public static final int SUB2_SSM_JOB_PERSIST = 0x80000000;

    /**
     *  \addtogroup lsb_submit_options3 lsb_submit_options3
     *  define statements used by \ref lsb_submit.
     */

    /**
     * < Application profile name. Equivalent to bsub -app.
     */
    public static final int SUB3_APP = 0x01;

    /**
     * < Job rerunable because of application profile
     */
    public static final int SUB3_APP_RERUNNABLE = 0x02;

    /**
     * < Job modified with absolute priority. Equivalent to bmod -aps.
     */
    public static final int SUB3_ABSOLUTE_PRIORITY = 0x04;

    /**
     * < Submit into a default job group. Equivalent to bsub -g.
     */
    public static final int SUB3_DEFAULT_JOBGROUP = 0x08;

    /**
     * < Run the specified post-execution command on the execution host after the job finishes. Equivalent to bsub -Ep.
     */
    public static final int SUB3_POST_EXEC = 0x10;
    /**
     * < Pass user shell limits to execution host. Equivalent to bsub -ul.
     */
    public static final int SUB3_USER_SHELL_LIMITS = 0x20;
    /**
     * < Current working directory specified on the command line with bsub -cwd
     */
    public static final int SUB3_CWD = 0x40;
    /**< Runtime estimate. Equivalent to bsub -We. Use in conjunction with SUB3_RUNTIME_ESTIMATION_ACC and SUB3_RUNTIME_ESTIMATION_PERC. */
    public static final int SUB3_RUNTIME_ESTIMATION = 0x80;

    /**
     * < Job is not rerunnable. Equivalent to bsub -rn.
     */
    public static final int SUB3_NOT_RERUNNABLE = 0x100;

    /**
     * < Job level requeue exit values.
     */
    public static final int SUB3_JOB_REQUEUE = 0x200;
    /**
     * < Initial checkpoint period. Equivalent to bsub -k initial_checkpoint_period.
     */
    public static final int SUB3_INIT_CHKPNT_PERIOD = 0x400;
    /**< Job migration threshold. Equivalent to bsub -mig migration_threshold. */
    public static final int SUB3_MIG_THRESHOLD = 0x800;

    /**
     * < Checkpoint dir was set by application profile
     */
    public static final int SUB3_APP_CHKPNT_DIR = 0x1000;
    /**
     * < Value of BSUB_CHK_RESREQ environment variable, used for select section resource requirement string syntax checking with bsub -R. bsub only checks the resreq syntax.
     */
    public static final int SUB3_BSUB_CHK_RESREQ = 0x2000;
    /**
     * < Runtime estimate that is the accumulated run time plus the runtime estimate. Equivalent to bmod -We+. Use in conjunction with SUB3_RUNTIME_ESTIMATION.
     */
    public static final int SUB3_RUNTIME_ESTIMATION_ACC = 0x4000;
    /**
     * < Runtime estimate in percentage of completion. Equivalent to bmod -Wep. Two digits after the decimal point are suported. The highest eight bits of runtimeEstimation in the submit structure are used for the integer; the remaining bits are used for the fraction. Use in conjunction with SUB3_RUNTIME_ESTIMATION.
     */
    public static final int SUB3_RUNTIME_ESTIMATION_PERC = 0x8000;

    /**
     * < Protects the sessions of interactive jobs with SSH encryption. Equivalent to bsub -IS|-ISp|-ISs.
     */
    public static final int SUB3_INTERACTIVE_SSH = 0x10000;
    /**< Protect the sessions of interactive x-window job with SSH encryption. Equivalent to bsub -IX.*/
    public static final int SUB3_XJOB_SSH = 0x20000;

    /**
     * < If set the submitted job is auto-resizable
     */
    public static final int SUB3_AUTO_RESIZE = 0x40000;

    /**
     * < If set, the resize notify cmd specified
     */
    public static final int SUB3_RESIZE_NOTIFY_CMD = 0x80000;


    /**
     * < Job broker bulk submit
     */
    public static final int SUB3_BULK_SUBMIT = 0x100000;

    /**
     * < tty mode for interactive job
     */
    public static final int SUB3_INTERACTIVE_TTY = 0x200000;

    /**
     * < Job submitted from floating client
     */
    public static final int SUB3_FLOATING_CLIENT = 0x400000;

    /**
     * < ssh X11 forwarding (bsub -XF)
     */
    public static final int SUB3_XFJOB = 0x800000;

    /**
     * < ssh X11 forwarding (bsub -XF) without bsub -I...
     */
    public static final int SUB3_XFJOB_EXCLUSIVE = 0x1000000;

    /**
     * < Job description.
     */
    public static final int SUB3_JOB_DESCRIPTION = 0x2000000;

    /**
     * < Job submitted from floating client
     */
    public static final int SUB3_SIMULATION = 0x4000000;

/* Check whether a job is symphony job. These macros should be used by all
*  components, including ("submit" actually):
*    - mbatchd: jData->submitReq
*    - sbatchd: jobCard->jobSpecs
*    - API: lsb_submit() and lsb_readjobinfo()
 */

    public static boolean IS_SSM_JOB(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SSM_JOB);
    }

    public static boolean IS_SSM_JOB_PERSIST(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SSM_JOB_PERSIST);
    }

    public static boolean IS_SYM_JOB(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SYM_JOB);
    }

    public static boolean IS_SYM_JOB_PARENT(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SYM_JOB_PARENT);
    }

    public static boolean IS_SYM_JOB_REALTIME(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SYM_JOB_REALTIME);
    }

    public static boolean IS_SYM_JOB_PERSIST_SRV(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SYM_JOB_PERSIST_SRV);
    }

    public static boolean IS_SRV_JOB(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SRV_JOB);
    }

    public static boolean IS_SYM_GRP(int option) {
        return JNAUtils.toBoolean((option) & SUB2_SYM_GRP);
    }

    public static boolean IS_SYM_JOB_OR_SYM_GRP (int option)  { return (IS_SYM_JOB(option) || IS_SYM_GRP(option)); }
/* symphony job for which resource usage should be collected */
    public static boolean IS_REAL_SYM_JOB (int option)  { return (IS_SYM_JOB(option) && !IS_SYM_JOB_PERSIST_SRV(option)); }

    public static boolean IS_WLM_JOB (int option)  { return (IS_SSM_JOB(option) || IS_SYM_JOB(option) || IS_SRV_JOB(option) || IS_SYM_GRP(option)); }
    public static boolean IS_BATCH_JOB (int option)  { return (!IS_WLM_JOB(option)); }
/* job for which resource usage should be collected */
    public static boolean IS_JOB_FOR_ACCT (int option)  { return (IS_REAL_SYM_JOB(option) || IS_BATCH_JOB(option)); }

    public static boolean IS_JOB_FOR_SYM (int option)  { return (IS_SYM_JOB(option) || IS_SRV_JOB(option) || IS_SYM_GRP(option)); }

/* Don't send IS_SYM_JOB/IS_SYM_GRP jobs to scheduler;
*  neither publish events nor brun the job allowed.
 */
    // NOTE: Don't know what this jp struct is.
    //public static boolean IS_SYM_JOB_OR_GRP (int jp)   { return (   (jp) != null && (jp)->shared != null && (  IS_SYM_JOB((jp)->shared->jobBill.options2) ||IS_SYM_GRP((jp)->shared->jobBill.options2))); }

/* name of the lost and find queue and host */
    public static final String LOST_AND_FOUND = "lost_and_found";

    public static final int DELETE_NUMBER = -2;
    public static final int DEL_NUMPRO = LibLsf.INFINIT_INT;
    public static final int DEFAULT_NUMPRO = LibLsf.INFINIT_INT - 1;
    /**
     *  \addtogroup calendar_command  calendar_command
     *  options  for user calendar commands
     */

    /**
     * < Add calenda
     */
    public static final int CALADD = 1;

    /**
     * < Modify calenda
     */
    public static final int CALMOD = 2;

    /**
     * < Delete calenda
     */
    public static final int CALDEL = 3;

    /**
     * < Undelete calenda
     */
    public static final int CALUNDEL = 4;

    /**
     * < Calenda occs
     */
    public static final int CALOCCS = 5;

/* for user event commands */
    public static final int EVEADD = 1;
    public static final int EVEMOD = 2;
    public static final int EVEDEL = 3;

    public static final int PLUGIN_REQUEUE = 126;
    public static final int PLUGIN_EXIT = 125;

    /**
     * \brief  xFile
     */
    public static class xFile extends Structure {
        public static class ByReference extends xFile implements Structure.ByReference {}
        public static class ByValue extends xFile implements Structure.ByValue {}
        public xFile() {}
        public xFile(Pointer p) { super(p); read(); }


        /**
         * < Pathname at submission host
         */
        public String subFn;

        /**
         * < Pathname at execution host
         */
        public String execFn;
        /**
         *  \addtogroup defs_lsb_XF_OP defs_lsb_XF_OP
         *  options  xFile operation
         */

        /**
         * < Transfer files from submit peer to  execution peer
         */
        public static final int XF_OP_SUB2EXEC = 0x1;

        /**
         * < Transfer files from execution peer to  submit peer
         */
        public static final int XF_OP_EXEC2SUB = 0x2;

        /**
         * < Transfer files from submit peer to  execution peer with appending mode
         */
        public static final int XF_OP_SUB2EXEC_APPEND = 0x4;

        /**
         * < Transfer files from execution peer to  submit peer with appending mode
         */
        public static final int XF_OP_EXEC2SUB_APPEND = 0x8;
        public static final int XF_OP_URL_SOURCE = 0x10;

        /**
         * < Defined in \ref defs_lsb_XF_OP
         */
        public int options;
    }



    /* For NQS */
    public static final int NQS_ROUTE = 0x1;
    public static final int NQS_SIG = 0x2;
    public static final int NQS_SERVER = 0x4;


    public static final int MAXNFA = 1024;
    public static final int MAXTAG = 10;

    public static final int OKP = 1;
    public static final int NOP = 0;

    public static final int CHR = 1;
    public static final int ANY = 2;
    public static final int CCL = 3;
    public static final int BOL = 4;
    public static final int EOL = 5;
    public static final int BOT = 6;
    public static final int EOT = 7;
    public static final int BOW = 8;
    public static final int EOW = 9;
    public static final int REF = 10;
    public static final int CLO = 11;

    public static final int END = 0;

    /**
     *  The following defines are not meant to be changeable.
     *  They are for readability only.
     */

    public static final int MAXCHR = 128;
    public static final int CHRBIT = 8;
    public static final int BITBLK = MAXCHR / CHRBIT;
    public static final int BLKIND = 0xAA;
    public static final int BITIND = 0x7;

    public static final int ASCIIB = 0x7F;

    /**
     *  byte classification table for word boundary operators BOW
     *  and EOW. the reason for not using ctype macros is that we can
     *  let the user add into our own table. see re_modw. This table
     *  is not in the bitset form, since we may wish to extend it in the
     *  future for other byte classifications.
     *
     *   TRUE for 0-9 A-Z a-z _
     */

    public static final byte[] chrtyp = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
            0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 1, 0, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 0, 0, 0, 0, 0
    };

    public static int inascii(int x) {
        return (0x7F & (x));
    }

    public static int iswordc(int x) {
        return chrtyp[inascii(x)];
    }

/*
*  skip values for CLO XXX to skip past the closure
 */


/* [CLO] ANY END ... */
    public static final int ANYSKIP = 2;

/* [CLO] CHR chr END ... */
    public static final int CHRSKIP = 3;

/* [CLO] CCL 16bytes END ... */
    public static final int CCLSKIP = 18;

/*  In LSF7.0.6, we introduce submit_ext structure to support
*   extended fields for furture added submit options.
*   Each new options should have a unique key defined here.
*   The new defined key should be bigger than 1000.
*   Keys below 1000 are used for internal use.
 */

/* submit_ext test */
    public static final int JDATA_EXT_TEST = 1001;

/* LSF simulator: simReq */
    public static final int JDATA_EXT_SIMREQ = 1002;

/* structure for lsb_submit() call */

    /**
     * \extend submit data structure
     */
    public static class submit_ext extends Structure {
        public static class ByReference extends submit_ext implements Structure.ByReference {}
        public static class ByValue extends submit_ext implements Structure.ByValue {}
        public submit_ext() {}
        public submit_ext(Pointer p) { super(p); read(); }


        /**
         * < number of key value pairs.
         */
        public int num;

        /**
         * < Array of keys of the extended fields.
         */
        public Pointer keys;

        /**
         * < Array of values of the extended fields
         */
        public Pointer values;
    }




    /**
     * \brief  submit request structure.
     */
    public static class submit extends Structure {
        public static class ByReference extends submit implements Structure.ByReference {}
        public static class ByValue extends submit implements Structure.ByValue {}
        public submit() {}
        public submit(Pointer p) { super(p); read(); }


        /**
         * <  <lsf/lsbatch.h> defines the flags in \ref lsb_submit_options constructed from bits. These flags correspond to some of the options of the bsub command line. Use the bitwise OR to set more than one flag.
         */
        public int options;


        /**
         * < Extended bitwise inclusive OR of some of the flags in \ref lsb_submit_options2.
         */
        public int options2;


        /**
         * < The job name. If jobName is null, command is used as the job name.
         */
        public String jobName;

        /**
         * < Submit the job to this queue. If queue is null, submit the job to a system default queue.
         */
        public String queue;

        /**
         * < The number of invoker specified candidate hosts for running the job. If numAskedHosts is 0, all qualified hosts will be considered.
         */
        public int numAskedHosts;

        /**
         * < The array of names of invoker specified candidate hosts.  The number of hosts is given by numAskedHosts.
         */
        public Pointer askedHosts;

        /**
         * < The resource requirements of the job. If resReq is null, the batch system will try to obtain resource requirements for command from the remote task lists (see \ref ls_task ). If the task does not appear in the remote task lists, then the default resource requirement is to run on host() of the same type.
         */
        public String resReq;

        /**
         * < Limits on the consumption of system resources by all processes belonging to this job. See getrlimit() for details. If an element of the array is -1, there is no limit for that resource. For the constants used to index the array, see \ref lsb_queueinfo .
         */
        public int[] rLimits = new int[LibLsf.LSF_RLIM_NLIMITS];

        /**
         * < Specify the host model to use for scaling rLimits[LSF_RLIMIT_CPU] and rLimits[LSF_RLIMIT_RUN]. (See \ref lsb_queueinfo). If hostSpec is null, the local host is assumed.
         */
        public String hostSpec;

        /**
         * <  The initial number of processors needed by a (parallel) job. The default is 1.
         */
        public int numProcessors;

        /**
         * < The job dependency condition.
         */
        public String dependCond;

        /**
         * <  Time event string
         */
        public String timeEvent;

        /**
         * <  Dispatch the job on or after beginTime, where beginTime is the number of seconds since 00:00:00 GMT, Jan. 1, 1970 (See time(), ctime()). If beginTime is 0, start the job as soon as possible.
         */
        public NativeLong beginTime;

        /**
         * <  The job termination deadline. If the job is still running at termTime, it will be sent a USR2 signal. If the job does not terminate within 10 minutes after being sent this signal, it will be ended. termTime has the same representation as beginTime. If termTime is 0, allow the job to run until it reaches a resource limit.
         */
        public NativeLong termTime;

        /**
         * < Applies to jobs submitted to a queue that has a run window (See \ref lsb_queueinfo). Send signal sigValue to the job 10 minutes before the run window is going to close. This allows the job to clean up or checkpoint itself, if desired. If the job does not terminate 10 minutes after being sent this signal, it will be suspended.
         */
        public int sigValue;

        /**
         * < The path name of the job's standard input file. If inFile is null, use /dev/null as the default.
         */
        public String inFile;

        /**
         * < The path name of the job's standard output file. If outFile is null, the job's output will be mailed to the submitter
         */
        public String outFile;

        /**
         * < The path name of the job's standard error output file. If errFile is null, the standard error output will be merged with the standard output of the job.
         */
        public String errFile;

        /**
         * < When submitting a job, the command line of the job.   When modifying a job, a mandatory parameter that  should be set to jobId in string format.
         */
        public String command;

        /**
         * < New command line for bmod.
         */
        public String newCommand;

        /**
         * < The job is checkpointable with a period of chkpntPeriod seconds. The value 0 disables periodic checkpointing.
         */
        public NativeLong chkpntPeriod;

        /**
         * < The directory where the chk directory for this job checkpoint files will be created. When a job is checkpointed, its checkpoint files are placed in chkpntDir/chk. chkpntDir can be a relative or absolute path name.
         */
        public String chkpntDir;

        /**
         * < The number of files to transfer.
         */
        public int nxf;

        /**
         * < The array of file transfer specifications. (The xFile structure is defined in <lsf/lsbatch.h>.)
         */
        public Pointer /* xFile.ByReference */ xf;

        /**
         * < The job pre-execution command.
         */
        public String preExecCmd;

        /**
         * < The user that results are mailed to.
         */
        public String mailUser;

        /**
         * < Delete options in options field.
         */
        public int delOptions;

        /**
         * < Extended delete options in options2 field.
         */
        public int delOptions2;

        /**
         * < The name of the project the job will be charged to.
         */
        public String projectName;

        /**
         * < Maximum number of processors required to run the job.
         */
        public int maxNumProcessors;

        /**
         * < Specified login shell used to initialize the execution environment for the job (see the -L option of bsub).
         */
        public String loginShell;

        /**
         * < The name of the LSF user group (see lsb.users) to which the job will belong. (see the -G option of bsub)
         */
        public String userGroup;

        /**
         * < Passes the exception handlers to mbatchd during a job. (see the -X option of bsub). Specifies execption handlers that tell the system how to respond to an exceptional condition for a job. An action is performed when any one of the following exceptions is detected: - \b missched - A job has not been scheduled within the time event specified in the -T option. - \b overrun - A job did not finish in its maximum time (maxtime). - \b underrun - A job finished before it reaches its minimum running time (mintime). - \b abend - A job terminated abnormally. Test an exit code that is one value, two or more comma separated values, or a range of values (two values separated by a `-' to indivate a range). If the job exits with one of the tested values, the abend condition is detected. - \b startfail - A job did not start due to insufficient system resources. - \b cantrun - A job did not start because a dependency condition (see the -w option of bsub) is invalid, or a startfail exception occurs 20 times in a row and the job is suspended. For jobs submitted with a time event (see the -T option of bsub), the cantrun exception condition can be detected once in each time event. - \b hostfail - The host running a job becomes unavailable. When one or more of the above exceptions is detected, you can specify one of the following actions to be taken: - \b alarm - Triggers an alarm incident (see balarms(1)). The alarm can be viewed, acknowledged and resolved. - \b setexcept - Causes the exception event event_name to be set. Other jobs waiting on the exception event event_name specified through the -w option can be triggered. event_name is an arbitrary string. - \b rerun - Causes the job to be rescheduled for execution. Any dependencies associated with the job must be satisfied before re-execution takes place. The rerun action can only be specified for the abend and hostfail exception conditions. The startfail exception condition automatically triggers the rerun action. - \b kill - Causes the current execution of the job to be terminated. This action can only be specified for the overrun exception condition.
         */
        public String exceptList;


        /**
         * < User priority for fairshare scheduling.
         */
        public int userPriority;

        /**
         * < Reservation ID for advance reservation.
         */
        public String rsvId;

        /**
         * < Job group under which the job runs.
         */
        public String jobGroup;

        /**
         * < SLA under which the job runs.
         */
        public String sla;

        /**
         * < External scheduler options.
         */
        public String extsched;

        /**
         * < Warning time period in seconds, -1 if unspecified.
         */
        public int warningTimePeriod;

        /**
         * < Warning action, SIGNAL | CHKPNT | command, null if unspecified.
         */
        public String warningAction;

        /**
         * < License Scheduler project name.
         */
        public String licenseProject;

        /**
         * < Extended bitwise inclusive OR of options flags in \ref lsb_submit_options3.
         */
        public int options3;

        /**
         * < Extended delete options in options3 field.
         */
        public int delOptions3;

        /**
         * < Application profile under which the job runs.
         */
        public String app;

        /**
         * < -1 if no -jsdl and -jsdl_strict options. - 0 -jsdl_strict option - 1 -jsdl option
         */
        public int jsdlFlag;

        /**
         * < JSDL filename
         */
        public String jsdlDoc;

        /**
         * < ARM correlator
         */
        public Pointer correlator;

        /**
         * <  Absolute priority scheduling string set by administrators to denote static system APS value or ADMIN factor APS value. This field is ignored by \ref lsb_submit.
         */
        public String apsString;

        /**
         * < Post-execution commands specified by -Ep option of bsub and bmod.
         */
        public String postExecCmd;

        /**
         * < Current working directory specified by -cwd option of bsub and bmod.
         */
        public String cwd;

        /**
         * < Runtime estimate specified by -We option of bsub and bmod.
         */
        public int runtimeEstimation;

        /**
         * < Job-level requeue exit values specified by -Q option of bsub and bmod.
         */
        public String requeueEValues;

        /**
         * < Initial checkpoint period specified by -k option of bsub and bmod.
         */
        public int initChkpntPeriod;

        /**
         * < Job migration threshold specified by -mig option of bsub and bmod.
         */
        public int migThreshold;

        /**
         * < Job resize notification command to be invoked on the first execution host when a resize request has been satisfied.
         */
        public String notifyCmd;

        /**
         * < Job description.
         */
        public String jobDescription;
/* #if defined(LSF_SIMULATOR)

/**< simulation related options */
        /*public String simReq;*/
        /* #endif */

        /**
         * < For new options in future
         */
        public submit_ext.ByReference submitExt;
    }




    /**
     * \brief submit reply.
     */
    public static class submitReply extends Structure {
        public static class ByReference extends submitReply implements Structure.ByReference {}
        public static class ByValue extends submitReply implements Structure.ByValue {}
        public submitReply() {}
        public submitReply(Pointer p) { super(p); read(); }


        /**
         * < The queue the job was submitted to.
         */
        public String queue;

        /**
         * < DependCond contained badJobId but badJobId does not exist in the system.
         */
        public long badJobId;

        /**
         * < DependCond contained badJobName but badJobName does not exist in the system. If the environment variable BSUB_CHK_RESREQ is set, the value of lsberrno is either LSBE_RESREQ_OK or LSBE_RESREQ_ERR, depending on the result of resource requirement string checking. The badJobName field contains the detailed error message.
         */
        public String badJobName;

        /**< If lsberrno is LSBE_BAD_HOST,
         *  (**askedHosts)[badReqIndx] is not a host known to the system.
         *  If lsberrno is LSBE_QUEUE_HOST, (**askedHosts)[badReqIndx]
         *  is not a host used by the specified queue. If lsberrno is
         *  LSBE_OVER_LIMIT, (*rLimits)[badReqIndx] exceeds the queue's
         *  limit for the resource. */
        public int badReqIndx;
    }



    /**
     * \brief  submit migration request.
     */
    public static class submig extends Structure {
        public static class ByReference extends submig implements Structure.ByReference {}
        public static class ByValue extends submig implements Structure.ByValue {}
        public submig() {}
        public submig(Pointer p) { super(p); read(); }


        /**
         * < The job ID of the job to be migrated.
         */
        public long jobId;

        /**
         * < Please refer to \ref lsb_submit_options.
         */
        public int options;

        /**
         * < The number of hosts supplied as candidates  for migration.
         */
        public int numAskedHosts;

        /**
         * < An array of pointers to the names of candidate hosts for migration.
         */
        public Pointer askedHosts;
    }



/* structure for lsb_addjgrp() call */

    public static class jgrpAdd extends Structure {
        public static class ByReference extends jgrpAdd implements Structure.ByReference {}
        public static class ByValue extends jgrpAdd implements Structure.ByValue {}
        public jgrpAdd() {}
        public jgrpAdd(Pointer p) { super(p); read(); }

        public String groupSpec;
        public String timeEvent;
        public String depCond;
        public String sla;
        public int maxJLimit;
    }



/* structure for lsb_modjgrp() call */

    public static class jgrpMod extends Structure {
        public static class ByReference extends jgrpMod implements Structure.ByReference {}
        public static class ByValue extends jgrpMod implements Structure.ByValue {}
        public jgrpMod() {}
        public jgrpMod(Pointer p) { super(p); read(); }

        public String destSpec;
        public jgrpAdd jgrp;
    }



/* structure for lsb_addjgrp() and lsb_modjgrp() call reply */

    public static class jgrpReply extends Structure {
        public static class ByReference extends jgrpReply implements Structure.ByReference {}
        public static class ByValue extends jgrpReply implements Structure.ByValue {}
        public jgrpReply() {}
        public jgrpReply(Pointer p) { super(p); read(); }

        public String badJgrpName;
        public int num;
        public Pointer delJgrpList;
    }



    /**
     * \brief Signal a group of jobs.
     */
    public static class signalBulkJobs extends Structure {
        public static class ByReference extends signalBulkJobs implements Structure.ByReference {}
        public static class ByValue extends signalBulkJobs implements Structure.ByValue {}
        public signalBulkJobs() {}
        public signalBulkJobs(Pointer p) { super(p); read(); }


        /**
         * < Signal type
         */
        public int signal;

        /**
         * < Number of jobs
         */
        public int njobs;

        /**
         * < Jobids list
         */
        public Pointer jobs;

        /**
         * < Flags
         */
        public int flags;
    }



/* structure for lsb_ctrljgrp() call */

    public static class jgrpCtrl extends Structure {
        public static class ByReference extends jgrpCtrl implements Structure.ByReference {}
        public static class ByValue extends jgrpCtrl implements Structure.ByValue {}
        public jgrpCtrl() {}
        public jgrpCtrl(Pointer p) { super(p); read(); }

        public String groupSpec;
        public String userSpec;
        public int options;

/* JGRP_RELEASE, JGRP_HOLD, JGRP_DEL */
        public int ctrlOp;
    }




/* Indicate no change in chkpnt period for lsb_chkpntjob() */
    public static final int LSB_CHKPERIOD_NOCHNG = -1;

    /**
     *  \addtogroup chkpnt_job_option  chkpnt_job_option
     *  checkpoint job options()
     */

    /**
     * < Kill process if successfully chkpnted
     */
    public static final int LSB_CHKPNT_KILL = 0x1;

    /**
     * < Force chkpnt even if non-chkpntable conditions exist.
     */
    public static final int LSB_CHKPNT_FORCE = 0x2;

    /**
     * < Copy all regular files in use by the  checkpointed process to the checkpoint directory.
     */
    public static final int LSB_CHKPNT_COPY = 0x3;

    /**
     * < Chkpnt for the purpose of migration
     */
    public static final int LSB_CHKPNT_MIG = 0x4;

    /**
     * < Stop  process if successfully chkpnted
     */
    public static final int LSB_CHKPNT_STOP = 0x8;

    /**
     *  \addtogroup kill_requeue  kill_requeue
     *  kill and requeue a job options()
     */

    /**
     * < Kill then re-queue a job
     */
    public static final int LSB_KILL_REQUEUE = 0x10;

/* options for lsb_openjobinfo() */
    /**
     *  \addtogroup defs_lsb_openjobinfo  defs_lsb_openjobinfo
     *  Information options about job.
     */

    /**
     * < Reserved user name
     */
    public static final String ALL_USERS = "all";
    /**
     * \defgroup defs_lsb_openjobinfo_a defs_lsb_openjobinfo_a
     * defs_lsb_openjobinfo_a is part of defs_lsb_openjobinfo
     */
    public static final int ALL_JOB = 0x0001;
    /**
     * < Information about all jobs, including unfinished jobs (pending, running or suspended) and recently finished jobs. LSF remembers jobs finished within the preceding period. This period is set by the parameter CLEAN_PERIOD in the lsb.params file. The default is 3600 seconds (1 hour). (See lsb.params). The command line equivalent is bjobs -a./
     * <p/>
     * /**< Information about recently finished jobs.
     */
    public static final int DONE_JOB = 0x0002;

    /**
     * < Information about pending jobs.
     */
    public static final int PEND_JOB = 0x0004;

    /**
     * < Information about suspended jobs.
     */
    public static final int SUSP_JOB = 0x0008;

    /**
     * < Information about all unfinished jobs.
     */
    public static final int CUR_JOB = 0x0010;

    /**
     * < Information about the last submitted job.
     */
    public static final int LAST_JOB = 0x0020;

    /**
     * < Information about all running jobs
     */
    public static final int RUN_JOB = 0x0040;

    /**
     * < Information about JobId only.
     */
    public static final int JOBID_ONLY = 0x0080;

    /**
     * < Internal use only.
     */
    public static final int HOST_NAME = 0x0100;

    /**
     * < Exclude pending jobs.
     */
    public static final int NO_PEND_REASONS = 0x0200;

    /**
     * < Return group info structures
     */
    public static final int JGRP_INFO = 0x0400;

    /**
     * < Recursively search job group tree
     */
    public static final int JGRP_RECURSIVE = 0x0800;

    /**
     * < Return job array info structures
     */
    public static final int JGRP_ARRAY_INFO = 0x1000;

    /**
     * < All jobs in the core
     */
    public static final int JOBID_ONLY_ALL = 0x02000;

    /**
     * < All zombie jobs
     */
    public static final int ZOMBIE_JOB = 0x04000;

    /**
     * < Display remote jobs by their submission jobid.
     */
    public static final int TRANSPARENT_MC = 0x08000;

    /**
     * < Exceptional jobs
     */
    public static final int EXCEPT_JOB = 0x10000;

    /**
     * < Display for murex jobs
     */
    public static final int MUREX_JOB = 0x20000;


    /**
     * < To symphony UA
     */
    public static final int TO_SYM_UA = 0x40000;

    /**
     * < Only show top-level symphony job
     */
    public static final int SYM_TOP_LEVEL_ONLY = 0x80000;

    /**
     * < For internal use only
     */
    public static final int JGRP_NAME = 0x100000;

    /**
     * < Condensed host group
     */
    public static final int COND_HOSTNAME = 0x200000;

    /**
     * < Called from command, for internal use only
     */
    public static final int FROM_BJOBSCMD = 0x400000;

    /**
     * < -l in command parameter, for internal use only
     */
    public static final int WITH_LOPTION = 0x800000;

    /**
     * < Jobs submitted to aps queue
     */
    public static final int APS_JOB = 0x1000000;

    /**
     * < Information about user group.
     */
    public static final int UGRP_INFO = 0x2000000;
    /** RFC#1531: -G option support*/

    /**
     * < -WL
     */
    public static final int TIME_LEFT = 0x4000000;
    /**
     * < Estimated time remaining based on the runtime estimate or runlimit.
     */

/* -WF*/
    public static final int FINISH_TIME = 0x8000000;
    /**
     * < Estimated finish time based on the runtime estimate or runlimit.
     */

/* -WP*/
    public static final int COM_PERCENTAGE = 0x10000000;
    /**
     * < Estimated completion percentage based on the runtime estimate or runlimit. If options is 0, default to CUR_JOB.
     */

/* -ss option */
    public static final int SSCHED_JOB = 0x20000000;

/* -G option */
    public static final int KILL_JGRP_RECURSIVE = 0x40000000;

    /**
     *  \addtogroup group_nodetypes group_nodetypes
     *  define statements group node types.
     */

    /**
     * <  Job
     */
    public static final int JGRP_NODE_JOB = 1;

    /**
     * <  Group
     */
    public static final int JGRP_NODE_GROUP = 2;

    /**
     * <  Array
     */
    public static final int JGRP_NODE_ARRAY = 3;

    /**
     * <  SLA
     */
    public static final int JGRP_NODE_SLA = 4;

/* jobId macros */
    public static final long LSB_MAX_ARRAY_JOBID = 0x0FFFFFFFFL;
    public static final long LSB_MAX_ARRAY_IDX = 0x07FFFFFFFL;
    public static final int LSB_MAX_SEDJOB_RUNID = (0x0F);
    public static long LSB_JOBID (int array_jobId, int array_idx)    { return (((long)array_idx << 32) | array_jobId); }
    public static int LSB_ARRAY_IDX (long jobId)   { return (((jobId) == -1) ? (0) : (int)(((long)jobId >> 32)  & LSB_MAX_ARRAY_IDX)); }
    public static int LSB_ARRAY_JOBID (long jobId)  { return (((jobId) == -1) ? (-1) : (int)(jobId)); }
    //public static int LSB_ARRAY_JOBID (long jobId)  { return (((jobId) == -1) ? (-1) : (int)(jobId & LSB_MAX_ARRAY_JOBID)); }

/* Status of a job group */

    public static final int JGRP_INACTIVE = 0;
    public static final int JGRP_ACTIVE = 1;
    public static final int JGRP_UNDEFINED = -1;

    /**
     *  \addtogroup jobgroup_controltypes jobgroup_controltypes
     *  define statements job group control types.
     */


    /**
     * < bgrelease
     */
    public static final int JGRP_RELEASE = 1;

    /**
     * < bghold
     */
    public static final int JGRP_HOLD = 2;

    /**
     * < bgdel
     */
    public static final int JGRP_DEL = 3;

    /**
     *  \addtogroup jobgroup_counterIndex jobgroup_counterIndex
     *   Following can be used to index  into 'counters' array.
     */

    /**
     * < Total jobs in the array
     */
    public static final int JGRP_COUNT_NJOBS = 0;

    /**
     * < Number of pending jobs in the array
     */
    public static final int JGRP_COUNT_PEND = 1;

    /**
     * < Number of held jobs in the array
     */
    public static final int JGRP_COUNT_NPSUSP = 2;

    /**
     * < Number of running jobs in the array
     */
    public static final int JGRP_COUNT_NRUN = 3;

    /**
     * < Number of jobs suspended by the system in the array
     */
    public static final int JGRP_COUNT_NSSUSP = 4;

    /**
     * < Number of jobs suspended by the user in the array
     */
    public static final int JGRP_COUNT_NUSUSP = 5;

    /**
     * < Number of exited jobs in the array
     */
    public static final int JGRP_COUNT_NEXIT = 6;

    /**
     * < Number of successfully completed jobs
     */
    public static final int JGRP_COUNT_NDONE = 7;

    /**
     * < Total slots in the array
     */
    public static final int JGRP_COUNT_NJOBS_SLOTS = 8;

    /**
     * < Number of pending slots in the array
     */
    public static final int JGRP_COUNT_PEND_SLOTS = 9;

    /**
     * < Number of running slots in the array
     */
    public static final int JGRP_COUNT_RUN_SLOTS = 10;

    /**
     * < Number of slots suspended by the system in the array
     */
    public static final int JGRP_COUNT_SSUSP_SLOTS = 11;

    /**
     * < Number of slots suspended by the user in the array
     */
    public static final int JGRP_COUNT_USUSP_SLOTS = 12;

    /**
     * < Number of reserverd slots in the array
     */
    public static final int JGRP_COUNT_RESV_SLOTS = 13;

/* job group modification types */
    public static final int JGRP_MOD_LIMIT = 0x1;

/*the number of counters of job group
* based on job level
*/
    public static final int NUM_JGRP_JOB_COUNTERS = 8;
/* the number of all counters of job group,
* including job level and slot level
*/
/* {njobs, npend, npsusp, nrun, nssusp nususp, nexit, ndone} */
    public static final int NUM_JGRP_COUNTERS = 14;

/* job group is created explicitly */
    public static final int JGRP_CREATE_EXP = 0x01;

/* job group is created implicitly */
    public static final int JGRP_CREATE_IMP = 0x02;
/* The LSF job group.
 */

    public static class jgrp extends Structure {
        public static class ByReference extends jgrp implements Structure.ByReference {}
        public static class ByValue extends jgrp implements Structure.ByValue {}
        public jgrp() {}
        public jgrp(Pointer p) { super(p); read(); }

        public String name;
        public String path;
        public String user;
        public String sla;
        public int[] counters = new int[NUM_JGRP_COUNTERS];
        public int maxJLimit;
    }



/* Structure for lsb_setjobattr() call */

    public static class jobAttrInfoEnt extends Structure {
        public static class ByReference extends jobAttrInfoEnt implements Structure.ByReference {}
        public static class ByValue extends jobAttrInfoEnt implements Structure.ByValue {}
        public jobAttrInfoEnt() {}
        public jobAttrInfoEnt(Pointer p) { super(p); read(); }


/* id of the job */
        public long jobId;

/* port number of the job */
        public short port;

/* first executing host of the job */
        public byte[] hostname = new byte[LibLsf.MAXHOSTNAMELEN];
    }



    /**
     * \brief  job attribute setting log.
     */
    public static class jobAttrSetLog extends Structure {
        public static class ByReference extends jobAttrSetLog implements Structure.ByReference {}
        public static class ByValue extends jobAttrSetLog implements Structure.ByValue {}
        public jobAttrSetLog() {}
        public jobAttrSetLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The user who requested the action
         */
        public int uid;

        /**
         * < Job attributes
         */
        public int port;

        /**
         * < Name of the host
         */
        public String hostname;
    }



    /**
     * \brief  job information head.
     */
    public static class jobInfoHead extends Structure {
        public static class ByReference extends jobInfoHead implements Structure.ByReference {}
        public static class ByValue extends jobInfoHead implements Structure.ByValue {}
        public jobInfoHead() {}
        public jobInfoHead(Pointer p) { super(p); read(); }


        /**
         * < The number of jobs in the connection
         */
        public int numJobs;

        /**
         * < An array of job identification numbers in the conection
         */
        public NativeLongByReference jobIds;

        /**
         * < The number of hosts in the connection
         */
        public int numHosts;

        /**
         * < An array of host names in the connection
         */
        public Pointer hostNames;

        /**
         * < The number of clusters in the connection
         */
        public int numClusters;

        /**
         * < An array of cluster names in the connection
         */
        public Pointer clusterNames;

        /**
         * < The number of remoteHosts in the connection
         */
        public IntByReference numRemoteHosts;

        /**
         * < An array of remoteHost names in the connection
         */
        public PointerByReference remoteHosts;
    }



    /**
     * \brief job Information head extent
     */
    public static class jobInfoHeadExt extends Structure {
        public static class ByReference extends jobInfoHeadExt implements Structure.ByReference {}
        public static class ByValue extends jobInfoHeadExt implements Structure.ByValue {}
        public jobInfoHeadExt() {}
        public jobInfoHeadExt(Pointer p) { super(p); read(); }


        /**
         * <  Job Information header
         */
        public jobInfoHead.ByReference jobInfoHead;

        /**
         * <  Group Information returned
         */
        public Pointer groupInfo;
    }



    /**
     * \brief structure reserveItem
     */
    public static class reserveItem extends Structure {
        public static class ByReference extends reserveItem implements Structure.ByReference {}
        public static class ByValue extends reserveItem implements Structure.ByValue {}
        public reserveItem() {}
        public reserveItem(Pointer p) { super(p); read(); }


        /**
         * < Name of the resource to reserve.
         */
        public String resName;

        /**
         * < The number of hosts to reserve this resource.
         */
        public int nHost;

        /**
         * < Amount of reservation is made on each host. Some hosts may reserve 0.
         */
        public FloatByReference value;

        /**
         * < Flag of shared or host-base resource
         */
        public int shared;
    }



    /**
     * \brief  job information entry.
     */
    public static class jobInfoEnt extends Structure {
        public static class ByReference extends jobInfoEnt implements Structure.ByReference {}
        public static class ByValue extends jobInfoEnt implements Structure.ByValue {}
        public jobInfoEnt() {}
        public jobInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < The job ID that the LSF system assigned to the job.
         */
        public long jobId;

        /**
         * < The name of the user who submitted the job.
         */
        public String user;

        /**
         * < The current status of the job.Possible values areshown in job_states
         */
        public int status;

        /**
         * < Pending or suspending reasons of the job
         */
        public IntByReference reasonTb;

        /**
         * < Length of reasonTb[]
         */
        public int numReasons;

        /**
         * < The reason a job is pending or suspended.
         */
        public int reasons;

        /**
         * < The reason a job is pending or suspended. If status is JOB_STAT_PEND, the values of reasons and subreasons are explained by \ref lsb_pendreason. If status is JOB_STAT_PSUSP, the values of reasons and subreasons are explained by \ref lsb_suspreason.   When reasons is PEND_HOST_LOAD or SUSP_LOAD_REASON,  subreasons indicates the load indices that are out of bounds. If reasons is PEND_HOST_LOAD, subreasons is the same as busySched in the hostInfoEnt structure; if reasons is SUSP_LOAD_REASON, subreasons is the same as busyStop in the hostInfoEnt structure. (See \ref lsb_hostinfo)
         */
        public int subreasons;

        /**
         * < The job process ID.
         */
        public int jobPid;

        /**
         * < The time the job was submitted, in seconds since 00:00:00 GMT, Jan. 1, 1970.
         */
        public NativeLong submitTime;

        /**
         * < Time when job slots are reserved
         */
        public NativeLong reserveTime;

        /**
         * < The time that the job started running, if it has been dispatched.
         */
        public NativeLong startTime;

        /**
         * < Job's predicted start time
         */
        public NativeLong predictedStartTime;

        /**
         * < The termination time of the job, if it has completed.
         */
        public NativeLong endTime;

        /**
         * < Last time event
         */
        public NativeLong lastEvent;

        /**
         * < Next time event
         */
        public NativeLong nextEvent;

        /**
         * < Duration time (minutes)
         */
        public int duration;

        /**
         * < CPU time consumed by the job
         */
        public float cpuTime;

        /**
         * < The file creation mask when the job was submitted.
         */
        public int umask;

        /**
         * < The current working directory when the job was submitted.
         */
        public String cwd;

        /**
         * < Home directory on submission host.
         */
        public String subHomeDir;

        /**
         * < The name of the host from which the job was  submitted.
         */
        public String fromHost;

        /**
         * < The array of names of hosts on which the job executes.
         */
        public Pointer exHosts;

        /**
         * < The number of hosts on which the job executes.
         */
        public int numExHosts;

        /**
         * < The CPU factor for normalizing CPU and wall clock time limits.
         */
        public float cpuFactor;

        /**
         * < The number of load indices in the loadSched and loadStop arrays.
         */
        public int nIdx;

        /**
         * < The values in the loadSched array specify the thresholds for the corresponding load indices. Only if the current values of all specified load indices of a host are within (below or above,  depending on the meaning of the load index) their corresponding thresholds may the suspended job be resumed on this host.  For an explanation of the entries in the loadSched, see \ref lsb_hostinfo.
         */
        public FloatByReference loadSched;

        /**
         * < The values in the loadStop array specify the thresholds for job suspension; if any of the current load index values of the host crosses its threshold, the job will be suspended.  For an explanation of the entries in the loadStop, see \ref lsb_hostinfo.
         */
        public FloatByReference loadStop;

        /**
         * < Structure for \ref lsb_submit call.
         */
        public submit submit;

        /**
         * < Job exit status.
         */
        public int exitStatus;

        /**
         * < Mapped UNIX user ID on the execution host.
         */
        public int execUid;

        /**
         * < Home directory for the job on the execution host.
         */
        public String execHome;

        /**
         * < Current working directory for the job on the execution host.
         */
        public String execCwd;

        /**
         * < Mapped user name on the execution host.
         */
        public String execUsername;

        /**
         * < Time of the last job resource usage update.
         */
        public NativeLong jRusageUpdateTime;

        /**
         * < Contains resource usage information for the job.
         */
        public LibLsf.jRusage runRusage;

        /**
         * < Job type.N_JOB, N_GROUP, N_HEAD
         */
        public int jType;

        /**
         * < The parent job group of a job or job group.
         */
        public String parentGroup;

        /**
         * < If jType is JGRP_NODE_GROUP, then it is the job group name. Otherwise, it is thejob name.
         */
        public String jName;

        /**
         * < Index into the counter array, only used for job arrays. Possible index values are shown in \ref jobgroup_counterIndex
         */
        public int[] counter = new int[NUM_JGRP_COUNTERS];

        /**
         * < Service port of the job.
         */
        public short port;

        /**
         * < Job dynamic priority
         */
        public int jobPriority;

        /**
         * < The number of external messages in the job.
         */
        public int numExternalMsg;

        /**
         * < This structure contains the information required to define an external message reply.
         */
        public Pointer externalMsg;

        /**
         * < MultiCluster cluster ID. If clusterId is greater than or equal to 0, the job is a pending remote job, and \ref lsb_readjobinfo checks for host_name\@cluster_name. If host name is needed, it should be found in  jInfoH->remoteHosts. If the remote host name is not available, the constant string remoteHost is used.
         */
        public int clusterId;

        /**
         * <  Detail reason field
         */
        public String detailReason;

        /**
         * < Idle factor for job exception handling. If the job idle factor is less than the specified threshold, LSF invokes LSF_SERVERDIR/eadmin to trigger the action for a job idle exception.
         */
        public float idleFactor;

        /**
         * < Job exception handling mask
         */
        public int exceptMask;


        /**
         * < Placement information of LSF HPC jobs.Placement information of LSF HPC jobs.Arbitrary information of a job stored as a string currently used by rms_rid  and rms_alloc
         */
        public String additionalInfo;

        /**
         * < Job termination reason. See lsbatch.h.
         */
        public int exitInfo;

        /**
         * < Job warning time period in seconds; -1 if unspecified.
         */
        public int warningTimePeriod;

        /**
         * < Warning action, SIGNAL | CHKPNT |command, null if unspecified
         */
        public String warningAction;

        /**
         * < SAAP charged for job
         */
        public String chargedSAAP;

        /**
         * < The rusage satisfied at job runtime
         */
        public String execRusage;

        /**
         * < The time when advance reservation expired or was deleted.
         */
        public NativeLong rsvInActive;

        /**
         * < The number of licenses reported from License Scheduler.
         */
        public int numLicense;

        /**
         * < License Scheduler license names.
         */
        public Pointer licenseNames;

        /**
         * < Absolute priority scheduling (APS) priority value.
         */
        public float aps;

        /**
         * < Absolute priority scheduling (APS) string set by administrators to denote static system APS value
         */
        public float adminAps;

        /**
         * < The real runtime on the execution host.
         */
        public int runTime;

        /**
         * < How many kinds of resource are reserved by this job
         */
        public int reserveCnt;

        /**
         * < Detail reservation information for each kind of resource
         */
        public Pointer /* reserveItem.ByReference */ items;

        /**
         * < Absolute priority scheduling (APS) string set by administrators to denote ADMIN factor APS value.
         */
        public float adminFactorVal;

        /**
         * < Pending resize min. 0, if no resize pending.
         */
        public int resizeMin;

        /**
         * < Pending resize max. 0, if no resize pending
         */
        public int resizeMax;

        /**
         * < Time when pending request was issued
         */
        public NativeLong resizeReqTime;

        /**
         * < Number of hosts when job starts
         */
        public int jStartNumExHosts;

        /**
         * < Host list when job starts
         */
        public Pointer jStartExHosts;

        /**
         * < Last time when job allocation changed
         */
        public NativeLong lastResizeTime;
    }


/* the bit set for jobInfoEnt->exceptMask */
    public static final int J_EXCEPT_OVERRUN = 0x02;
    public static final int J_EXCEPT_UNDERUN = 0x04;
    public static final int J_EXCEPT_IDLE = 0x80;
    public static final int J_EXCEPT_RUNTIME_EST_EXCEEDED = 0x100;

/* exception showed by bjobs -l and bacct -l*/
    public static final String OVERRUN = "overrun";
    public static final String UNDERRUN = "underrun";
    public static final String IDLE = "idle";
    public static final String SPACE = "  ";
    public static final String RUNTIME_EST_EXCEEDED = "runtime_est_exceeded";

/* LSF7.0 moved jobInfoReq structure definition from
*  daemonout.h to lsbatch.h. This structure will work
*  with new API \ref lsb_openjobinfo_req
 */

    /**
     * \brief  job Information Request
     */
    public static class jobInfoReq extends Structure {
        public static class ByReference extends jobInfoReq implements Structure.ByReference {}
        public static class ByValue extends jobInfoReq implements Structure.ByValue {}
        public jobInfoReq() {}
        public jobInfoReq(Pointer p) { super(p); read(); }


        /**
         * < Options defined in \ref defs_lsb_openjobinfo
         */
        public int options;

        /**
         * < Name of user whose jobs to be checked
         */
        public String userName;

        /**
         * < Job id, 0 means all jobs
         */
        public long jobId;

        /**
         * < Job name
         */
        public String jobName;

        /**
         * < Queue name
         */
        public String queue;

        /**
         * < Check jobs running on this host
         */
        public String host;

        /**
         * < Job application
         */
        public String app;

        /**
         * < Job description
         */
        public String jobDescription;

        /**
         * < For new options in future
         */
        public submit_ext.ByReference submitExt;
    }



    /**
     * \brief  user information entry.
     */
    public static class userInfoEnt extends Structure {
        public static class ByReference extends userInfoEnt implements Structure.ByReference {}
        public static class ByValue extends userInfoEnt implements Structure.ByValue {}
        public userInfoEnt() {}
        public userInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < Name of the user or user group
         */
        public String user;

        /**
         * < The maximum number of job slots the user or user group can use on each processor. The job slots can be used by started jobs or reserved for PEND jobs.
         */
        public float procJobLimit;

        /**
         * < The maximum number of job slots that the user or user group can use simultaneously in the local LSF cluster. The job slots can be used by started jobs or reserved for PEND jobs.
         */
        public int maxJobs;

        /**
         * < The current number of job slots used by running and suspended jobs belonging to the user or user group.
         */
        public int numStartJobs;

        /**
         * < The total number of job slots in the LSF cluster for the jobs submitted by the user or user group.
         */
        public int numJobs;

        /**
         * < The number of job slots the user or user group has for pending jobs.
         */
        public int numPEND;

        /**
         * < The number of job slots the user or user group has for running jobs.
         */
        public int numRUN;

        /**
         * < The number of job slots for the jobs belonging to the user or user group that have been suspended by the system.
         */
        public int numSSUSP;

        /**
         * < The number of job slots for the jobs belonging to the user or user group that have been suspended by the user or the LSF system administrator.
         */
        public int numUSUSP;

        /**
         * < The number of job slots reserved for the pending jobs belonging to the user or user group.
         */
        public int numRESERVE;

        /**
         * < The maximum number of pending jobs allowed.
         */
        public int maxPendJobs;
    }



/* UserEquivalent info */

    public static class userEquivalentInfoEnt extends Structure {
        public static class ByReference extends userEquivalentInfoEnt implements Structure.ByReference {}
        public static class ByValue extends userEquivalentInfoEnt implements Structure.ByValue {}
        public userEquivalentInfoEnt() {}
        public userEquivalentInfoEnt(Pointer p) { super(p); read(); }

        public String equivalentUsers;
    }



/* UserMapping info */

    public static class userMappingInfoEnt extends Structure {
        public static class ByReference extends userMappingInfoEnt implements Structure.ByReference {}
        public static class ByValue extends userMappingInfoEnt implements Structure.ByValue {}
        public userMappingInfoEnt() {}
        public userMappingInfoEnt(Pointer p) { super(p); read(); }


/* Users in the local cluster */
        public String localUsers;

/* Users in remote clusters */
        public String remoteUsers;

/* "export" or "import" */
        public String direction;
    }




/* APS structures used for mapping between factors */

    /**
     * \brief  APS structures used for mapping between factors
     */
    public static class apsFactorMap extends Structure {
        public static class ByReference extends apsFactorMap implements Structure.ByReference {}
        public static class ByValue extends apsFactorMap implements Structure.ByValue {}
        public apsFactorMap() {}
        public apsFactorMap(Pointer p) { super(p); read(); }


        /**
         * < Name of factor.
         */
        public String factorName;

        /**
         * < SubFactor names.
         */
        public String subFactorNames;
    }



    /**
     * \brief  APS structures used for mapping between factors
     */
    public static class apsLongNameMap extends Structure {
        public static class ByReference extends apsLongNameMap implements Structure.ByReference {}
        public static class ByValue extends apsLongNameMap implements Structure.ByValue {}
        public apsLongNameMap() {}
        public apsLongNameMap(Pointer p) { super(p); read(); }


        /**
         * < Short name
         */
        public String shortName;

        /**
         * < Long name
         */
        public String longName;
    }




/* options for lsb_queueinfo() , some values should not
*  conflict with the option values for lsb_usergrpinfo() and lsb_hostinfo_ex()
*  since they share the same xdr_infoReq()
*/

/* for compatibility for 2.0 */
    public static final int ALL_QUEUE = 0x01;

/* for compatibility for 2.0 */
    public static final int DFT_QUEUE = 0x02;
    public static final int CHECK_HOST = 0x80;
    public static final int CHECK_USER = 0x100;
    public static final int SORT_HOST = 0x200;

/* not bqueues -l or -r */
    public static final int QUEUE_SHORT_FORMAT = 0x400;
/* expand hostname into official hostname in lsb_queueinfo */
    public static final int EXPAND_HOSTNAME = 0x800;

/* only retrieve batch partitions */
    public static final int RETRIEVE_BATCH = 0x1000;

/* Signal number in each version LSB_SIG_NUM must be equal to
*  signal number in the latest version.
 */
    public static final int LSB_SIG_NUM_40 = 25;
    public static final int LSB_SIG_NUM_41 = 26;

/* Solutions #38347 */
    public static final int LSB_SIG_NUM_51 = 30;
    public static final int LSB_SIG_NUM_60 = 30;
    public static final int LSB_SIG_NUM = 30;

/* Dynamic CPU provision
*  to indicate whether a SP can lend or borrow hosts
 */
    public static final int DCP_LEND_HOSTS = 0x0001;
    public static final int DCP_BORROW_HOSTS = 0x0002;

/* status to indicate the current situation of Dynamic CPU provision
*  DCP_UNDER_ALLOC_AND_STARVING means a partition is under allocation
*  of dynamic cpu and its pending jobs are starving for more cpus.
 */
    public static final int DCP_ALLOC_CPU_OK = 0x0;
    public static final int DCP_UNDER_ALLOC_CPU = 0x0001;
    public static final int DCP_JOB_WAIT_FOR_CPU = 0x0002;
    public static final int DCP_ALLOC_CPU_BUSY = 0x0004;

/* Structure for lsb_queueinfo() call */
/* !!! IMPORTANT !!!
*  If you change queueInfoEnt, you have to change Intlib/ade.lsbatch.h too!
 */

    /**
     * queueInfoEnt  queue information entry.
     */
    public static class queueInfoEnt extends Structure {
        public static class ByReference extends queueInfoEnt implements Structure.ByReference {}
        public static class ByValue extends queueInfoEnt implements Structure.ByValue {}
        public queueInfoEnt() {}
        public queueInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < The name of the queue.
         */
        public String queue;

        /**
         * < Describes the typical use of the queue.
         */
        public String description;

        /**
         * < Defines the priority of the queue. This determines the order in which the job queues are searched at job dispatch time: queues with higher priority values are searched first. (This is contrary to UNIX process priority ordering.)
         */
        public int priority;

        /**
         * < Defines the nice value at which jobs in this queue will be run.
         */
        public short nice;

        /**
         * < A blank-separated list of names of users allowed to submit jobs to this queue.
         */
        public String userList;

        /**
         * < A blank-separated list of names of hosts to which jobs in this queue may be dispatched.
         */
        public String hostList;

        /**
         * < Original HOSTS string in case "-" is used.
         */
        public String hostStr;

        /**
         * < The number of load indices in the loadSched and loadStop arrays.
         */
        public int nIdx;

        /**
         * < The queue and host loadSched and loadStop arrays control batch job dispatch, suspension, and resumption. The values in the loadSched array specify thresholds for the corresponding load indices. Only if the current values of all specified load indices of a host are within (below or above, depending on the meaning of the load index) the corresponding thresholds of this queue, will jobs in this queue be dispatched to the host. The same conditions are used to resume jobs dispatched from this queue that have been suspended on the host.
         */
        public FloatByReference loadSched;

        /**
         * < The values in the loadStop array specify the thresholds for job suspension. If any of the current load index values of a host goes beyond a queue's threshold, jobs from the queue will be suspended. For an explanation of the fields in the loadSched and loadStop arrays, see \ref lsb_hostinfo.
         */
        public FloatByReference loadStop;

        /**
         * < Per-user limit on the number of jobs that can be dispatched from this queue and executed concurrently.
         */
        public int userJobLimit;

        /**
         * < Per-processor limit on the number of jobs that can be dispatched from this queue and executed concurrently.
         */
        public float procJobLimit;

        /**
         * < A blank-separated list of time windows describing the run window of the queue. When a queue's run window is closed, no job from this queue will be dispatched. When the run window closes, any running jobs from this queue will be suspended until the run window reopens, when they will be resumed. The default is no restriction, or always open (i.e., 24 hours a day, seven days a week). A time window has the format begin_time-end_time. Time is specified in the format [day:]hour[:minute], where all fields are numbers in their respective legal ranges: 0(Sunday)-6 for day, 0-23 for hour, and 0-59 for minute. The default value for minute is 0 (on the hour); the default value for day is every day of the week. The begin_time and end_time of a window are separated by `-', with no white space (i.e., blank or TAB) in between. Both begin_time and end_time must be present for a time window. Note that this run window only applies to batch jobs;interactive jobs scheduled by the LSF Load Information Manager (LIM) are controlled by another set of run windows.
         */
        public String windows;

        /**
         * < The per-process UNIX hard resource limits for all jobs submitted to this queue (see getrlimit() and lsb.queues). The default values for the resource limits are unlimited, indicated by -1. The constants used to index the rLimits array and the corresponding resource limits are listed below. <br> LSF_RLIMIT_CPU (CPULIMIT) <br> LSF_RLIMIT_FSIZE (FILELIMIT) <br> LSF_RLIMIT_DATA (DATALIMIT) <br> LSF_RLIMIT_STACK    (STACKLIMIT) <br> LSF_RLIMIT_CORE     (CORELIMIT) <br> LSF_RLIMIT_RSS      (MEMLIMIT) <br> LSF_RLIMIT_RUN      (RUNLIMIT) <br> LSF_RLIMIT_PROCESS     (PROCESSLIMIT) <br> LSF_RLIMIT_SWAP     (SWAPLIMIT) <br> LSF_RLIMIT_THREAD <br> LSF_RLIMIT_NOFILE <br> LSF_RLIMIT_OPENMAX <br> LSF_RLIMIT_VMEM
         */
        public int[] rLimits = new int[LibLsf.LSF_RLIM_NLIMITS];

        /**
         * < A host name or host model name. If the queue CPULIMIT or RUNLIMIT gives a host specification, hostSpec will be that specification. Otherwise, if defaultHostSpec (see below) is not null, hostSpec will be defaultHostSpec. Otherwise, if DEFAULT_HOST_SPEC is defined in the lsb.params file, (see lsb.params), hostSpec will be this value. Otherwise, hostSpec will be the name of the host with the largest CPU factor in the cluster.
         */
        public String hostSpec;

        /**
         * < The attributes of the queue.
         */
        public int qAttrib;

        /**
         * < The status of the queue.
         */
        public int qStatus;

        /**
         * < The maximum number of jobs dispatched by the queue and not yet finished.
         */
        public int maxJobs;

        /**
         * < Number of jobs in the queue, including pending, running, and suspended jobs.
         */
        public int numJobs;

        /**
         * < Number of pending jobs in the queue.
         */
        public int numPEND;

        /**
         * < Number of running jobs in the queue.
         */
        public int numRUN;

        /**
         * < Number of system suspended jobs in the queue.
         */
        public int numSSUSP;

        /**
         * < Number of user suspended jobs in the queue.
         */
        public int numUSUSP;

        /**
         * < The queue migration threshold in minutes.
         */
        public int mig;

        /**
         * < The number of seconds that a new job waits, before being scheduled. A value of zero (0) means the job is scheduled without any delay.
         */
        public int schedDelay;

        /**
         * < The number of seconds for a host to wait after dispatching a job to a host, before accepting a second job to dispatch to the same host.
         */
        public int acceptIntvl;

        /**
         * < A blank-separated list of time windows describing the dispatch window of the queue. When a queue's dispatch window is closed, no job from this queue will be dispatched.The default is no restriction, or always open (i.e., 24 hours a day, seven days a week). For the time window format, see windows (above).
         */
        public String windowsD;

        /**
         * < A blank-separated list of queue specifiers. Each queue specifier is of the form queue\@host where host is an NQS host name and queue is the name of a queue on that host.
         */
        public String nqsQueues;

        /**
         * < A blank-separated list of user shares. Each share is of the form [user, share] where user is a user name, a user group name, the reserved word default or the reserved word others, and share is the number of shares the user gets.
         */
        public String userShares;

        /**
         * < The value of DEFAULT_HOST_SPEC in the Queue section for this queue in the lsb.queues file.
         */
        public String defaultHostSpec;

        /**
         * < An LSF resource limit used to limit the number of job slots (processors) a (parallel) job in the queue will use. A job submitted to this queue must specify a number of processors not greater than this limit.
         */
        public int procLimit;

        /**
         * < A list of administrators of the queue. The users whose names are here are allowed to operate on the jobs in the queue and on the queue itself.
         */
        public String admins;

        /**
         * < Queue's pre-exec command. The command is executed before the real batch job is run on the execution host (or on the first host selected for a parallel batch job).
         */
        public String preCmd;

        /**
         * < Queue's post-exec command. The command is run when a job terminates.
         */
        public String postCmd;

        /**
         * < Jobs that exit with these values are automatically requeued.
         */
        public String requeueEValues;

        /**
         * < The maximum number of job slots a host can process from this queue, including job slots of dispatched jobs which have not finished yet and reserved slots for some PEND jobs. This limit controls the number of jobs sent to each host, regardless of a uniprocessor host or multiprocessor host. Default value for this limit is infinity.
         */
        public int hostJobLimit;

        /**
         * < Resource requirement string used to determine eligible hosts for a job.
         */
        public String resReq;

        /**
         * < Number of reserved job slots for pending jobs.
         */
        public int numRESERVE;

        /**
         * < The time used to hold the reserved job slots for a PEND job in this queue.
         */
        public int slotHoldTime;

        /**
         * < Remote MultiCluster send-jobs queues to forward jobs to.
         */
        public String sndJobsTo;

        /**
         * < Remote MultiCluster receive-jobs queues that can forward to this queue.
         */
        public String rcvJobsFrom;

        /**
         * < Resume threshold conditions for a suspended job in this queue.
         */
        public String resumeCond;

        /**
         * < Stop threshold conditions for a running job in this queue.
         */
        public String stopCond;

        /**
         * < Job starter command for a running job in this queue
         */
        public String jobStarter;

        /**
         * < Command configured for the SUSPEND action.
         */
        public String suspendActCmd;

        /**
         * < Command configured for the RESUME action.
         */
        public String resumeActCmd;

        /**
         * < Command configured for the TERMINATE action.
         */
        public String terminateActCmd;

        /**
         * < Configurable signal mapping
         */
        public int[] sigMap = new int[LSB_SIG_NUM];

        /**
         * < Preemptive scheduling and preemption policy specified for the queue.
         */
        public String preemption;

        /**
         * < Time period for a remote cluster to schedule a job. MultiCluster job forwarding model only. Determines how long a MultiCluster job stays pending in the execution cluster before returning to the submission cluster. The remote timeout limit in seconds is: \li MAX_RSCHED_TIME.ByReference  MBD_SLEEP_TIME=timeout
         */
        public int maxRschedTime;


        /**
         * < Number of share accounts in the queue.
         */
        public int numOfSAccts;

        /**
         * < (Only used for queues with fairshare policy) a share account vector capturing the fairshare information of the users using the queue. The storage for the array of queueInfoEnt structures will be reused by the next call.
         */
        public Pointer /* shareAcctInfoEnt.ByReference */ shareAccts;

        /**
         * < The directory where the checkpoint files are created.
         */
        public String chkpntDir;

        /**
         * < The checkpoint period in minutes.
         */
        public int chkpntPeriod;

        /**
         * < MultiCluster job forwarding model only. Specifies the MultiCluster pending job limit for a receive-jobs queue. This represents the maximum number of MultiCluster import jobs that can be pending in the queue; once the limit has been reached, the queue stops accepting jobs from remote clusters.
         */
        public int imptJobBklg;

        /**
         * < The default (soft) resource limits for all jobs submitted to this queue (see getrlimit() and lsb.queues).
         */
        public int[] defLimits = new int[LibLsf.LSF_RLIM_NLIMITS];

        /**
         * < The maximum number of jobs allowed to be dispatched together in one job chunk. Must be a positive integer greater than 1.
         */
        public int chunkJobSize;

        /**
         * < The minimum number of job slots (processors) that a job in the queue will use.
         */
        public int minProcLimit;

        /**
         * < The default (soft) limit on the number of job slots (processors) that a job in the queue will use.
         */
        public int defProcLimit;

        /**
         * < The list of queues for cross-queue fairshare.
         */
        public String fairshareQueues;

        /**
         * < Default external scheduling for the queue.
         */
        public String defExtSched;

        /**
         * < Mandatory external scheduling options for the queue.
         */
        public String mandExtSched;

        /**
         * < Share of job slots for queue-based fairshare. Represents the percentage of running jobs (job slots) in use from the queue. SLOT_SHARE must be greater than zero (0) and less than or equal to 100. The sum of SLOT_SHARE for all queues in the pool does not need to be 100%. It can be more or less, depending on your needs.
         */
        public int slotShare;

        /**
         * < Name of the pool of job slots the queue belongs to for queue-based fairshare. A queue can only belong to one pool. All queues in the pool must share the same set of hosts. Specify any ASCII string up to 60 chars long. You can use letters, digits, underscores (_) or dashes (-). You cannot use blank spaces.
         */
        public String slotPool;

        /**
         * < Specifies a threshold for job underrun exception handling. If a job exits before the specified number of minutes, LSF invokes LSF_SERVERDIR/eadmin to trigger the action for a job underrun exception.
         */
        public int underRCond;

        /**
         * < Specifies a threshold for job overrun exception handling. If a job runs longer than the specified run time, LSF invokes LSF_SERVERDIR/eadmin to trigger the action for a job overrun exception.
         */
        public int overRCond;

        /**
         * < Specifies a threshold for idle job exception handling. The value should be a number between 0.0 and 1.0 representing CPU time/runtime. If the job idle factor is less than the specified threshold, LSF invokes LSF_SERVERDIR/eadmin to trigger the action for a job idle exception.
         */
        public float idleCond;

        /**
         * < The number of underrun jobs in the queue.
         */
        public int underRJobs;

        /**
         * < The number of overrun jobs in the queue.
         */
        public int overRJobs;

        /**
         * < The number of idle jobs in the queue.
         */
        public int idleJobs;

        /**
         * < Specifies the amount of time before a job control action occurs that a job warning action is to be taken. For example, 2 minutes before the job reaches run time limit or termination deadline, or the queue's run window is closed, an URG signal is sent to the job. Job action warning time is not normalized. A job action warning time must be specified with a job warning action in order for job warning to take effect.
         */
        public int warningTimePeriod;

        /**
         * < Specifies the job action to be taken before a job control action occurs. For example, 2 minutes before the job reaches run time limit or termination deadline, or the queue's run window is closed, an URG signal is sent to the job. A job warning action must be specified with a job action warning time in order for job warning to take effect. If specified, LSF sends the warning action to the job before the actual control action is taken. This allows the job time to save its result before being terminated by the job control action. You can specify actions similar to the JOB_CONTROLS queue level parameter: send a signal, invoke a command, or checkpoint the job.
         */
        public String warningAction;

        /**
         * < AdminAction - queue control message
         */
        public String qCtrlMsg;

        /**
         * < Acept resource request.
         */
        public String acResReq;

        /**
         * < Limit of running session scheduler jobs.
         */
        public int symJobLimit;

        /**
         * < cpu_req for service partition of session scheduler
         */
        public String cpuReq;

        /**
         * < Indicate whether it would be willing to donate/borrow.
         */
        public int proAttr;

        /**
         * < The maximum number of hosts to lend.
         */
        public int lendLimit;

        /**
         * < The grace period to lend/return idle hosts.
         */
        public int hostReallocInterval;

        /**
         * < Number of CPUs required by CPU provision.
         */
        public int numCPURequired;

        /**
         * < Number of CPUs actually allocated.
         */
        public int numCPUAllocated;

        /**
         * < Number of CPUs borrowed.
         */
        public int numCPUBorrowed;

        /**
         * < Number of CPUs lent.
         */
        public int numCPULent;
        /* the number of reserved cpu(numCPUReserved) = numCPUAllocated - numCPUBorrowed + numCPULent */


        /* the following fields are for real-time app(ex. murex) of symphony */

        /**
         * < Scheduling granularity. in milliseconds.
         */
        public int schGranularity;

        /**
         * < The grace period for stopping session scheduler tasks.
         */
        public int symTaskGracePeriod;

        /**
         * < Minimum number of SSMs.
         */
        public int minOfSsm;

        /**
         * < Maximum number of SSMs.
         */
        public int maxOfSsm;

        /**
         * < Number of allocated slots.
         */
        public int numOfAllocSlots;

        /**
         * < Service preemptin policy.
         */
        public String servicePreemption;


        /**
         * < Dynamic cpu provision status.
         */
        public int provisionStatus;

        /**
         * < The minimum time for preemption and backfill, in seconds.
         */
        public int minTimeSlice;

        /**
         * < List of queues defined in a queue group for absolute priority scheduling (APS) across multiple queues.
         */
        public String queueGroup;

        /**
         * < The number of calculation factors for absolute priority scheduling (APS).
         */
        public int numApsFactors;

        /**
         * < List of calculation factors for absolute priority scheduling (APS)
         */
        public Pointer /* apsFactorInfo.ByReference */ apsFactorInfoList;

        /**
         * < The mapping of factors to subfactors for absolute priority scheduling (APS).
         */
        public Pointer /* apsFactorMap.ByReference */ apsFactorMaps;

        /**
         * < The mapping of factors to their long names for absolute priority scheduling (APS).
         */
        public Pointer /* apsLongNameMap.ByReference */ apsLongNames;

        /**
         * < Maximum number of job preempted times.
         */
        public int maxJobPreempt;

        /**
         * < Maximum number of pre-exec retry times.
         */
        public int maxPreExecRetry;

        /**
         * < Maximum number of pre-exec retry times for local cluster
         */
        public int localMaxPreExecRetry;

        /**
         * < Maximum number of job re-queue times.
         */
        public int maxJobRequeue;

        /**
         * < Use Linux-PAM
         */
        public int usePam;
        /* compute unit exclusive */

        /**
         * < Compute unit type
         */
        public int cu_type_exclusive;

        /**
         * < A string specified in EXCLUSIVE=CU[\<string>]
         */
        public String cu_str_exclusive;

        /**
         * < Resource reservation limit
         */
        public String resRsvLimit;

    }



    /**
     *  \addtogroup signal_action signal_action
     *  define status for signal action
     */

    /**
     * <  No action
     */
    public static final int ACT_NO = 0;

    /**
     * <  Start
     */
    public static final int ACT_START = 1;

    /**
     * <  Preempt
     */
    public static final int ACT_PREEMPT = 2;

    /**
     * <  Done
     */
    public static final int ACT_DONE = 3;

    /**
     * <  Fail
     */
    public static final int ACT_FAIL = 4;

    /**
     * \brief  host information entry.
     */
    public static class hostInfoEnt extends Structure {
        public static class ByReference extends hostInfoEnt implements Structure.ByReference {}
        public static class ByValue extends hostInfoEnt implements Structure.ByValue {}
        public hostInfoEnt() {}
        public hostInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < The name of the host.
         */
        public String host;

        /**
         * < The status of the host. It is the bitwise  inclusive OR.  see \ref host_status
         */
        public int hStatus;

        /**
         * < Indicate host loadSched busy reason
         */
        public IntByReference busySched;

        /**
         * < Indicate host loadStop  busy reason.
         */
        public IntByReference busyStop;

        /**
         * < The host CPU factor used to scale CPU load values to account for differences in CPU speeds. The faster the CPU, the larger the CPU factor.
         */
        public float cpuFactor;

        /**
         * < The number of load indices in the load, loadSched and loadStop arrays.
         */
        public int nIdx;

        /**
         * < Load information array on a host. This array gives the load information that is used for  scheduling batch jobs. This load information  is the effective load information from \ref ls_loadofhosts (see \ref ls_loadofhosts) plus the load reserved for running jobs (see lsb.queues for details on resource reservation). The load array is indexed the same as loadSched and loadStop  (see loadSched and loadStop below).
         */
        public FloatByReference load;

        /**
         * < Stop scheduling new jobs if over
         */
        public FloatByReference loadSched;

        /**
         * < Stop jobs if over this load. The loadSched and loadStop arrays control batch job scheduling, suspension, and resumption. The values in the loadSched array specify the scheduling thresholds for the corresponding load indices. Only if the current values of all specified load indices of this host are within (below or above, depending on the meaning of the load index) the corresponding thresholds of this host, will jobs be scheduled to run on this host. Similarly, the values in the loadStop array specify the stop thresholds for the corresponding load indices. If any of the load index values of the host goes beyond its stop threshold, the job will be suspended. The loadSched and loadStop arrays are indexed by the following constants:\n R15S\n 15-second average CPU run queue length.\n R1M\n 1-minute average CPU run queue length.\n R15M\n 15-minute average CPU run queue length.\n UT\n CPU utilization over the last minute.\n PG\n Average memory paging rate, in pages per second.\n IO\n Average disk I/O rate, in KB per second.\n LS\n Number of current login users.\n IT\n Idle time of the host in minutes.\n TMP\n The amount of free disk space in the file system containing /tmp, in MB.\n SWP\n The amount of swap space available, in MB.\n MEM\n The amount of available user memory on this host, in MB.
         */
        public FloatByReference loadStop;

        /**
         * < ASCII desp of run windows.One or more time windows in a week during which batch jobs may be dispatched to run on this host . The default is no restriction, or always open (i.e., 24 hours a day seven days a week). These windows are similar to the dispatch windows of batch job queues. See \ref lsb_queueinfo.
         */
        public String windows;

        /**
         * < The maximum number of job slots any user is allowed to use on this host.
         */
        public int userJobLimit;

        /**
         * < The maximum number of job slots that the host can process concurrently.
         */
        public int maxJobs;

        /**
         * < The number of job slots running or suspended on the host.
         */
        public int numJobs;

        /**
         * < The number of job slots running on the host.
         */
        public int numRUN;

        /**
         * < The number of job slots suspended by the batch daemon on the host.
         */
        public int numSSUSP;

        /**
         * < The number of job slots suspended by the job submitter or the LSF system administrator.
         */
        public int numUSUSP;

        /**
         * < The migration threshold in minutes after which a suspended job will be considered for migration.
         */
        public int mig;


        /**
         * < The host attributes; the bitwise inclusive OR of some of \ref host_attributes
         */
        public int attr;
        /**
         *  \addtogroup host_attributes host_attributes
         *  The host attributes
         */

        /**
         * < This host can checkpoint jobs
         */
        public static final int H_ATTR_CHKPNTABLE = 0x1;

        /**
         * < This host provides kernel support for checkpoint copy.
         */
        public static final int H_ATTR_CHKPNT_COPY = 0x2;

        /**
         * < The effective load of the host.
         */
        public FloatByReference realLoad;

        /**
         * < The number of job slots reserved by LSF for the PEND jobs.
         */
        public int numRESERVE;

        /**
         * < If attr has an H_ATTR_CHKPNT_COPY attribute, chkSig is set to the signal which triggers  checkpoint and copy operation. Otherwise,  chkSig is set to the signal which triggers  checkpoint operation on the host
         */
        public int chkSig;


        /**
         * < Num of resource used by the consumer
         */
        public float cnsmrUsage;

        /**
         * < Num of resource used by the provider
         */
        public float prvdrUsage;

        /**
         * < Num of resource available for the consumer to use
         */
        public float cnsmrAvail;

        /**
         * < Num of resource available for the provider to use
         */
        public float prvdrAvail;

        /**
         * < Num maximum of resource available in total
         */
        public float maxAvail;

        /**
         * < The job exit rate threshold on the host
         */
        public float maxExitRate;

        /**
         * < Number of job exit rate on the host
         */
        public float numExitRate;

        /**
         * < AdminAction - host control message
         */
        public String hCtrlMsg;

    }



    /**
     * \brief  Host information condition entry.
     */
    public static class condHostInfoEnt extends Structure {
        public static class ByReference extends condHostInfoEnt implements Structure.ByReference {}
        public static class ByValue extends condHostInfoEnt implements Structure.ByValue {}
        public condHostInfoEnt() {}
        public condHostInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < Host name
         */
        public String name;


        /**
         * < How many hosts are in the ok status
         */
        public int howManyOk;

        /**
         * < How many hosts are in the busy status
         */
        public int howManyBusy;

        /**
         * < How many hosts are in the closed status
         */
        public int howManyClosed;

        /**
         * < How many hosts are in the full status
         */
        public int howManyFull;

        /**
         * < How many hosts are in the unreach status
         */
        public int howManyUnreach;

        /**
         * < How many hosts are in the unavail status
         */
        public int howManyUnavail;


        /**
         * < The status of each host in the host group
         */
        public Pointer /* hostInfoEnt.ByReference */ hostInfo;

    }



    public static class adjustParam extends Structure {
        public static class ByReference extends adjustParam implements Structure.ByReference {}
        public static class ByValue extends adjustParam implements Structure.ByValue {}
        public adjustParam() {}
        public adjustParam(Pointer p) { super(p); read(); }


/* key name of share adjustment */
        public String key;

/* value of the key */
        public float value;
    }




/* cpu time factor */
    public static final int FAIR_ADJUST_CPU_TIME_FACTOR = 0;

/* run time factor */
    public static final int FAIR_ADJUST_RUN_TIME_FACTOR = 1;

/* run job factor */
    public static final int FAIR_ADJUST_RUN_JOB_FACTOR = 2;

/* committed run time factor */
    public static final int FAIR_ADJUST_COMMITTED_RUN_TIME_FACTOR = 3;

/* enable hist run time */
    public static final int FAIR_ADJUST_ENABLE_HIST_RUN_TIME = 4;

/* cpu time of finished jobs with decay */
    public static final int FAIR_ADJUST_HIST_CPU_TIME = 5;

/* cpu time of finished jobs within decay */
    public static final int FAIR_ADJUST_NEW_USED_CPU_TIME = 6;

/* total time that job spend in RUN state */
    public static final int FAIR_ADJUST_RUN_TIME = 7;

/* historical run time of finished jobs */
    public static final int FAIR_ADJUST_HIST_RUN_TIME = 8;

/* committed run time of started jobs */
    public static final int FAIR_ADJUST_COMMITTED_RUN_TIME = 9;

/* number of job slots used by started jobs */
    public static final int FAIR_ADJUST_NUM_START_JOBS = 10;

/* number of reserved slots used by pending jobs */
    public static final int FAIR_ADJUST_NUM_RESERVE_JOBS = 11;

/* total amount of memory used by started jobs */
    public static final int FAIR_ADJUST_MEM_USED = 12;

/* average memory allocated per slot */
    public static final int FAIR_ADJUST_MEM_ALLOCATED = 13;

/* total number of fairshare adjustment key value pairs */
    public static final int FAIR_ADJUST_KVPS_SUM = 14;

    //public String[] FairAdjustPairArrayName = new String[FAIR_ADJUST_KVPS_SUM];

    public static class shareAdjustPair extends Structure {
        public static class ByReference extends shareAdjustPair implements Structure.ByReference {}
        public static class ByValue extends shareAdjustPair implements Structure.ByValue {}
        public shareAdjustPair() {}
        public shareAdjustPair(Pointer p) { super(p); read(); }


/* queue share account */
        public static int SHAREACCTTYPEQUEUE = 0x01;

/* host partition share account */
        public static final int SHAREACCTTYPEHP = 0x02;

/* SLA share account */
        public static final int SHAREACCTTYPESLA = 0x04;

/* type of share account*/
        public int shareAcctType;

/* name of the share holder that use the share */
        public String holderName;

/* name of the provider policy name(name of queue, host partition or SLA) */
        public String providerName;

/* number of share adjustment key value pair */
        public int numPair;

/* share adjustment key value pair */
        public Pointer /* adjustParam.ByReference */ adjustParam;
    }



    // NOTE: Not in libbat
    //public static native float fairshare_adjustment(shareAdjustPair shareAdjustPair1);

/* For lsb_hostpartinfo() call */

    /**
     * \brief   gets user information about host partitions.
     */
    public static class hostPartUserInfo extends Structure {
        public static class ByReference extends hostPartUserInfo implements Structure.ByReference {}
        public static class ByValue extends hostPartUserInfo implements Structure.ByValue {}
        public hostPartUserInfo() {}
        public hostPartUserInfo(Pointer p) { super(p); read(); }


        /**
         * < The user name or user group name.  See \ref lsb_userinfo  and \ref lsb_usergrpinfo
         */
        public String user;

        /**
         * < The number of shares assigned to the user or user group, as configured in the file lsb.hosts. (See lsb.hosts.)
         */
        public int shares;

        /**
         * < The priority of the user or user group to use the host partition. Bigger values represent higher priorities. Jobs belonging to the user or user group with the highest priority are considered first for dispatch when resources in the host partition are being contended for. In general, a user or user group with more shares, fewer numStartJobs and less histCpuTime has higher priority. The storage for the array of hostPartInfoEnt structures will be reused by the next call.
         */
        public float priority;

        /**
         * < The number of job slots belonging to the user or user group that are running or suspended in the host partition.
         */
        public int numStartJobs;

        /**
         * < The normalized CPU time accumulated in the host partition during the recent period by finished jobs belonging to the user or user group. The period may be configured in the file lsb.params (see lsb.params), with a default value of five (5) hours.
         */
        public float histCpuTime;

        /**
         * < The number of job slots that are reserved for the PEND jobs belonging to the user or user group in the host partition.
         */
        public int numReserveJobs;

        /**
         * < The time unfinished jobs spend  in RUN state
         */
        public int runTime;

        /**
         * < The fairshare adjustment value from the fairshare plugin  (libfairshareadjust.ByReference ). The adjustment is enabled and weighted by setting the value of FAIRSHARE_ADJUSTMENT_FACTOR in lsb.params.
         */
        public float shareAdjustment;
    }



/* For lsb_hostpartinfo() call */

    /**
     * \brief  gets information entry about host partitions.
     */
    public static class hostPartInfoEnt extends Structure {
        public static class ByReference extends hostPartInfoEnt implements Structure.ByReference {}
        public static class ByValue extends hostPartInfoEnt implements Structure.ByValue {}
        public hostPartInfoEnt() {}
        public hostPartInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < The name of the host partition
         */
        public byte[] hostPart = new byte[MAX_LSB_NAME_LEN];

        /**
         * < A blank-separated list of names of hosts and host groups which are members of the host partition. The name of a host group has a '/' appended. see \ref lsb_hostgrpinfo.
         */
        public String hostList;

        /**
         * < The number of users in this host partition. i.e., the number of hostPartUserInfo structures.
         */
        public int numUsers;

        /**
         * < An array of hostPartUserInfo structures which hold information on users in this host partition.
         */
        public Pointer /* hostPartUserInfo.ByReference */ users;
    }



/* Library rappresentation of the share account */

    /**
     * \brief Library rappresentation of the share account
     */
    public static class shareAcctInfoEnt extends Structure {
        public static class ByReference extends shareAcctInfoEnt implements Structure.ByReference {}
        public static class ByValue extends shareAcctInfoEnt implements Structure.ByValue {}
        public shareAcctInfoEnt() {}
        public shareAcctInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < The user name or user group name. (See \ref lsb_userinfo and \ref lsb_usergrpinfo.)
         */
        public String shareAcctPath;

        /**
         * < The number of shares assigned to  the user or user group, as configured in the file lsb.queues.
         */
        public int shares;

        /**
         * < The priority of the user or user group in the fairshare queue. Larger values represent higher priorities. Job belonging to the user or user group with the highest priority are considered first for dispatch in the fairshare queue. In general, a user or user group with more shares, fewer numStartJobs and less histCpuTime has higher priority.
         */
        public float priority;

        /**
         * < The number of job slots (belonging to the user or user group) that are running or suspended in the fairshare queue.
         */
        public int numStartJobs;

        /**
         * < The normalized CPU time accumulated in the fairshare queue by jobs belonging to the user or user group, over the time period configured in the file lsb.params. The default time period is 5 hours.
         */
        public float histCpuTime;

        /**
         * < The number of job slots that are reserved for the PEND jobs belonging to the user or user group in the host partition.
         */
        public int numReserveJobs;

        /**
         * < The time unfinished jobs spend in the RUN state.
         */
        public int runTime;

        /**
         * < The fairshare adjustment value from the fairshare plugin  (libfairshareadjust.SOEXT). The adjustment is enabled and weighted  by setting the value of FAIRSHARE_ADJUSTMENT_FACTOR in lsb.params.
         */
        public float shareAdjustment;
    }



/* boundaries and default value used by mbatchd for the maxJobId */
    public static final int DEF_MAX_JOBID = 999999;
    public static final int MAX_JOBID_LOW = 999999;
    public static final int MAX_JOBID_HIGH = (LibLsf.INFINIT_INT - 1);


/* default preemption wait time */
    public static final int DEF_PREEMPTION_WAIT_TIME = 300;

/* default number of hosts specified by -m */
    public static final int DEF_MAX_ASKED_HOSTS = 512;

/* For lsb_parameterinfo() call */

    /**
     * \brief The parameterInfo structure contains the following fields:
     */
    public static class parameterInfo extends Structure {
        public static class ByReference extends parameterInfo implements Structure.ByReference {}
        public static class ByValue extends parameterInfo implements Structure.ByValue {}
        public parameterInfo() {}
        public parameterInfo(Pointer p) { super(p); read(); }


        /**
         * < DEFAULT_QUEUE: A blank_separated list of queue names for automatic queue selection.
         */
        public String defaultQueues;

        /**
         * < DEFAULT_HOST_SPEC: The host name or host model name used as the system default for scaling CPULIMIT and RUNLIMIT.
         */
        public String defaultHostSpec;

        /**
         * < MBD_SLEEP_TIME: The interval in seconds at which the mbatchd dispatches jobs.
         */
        public int mbatchdInterval;

        /**
         * < SBD_SLEEP_TIME: The interval in seconds at which the sbatchd suspends or resumes jobs.
         */
        public int sbatchdInterval;

        /**
         * < JOB_ACCEPT_INTERVAL: The interval at which  a host accepts two successive jobs. (In units of SBD_SLEEP_TIME.)
         */
        public int jobAcceptInterval;

        /**
         * < MAX_RETRY: The maximum number of retries for dispatching a job.
         */
        public int maxDispRetries;

        /**
         * < MAX_SBD_FAIL: The maximum number of retries for reaching an sbatchd.
         */
        public int maxSbdRetries;

        /**
         * < PREEM_PERIOD: The interval in seconds for preempting jobs running on the same host.
         */
        public int preemptPeriod;

        /**
         * < CLEAN_PERIOD: The interval in seconds during which finished jobs are kept in core.
         */
        public int cleanPeriod;

        /**
         * < MAX_JOB_NUM: The maximum number of finished jobs that are logged in the current event file.
         */
        public int maxNumJobs;

        /**
         * < HIST_HOURS: The number of hours of resource consumption history used for fair share scheduling and scheduling within a host partition.
         */
        public float historyHours;

        /**
         * < PG_SUSP_IT: The interval a host must be idle before resuming a job suspended for excessive paging.
         */
        public int pgSuspendIt;

        /**
         * < The default project assigned to jobs.
         */
        public String defaultProject;

        /**
         * < Job submission retry interval
         */
        public int retryIntvl;

        /**
         * < For Cray NQS compatiblilty only. Used by LSF to get the NQS queue information
         */
        public int nqsQueuesFlags;

        /**
         * < nqsRequestsFlags
         */
        public int nqsRequestsFlags;

        /**
         * < The maximum number of times to attempt the preexecution command of a job from a remote cluster ( MultiCluster only)
         */
        public int maxPreExecRetry;

        /**
         * < Maximum number of pre-exec retry times for local cluster
         */
        public int localMaxPreExecRetry;

        /**
         * < Event watching Interval in seconds
         */
        public int eventWatchTime;

        /**
         * < Run time weighting factor for fairshare scheduling
         */
        public float runTimeFactor;

        /**
         * < Used for calcultion of the fairshare scheduling formula
         */
        public float waitTimeFactor;

        /**
         * < Job slots weighting factor for fairshare scheduling
         */
        public float runJobFactor;

        /**
         * < Default check interval
         */
        public int eEventCheckIntvl;

        /**
         * < sbatchd report every sbd_sleep_time
         */
        public int rusageUpdateRate;

        /**
         * < sbatchd updates jobs jRusage in mbatchd if more than 10% changes
         */
        public int rusageUpdatePercent;

        /**
         * < Time period to check for reconfig
         */
        public int condCheckTime;

        /**
         * < The maximum number of connections between master and slave batch daemons
         */
        public int maxSbdConnections;

        /**
         * < The interval for rescheduling jobs
         */
        public int rschedInterval;

        /**
         * < Max time mbatchd stays in scheduling routine, after which take a breather
         */
        public int maxSchedStay;

        /**
         * < During which load remains fresh
         */
        public int freshPeriod;

        /**
         * < The preemption behavior, GROUP_MAX, GROUP_JLP, USER_JLP, HOST_JLU,MINI_JOB, LEAST_RUN_TIME
         */
        public int preemptFor;

        /**
         * < Flags whether users can resume their jobs when suspended by the LSF administrator
         */
        public int adminSuspend;

        /**
         * < Flags to enable/disable normal user to create advance reservation
         */
        public int userReservation;

        /**
         * < CPU time weighting factor for fairshare scheduling
         */
        public float cpuTimeFactor;

        /**
         * < The starting month for a fiscal year
         */
        public int fyStart;

        /**
         * < The maximum number of jobs in a job array
         */
        public int maxJobArraySize;

        /**
         * < Replay period for exceptions, in seconds
         */
        public NativeLong exceptReplayPeriod;

        /**
         * < The interval to terminate a job
         */
        public int jobTerminateInterval;

        /**
         * <  User level account mapping for remote jobs is disabled
         */
        public int disableUAcctMap;

        /**
         * < If set to TRUE, Project name for a job will be considerred when doing fairshare scheduling, i.e., as if user has submitted jobs using -G
         */
        public int enforceFSProj;

        /**
         * < Enforces the check to see if the invoker of bsub is in the specifed group when the -P option is used
         */
        public int enforceProjCheck;

        /**
         * < Run time for a job
         */
        public int jobRunTimes;

        /**
         * < Event table Job default interval
         */
        public int dbDefaultIntval;

        /**
         * < Event table Job Host Count
         */
        public int dbHjobCountIntval;

        /**
         * < Event table Job Queue Count
         */
        public int dbQjobCountIntval;

        /**
         * < Event table Job User Count
         */
        public int dbUjobCountIntval;

        /**
         * < Event table Job Resource Interval
         */
        public int dbJobResUsageIntval;

        /**
         * < Event table Resource Load Interval
         */
        public int dbLoadIntval;

        /**
         * < Event table Job Info
         */
        public int dbJobInfoIntval;

        /**
         * < Used with job dependency scheduling
         */
        public int jobDepLastSub;

        /**
         * < Used with job dependency scheduling,  deprecated
         */
        public int maxJobNameDep;

        /**
         * < Select resources to be logged
         */
        public String dbSelectLoad;

        /**
         * < Job synchronizes its group status
         */
        public int jobSynJgrp;

        /**
         * < The batch jobs' temporary output directory
         */
        public String pjobSpoolDir;


        /**
         * < Maximal job priority defined for all users
         */
        public int maxUserPriority;

        /**
         * < Job priority is increased by the system dynamically based on waiting time
         */
        public int jobPriorityValue;

        /**
         * < Waiting time to increase Job priority by the system dynamically
         */
        public int jobPriorityTime;

        /**
         * < Enable internal statistical adjustment
         */
        public int enableAutoAdjust;

        /**
         * < Start to autoadjust when the user has  this number of pending jobs
         */
        public int autoAdjustAtNumPend;

        /**
         * < If this number of jobs has been visited skip the user
         */
        public float autoAdjustAtPercent;

        /**
         * <  Static shared resource update interval for the cluster actor
         */
        public int sharedResourceUpdFactor;

        /**
         * < Schedule job based on raw load info
         */
        public int scheRawLoad;

        /**
         * <  The batch jobs' external storage for attached data
         */
        public String jobAttaDir;

        /**
         * < Maximum message number for each job
         */
        public int maxJobMsgNum;

        /**
         * < Maximum attached data size to be transferred for each message
         */
        public int maxJobAttaSize;

        /**
         * < The life time of a child MBD to serve queries in the MT way
         */
        public int mbdRefreshTime;

        /**
         * < The interval of the execution cluster updating the job's resource usage
         */
        public int updJobRusageInterval;

        /**
         * < The account to which all windows workgroup users are to be mapped
         */
        public String sysMapAcct;

        /**
         * < Dispatch delay internal
         */
        public int preExecDelay;

        /**
         * < Update duplicate event interval
         */
        public int updEventUpdateInterval;

        /**
         * < Resources are reserved for parallel jobs on a per-slot basis
         */
        public int resourceReservePerSlot;

        /**
         * < Maximum job id --- read from the lsb.params
         */
        public int maxJobId;

        /**
         * < Define a list of preemptable resource  names
         */
        public String preemptResourceList;

        /**
         * < The preemption wait time
         */
        public int preemptionWaitTime;

        /**
         * < Maximum number of rollover lsb.acct files kept by mbatchd.
         */
        public int maxAcctArchiveNum;

        /**
         * < mbatchd Archive Interval
         */
        public int acctArchiveInDays;

        /**
         * < mbatchd Archive threshold
         */
        public int acctArchiveInSize;

        /**
         * < Committed run time weighting factor
         */
        public float committedRunTimeFactor;

        /**
         * < Enable the use of historical run time in the calculation of fairshare scheduling priority, Disable the use of historical run time in the calculation of fairshare scheduling priority
         */
        public int enableHistRunTime;

/*#ifdef PS_SXNQS */
/**< NQS resource usage update interval */
/*    public int   nqsUpdateInterval;*/
/*#endif */

        /**
         * < Open lease reclaim time
         */
        public int mcbOlmReclaimTimeDelay;

        /**
         * < Enable chunk job dispatch for jobs with CPU limit or run limits
         */
        public int chunkJobDuration;

        /**
         * < The interval for scheduling jobs by scheduler daemon
         */
        public int sessionInterval;

        /**
         * < The number of jobs per user per queue whose pending reason is published at the PEND_REASON_UPDATE_INTERVAL interval
         */
        public int publishReasonJobNum;

        /**
         * < The interval for publishing job pending reason by scheduler daemon
         */
        public int publishReasonInterval;

        /**
         * < Interval(in seconds) of pending reason  publish for all jobs
         */
        public int publishReason4AllJobInterval;

        /**
         * < MC pending reason update interval (0 means no updates)
         */
        public int mcUpdPendingReasonInterval;

        /**
         * < MC pending reason update package size (0 means no limit)
         */
        public int mcUpdPendingReasonPkgSize;

        /**
         * < No preemption if the run time is greater  than the value defined in here
         */
        public int noPreemptRunTime;

        /**
         * < No preemption if the finish time is less than the value defined in here
         */
        public int noPreemptFinishTime;

        /**
         * < mbatchd Archive Time
         */
        public String acctArchiveAt;

        /**
         * < Absolute run limit for job
         */
        public int absoluteRunLimit;

        /**
         * < The job exit rate duration
         */
        public int lsbExitRateDuration;

        /**
         * <  The duration to trigger eadmin
         */
        public int lsbTriggerDuration;

        /**
         * < Maximum time for job information query commands (for example,with bjobs) to wait
         */
        public int maxJobinfoQueryPeriod;

        /**
         * < Job submission retrial interval for client
         */
        public int jobSubRetryInterval;

        /**
         * < System wide max pending jobs
         */
        public int pendingJobThreshold;


        /**
         * < Max number of concurrent query
         */
        public int maxConcurrentJobQuery;

        /**
         * < Min event switch time period
         */
        public int minSwitchPeriod;


        /**
         * < Condense pending reasons enabled
         */
        public int condensePendingReasons;

        /**
         * < Schedule Parallel jobs based on slots instead of CPUs
         */
        public int slotBasedParallelSched;

        /**
         * < Disable user job movement operations, like btop/bbot.
         */
        public int disableUserJobMovement;

        /**
         * < Detect and report idle jobs only after specified minutes.
         */
        public int detectIdleJobAfter;
        public int useSymbolPriority;
        /**
         * < Use symbolic when specifing priority of symphony jobs/
         * <p/>
         * /**< Priority rounding for symphony jobs
         */
        public int JobPriorityRound;

        /**
         * < The mapping of the symbolic priority  for symphony jobs
         */
        public String priorityMapping;

        /**
         * < Maximum number of subdirectories under LSB_SHAREDIR/cluster/logdir/info
         */
        public int maxInfoDirs;

        /**
         * < The minimum period of a child MBD to serve queries in the MT way
         */
        public int minMbdRefreshTime;

        /**
         * < Stop asking license to LS not due to lack license
         */
        public int enableStopAskingLicenses2LS;

        /**
         * < Expire time for finished job which will not taken into account when calculating queue fairshare priority
         */
        public int expiredTime;

        /**
         * < MBD child query processes will only run on the following CPUs
         */
        public String mbdQueryCPUs;

        /**
         * < The default application profile assigned to jobs
         */
        public String defaultApp;

        /**
         * < Enable or disable data streaming
         */
        public int enableStream;

        /**
         * < File to which lsbatch data is streamed
         */
        public String streamFile;

        /**
         * < File size in MB to which lsbatch data is streamed
         */
        public int streamSize;

        /**
         * < Sync up host status with master LIM is enabled
         */
        public int syncUpHostStatusWithLIM;

        /**
         * < Project schedulign default SLA
         */
        public String defaultSLA;

        /**
         * < EGO Enabled SLA scheduling timer period
         */
        public int slaTimer;

        /**
         * < EGO Enabled SLA scheduling time to live
         */
        public int mbdEgoTtl;

        /**
         * < EGO Enabled SLA scheduling connection timeout
         */
        public int mbdEgoConnTimeout;

        /**
         * < EGO Enabled SLA scheduling read timeout
         */
        public int mbdEgoReadTimeout;

        /**
         * < EGO Enabled SLA scheduling use MXJ flag
         */
        public int mbdUseEgoMXJ;

        /**
         * < EGO Enabled SLA scheduling reclaim by queue
         */
        public int mbdEgoReclaimByQueue;

        /**
         * < EGO Enabled SLA scheduling default velocity
         */
        public int defaultSLAvelocity;

        /**
         * < Type of host exit rate exception handling types: EXIT_RATE_TYPE
         */
        public String exitRateTypes;

        /**
         * < Type of host exit rate exception handling types: GLOBAL_EXIT_RATE
         */
        public float globalJobExitRate;

        /**
         * < Type of host exit rate exception handling types ENABLE_EXIT_RATE_PER_SLOT
         */
        public int enableJobExitRatePerSlot;

        /**
         * < Performance metrics monitor is enabled  flag
         */
        public int enableMetric;

        /**
         * < Performance metrics monitor sample period flag
         */
        public int schMetricsSample;

        /**
         * < Used to bound: (1) factors, (2) weights, and (3) APS values
         */
        public float maxApsValue;

        /**
         * < Child mbatchd gets updated information about new jobs from the parent mbatchd
         */
        public int newjobRefresh;

        /**
         * < Job type to preempt, PREEMPT_JOBTYPE_BACKFILL, PREEMPT_JOBTYPE_EXCLUSIVE
         */
        public int preemptJobType;

        /**
         * < The default job group assigned to jobs
         */
        public String defaultJgrp;

        /**
         * < Max ratio between run limit and runtime estimation
         */
        public int jobRunlimitRatio;

        /**
         * < Enable the post-execution processing of the job to be included as part of the job flag
         */
        public int jobIncludePostproc;

        /**
         * < Timeout of post-execution processing
         */
        public int jobPostprocTimeout;

        /**
         * < The interval, in seconds, for updating the session scheduler status summary
         */
        public int sschedUpdateSummaryInterval;

        /**
         * < The number of completed tasks for updating the session scheduler status summary
         */
        public int sschedUpdateSummaryByTask;

        /**
         * < The maximum number of times a task can be requeued via requeue exit values
         */
        public int sschedRequeueLimit;

        /**
         * < The maximum number of times a task can be retried after a dispatch error
         */
        public int sschedRetryLimit;

        /**
         * < The maximum number of tasks that can be submitted in one session
         */
        public int sschedMaxTasks;

        /**
         * < The maximum run time of a single task
         */
        public int sschedMaxRuntime;

        /**
         * < The output directory for task accounting files
         */
        public String sschedAcctDir;

        /**
         * < If TRUE enable the job group automatic deletion functionality (default is FALSE).
         */
        public int jgrpAutoDel;

        /**
         * < Maximum number of job preempted times
         */
        public int maxJobPreempt;

        /**
         * < Maximum number of job re-queue times
         */
        public int maxJobRequeue;

        /**
         * < No preempt run time percent
         */
        public int noPreemptRunTimePercent;

        /**
         * < No preempt finish time percent
         */
        public int noPreemptFinishTimePercent;


        /**
         * < The reservation request being within JL/U.
         */
        public int slotReserveQueueLimit;

        /**
         * < Job accept limit percentage.
         */
        public int maxJobPercentagePerSession;

        /**
         * < The low priority job will use the slots freed by preempted jobs.
         */
        public int useSuspSlots;


        /**
         * < Maximum number of the backup stream.utc files
         */
        public int maxStreamFileNum;

        /**
         * < If enforced only admin can use bkill -r option
         */
        public int privilegedUserForceBkill;

        /**
         * < It controls the remote queue selection flow.
         */
        public int mcSchedulingEnhance;

        /**
         * < It controls update interval of the counters  and other original data in MC implementation
         */
        public int mcUpdateInterval;

        /**
         * < Jobs run on only on hosts belonging to the intersection of the queue the job was submitted to, advance reservation hosts, and any hosts specified by bsub -m at the time of submission.
         */
        public int intersectCandidateHosts;

        /**
         * < Enforces the limitations of a single specified user group.
         */
        public int enforceOneUGLimit;

        /**
         * < Enable or disable logging runtime estimation exceeded event
         */
        public int logRuntimeESTExceeded;

        /**
         * < Compute unit types.
         */
        public String computeUnitTypes;

        /**
         * < Fairshare adjustment weighting factor
         */
        public float fairAdjustFactor;

        /**
         * < abs runtime and cputime for LSF simulator
         */
        public int simAbsoluteTime;

        /**
         * < switch for job exception enhancement
         */
        public int extendJobException;
    }

    /* parameterInfo */


/* Bits for preemptFor parameter */
    public static final int GROUP_MAX = 0x0001;
    public static final int GROUP_JLP = 0x0002;
    public static final int USER_JLP = 0x0004;
    public static final int HOST_JLU = 0x0008;

/* minimum of job */
    public static final int MINI_JOB = 0x0010;

/* least run time */
    public static final int LEAST_RUN_TIME = 0x0020;

/* optimal mini job */
    public static final int OPTIMAL_MINI_JOB = 0x0040;

/* Bits for mcSchedulingEnhance parameter */
    public static final int RESOURCE_ONLY = 0x0001;
    public static final int COUNT_PREEMPTABLE = 0x0002;
    public static final int HIGH_QUEUE_PRIORITY = 0x0004;
    public static final int PREEMPTABLE_QUEUE_PRIORITY = 0x0008;
    public static final int PENDING_WHEN_NOSLOTS = 0x0010;

/* options for bcaladd, bcalmod, bcaldel */
    public static final int CAL_FORCE = 0x0001;

/* Bits for preemptJobType parameter,
*  used to enable backfill and exclusive
*  preemption */
    public static final int PREEMPT_JOBTYPE_EXCLUSIVE = 0x0001;
    public static final int PREEMPT_JOBTYPE_BACKFILL = 0x0002;

/* For lsb_calendarinfo() call */

    /**
     * \brief  calendar Information Entry.
     */
    public static class calendarInfoEnt extends Structure {
        public static class ByReference extends calendarInfoEnt implements Structure.ByReference {}
        public static class ByValue extends calendarInfoEnt implements Structure.ByValue {}
        public calendarInfoEnt() {}
        public calendarInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < A pointer to the name of the calendar.
         */
        public String name;

        /**
         * < A description string associated with the calendar.
         */
        public String desc;

        /**
         * < Calendar Expression
         */
        public String calExpr;

        /**
         * < User name
         */
        public String userName;

        /**
         * < Calendar status
         */
        public int status;

        /**
         * < For future use
         */
        public int options;

        /**
         * < Last time event of the calendar
         */
        public int lastDay;

        /**
         * < Next time event of the calendar
         */
        public int nextDay;

        /**
         * < Create Time
         */
        public NativeLong creatTime;

        /**
         * < Last Modify Time
         */
        public NativeLong lastModifyTime;

        /**
         * < Type of calendar, etc.
         */
        public int flags;
    }



    public static final int ALL_CALENDARS = 0x1;

    public static final int EVE_HIST = 0x1;
    public static final int EVENT_ACTIVE = 1;
    public static final int EVENT_INACTIVE = 2;
    public static final int EVENT_REJECT = 3;

    public static final int EVENT_TYPE_UNKNOWN = 0;
    public static final int EVENT_TYPE_LATCHED = 1;
    public static final int EVENT_TYPE_PULSEALL = 2;
    public static final int EVENT_TYPE_PULSE = 3;
    public static final int EVENT_TYPE_EXCLUSIVE = 4;

/* define event types */
    public static final int EV_UNDEF = 0;
    public static final int EV_FILE = 1;
    public static final int EV_EXCEPT = 2;
    public static final int EV_USER = 3;

    public static class loadInfoEnt extends Structure {
        public static class ByReference extends loadInfoEnt implements Structure.ByReference {}
        public static class ByValue extends loadInfoEnt implements Structure.ByValue {}
        public loadInfoEnt() {}
        public loadInfoEnt(Pointer p) { super(p); read(); }

        public String hostName;
        public int status;
        public FloatByReference load;
    }



    public static class queuePairEnt extends Structure {
        public static class ByReference extends queuePairEnt implements Structure.ByReference {}
        public static class ByValue extends queuePairEnt implements Structure.ByValue {}
        public queuePairEnt() {}
        public queuePairEnt(Pointer p) { super(p); read(); }

        public String local;
        public String remote;
        public int send;
        public int status;
    }



    public static class rmbCluAppEnt extends Structure {
        public static class ByReference extends rmbCluAppEnt implements Structure.ByReference {}
        public static class ByValue extends rmbCluAppEnt implements Structure.ByValue {}
        public rmbCluAppEnt() {}
        public rmbCluAppEnt(Pointer p) { super(p); read(); }

        public String name;
        public String description;
    }



/* define 'cluster status' in lease model
*  for bclusters command
 */


/* disconnection */
    public static final int LEASE_CLU_STAT_DISC = 1;

/* policy is exchanged but no lease is signed */
    public static final int LEASE_CLU_STAT_CONN = 2;

/* there are leases signed between two clusters */
    public static final int LEASE_CLU_STAT_OK = 3;
    public static final int LEASE_CLU_STAT_NUMBER = 3;
/* consumer cluster status in lease model */

    public static class consumerCluEnt extends Structure {
        public static class ByReference extends consumerCluEnt implements Structure.ByReference {}
        public static class ByValue extends consumerCluEnt implements Structure.ByValue {}
        public consumerCluEnt() {}
        public consumerCluEnt(Pointer p) { super(p); read(); }


/* consumer cluster name */
        public String cluName;

/* cluster status, Ref- 'cluster status' definitions */
        public int status;
    }


/* provider cluster status in lease model */

    public static class providerCluEnt extends Structure {
        public static class ByReference extends providerCluEnt implements Structure.ByReference {}
        public static class ByValue extends providerCluEnt implements Structure.ByValue {}
        public providerCluEnt() {}
        public providerCluEnt(Pointer p) { super(p); read(); }


/* provider cluster name */
        public String cluName;

/* cluster status, Ref- 'cluster status' definitions */
        public int status;
    }


/* for remote batch model, its definition is same as  clusterInfoEnt*/

    public static class rmbCluInfoEnt extends Structure {
        public static class ByReference extends rmbCluInfoEnt implements Structure.ByReference {}
        public static class ByValue extends rmbCluInfoEnt implements Structure.ByValue {}
        public rmbCluInfoEnt() {}
        public rmbCluInfoEnt(Pointer p) { super(p); read(); }

        public String cluster;
        public int numPairs;
        public Pointer /* queuePairEnt.ByReference */ queues;
        public int numApps;
        public Pointer /* rmbCluAppEnt.ByReference */ apps;
    }



/* for leasing model */

    public static class leaseCluInfoEnt extends Structure {
        public static class ByReference extends leaseCluInfoEnt implements Structure.ByReference {}
        public static class ByValue extends leaseCluInfoEnt implements Structure.ByValue {}
        public leaseCluInfoEnt() {}
        public leaseCluInfoEnt(Pointer p) { super(p); read(); }


/* 1, import from all if "allremote" defined in lease queue*/
        public int flags;

/* the array size of consumer cluster array */
        public int numConsumer;

/* the consumer cluster array */
        public Pointer /* consumerCluEnt.ByReference */ consumerClus;

/* the array size of provider cluster array */
        public int numProvider;

/* the provider cluster array */
        public Pointer /* providerCluEnt.ByReference */ providerClus;
    }



/* This is the old data structure, we
*  leave it here to keep backward compatibility.
*  It's definition is same as structure rmbCluInfoEnt.
*  It is to transfer cluster status between mbatchd with
*  old(4.x) bclusters command and old API-lsb_clusterinfo()
 */

    public static class clusterInfoEnt extends Structure {
        public static class ByReference extends clusterInfoEnt implements Structure.ByReference {}
        public static class ByValue extends clusterInfoEnt implements Structure.ByValue {}
        public clusterInfoEnt() {}
        public clusterInfoEnt(Pointer p) { super(p); read(); }

        public String cluster;
        public int numPairs;
        public Pointer /* queuePairEnt.ByReference */ queues;
        public int numApps;
        public Pointer /* rmbCluAppEnt.ByReference */ apps;
    }


/* the new data structure to transfer cluster status between mbatchd with
*  new(5.0) bclusters command and new API-lsb_clusterinfoEx()
 */

    public static class clusterInfoEntEx extends Structure {
        public static class ByReference extends clusterInfoEntEx implements Structure.ByReference {}
        public static class ByValue extends clusterInfoEntEx implements Structure.ByValue {}
        public clusterInfoEntEx() {}
        public clusterInfoEntEx(Pointer p) { super(p); read(); }


/* cluster status related to remote batch*/
        public rmbCluInfoEnt.ByReference rmbCluInfo;

/* cluster status related to resource lease*/
        public leaseCluInfoEnt leaseCluInfo;
    }



    public static class eventInfoEnt extends Structure {
        public static class ByReference extends eventInfoEnt implements Structure.ByReference {}
        public static class ByValue extends eventInfoEnt implements Structure.ByValue {}
        public eventInfoEnt() {}
        public eventInfoEnt(Pointer p) { super(p); read(); }


/* name of event */
        public String name;

/* one of ACTIVE or INACTIVE */
        public int status;

/* one of LATCHED, PULSE and EXCLUSIVE */
        public int type;

/* one of FILE, ALARM, USER */
        public int eType;

/* user who created the event */
        public String userName;

/* event's attributes sent back from eeventd */
        public String attributes;

/* number of expression dependent on the event */
        public int numDependents;

/* last time when eeventd sent back message */
        public NativeLong updateTime;

/* last dispatched job dependent on the event */
        public long lastDisJob;

/* the time when the last job was dispatched */
        public NativeLong lastDisTime;
    }


    public static final int ALL_EVENTS = 0x01;

    /**
     *  \addtogroup groupinfo_define groupinfo_define
     *  define options for \ref lsb_usergrpinfo and \ref lsb_hostgrpinfo calls
     */

    /**
     * < User group
     */
    public static final int USER_GRP = 0x1;

    /**
     * < Host group
     */
    public static final int HOST_GRP = 0x2;

    /**
     * < Host part group
     */
    public static final int HPART_HGRP = 0x4;
    /**
     *  \defgroup group_membership_option group_membership_option
     *  \ingroup groupinfo_define
     *  group membership options
     */

    /**
     * < Expand the group membership recursively. That is, if a member of a group is itself a group, give the names of its members recursively, rather than its name, which is the default.
     */
    public static final int GRP_RECURSIVE = 0x8;

    /**
     * < Get membership of all groups.
     */
    public static final int GRP_ALL = 0x10;

    /**
     * < NQSQ_GRP
     */
    public static final int NQSQ_GRP = 0x20;

    /**
     * < Group shares
     */
    public static final int GRP_SHARES = 0x40;

    /**
     * < Dynamic group
     */
    public static final int DYNAMIC_GRP = 0x800;

    /**
     * < Group cu
     */
    public static final int GRP_CU = 0x1000;

    /**
     * \brief Structure for representing the shares assigned to a user group.
     */
    public static class userShares extends Structure {
        public static class ByReference extends userShares implements Structure.ByReference {}
        public static class ByValue extends userShares implements Structure.ByValue {}
        public userShares() {}
        public userShares(Pointer p) { super(p); read(); }


        /**
         * < This can be a user or a keyword "default" or others
         */
        public String user;

        /**
         * < The number of shares assigned to the user
         */
        public int shares;
    }




    /**
     * \brief  group information entry.
     */
    public static class groupInfoEnt extends Structure {
        public static class ByReference extends groupInfoEnt implements Structure.ByReference {}
        public static class ByValue extends groupInfoEnt implements Structure.ByValue {}
        public groupInfoEnt() {}
        public groupInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < Group name
         */
        public String group;

        /**
         * < ASCII list of member names
         */
        public String memberList;

        /**
         * < ASCII list of admin member names
         */
        public String adminMemberList;

        /**
         * < The number of users with shares
         */
        public int numUserShares;

        /**
         * < The user shares rappresentation
         */
        public Pointer /* userShares.ByReference */ userShares;

        /**
         *  \addtogroup group_define group_define
         *   group define statements
         */

        /**
         * < Group output is in regular (uncondensed) format.
         */
        public static final int GRP_NO_CONDENSE_OUTPUT = 0x01;

        /**
         * < Group output is in condensed format.
         */
        public static final int GRP_CONDENSE_OUTPUT = 0x02;

        /**
         * < Group have regular expresion
         */
        public static final int GRP_HAVE_REG_EXP = 0x04;

        /**
         * < Group is a service class.
         */
        public static final int GRP_SERVICE_CLASS = 0x08;

        /**
         * < Group is a compute unit.
         */
        public static final int GRP_IS_CU = 0x10;

        /**
         * < Options.see \ref group_define
         */
        public int options;

        /**
         * < Host membership pattern
         */
        public String pattern;

        /**
         * < Negation membership pattern
         */
        public String neg_pattern;

        /**
         * < Compute unit type
         */
        public int cu_type;
    }



    /**
     * \brief  run job request.
     */
    public static class runJobRequest extends Structure {
        public static class ByReference extends runJobRequest implements Structure.ByReference {}
        public static class ByValue extends runJobRequest implements Structure.ByValue {}
        public runJobRequest() {}
        public runJobRequest(Pointer p) { super(p); read(); }


        /**
         * < Jobid of the requested job
         */
        public long jobId;

        /**
         * < The number of hosts
         */
        public int numHosts;

        /**
         * < Vector of hostnames
         */
        public Pointer hostname;
        /**
         *  \addtogroup runjob_option runjob_option
         *  Options used for lsb_runjob:
         */


        /**
         * < Normal jobs
         */
        public static final int RUNJOB_OPT_NORMAL = 0x01;

        /**
         * < Nostop jobs
         */
        public static final int RUNJOB_OPT_NOSTOP = 0x02;

        /**
         * < Pending jobs only, no finished jobs
         */
        public static final int RUNJOB_OPT_PENDONLY = 0x04;

        /**
         * < Check point job only, from beginning
         */
        public static final int RUNJOB_OPT_FROM_BEGIN = 0x08;

        /**
         * < brun to use free CPUs only
         */
        public static final int RUNJOB_OPT_FREE = 0x10;

        /**
         * < brun ignoring rusage
         */
        public static final int RUNJOB_OPT_IGNORE_RUSAGE = 0x20;

        /**
         * < Run job request options, see \ref runjob_option
         */
        public int options;

        /**
         * < Vector of number of slots per host
         */
        public IntByReference slots;
    }



    /**
     *  \addtogroup external_msg_processing external_msg_processing
     *  options for \ref lsb_readjobmsg call
     */

    /**
     *  \defgroup external_msg_post external_msg_post
     *  options specifying if the message has an attachment to be posted
     */

    /**
     * < Post the external job message. There  is no attached data file.
     */
    public static final int EXT_MSG_POST = 0x01;

    /**
     * < Post the external job message and data file posted to the job.
     */
    public static final int EXT_ATTA_POST = 0x02;

    /**
     * <Read the external job message. There is no attached data file.
     */
    public static final int EXT_MSG_READ = 0x04;

    /**
     * < Read the external job message and data file posted to the job.If there is no data file attached, the error message "The attached data of the message is not available" is displayed, and the external job message is displayed.
     */
    public static final int EXT_ATTA_READ = 0x08;

    /**
     * < Replay the external message
     */
    public static final int EXT_MSG_REPLAY = 0x10;

    /**
     * < Post the external job noevent message
     */
    public static final int EXT_MSG_POST_NOEVENT = 0x20;


    /**
     * \brief structure jobExternalMsgReq contains the information required to
     * define an external message of a job.
     */
    public static class jobExternalMsgReq extends Structure {
        public static class ByReference extends jobExternalMsgReq implements Structure.ByReference {}
        public static class ByValue extends jobExternalMsgReq implements Structure.ByValue {}
        public jobExternalMsgReq() {}
        public jobExternalMsgReq(Pointer p) { super(p); read(); }


        /**
         * < Specifies if the message has an attachment to be read.<lsf/lsbatch.h> defines the following flags constructed from bits. These flags correspond to options.\n EXT_MSG_READ\n Read the external job message. There is no attached data file.\n EXT_ATTA_READ\n Read the external job message and data file posted to the job.\n If there is no data file attached, the error message "The attached data of the message is not available" is displayed, and the external job  message is displayed.
         */
        public int options;

        /**
         * < The system generated job Id of the job.
         */
        public long jobId;

        /**
         * < The name of the job if jobId is undefined (<=0)
         */
        public String jobName;

        /**
         * < The message index. A job can have more than one message. Use msgIdx in an array to index messages.
         */
        public int msgIdx;

        /**
         * < Text description of the msg
         */
        public String desc;

        /**
         * < The userId of the author of the message.
         */
        public int userId;

        /**
         * < The size of the data file. If no data file is attached, the size is 0.
         */
        public NativeLong dataSize;

        /**
         * < The time the author posted the message.
         */
        public NativeLong postTime;

        /**
         * < The author of the message.
         */
        public String userName;
    }



    /**
     *  \addtogroup ext_data_status ext_data_status
     */

    /**
     * < Transferring the message's data file.
     */
    public static final int EXT_DATA_UNKNOWN = 0;

    /**
     * < The message does not have an attached  data file.
     */
    public static final int EXT_DATA_NOEXIST = 1;

    /**
     * < The message's data file is available.
     */
    public static final int EXT_DATA_AVAIL = 2;

    /**
     * < The message's data file is corrupt.
     */
    public static final int EXT_DATA_UNAVAIL = 3;

    /**
     * \brief structure jobExternalMsgReply contains the information required to
     * define an external message reply.
     */
    public static class jobExternalMsgReply extends Structure {
        public static class ByReference extends jobExternalMsgReply implements Structure.ByReference {}
        public static class ByValue extends jobExternalMsgReply implements Structure.ByValue {}
        public jobExternalMsgReply() {}
        public jobExternalMsgReply(Pointer p) { super(p); read(); }


        /**
         * < The system generated job Id of the job associated with the message.
         */
        public long jobId;

        /**
         * < The message index. A job can have more than one message. Use msgIdx in an array to index messages.
         */
        public int msgIdx;

        /**
         * < The message you want to read.
         */
        public String desc;

        /**
         * < The user Id of the author of the message.
         */
        public int userId;

        /**
         * < The size of the data file attached. If no data file is attached, the size is 0.
         */
        public NativeLong dataSize;

        /**
         * < The time the message was posted.
         */
        public NativeLong postTime;

        /**
         * < The status of the attached data file.  The status of the data file can be one of the following:\n EXT_DATA_UNKNOWN\n Transferring the message's data file.\n EXT_DATA_NOEXIST\n The message does not have an attached data file.\n EXT_DATA_AVAIL\n The message's data file is available. \n EXT_DATA_UNAVAIL\n The message's data file is corrupt.
         */
        public int dataStatus;

        /**
         * < The author of the msg
         */
        public String userName;
    }




    /**
     * Data structures representing the symphony job status update request.
     */
    public static class symJobInfo extends Structure {
        public static class ByReference extends symJobInfo implements Structure.ByReference {}
        public static class ByValue extends symJobInfo implements Structure.ByValue {}
        public symJobInfo() {}
        public symJobInfo(Pointer p) { super(p); read(); }


/* the service parititon that SSM works for */
        public String partition;

/* the priority of the symphony job */
        public int priority;

/* the full name that indicates the job relationship */
        public String jobFullName;

/* the auxiliary description to help updating command info */
        public String auxCmdDesc;

/* the auxiliary description to help updating job description info */
        public String auxJobDesc;
    }



    public static class symJobStatus extends Structure {
        public static class ByReference extends symJobStatus implements Structure.ByReference {}
        public static class ByValue extends symJobStatus implements Structure.ByValue {}
        public symJobStatus() {}
        public symJobStatus(Pointer p) { super(p); read(); }


/* text description of the symphony job status */
        public String desc;
    }



    public static class symJobProgress extends Structure {
        public static class ByReference extends symJobProgress implements Structure.ByReference {}
        public static class ByValue extends symJobProgress implements Structure.ByValue {}
        public symJobProgress() {}
        public symJobProgress(Pointer p) { super(p); read(); }


/* text description of the symphony job progress */
        public String desc;
    }




    public static class symJobStatusUpdateReq extends Structure {
        public static class ByReference extends symJobStatusUpdateReq implements Structure.ByReference {}
        public static class ByValue extends symJobStatusUpdateReq implements Structure.ByValue {}
        public symJobStatusUpdateReq() {}
        public symJobStatusUpdateReq(Pointer p) { super(p); read(); }


/* the job to be update info into MBD */
        public long jobId;

        public static final int SYM_JOB_UPDATE_NONE = 0x0;
        public static final int SYM_JOB_UPDATE_INFO = 0x1;
        public static final int SYM_JOB_UPDATE_STATUS = 0x2;
        public static final int SYM_JOB_UPDATE_PROGRESS = 0x4;

/* the option to update the info */
        public int bitOption;
        public symJobInfo info;
        public int numOfJobStatus;
        public Pointer /* symJobStatus.ByReference */ status;
        public symJobProgress progress;
    }



    public static class symJobStatusUpdateReqArray extends Structure {
        public static class ByReference extends symJobStatusUpdateReqArray implements Structure.ByReference {}
        public static class ByValue extends symJobStatusUpdateReqArray implements Structure.ByValue {}
        public symJobStatusUpdateReqArray() {}
        public symJobStatusUpdateReqArray(Pointer p) { super(p); read(); }

        public int numOfJobReq;
        public Pointer /* symJobStatusUpdateReq.ByReference */ symJobReqs;
    }




    /**
     * Data structures representing the symphony job status update reply.
     */

    public static class symJobUpdateAck extends Structure {
        public static class ByReference extends symJobUpdateAck implements Structure.ByReference {}
        public static class ByValue extends symJobUpdateAck implements Structure.ByValue {}
        public symJobUpdateAck() {}
        public symJobUpdateAck(Pointer p) { super(p); read(); }

        public static int SYM_UPDATE_ACK_OK = 0;
        public static final int SYM_UPDATE_ACK_ERR = 1;
        public int ackCode;

/* text description of job info update acknowledgement */
        public String desc;
    }



    public static class symJobStatusUpdateReply extends Structure {
        public static class ByReference extends symJobStatusUpdateReply implements Structure.ByReference {}
        public static class ByValue extends symJobStatusUpdateReply implements Structure.ByValue {}
        public symJobStatusUpdateReply() {}
        public symJobStatusUpdateReply(Pointer p) { super(p); read(); }


/* the job to be update info into MBD */
        public long jobId;
        public static final int SYM_UPDATE_INFO_IDX = 0;
        public static final int SYM_UPDATE_STATUS_IDX = 1;
        public static final int SYM_UPDATE_PROGRESS_IDX = 2;
        public static final int NUM_SYM_UPDATE_ACK = 3;
        public symJobUpdateAck[] acks = new symJobUpdateAck[NUM_SYM_UPDATE_ACK];
    }



    public static class symJobStatusUpdateReplyArray extends Structure {
        public static class ByReference extends symJobStatusUpdateReplyArray implements Structure.ByReference {}
        public static class ByValue extends symJobStatusUpdateReplyArray implements Structure.ByValue {}
        public symJobStatusUpdateReplyArray() {}
        public symJobStatusUpdateReplyArray(Pointer p) { super(p); read(); }

        public int numOfJobReply;
        public Pointer /* symJobStatusUpdateReply.ByReference */ symJobReplys;
    }




/* Data structure representing the job array requeue operation.
*  o jobId is the Lsbatch id of the job array to be requeued
*  o status is the desired requeue status of the job, by default
*    it is JOB_STAT_PEND, or user specified JOB_STAT_PSUSP
*  o options specifies the status of the array elements that have
*    to be requeued.
*
*  The function that operates on the data is lsb_requeuejob()
 */

    /**
     *  \addtogroup requeuejob_options requeuejob_options
     *  define statements used by \ref lsb_requeuejob.
     */

    /**
     * < Requeues jobs that have finished running. Jobs that have exited are not re-run. Equivalent to brequeue -d command line option.
     */
    public static final int REQUEUE_DONE = 0x1;

    /**
     * < Requeues jobs that have exited. Finished jobs are not re-run. Equivalent to brequeue -e command line option.
     */
    public static final int REQUEUE_EXIT = 0x2;

    /**
     * < Requeues running jobs and puts them in PEND state. Equivalent to brequeue -r command line option.
     */
    public static final int REQUEUE_RUN = 0x4;

    /**
     * \brief  requeued job
     */
    public static class jobrequeue extends Structure {
        public static class ByReference extends jobrequeue implements Structure.ByReference {}
        public static class ByValue extends jobrequeue implements Structure.ByValue {}
        public jobrequeue() {}
        public jobrequeue(Pointer p) { super(p); read(); }


        /**
         * < Specifies the jobid of a single job or an array of jobs.
         */
        public long jobId;

        /**
         * < Specifies the lsbatch status of the requeued job after it has been requeued. The job status can be JOB_STAT_PEND or JOB_STATE_PSUSP. The default status is JOB_STAT_PEND.
         */
        public int status;

        /**
         * < Specifies the array elements to be requeued.  see \ref requeuejob_options
         */
        public int options;
    }



    public static class requeueEStruct extends Structure {
        public static class ByReference extends requeueEStruct implements Structure.ByReference {}
        public static class ByValue extends requeueEStruct implements Structure.ByValue {}
        public requeueEStruct() {}
        public requeueEStruct(Pointer p) { super(p); read(); }


/* requeue type: normal, exclude, other, prefer_other, etc. */
        public int type;

/* requeue type: normal - as in 2.2 */
        public static final int RQE_NORMAL = 0;

/* requeue type: exclude */
        public static final int RQE_EXCLUDE = 1;

/* indicate the end of the list */
        public static final int RQE_END = 255;

/* requeue exit value */
        public int value;

/* requeue interval */
        public int interval;
    }



    public static class requeue extends Structure {
        public static class ByReference extends requeue implements Structure.ByReference {}
        public static class ByValue extends requeue implements Structure.ByValue {}
        public requeue() {}
        public requeue(Pointer p) { super(p); read(); }

        public int numReqValues;
        public Pointer /* requeueEStruct.ByReference */ reqValues;
    }



/* The Service Level Agreement in LSF
 */


/* This is the library representation of the
*  service class.
 */

    public static class serviceClass extends Structure {
        public static class ByReference extends serviceClass implements Structure.ByReference {}
        public static class ByValue extends serviceClass implements Structure.ByValue {}
        public serviceClass() {}
        public serviceClass(Pointer p) { super(p); read(); }


/* SLA name */
        public String name;

/* SLA priority */
        public float priority;

/* The number of goals */
        public int ngoals;

/* The array of goals */
        public Pointer /* objective.ByReference */ goals;

/* Users allowed to use the SLA */
        public String userGroups;

/* SLA description */
        public String description;

/* SLA control action */
        public String controlAction;

/* Finished jobs per CLEAN_PERIOD */
        public float throughput;

/* Job counters */
        public int[] counters = new int[NUM_JGRP_COUNTERS + 1];

/* project scheduling enabled sla */
        public String consumer;

/* SLA EGO control parameters */
        public slaControl.ByReference ctrl;

/* SLA EGO control parameters */
        public slaControlExt.ByReference ctrlExt;
    }



/* This is the library representation of the
*  Service Level Objective.
 */

    public static final int GOAL_WINDOW_OPEN = 0x1;
    public static final int GOAL_WINDOW_CLOSED = 0x2;
    public static final int GOAL_ONTIME = 0x4;
    public static final int GOAL_DELAYED = 0x8;
    public static final int GOAL_DISABLED = 0x10;

/* Enumerate all the possible performance goals
*  for a service class.
 */

    public static interface objectives {
        public static int GOAL_DEADLINE = 0;
        public static int GOAL_VELOCITY = 1;
        public static int GOAL_THROUGHPUT = 2;
    }



/* The objective of a goal, also called SLO, is represented
*  by this data structure.
 */

    public static class objective extends Structure {
        public static class ByReference extends objective implements Structure.ByReference {}
        public static class ByValue extends objective implements Structure.ByValue {}
        public objective() {}
        public objective(Pointer p) { super(p); read(); }


/* goal specs from lsb.serviceclasses */
        public String spec;

/* goal type */
        public int type;

/* the state of the goal OnTime || Delayed */
        public int state;

/* the configured value */
        public int goal;

/* the actual value */
        public int actual;

/* the optimum value */
        public int optimum;

/* the minimum value */
        public int minimum;
    }



/* Control parameters for SLA management of hosts belonging
*  to the EGO cluster. The control parameters are for each
*  SLA that gets its hosts from EGO.
 */

    public static class slaControl extends Structure {
        public static class ByReference extends slaControl implements Structure.ByReference {}
        public static class ByValue extends slaControl implements Structure.ByValue {}
        public slaControl() {}
        public slaControl(Pointer p) { super(p); read(); }


/* sla name */
        public String sla;

/* EGO consumer the sla is mapped to */
        public String consumer;

/* timeout for returning hosts to EGO */
        public int maxHostIdleTime;

/* timeout left before EGO forcefully reclaims */
        public int recallTimeout;

/* number of hosts beign recalled */
        public int numHostRecalled;

/* EGO resource requirement */
        public String egoResReq;
    }



    public static class slaControlExt extends Structure {
        public static class ByReference extends slaControlExt implements Structure.ByReference {}
        public static class ByValue extends slaControlExt implements Structure.ByValue {}
        public slaControlExt() {}
        public slaControlExt(Pointer p) { super(p); read(); }


/* whether exclusive allocation */
        public int allocflags;

/* tile parameter */
        public int tile;
    }



/* Application Encapsulation in LSF
*
*  This is the library representation of the
*  application.
 */

    public static class appInfoEnt extends Structure {
        public static class ByReference extends appInfoEnt implements Structure.ByReference {}
        public static class ByValue extends appInfoEnt implements Structure.ByValue {}
        public appInfoEnt() {}
        public appInfoEnt(Pointer p) { super(p); read(); }


/* app name */
        public String name;

/* app description */
        public String description;

/* num of total jobs */
        public int numJobs;

/* num of pending slots */
        public int numPEND;

/* num of running slots */
        public int numRUN;

/* num of suspend slots */
        public int numSSUSP;

/* num of ususp slots */
        public int numUSUSP;

/* reserved job slots */
        public int numRESERVE;

/* app attributes */
        public int aAttrib;

/* number of jobs in one chunk */
        public int chunkJobSize;

/* requeue exit values */
        public String requeueEValues;

/* success exit values */
        public String successEValues;

/* app pre execution */
        public String preCmd;

/* app post execution */
        public String postCmd;

/* Job starter command(s) */
        public String jobStarter;

/* suspend action command */
        public String suspendActCmd;

/* resume action command */
        public String resumeActCmd;

/* terimate action command */
        public String terminateActCmd;

/*memory limit level type */
        public int memLimitType;

/* LSF resource limits (soft)*/
        public int[] defLimits = new int[LibLsf.LSF_RLIM_NLIMITS];

/* host spec from CPULIMIT or  RUNLIMIT */
        public String hostSpec;

/* resource requirement string */
        public String resReq;

/* maximal processor limit */
        public int maxProcLimit;

/* default processor limit */
        public int defProcLimit;

/* minimal processor limit */
        public int minProcLimit;

/* estimated run time */
        public int runTime;

/* include postproc as part of job */
        public int jobIncludePostProc;

/* time window for postproc */
        public int jobPostProcTimeOut;

/* remote task gone action */
        public String rTaskGoneAction;

/* pathname of pjob env script */
        public String djobEnvScript;

/* DJOB rusage interval */
        public int djobRuInterval;

/* DJOB heartbeat interval */
        public int djobHbInterval;

/* DJOB communication fail action */
        public String djobCommfailAction;

/* disable Distributed Application Framework */
        public int djobDisabled;

/* grace period (in seconds) before terminating tasks when a job shrinks*/
        public int djobResizeGracePeriod;

/* chkpnt directory */
        public String chkpntDir;

/* chlpnt method */
        public String chkpntMethod;

/* chkpnt period */
        public int chkpntPeriod;

/* initial chkpnt period */
        public int initChkpntPeriod;

/* migration  threshold */
        public int migThreshold;

/* maximum number of job preempted times */
        public int maxJobPreempt;

/* maximum number of pre-exec retry times */
        public int maxPreExecRetry;

/* maximum number of pre-exec retry times for local cluster */
        public int localMaxPreExecRetry;

/* maximum number of job re-queue times */
        public int maxJobRequeue;

/* no preempt run time */
        public int noPreemptRunTime;

/* no preempt finish time */
        public int noPreemptFinishTime;

/* no preempt run time percent */
        public int noPreemptRunTimePercent;

/* no preempt finish time percent */
        public int noPreemptFinishTimePercent;

/* use Linux-PAM */
        public int usePam;

/* processor binding options */
        public int bindingOption;

/* persistent same hosts and same order */
        public int persistHostOrder;

/* job resize notification cmd */
        public String resizeNotifyCmd;
    }



/* application attributes
 */

/* rerunnable application */
    public static final int A_ATTRIB_RERUNNABLE = 0x01;

/* non rerunnable application */
    public static final int A_ATTRIB_NONRERUNNABLE = 0x02;

/* default application */
    public static final int A_ATTRIB_DEFAULT = 0x04;

/* runtime is absolute */
    public static final int A_ATTRIB_ABS_RUNLIMIT = 0x08;

/* process binding application */
    public static final int A_ATTRIB_JOBBINDING = 0x10;

/* process binding application */
    public static final int A_ATTRIB_NONJOBBINDING = 0x20;

/* checkpointable application */
    public static final int A_ATTRIB_CHKPNT = 0x40;

/* Job can be resizable manually */
    public static final int A_ATTRIB_RESIZABLE = 0x80;

/* Job can be resized automatically */
    public static final int A_ATTRIB_AUTO_RESIZABLE = 0x100;


/* processor binding options */
    public static final int BINDING_OPTION_BALANCE = 0x1;
    public static final int BINDING_OPTION_PACK = 0x2;
    public static final int BINDING_OPTION_ANY = 0x4;
    public static final int BINDING_OPTION_USER = 0x8;
    public static final int BINDING_OPTION_USER_CPU_LIST = 0x10;
    public static final int BINDING_OPTION_NONE = 0x20;

    /**
     *  \addtogroup movejob_options movejob_options
     *  options for \ref lsb_movejob call
     */

    /**
     * <  To top
     */
    public static final int TO_TOP = 1;

    /**
     * <  To bottom
     */
    public static final int TO_BOTTOM = 2;

    /**
     *  \addtogroup queue_ctrl_option queue_ctrl_option
     *  options for \ref lsb_queuecontrol call
     */

    /**
     * < Open the queue to accept jobs.
     */
    public static final int QUEUE_OPEN = 1;

    /**
     * < Close the queue so it will not accept jobs.
     */
    public static final int QUEUE_CLOSED = 2;

    /**
     * < Activate the queue to dispatch jobs.
     */
    public static final int QUEUE_ACTIVATE = 3;

    /**
     * < Inactivate the queue so it will not dispatch jobs.
     */
    public static final int QUEUE_INACTIVATE = 4;

    /**
     * < Clean the queue
     */
    public static final int QUEUE_CLEAN = 5;

    /**
     * \brief The structure of queueCtrlReq
     */
    public static class queueCtrlReq extends Structure {
        public static class ByReference extends queueCtrlReq implements Structure.ByReference {}
        public static class ByValue extends queueCtrlReq implements Structure.ByValue {}
        public queueCtrlReq() {}
        public queueCtrlReq(Pointer p) { super(p); read(); }


        /**
         * < The name of the queue to be controlled.
         */
        public String queue;

        /**
         * < Operations to be applied, for example, QUEUE_OPEN. You can refer to \ref queue_ctrl_option for more options.
         */
        public int opCode;

        /**
         * < The message attached by the admin
         */
        public String message;
    }



/* options for lsb_hostcontrol() call */
    /**
     *  \addtogroup host_ctrl_option host_ctrl_option
     *  options operations to be applied
     */

    /**
     * < Opens the host to accept jobs.
     */
    public static final int HOST_OPEN = 1;

    /**
     * < Closes the host so that no jobs can be dispatched to it.
     */
    public static final int HOST_CLOSE = 2;

    /**
     * < Restarts sbatchd on the host. sbatchd will receive a request from mbatchd and re-execute. This permits the sbatchd binary to be updated. This operation fails if no sbatchd is running on the specified host.
     */
    public static final int HOST_REBOOT = 3;

    /**
     * < The sbatchd on the host will exit.
     */
    public static final int HOST_SHUTDOWN = 4;

    /**
     * < Used for closing leased host on the submission cluster
     */
    public static final int HOST_CLOSE_REMOTE = 5;

    /**
     * \brief  Host control request.
     */
    public static class hostCtrlReq extends Structure {
        public static class ByReference extends hostCtrlReq implements Structure.ByReference {}
        public static class ByValue extends hostCtrlReq implements Structure.ByValue {}
        public hostCtrlReq() {}
        public hostCtrlReq(Pointer p) { super(p); read(); }


        /**
         * < The host to be controlled. If host is null, the local host is assumed.
         */
        public String host;

        /**
         * < Operations to be applied in \ref host_ctrl_option.
         */
        public int opCode;

        /**
         * < Message attached by the administrator.
         */
        public String message;
    }



/* options for lsb_hgcontrol() call */
    public static final int HGHOST_ADD = 1;
    public static final int HGHOST_DEL = 2;

    public static class hgCtrlReq extends Structure {
        public static class ByReference extends hgCtrlReq implements Structure.ByReference {}
        public static class ByValue extends hgCtrlReq implements Structure.ByValue {}
        public hgCtrlReq() {}
        public hgCtrlReq(Pointer p) { super(p); read(); }

        public int opCode;
        public String grpname;
        public int numhosts;
        public Pointer hosts;
        public String message;
    }



    public static class hgCtrlReply extends Structure {
        public static class ByReference extends hgCtrlReply implements Structure.ByReference {}
        public static class ByValue extends hgCtrlReply implements Structure.ByValue {}
        public hgCtrlReply() {}
        public hgCtrlReply(Pointer p) { super(p); read(); }

        public int numsucc;
        public int numfail;
        public Pointer succHosts;
        public Pointer failHosts;
        public IntByReference failReasons;
    }



/* options for lsb_reconfig() call */
    /**
     *  \addtogroup mbd_operation mbd_operation
     *   options for \ref lsb_reconfig call
     */

    /**
     * < mbatchd restart
     */
    public static final int MBD_RESTART = 0;

    /**
     * < mbatchd reread configuration files
     */
    public static final int MBD_RECONFIG = 1;

    /**
     * < mbatchd check validity of configuration files
     */
    public static final int MBD_CKCONFIG = 2;

    /**
     * \brief  mbatchd control request.
     */
    public static class mbdCtrlReq extends Structure {
        public static class ByReference extends mbdCtrlReq implements Structure.ByReference {}
        public static class ByValue extends mbdCtrlReq implements Structure.ByValue {}
        public mbdCtrlReq() {}
        public mbdCtrlReq(Pointer p) { super(p); read(); }


        /**
         * < Operation applied, defined in \ref mbd_operation
         */
        public int opCode;

        /**
         * < Not used so far
         */
        public String name;

        /**
         * < The message attached by the admin
         */
        public String message;
    }



/* opcode for turn on or off the perfmon monitor */
    public static final int PERFMON_START = 1;
    public static final int PERFMON_STOP = 2;
    public static final int PERFMON_SET_PERIOD = 3;


/* defualt sample period 60 */
    public static final int DEF_PERFMON_PERIOD = 60;


    public static class perfmonMetricsEnt extends Structure {
        public static class ByReference extends perfmonMetricsEnt implements Structure.ByReference {}
        public static class ByValue extends perfmonMetricsEnt implements Structure.ByValue {}
        public perfmonMetricsEnt() {}
        public perfmonMetricsEnt(Pointer p) { super(p); read(); }

/* metrice name */
        public String name;

/* last period counters */
        public NativeLong current;

/* max of (counter/interval)*sample period for one period */
        public NativeLong max;

/* min of (counter/interval)*sample period for one period */
        public NativeLong min;

/* avg of (total/interval)*sample period for one period */
        public NativeLong avg;

/* total counters from performance monitor turn on */
        public String total;
    }



/*performance monitor info*/

    public static class perfmonInfo extends Structure {
        public static class ByReference extends perfmonInfo implements Structure.ByReference {}
        public static class ByValue extends perfmonInfo implements Structure.ByValue {}
        public perfmonInfo() {}
        public perfmonInfo(Pointer p) { super(p); read(); }

/* number of metrics*/
        public int num;

/* array of metrics counter */
        public Pointer /* perfmonMetricsEnt.ByReference */ record;

/* sample period */
        public int period;

/* time when the performance moniter turn on */
        public NativeLong start;

/* time when the performance moniter turn off */
        public NativeLong end;
    }



/* options for lsb_reljgrp() call */
    public static final int JGRP_RELEASE_PARENTONLY = 0x01;


    /**
     * \brief Records of logged events
     */
    public static class logSwitchLog extends Structure {
        public static class ByReference extends logSwitchLog implements Structure.ByReference {}
        public static class ByValue extends logSwitchLog implements Structure.ByValue {}
        public logSwitchLog() {}
        public logSwitchLog(Pointer p) { super(p); read(); }


        /**
         * < The last jobId so far
         */
        public int lastJobId;
/*#if defined(LSF_SIMULATOR)*/

/**< last trace record time */
/*    public NativeLong lastTraceTime;*/

        /**< last trace record type */
/*public int    lastTraceType;*

    /**< last trace record info */
/*public String lastTraceInfo;*/
        /*#endif*/
    }



    /**
     * \brief Records of job CPU data logged event
     */
    public static class dataLoggingLog extends Structure {
        public static class ByReference extends dataLoggingLog implements Structure.ByReference {}
        public static class ByValue extends dataLoggingLog implements Structure.ByValue {}
        public dataLoggingLog() {}
        public dataLoggingLog(Pointer p) { super(p); read(); }


        /**
         * < The time of last job cpu data logging
         */
        public NativeLong loggingTime;
    }



    /**
     * \brief  new job group log.
     */
    public static class jgrpNewLog extends Structure {
        public static class ByReference extends jgrpNewLog implements Structure.ByReference {}
        public static class ByValue extends jgrpNewLog implements Structure.ByValue {}
        public jgrpNewLog() {}
        public jgrpNewLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The job submission time
         */
        public NativeLong submitTime;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < The job dependency condition
         */
        public String depCond;

        /**
         * < Time event string
         */
        public String timeEvent;

        /**
         * < Job group name
         */
        public String groupSpec;

        /**
         * < New job group name
         */
        public String destSpec;

        /**
         * < Delete options in options field
         */
        public int delOptions;

        /**
         * < Extended Delete options in options2 field
         */
        public int delOptions2;

        /**
         * < Platform type: such as Unix, Windows
         */
        public int fromPlatform;

        /**
         * < SLA service class name under which the job runs
         */
        public String sla;

        /**
         * < Max job group slots limit
         */
        public int maxJLimit;

        /**
         * < Job group creation method: implicit or explicit
         */
        public int options;
    }



    /**
     * \brief  job group control log.
     */
    public static class jgrpCtrlLog extends Structure {
        public static class ByReference extends jgrpCtrlLog implements Structure.ByReference {}
        public static class ByValue extends jgrpCtrlLog implements Structure.ByValue {}
        public jgrpCtrlLog() {}
        public jgrpCtrlLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Job group name
         */
        public String groupSpec;

        /**
         * < Options
         */
        public int options;

        /**
         * < Job control JGRP_RELEASE, JGRP_HOLD, JGRP_DEL
         */
        public int ctrlOp;
    }



    /**
     * \brief  job group status log.
     */
    public static class jgrpStatusLog extends Structure {
        public static class ByReference extends jgrpStatusLog implements Structure.ByReference {}
        public static class ByValue extends jgrpStatusLog implements Structure.ByValue {}
        public jgrpStatusLog() {}
        public jgrpStatusLog(Pointer p) { super(p); read(); }


        /**
         * < The full group path name for the job group
         */
        public String groupSpec;

        /**
         * < Job group status
         */
        public int status;

        /**
         * < Prior status
         */
        public int oldStatus;
    }



    /**
     * \brief jobNewLog logged in lsb.events when a job is submitted.
     */
    public static class jobNewLog extends Structure {
        public static class ByReference extends jobNewLog implements Structure.ByReference {}
        public static class ByValue extends jobNewLog implements Structure.ByValue {}
        public jobNewLog() {}
        public jobNewLog(Pointer p) { super(p); read(); }


        /**
         * < The job ID that the LSF assigned to the job
         */
        public int jobId;

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Job submission options. see \ref lsb_submit.
         */
        public int options;

        /**
         * < Job submission options. see \ref lsb_submit.
         */
        public int options2;

        /**
         * < The number of processors requested for execution
         */
        public int numProcessors;

        /**
         * < The job submission time
         */
        public NativeLong submitTime;

        /**
         * < The job should be started on or after this time
         */
        public NativeLong beginTime;

        /**
         * < If the job has not finished by this time, it will be killed
         */
        public NativeLong termTime;

        /**
         * < The signal value sent to the job 10 minutes before its run window closes
         */
        public int sigValue;

        /**
         * < The checkpointing period
         */
        public int chkpntPeriod;

        /**
         * < The process ID assigned to the job when it was restarted
         */
        public int restartPid;

        /**
         * < The user's resource limits
         */
        public int[] rLimits = new int[LibLsf.LSF_RLIM_NLIMITS];

        /**
         * < The model, host name or host type for scaling CPULIMIT and RUNLIMIT
         */
        public byte[] hostSpec = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The CPU factor for the above model, host name or host type
         */
        public float hostFactor;

        /**
         * < The file creation mask for this job
         */
        public int umask;

        /**
         * < The name of the queue to which this job was submitted
         */
        public byte[] queue = new byte[MAX_LSB_NAME_LEN];

        /**
         * < The resource requirements of the job
         */
        public String resReq;

        /**
         * < The submission host name
         */
        public byte[] fromHost = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The current working directory
         */
        public String cwd;

        /**
         * < The checkpoint directory
         */
        public String chkpntDir;

        /**
         * < The input file name
         */
        public String inFile;

        /**
         * < The output file name
         */
        public String outFile;

        /**
         * < The error output file name
         */
        public String errFile;

        /**
         * < Job spool input file
         */
        public String inFileSpool;

        /**
         * < Job spool command file
         */
        public String commandSpool;

        /**
         * < Job spool directory
         */
        public String jobSpoolDir;

        /**
         * < The home directory of the submitter
         */
        public String subHomeDir;

        /**
         * < The job file name
         */
        public String jobFile;

        /**
         * < The number of hosts considered for dispatching this job
         */
        public int numAskedHosts;

        /**
         * < The array of names of hosts considered for dispatching this job
         */
        public Pointer askedHosts;

        /**
         * < The job dependency condition
         */
        public String dependCond;

        /**
         * < Time event string
         */
        public String timeEvent;

        /**
         * < The job name
         */
        public String jobName;

        /**
         * < The job command
         */
        public String command;

        /**
         * < The number of files to transfer
         */
        public int nxf;

        /**
         * < The array of file transfer specifications. (The xFile structure is defined in <lsf/lsbatch.h>)
         */
        public Pointer /* xFile.ByReference */ xf;

        /**
         * < The command string to be pre_executed
         */
        public String preExecCmd;

        /**
         * < User option mail string
         */
        public String mailUser;

        /**
         * < The project name for this job, used for accounting purposes
         */
        public String projectName;

        /**
         * < Port to be used for interactive jobs
         */
        public int niosPort;

        /**
         * < Maximum number of processors
         */
        public int maxNumProcessors;

        /**
         * < Execution host type
         */
        public String schedHostType;

        /**
         * < Login shell specified by user
         */
        public String loginShell;

        /**
         * < The user group name for this job
         */
        public String userGroup;

        /**
         * < List of alarm conditions for job
         */
        public String exceptList;

        /**
         * < Array idx, must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < User priority
         */
        public int userPriority;

        /**
         * < Advance reservation ID
         */
        public String rsvId;

        /**
         * < The job group under which the job runs.
         */
        public String jobGroup;

        /**
         * < External scheduling options
         */
        public String extsched;

        /**
         * < Warning time period in seconds, -1 if unspecified
         */
        public int warningTimePeriod;

        /**
         * < Warning action, SIGNAL | CHKPNT | command, null if unspecified
         */
        public String warningAction;

        /**
         * < The service class under which the job runs.
         */
        public String sla;

        /**
         * < The absolute run limit of the job
         */
        public int SLArunLimit;

        /**
         * < License Project
         */
        public String licenseProject;

        /**
         * < Extended bitwise inclusive OR of options flags. See \ref lsb_submit.
         */
        public int options3;

        /**
         * < Application profile under which the job runs.
         */
        public String app;

        /**
         * < Post-execution commands.
         */
        public String postExecCmd;

        /**
         * < Runtime estimate specified.
         */
        public int runtimeEstimation;

        /**
         * < Job-level requeue exit values.
         */
        public String requeueEValues;

        /**
         * < Initial checkpoint period
         */
        public int initChkpntPeriod;

        /**
         * < Job migration threshold.
         */
        public int migThreshold;

        /**
         * < Resize notify command
         */
        public String notifyCmd;

        /**
         * < Job description.
         */
        public String jobDescription;

        /**
         * < For new options in future
         */
        public submit_ext.ByReference submitExt;

/*#if defined(LSF_SIMULATOR)*/

/**< maximum memory */
        /*public int    maxmem;*/

        /**< exit status */
        /*public int    exitstatus;*/

        /**< job run time */
        /*public int    runtime;*/

        /**< system cpu time */
        /*public int    cputime;*/

        /**< allocated slots */
        /*public int    slots;*/

        /**< cpu factor */
        /*public float  cpufactor;*/

        /*#endif*/
    }



/*
#if defined(LSF_SIMULATOR)
public static class jobArrayElementLog extends Structure {
public static class ByReference extends jobArrayElementLog implements Structure.ByReference {}
public static class ByValue extends jobArrayElementLog implements Structure.ByValue {}

    public int jobId;
*/
/* Copy LSF simulator related fields from jobNewLog */
/*
    public int idx;
        public int maxmem;
        public int exitstatus;
        public int runtime;
        public int cputime;
        public int slots;
        public float cpufactor;
    };
    #endif
    */

    /**
     * \brief  job modified log.
     */
    public static class jobModLog extends Structure {
        public static class ByReference extends jobModLog implements Structure.ByReference {}
        public static class ByValue extends jobModLog implements Structure.ByValue {}
        public jobModLog() {}
        public jobModLog(Pointer p) { super(p); read(); }


        /**
         * < JobId or jobName in String/
         * public String jobIdStr;
         * <p/>
         * /**< Job submission options(See \ref lsb_submit)
         */
        public int options;

        /**
         * < Job submission options(See \ref lsb_submit)
         */
        public int options2;

        /**
         * < Delete options in options field
         */
        public int delOptions;

        /**
         * < Extended delete options in options2 field .
         */
        public int delOptions2;

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the submitter
         */
        public String userName;

        /**
         * < The job submission time
         */
        public int submitTime;

        /**
         * < The file creation mask for this job
         */
        public int umask;

        /**
         * < The number of processors requested for execution
         */
        public int numProcessors;

        /**
         * < The job should be started on or after this time
         */
        public NativeLong beginTime;

        /**
         * < If the job has not finished by this time,  it will be killed
         */
        public NativeLong termTime;

        /**
         * < The signal value sent to the job 10 minutes before its run window closes
         */
        public int sigValue;

        /**
         * < The process ID assigned to the job when it was restarted
         */
        public int restartPid;


        /**
         * < The job name
         */
        public String jobName;

        /**
         * < The name of the queue to which this job was submitted
         */
        public String queue;


        /**
         * < The number of hosts considered for dispatching this job
         */
        public int numAskedHosts;

        /**
         * < List of asked hosts
         */
        public Pointer askedHosts;


        /**
         * < The resource requirements of the job
         */
        public String resReq;

        /**
         * < User's resource limits (soft)
         */
        public int[] rLimits = new int[LibLsf.LSF_RLIM_NLIMITS];

        /**
         * < The model, host name or host type for scaling CPULIMIT and RUNLIMIT
         */
        public String hostSpec;


        /**
         * < The job dependency condition
         */
        public String dependCond;

        /**
         * < Time event string.
         */
        public String timeEvent;


        /**
         * < The home directory of the submitter
         */
        public String subHomeDir;

        /**
         * < The input file name
         */
        public String inFile;

        /**
         * < The output file name
         */
        public String outFile;

        /**
         * < The error output file name
         */
        public String errFile;

        /**
         * < Command description - this is really a job description field
         */
        public String command;

        /**
         * < Job spool input file
         */
        public String inFileSpool;

        /**
         * < Job spool command file
         */
        public String commandSpool;

        /**
         * < The checkpointing period
         */
        public int chkpntPeriod;

        /**
         * < The checkpoint directory
         */
        public String chkpntDir;

        /**
         * < The number of files to transfer
         */
        public int nxf;

        /**
         * < The array of file transfer specifications.  (The xFile structure is defined in <lsf/lsbatch.h>)
         */
        public Pointer /* xFile.ByReference */ xf;


        /**
         * < The job file name: If == '\\0', indicate let mbatchd make up name, otherwise, mbatchd will use given name.  It is '\\0' if it is a regular job,non-nil means it is a restart job.
         */
        public String jobFile;

        /**
         * < The submission host name
         */
        public String fromHost;

        /**
         * < The current working directory
         */
        public String cwd;


        /**
         * < The pre-execution command
         */
        public String preExecCmd;

        /**
         * < User option mail string
         */
        public String mailUser;

        /**
         * < Project name for the job; used for accounting purposes
         */
        public String projectName;


        /**
         * < NIOS callback port to be used for interactive jobs
         */
        public int niosPort;

        /**
         * < Maximum number of processors
         */
        public int maxNumProcessors;


        /**
         * < The login shell specified by user
         */
        public String loginShell;

        /**
         * < Restart job's submission host type
         */
        public String schedHostType;

        /**
         * < The user group name for this job
         */
        public String userGroup;

        /**
         * < List of job exception conditions
         */
        public String exceptList;

        /**
         * < User priority
         */
        public int userPriority;

        /**
         * < Advance reservation ID
         */
        public String rsvId;

        /**
         * < External scheduling options
         */
        public String extsched;

        /**
         * < Job warning time period in seconds; -1 if unspecified
         */
        public int warningTimePeriod;

        /**
         * < Job warning action: SIGNAL | CHKPNT | command; null if unspecified
         */
        public String warningAction;

        /**
         * < The job group under which the job runs
         */
        public String jobGroup;

        /**
         * < SLA service class name under which the job runs
         */
        public String sla;

        /**
         * < LSF License Scheduler project name
         */
        public String licenseProject;

        /**
         * < Extended bitwise inclusive OR of options flags. see \ref lsb_submit.
         */
        public int options3;

        /**
         * < Extended delete options in options3 field.
         */
        public int delOptions3;

        /**
         * < Application profile under which the job runs.
         */
        public String app;

        /**
         * < Absolute priority scheduling string set by administrators to denote static  system APS value or ADMIN factor APS value.
         */
        public String apsString;

        /**
         * < Post-execution commands.
         */
        public String postExecCmd;

        /**
         * < Runtime estimate.
         */
        public int runtimeEstimation;

        /**
         * < Job-level requeue exit values.
         */
        public String requeueEValues;

        /**
         * < Initial checkpoint period
         */
        public int initChkpntPeriod;

        /**
         * < Job migration threshold.
         */
        public int migThreshold;

        /**
         * < Resize notify command
         */
        public String notifyCmd;

        /**
         * < Job description.
         */
        public String jobDescription;

        /**
         * < For new options in future
         */
        public submit_ext.ByReference submitExt;
    }



    /**
     * \brief  logged in lsb.events when a job is started.
     */
    public static class jobStartLog extends Structure {
        public static class ByReference extends jobStartLog implements Structure.ByReference {}
        public static class ByValue extends jobStartLog implements Structure.ByValue {}
        public jobStartLog() {}
        public jobStartLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < The status of the job (see  \ref lsb_readjobinfo )
         */
        public int jStatus;

        /**
         * < The job process ID
         */
        public int jobPid;

        /**
         * < The job process group ID
         */
        public int jobPGid;

        /**
         * < The CPU factor of the first execution host
         */
        public float hostFactor;

        /**
         * < The number of processors used for execution
         */
        public int numExHosts;

        /**
         * < The array of execution host names
         */
        public Pointer execHosts;

        /**
         * < Pre-execution command defined in the queue
         */
        public String queuePreCmd;

        /**
         * < Post-execution command defined in the queue
         */
        public String queuePostCmd;

        /**
         * < Job processing flags
         */
        public int jFlags;

        /**
         * < The user group name for this job
         */
        public String userGroup;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < Placement information of LSF HPC jobs
         */
        public String additionalInfo;

        /**
         * < How long a backfilled job can run; used for preemption backfill jobs
         */
        public int duration4PreemptBackfill;

        /**
         * <  Job Flags2
         */
        public int jFlags2;
    }



    /**
     * \brief logged in lsb.events when a job start request is accepted.
     */
    public static class jobStartAcceptLog extends Structure {
        public static class ByReference extends jobStartAcceptLog implements Structure.ByReference {}
        public static class ByValue extends jobStartAcceptLog implements Structure.ByValue {}
        public jobStartAcceptLog() {}
        public jobStartAcceptLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < The job process ID
         */
        public int jobPid;

        /**
         * < The job process group ID
         */
        public int jobPGid;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief logged in lsb.events when a job is executed.
     */
    public static class jobExecuteLog extends Structure {
        public static class ByReference extends jobExecuteLog implements Structure.ByReference {}
        public static class ByValue extends jobExecuteLog implements Structure.ByValue {}
        public jobExecuteLog() {}
        public jobExecuteLog(Pointer p) { super(p); read(); }

        /* logged in lsb.events when a job is executed */

        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < User ID under which the job is running
         */
        public int execUid;

        /**
         * < Home directory of the user denoted by execUid
         */
        public String execHome;

        /**
         * < Current working directory where job is running
         */
        public String execCwd;

        /**
         * < The job process group ID
         */
        public int jobPGid;

        /**
         * < User name under which the job is running
         */
        public String execUsername;

        /**
         * < The job process ID
         */
        public int jobPid;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < Placement information of LSF HPC jobs
         */
        public String additionalInfo;

        /**
         * < The run limit scaled by the exec host
         */
        public int SLAscaledRunLimit;

        /**
         * < The position of the job
         */
        public int position;

        /**
         * < The rusage satisfied at job runtime
         */
        public String execRusage;

        /**
         * < The duration for preemptive backfill class in seconds
         */
        public int duration4PreemptBackfill;
    }




    /**
     * \brief logged when a job's status is changed.
     */
    public static class jobStatusLog extends Structure {
        public static class ByReference extends jobStatusLog implements Structure.ByReference {}
        public static class ByValue extends jobStatusLog implements Structure.ByValue {}
        public jobStatusLog() {}
        public jobStatusLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < The job status (see \ref lsb_readjobinfo )
         */
        public int jStatus;

        /**
         * < The reason the job is pending or suspended  (see \ref lsb_pendreason and \ref lsb_suspreason )
         */
        public int reason;

        /**
         * < The load indices that have overloaded the host (see \ref lsb_pendreason  and \ref lsb_suspreason )
         */
        public int subreasons;

        /**
         * < The CPU time consumed before this event occurred
         */
        public float cpuTime;

        /**
         * < The job completion time
         */
        public NativeLong endTime;

        /**
         * < Boolean indicating lsfRusage is logged
         */
        public int ru;

        /**
         * < Resource usage statisticsThe lsfRusage structure is defined in <lsf/lsf.h>. Note that the availability of certain fields depends on the platform on which the sbatchd runs. The fields that do not make sense on the platform will be logged as -1.0.
         */
        public LibLsf.lsfRusage lsfRusage;

        /**
         * < Job exit status
         */
        public int jFlags;

        /**
         * < Job's exit status
         */
        public int exitStatus;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < Job termination reason, see <lsf/lsbatch.h>
         */
        public int exitInfo;
    }




    /**
     * \brief logged when a job's status is changed
     */
    public static class sbdJobStatusLog extends Structure {
        public static class ByReference extends sbdJobStatusLog implements Structure.ByReference {}
        public static class ByValue extends sbdJobStatusLog implements Structure.ByValue {}
        public sbdJobStatusLog() {}
        public sbdJobStatusLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < The status of the job (see \ref lsb_readjobinfo)
         */
        public int jStatus;

        /**
         * < The reason the job is pending or suspended (See \ref lsb_pendreason and \ref lsb_suspreason)
         */
        public int reasons;

        /**
         * < The load indices that have overloaded the host (See \ref lsb_pendreason and \ref lsb_suspreason)
         */
        public int subreasons;

        /**
         * < Action process ID
         */
        public int actPid;

        /**
         * < Action Value SIG_CHKPNT | SIG_CHKPNT_COPY |  SIG_WARNING
         */
        public int actValue;

        /**
         * < Action period
         */
        public NativeLong actPeriod;

        /**
         * < Action flag
         */
        public int actFlags;

        /**
         * < Action logging status
         */
        public int actStatus;

        /**
         * < Action Reason SUSP_MBD_LOCK | SUSP_USER_STOP | SUSP_USER_RESUME | SUSP_SBD_STARTUP
         */
        public int actReasons;

        /**
         * < Sub Reason SUB_REASON_RUNLIMIT | SUB_REASON_DEADLINE |SUB_REASON_PROCESSLIMIT | SUB_REASON_MEMLIMIT |SUB_REASON_CPULIMIT
         */
        public int actSubReasons;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The signal value
         */
        public int sigValue;

        /**
         * < The termination reason of a job
         */
        public int exitInfo;
    }



    /**
     * \brief job status that we could send to MBD
     */
    public static class sbdUnreportedStatusLog extends Structure {
        public static class ByReference extends sbdUnreportedStatusLog implements Structure.ByReference {}
        public static class ByValue extends sbdUnreportedStatusLog implements Structure.ByValue {}
        public sbdUnreportedStatusLog() {}
        public sbdUnreportedStatusLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Action process ID
         */
        public int actPid;

        /**
         * < The job process ID
         */
        public int jobPid;

        /**
         * < The job process group ID
         */
        public int jobPGid;

        /**
         * < New status of the job
         */
        public int newStatus;

        /**
         * < Pending or suspending reason code
         */
        public int reason;

        /**
         * < Pending or suspending subreason code
         */
        public int subreasons;

        /**
         * < Resource usage information for the job  (see jobFinishLog)
         */
        public LibLsf.lsfRusage lsfRusage;

        /**
         * < User ID under which the job is running
         */
        public int execUid;

        /**
         * < Job exit status
         */
        public int exitStatus;

        /**
         * < Current working directory where job is running
         */
        public String execCwd;

        /**
         * < Home directory of the user denoted by execUid
         */
        public String execHome;

        /**
         * < User name under which the job is running
         */
        public String execUsername;

        /**
         * < Message index
         */
        public int msgId;

        /**
         * < Job's resource usage
         */
        public LibLsf.jRusage runRusage;

        /**
         * < Signal value
         */
        public int sigValue;

        /**
         * < Action logging status
         */
        public int actStatus;

        /**
         * < Sequence status of the job
         */
        public int seq;

        /**
         * < Job array index
         */
        public int idx;

        /**
         * < The termination reason of a job
         */
        public int exitInfo;
    }



    /**
     * \brief logged when a job is switched to another queue
     */
    public static class jobSwitchLog extends Structure {
        public static class ByReference extends jobSwitchLog implements Structure.ByReference {}
        public static class ByValue extends jobSwitchLog implements Structure.ByValue {}
        public jobSwitchLog() {}
        public jobSwitchLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < The name of the queue the job has been switched to
         */
        public byte[] queue = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief logged when a job is moved to another position
     */
    public static class jobMoveLog extends Structure {
        public static class ByReference extends jobMoveLog implements Structure.ByReference {}
        public static class ByValue extends jobMoveLog implements Structure.ByValue {}
        public jobMoveLog() {}
        public jobMoveLog(Pointer p) { super(p); read(); }

        /* logged when a job is moved to another position */

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < The new position of the job
         */
        public int position;

        /**
         * < The operation code for the move (see  \ref lsb_movejob)
         */
        public int base;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief  check point log.
     */
    public static class chkpntLog extends Structure {
        public static class ByReference extends chkpntLog implements Structure.ByReference {}
        public static class ByValue extends chkpntLog implements Structure.ByValue {}
        public chkpntLog() {}
        public chkpntLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < The new checkpointing period
         */
        public NativeLong period;

        /**
         * < The process ID of the checkpointing process (a child sbatchd)
         */
        public int pid;

        /**
         * < 0: checkpoint started; 1: checkpoint succeeded
         */
        public int ok;

        /**
         * < One of the following: \n LSB_CHKPNT_KILL : Kill process if checkpoint successful \n LSB_CHKPNT_FORCE : Force checkpoint even if non-checkpointable conditions exist \n LSB_CHKPNT_MIG : Checkpoint for the purpose of migration
         */
        public int flags;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  job requeue log.
     */
    public static class jobRequeueLog extends Structure {
        public static class ByReference extends jobRequeueLog implements Structure.ByReference {}
        public static class ByValue extends jobRequeueLog implements Structure.ByValue {}
        public jobRequeueLog() {}
        public jobRequeueLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  job clean log.
     */
    public static class jobCleanLog extends Structure {
        public static class ByReference extends jobCleanLog implements Structure.ByReference {}
        public static class ByValue extends jobCleanLog implements Structure.ByValue {}
        public jobCleanLog() {}
        public jobCleanLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  job exception log.
     */
    public static class jobExceptionLog extends Structure {
        public static class ByReference extends jobExceptionLog implements Structure.ByReference {}
        public static class ByValue extends jobExceptionLog implements Structure.ByValue {}
        public jobExceptionLog() {}
        public jobExceptionLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Job exception handling mask
         */
        public int exceptMask;

        /**
         * < Action Id (kill | alarm | rerun | setexcept)
         */
        public int actMask;

        /**
         * < Time event string
         */
        public NativeLong timeEvent;

        /**
         * < Except Info, pending reason for missched or cantrun exception, the exit code of thejob for the abend exception, otherwise 0.
         */
        public int exceptInfo;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  signal action log.
     */
    public static class sigactLog extends Structure {
        public static class ByReference extends sigactLog implements Structure.ByReference {}
        public static class ByValue extends sigactLog implements Structure.ByValue {}
        public sigactLog() {}
        public sigactLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < Action period
         */
        public NativeLong period;

        /**
         * < Action process ID
         */
        public int pid;

        /**
         * < Job status
         */
        public int jStatus;

        /**
         * < Pending reasons
         */
        public int reasons;

        /**
         * < Action flag
         */
        public int flags;

        /**
         * < Signal symbol from the set: DELETEJOB |  KILL | KILLREQUEUE |REQUEUE_DONE | REQUEUE_EXIT | REQUEUE_PEND |REQUEUE_PSUSP_ADMIN | REQUEUE_PSUSP_USER | SIG_CHKPNT |  SIG_CHKPNT_COPY
         */
        public String signalSymbol;

        /**
         * < Action logging status (ACT_NO | ACT_START | ACT_PREEMPT | ACT_DONE |  ACT_FAIL) .Shown in signal_action
         */
        public int actStatus;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  migration log.
     */
    public static class migLog extends Structure {
        public static class ByReference extends migLog implements Structure.ByReference {}
        public static class ByValue extends migLog implements Structure.ByValue {}
        public migLog() {}
        public migLog(Pointer p) { super(p); read(); }


        /**
         * < The job to be migrated
         */
        public int jobId;

        /**
         * < The number of candidate hosts for migration
         */
        public int numAskedHosts;

        /**
         * < The array of candidate host names
         */
        public Pointer askedHosts;

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The user name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief  signal log.
     */
    public static class signalLog extends Structure {
        public static class ByReference extends signalLog implements Structure.ByReference {}
        public static class ByValue extends signalLog implements Structure.ByValue {}
        public signalLog() {}
        public signalLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < Signal symbol from the set: DELETEJOB | KILL | KILLREQUEUE |REQUEUE_DONE | REQUEUE_EXIT | REQUEUE_PEND |REQUEUE_PSUSP_ADMIN | REQUEUE_PSUSP_USER | SIG_CHKPNT | SIG_CHKPNT_COPY
         */
        public String signalSymbol;

        /**
         * < The number of running times
         */
        public int runCount;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief logged when bqc command is invoked.
     */
    public static class queueCtrlLog extends Structure {
        public static class ByReference extends queueCtrlLog implements Structure.ByReference {}
        public static class ByValue extends queueCtrlLog implements Structure.ByValue {}
        public queueCtrlLog() {}
        public queueCtrlLog(Pointer p) { super(p); read(); }


        /**
         * < The queue control operation (see \ref lsb_queuecontrol)
         */
        public int opCode;

        /**
         * < The name of the queue
         */
        public byte[] queue = new byte[MAX_LSB_NAME_LEN];

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Queue control message
         */
        public byte[] message = new byte[LibLsf.MAXLINELEN];
    }



/*
*  \brief  new debug log.
 */

    public static class newDebugLog extends Structure {
        public static class ByReference extends newDebugLog implements Structure.ByReference {}
        public static class ByValue extends newDebugLog implements Structure.ByValue {}
        public newDebugLog() {}
        public newDebugLog(Pointer p) { super(p); read(); }


        /**
         * < The queue control operation
         */
        public int opCode;

        /**
         * < Debug level
         */
        public int level;

        /**
         * < Class of log
         */
        public int _logclass;

        /**
         * < Log enabled, disabled
         */
        public int turnOff;

        /**
         * < Name of log file
         */
        public byte[] logFileName = new byte[LibLsf.MAXLSFNAMELEN];

        /**
         * < The user ID of the submitter
         */
        public int userId;
    }



    /**
     * \brief log the host control information.
     */
    public static class hostCtrlLog extends Structure {
        public static class ByReference extends hostCtrlLog implements Structure.ByReference {}
        public static class ByValue extends hostCtrlLog implements Structure.ByValue {}
        public hostCtrlLog() {}
        public hostCtrlLog(Pointer p) { super(p); read(); }


        /**
         * < The host control operation (See  \ref lsb_hostcontrol)
         */
        public int opCode;

        /**
         * < The name of the host
         */
        public byte[] host = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Host control message
         */
        public byte[] message = new byte[LibLsf.MAXLINELEN];
    }



    /**
     * \brief logged when dynamic hosts are added to group.
     */
    public static class hgCtrlLog extends Structure {
        public static class ByReference extends hgCtrlLog implements Structure.ByReference {}
        public static class ByValue extends hgCtrlLog implements Structure.ByValue {}
        public hgCtrlLog() {}
        public hgCtrlLog(Pointer p) { super(p); read(); }


        /**
         * < The host control operation  (see \ref lsb_hostcontrol)
         */
        public int opCode;

        /**
         * < The name of the host
         */
        public byte[] host = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The name of the host group
         */
        public byte[] grpname = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Host group control message
         */
        public byte[] message = new byte[LibLsf.MAXLINELEN];
    }




/* simulator is ready to schedule jobs */
    public static final int SIMU_STATUS_READYSCHEDULE = 0x01;

    /**
     * \brief  mbatchd start log.
     */
    public static class mbdStartLog extends Structure {
        public static class ByReference extends mbdStartLog implements Structure.ByReference {}
        public static class ByValue extends mbdStartLog implements Structure.ByValue {}
        public mbdStartLog() {}
        public mbdStartLog(Pointer p) { super(p); read(); }


        /**
         * < The master host name
         */
        public byte[] master = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The cluster name
         */
        public byte[] cluster = new byte[LibLsf.MAXLSFNAMELEN];

        /**
         * < The number of hosts in the cluster
         */
        public int numHosts;

        /**
         * < The number of queues in the cluster
         */
        public int numQueues;
/*
    public int    simDiffTime;
    public int    pendJobsThreshold;
    public int    simStatus;
*/
    }



    public static class mbdSimStatusLog extends Structure {
        public static class ByReference extends mbdSimStatusLog implements Structure.ByReference {}
        public static class ByValue extends mbdSimStatusLog implements Structure.ByValue {}
        public mbdSimStatusLog() {}
        public mbdSimStatusLog(Pointer p) { super(p); read(); }


/* simulator status */
        public int simStatus;
    }



    /**
     * \brief  mbatchd die log.
     */
    public static class mbdDieLog extends Structure {
        public static class ByReference extends mbdDieLog implements Structure.ByReference {}
        public static class ByValue extends mbdDieLog implements Structure.ByValue {}
        public mbdDieLog() {}
        public mbdDieLog(Pointer p) { super(p); read(); }


        /**
         * < The master host name
         */
        public byte[] master = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < The number of finished jobs that have been removed from the system and logged in the current event file
         */
        public int numRemoveJobs;

        /**
         * < The exit code from the master batch daemon
         */
        public int exitCode;

        /**
         * < mbatchd administrator control message
         */
        public byte[] message = new byte[LibLsf.MAXLINELEN];
    }



    /**
     * \brief logged before mbatchd dies.
     */
    public static class unfulfillLog extends Structure {
        public static class ByReference extends unfulfillLog implements Structure.ByReference {}
        public static class ByValue extends unfulfillLog implements Structure.ByValue {}
        public unfulfillLog() {}
        public unfulfillLog(Pointer p) { super(p); read(); }


        /**
         * < The job ID.
         */
        public int jobId;

        /**
         * < The mbatchd has switched the job to a new queue but the sbatchd has not been informed of the switch
         */
        public int notSwitched;

        /**
         * < This signal was not sent to the job
         */
        public int sig;

        /**
         * < The job was not signaled to checkpoint itself
         */
        public int sig1;

        /**
         * < Checkpoint flags. see the chkpntLog structure below.
         */
        public int sig1Flags;

        /**
         * < The new checkpoint period for the job
         */
        public NativeLong chkPeriod;

        /**
         * < Flag for bmod running job's parameters
         */
        public int notModified;

        /**
         * < Job array index
         */
        public int idx;

        /**
         * < Option flags for pending job signals
         */
        public int miscOpts4PendSig;
    }



    public static final int TERM_UNKNOWN = 0;
    public static final int TERM_PREEMPT = 1;
    public static final int TERM_WINDOW = 2;
    public static final int TERM_LOAD = 3;
    public static final int TERM_OTHER = 4;
    public static final int TERM_RUNLIMIT = 5;
    public static final int TERM_DEADLINE = 6;
    public static final int TERM_PROCESSLIMIT = 7;
    public static final int TERM_FORCE_OWNER = 8;
    public static final int TERM_FORCE_ADMIN = 9;
    public static final int TERM_REQUEUE_OWNER = 10;
    public static final int TERM_REQUEUE_ADMIN = 11;
    public static final int TERM_CPULIMIT = 12;
    public static final int TERM_CHKPNT = 13;
    public static final int TERM_OWNER = 14;
    public static final int TERM_ADMIN = 15;
    public static final int TERM_MEMLIMIT = 16;
    public static final int TERM_EXTERNAL_SIGNAL = 17;
    public static final int TERM_RMS = 18;
    public static final int TERM_ZOMBIE = 19;
    public static final int TERM_SWAP = 20;
    public static final int TERM_THREADLIMIT = 21;
    public static final int TERM_SLURM = 22;
    public static final int TERM_BUCKET_KILL = 23;
    public static final int TERM_CTRL_PID = 24;
    public static final int TERM_CWD_NOTEXIST = 25;

    /**
     * \brief logged in lsb.acct when a job finished.
     */
    public static class jobFinishLog extends Structure {
        public static class ByReference extends jobFinishLog implements Structure.ByReference {}
        public static class ByValue extends jobFinishLog implements Structure.ByValue {}
        public jobFinishLog() {}
        public jobFinishLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The user name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Job submission options (see  \ref lsb_submit)
         */
        public int options;

        /**
         * < The number of processors requested for execution
         */
        public int numProcessors;

        /**
         * < The status of the job (See \ref lsb_readjobinfo)
         */
        public int jStatus;

        /**
         * < Job submission time
         */
        public NativeLong submitTime;

        /**
         * < The job started at or after this time
         */
        public NativeLong beginTime;

        /**
         * < If the job was not finished by this time, it was killed
         */
        public NativeLong termTime;

        /**
         * < Job dispatch time
         */
        public NativeLong startTime;

        /**
         * < The time the job finished
         */
        public NativeLong endTime;

        /**
         * < The name of the queue to which this job was submitted
         */
        public byte[] queue = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Resource requirements
         */
        public String resReq;

        /**
         * < Submission host name
         */
        public byte[] fromHost = new byte[LibLsf.MAXHOSTNAMELEN];

        /**
         * < Current working directory
         */
        public String cwd;

        /**
         * < Input file name
         */
        public String inFile;

        /**
         * < Output file name
         */
        public String outFile;

        /**
         * < Error output file name
         */
        public String errFile;

        /**
         * < Job spool input file
         */
        public String inFileSpool;

        /**
         * < Job spool command file
         */
        public String commandSpool;

        /**
         * < Job file name
         */
        public String jobFile;

        /**
         * < The number of hosts considered for dispatching this job
         */
        public int numAskedHosts;

        /**
         * < The array of names of hosts considered for dispatching this job
         */
        public Pointer askedHosts;

        /**
         * < The CPU factor of the first execution host
         */
        public float hostFactor;

        /**
         * < The number of processors used for execution
         */
        public int numExHosts;

        /**
         * < The array of names of execution hosts
         */
        public Pointer execHosts;

        /**
         * < The total CPU time consumed by the job
         */
        public float cpuTime;

        /**
         * < Job name
         */
        public String jobName;

        /**
         * < Job command
         */
        public String command;

        /**
         * < Resource usage statistics.The lsfRusage structure is defined in <lsf/lsf.h>. Note that the availability of certain fields depends on the platform on which the sbatchd runs. The fields that do not make sense on this platform will be logged as -1.0.
         */
        public LibLsf.lsfRusage lsfRusage;

        /**
         * < The job dependency condition
         */
        public String dependCond;

        /**
         * < Time event string
         */
        public String timeEvent;

        /**
         * < The pre-execution command
         */
        public String preExecCmd;

        /**
         * < Name of the user to whom job related mail was sent
         */
        public String mailUser;

        /**
         * < The project name, used for accounting purposes.
         */
        public String projectName;

        /**
         * < Job's exit status
         */
        public int exitStatus;

        /**
         * < Maximum number of processors specified for the job
         */
        public int maxNumProcessors;

        /**
         * < Login shell specified by user
         */
        public String loginShell;

        /**
         * < Job array index
         */
        public int idx;

        /**
         * < Maximum memory used by job
         */
        public int maxRMem;

        /**
         * < Maximum swap used by job
         */
        public int maxRSwap;

        /**
         * < Advanced reservation ID
         */
        public String rsvId;

        /**
         * < Service class of the job
         */
        public String sla;

        /**
         * < Job exception handling mask
         */
        public int exceptMask;

        /**
         * < Placement information of LSF HPC jobs
         */
        public String additionalInfo;

        /**
         * < Job termination reason, see <lsf/lsbatch.h>
         */
        public int exitInfo;

        /**
         * < Job warning time period in seconds; -1 if unspecified
         */
        public int warningTimePeriod;

        /**
         * < Warning action, SIGNAL | CHKPNT | command, null if unspecified
         */
        public String warningAction;

        /**
         * < SAAP charged for job
         */
        public String chargedSAAP;

        /**
         * < LSF License Scheduler project name
         */
        public String licenseProject;

        /**
         * < Application profile under which the job runs.
         */
        public String app;

        /**
         * < Post-execution commands.
         */
        public String postExecCmd;

        /**
         * < Runtime estimate specified.
         */
        public int runtimeEstimation;

        /**
         * < Job group name
         */
        public String jgroup;

        /**
         * < Option2
         */
        public int options2;

        /**
         * < Job requeue exit values
         */
        public String requeueEValues;

        /**
         * < Resize notify command
         */
        public String notifyCmd;

        /**
         * < Last resize start time
         */
        public NativeLong lastResizeTime;

        /**
         * < Job description.
         */
        public String jobDescription;

        /**
         * < For new options in future
         */
        public submit_ext.ByReference submitExt;
    }




    /**
     * \brief  load index log.
     */

    public static class loadIndexLog extends Structure {
        public static class ByReference extends loadIndexLog implements Structure.ByReference {}
        public static class ByValue extends loadIndexLog implements Structure.ByValue {}
        public loadIndexLog() {}
        public loadIndexLog(Pointer p) { super(p); read(); }


        /**
         * < The number of load indices
         */
        public int nIdx;

        /**
         * < The array of load index names
         */
        public Pointer name;
    }



    /**
     * \brief  calendar log.
     */
    public static class calendarLog extends Structure {
        public static class ByReference extends calendarLog implements Structure.ByReference {}
        public static class ByValue extends calendarLog implements Structure.ByValue {}
        public calendarLog() {}
        public calendarLog(Pointer p) { super(p); read(); }


        /**
         * < Reserved for future use
         */
        public int options;

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The name of the calendar
         */
        public String name;

        /**
         * < Description
         */
        public String desc;

        /**
         * < Calendar expression
         */
        public String calExpr;
    }



    /**
     * \brief  job forward log.
     */
    public static class jobForwardLog extends Structure {
        public static class ByReference extends jobForwardLog implements Structure.ByReference {}
        public static class ByValue extends jobForwardLog implements Structure.ByValue {}
        public jobForwardLog() {}
        public jobForwardLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < The cluster name
         */
        public String cluster;

        /**
         * < Number of Reserved Hosts
         */
        public int numReserHosts;

        /**
         * < Reserved Host Names
         */
        public Pointer reserHosts;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < Remote job attributes from: \n JOB_FORWARD Remote batch job on submission side \n JOB_LEASE Lease job on submission side \n JOB_REMOTE_BATCH Remote batch job on execution side \n JOB_REMOTE_LEASE Lease job on execution side \n JOB_LEASE_RESYNC Lease job resync during restart \n JOB_REMOTE_RERUNNABLE Remote batch job rerunnable on execution cluster
         */
        public int jobRmtAttr;
    }



    /**
     * \brief  job accept log.
     */
    public static class jobAcceptLog extends Structure {
        public static class ByReference extends jobAcceptLog implements Structure.ByReference {}
        public static class ByValue extends jobAcceptLog implements Structure.ByValue {}
        public jobAcceptLog() {}
        public jobAcceptLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < The unique ID of the remote job
         */
        public long remoteJid;

        /**
         * < The cluster name
         */
        public String cluster;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < Remote job attributes from: \n JOB_FORWARD Remote batch job on submission side \n JOB_LEASE Lease job on submission side \n JOB_REMOTE_BATCH Remote batch job on execution side \n JOB_REMOTE_LEASE Lease job on execution side \n JOB_LEASE_RESYNC Lease job resync during restart \n JOB_REMOTE_RERUNNABLE Remote batch job rerunnable on execution cluster
         */
        public int jobRmtAttr;
    }



    /**
     * \brief  status Ack log.
     */
    public static class statusAckLog extends Structure {
        public static class ByReference extends statusAckLog implements Structure.ByReference {}
        public static class ByValue extends statusAckLog implements Structure.ByValue {}
        public statusAckLog() {}
        public statusAckLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobId;

        /**
         * < Line number of Status
         */
        public int statusNum;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  job message log.
     */
    public static class jobMsgLog extends Structure {
        public static class ByReference extends jobMsgLog implements Structure.ByReference {}
        public static class ByValue extends jobMsgLog implements Structure.ByValue {}
        public jobMsgLog() {}
        public jobMsgLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int usrId;

        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Message index
         */
        public int msgId;

        /**
         * < Message type
         */
        public int type;

        /**
         * < Message source
         */
        public String src;

        /**
         * < Message destination
         */
        public String dest;

        /**
         * < Message
         */
        public String msg;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  job message ack log.
     */
    public static class jobMsgAckLog extends Structure {
        public static class ByReference extends jobMsgAckLog implements Structure.ByReference {}
        public static class ByValue extends jobMsgAckLog implements Structure.ByValue {}
        public jobMsgAckLog() {}
        public jobMsgAckLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int usrId;

        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Message index
         */
        public int msgId;

        /**
         * < Message type
         */
        public int type;

        /**
         * < Message source
         */
        public String src;

        /**
         * < Message destination
         */
        public String dest;

        /**
         * < Message
         */
        public String msg;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;
    }



    /**
     * \brief  job occupy request log. jobOccupyReqLog is for future use.
     */
    public static class jobOccupyReqLog extends Structure {
        public static class ByReference extends jobOccupyReqLog implements Structure.ByReference {}
        public static class ByValue extends jobOccupyReqLog implements Structure.ByValue {}
        public jobOccupyReqLog() {}
        public jobOccupyReqLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Number of Jobs Slots desired
         */
        public int numOccupyRequests;

        /**
         * < List of slots occupied
         */
        public Pointer occupyReqList;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief  job vacate log.jobVacatedLog is for future use.
     */
    public static class jobVacatedLog extends Structure {
        public static class ByReference extends jobVacatedLog implements Structure.ByReference {}
        public static class ByValue extends jobVacatedLog implements Structure.ByValue {}
        public jobVacatedLog() {}
        public jobVacatedLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief  job force request log.
     */
    public static class jobForceRequestLog extends Structure {
        public static class ByReference extends jobForceRequestLog implements Structure.ByReference {}
        public static class ByValue extends jobForceRequestLog implements Structure.ByValue {}
        public jobForceRequestLog() {}
        public jobForceRequestLog(Pointer p) { super(p); read(); }


        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < >1 for local/lease jobs; 0 for remote batch model
         */
        public int numExecHosts;

        /**
         * < The array of execution host names
         */
        public Pointer execHosts;

        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < Job run options (RUNJOB_OPT_NOSTOP | JFLAG_URGENT_NOSTOP |JFLAG_URGENT)
         */
        public int options;

        /**
         * < The name of the submitter
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < The name of the queue to which this job was submitted
         */
        public String queue;
    }



    /**
     * \brief  job chunck log.
     */
    public static class jobChunkLog extends Structure {
        public static class ByReference extends jobChunkLog implements Structure.ByReference {}
        public static class ByValue extends jobChunkLog implements Structure.ByValue {}
        public jobChunkLog() {}
        public jobChunkLog(Pointer p) { super(p); read(); }


        /**
         * < Size of array membJobId
         */
        public NativeLong membSize;

        /**
         * < Job ids of jobs in the chunk
         */
        public LongByReference membJobId;

        /**
         * < The number of processors used for execution
         */
        public NativeLong numExHosts;

        /**
         * < The array of names of execution hosts
         */
        public Pointer execHosts;
    }



    /**
     * \brief  job external message log.
     */
    public static class jobExternalMsgLog extends Structure {
        public static class ByReference extends jobExternalMsgLog implements Structure.ByReference {}
        public static class ByValue extends jobExternalMsgLog implements Structure.ByValue {}
        public jobExternalMsgLog() {}
        public jobExternalMsgLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID for the job
         */
        public int jobId;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < The message index
         */
        public int msgIdx;

        /**
         * < Message description
         */
        public String desc;

        /**
         * < The user ID of the submitter
         */
        public int userId;

        /**
         * < Size of the message
         */
        public NativeLong dataSize;

        /**
         * < The time the author posted the message.
         */
        public NativeLong postTime;

        /**
         * < The status of the message
         */
        public int dataStatus;

        /**
         * < Name of attached data file. If no file is attached, use null.
         */
        public String fileName;

        /**
         * < The author of the message
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];
    }



    /**
     * \brief  reservation request.
     */
    public static class rsvRes extends Structure {
        public static class ByReference extends rsvRes implements Structure.ByReference {}
        public static class ByValue extends rsvRes implements Structure.ByValue {}
        public rsvRes() {}
        public rsvRes(Pointer p) { super(p); read(); }


        /**
         * < Name of the resource (currently: host)
         */
        public String resName;

        /**
         * < Reserved counter (currently: cpu number)
         */
        public int count;

        /**
         * < Used of the reserved counter (not used)
         */
        public int usedAmt;
    }



    /**
     * \brief for advanced reservation.
     */
    public static class rsvFinishLog extends Structure {
        public static class ByReference extends rsvFinishLog implements Structure.ByReference {}
        public static class ByValue extends rsvFinishLog implements Structure.ByValue {}
        public rsvFinishLog() {}
        public rsvFinishLog(Pointer p) { super(p); read(); }


        /**
         * < Time when the reservation is required
         */
        public NativeLong rsvReqTime;

        /**
         * < Same as the options field in the addRsvRequest(lsbatch.h)
         */
        public int options;

        /**
         * < The user who creat the reservation
         */
        public int uid;

        /**
         * < Reservation ID
         */
        public String rsvId;

        /**
         * < Client of the reservation
         */
        public String name;

        /**
         * < Number of resources reserved
         */
        public int numReses;

        /**
         * < Allocation vector
         */
        public Pointer /* rsvRes.ByReference */ alloc;

        /**
         * < Time window within which the reservation is active \n Two forms: int1-int2 or [day1]:hour1:0-[day2]:hour2:0
         */
        public String timeWindow;

        /**
         * < Duration in seconds. duration = to - from : when the reservation expired
         */
        public NativeLong duration;

        /**
         * < Creator of the reservation
         */
        public String creator;
    }



    /**
     * \brief  CPU Profile Log
     */
    public static class cpuProfileLog extends Structure {
        public static class ByReference extends cpuProfileLog implements Structure.ByReference {}
        public static class ByValue extends cpuProfileLog implements Structure.ByValue {}
        public cpuProfileLog() {}
        public cpuProfileLog(Pointer p) { super(p); read(); }


        /**
         * < Queue name
         */
        public byte[] servicePartition = new byte[MAX_LSB_NAME_LEN];

        /**
         * < The number of CPU required
         */
        public int slotsRequired;

        /**
         * < The number of CPU actually allocated
         */
        public int slotsAllocated;

        /**
         * < The number of CPU borrowed
         */
        public int slotsBorrowed;

        /**
         * < The number of CPU lent
         */
        public int slotsLent;
        /** note:  the number of CPU reserved = slotsAllocated - slotsBorrowed + slotsLent */
    }



    /**
     * \brief  job resize start notify log.
     */
    public static class jobResizeNotifyStartLog extends Structure {
        public static class ByReference extends jobResizeNotifyStartLog implements Structure.ByReference {}
        public static class ByValue extends jobResizeNotifyStartLog implements Structure.ByValue {}
        public jobResizeNotifyStartLog() {}
        public jobResizeNotifyStartLog(Pointer p) { super(p); read(); }


        /**
         * <  JobId
         */
        public int jobId;

        /**
         * <  Index
         */
        public int idx;

        /**
         * <  Notify Id
         */
        public int notifyId;

        /**
         * <  Number of resized hosts.
         */
        public int numResizeHosts;

        /**
         * <  Resize Hosts
         */
        public Pointer resizeHosts;
    }



    /**
     * \brief  job resize accept notify log.
     */
    public static class jobResizeNotifyAcceptLog extends Structure {
        public static class ByReference extends jobResizeNotifyAcceptLog implements Structure.ByReference {}
        public static class ByValue extends jobResizeNotifyAcceptLog implements Structure.ByValue {}
        public jobResizeNotifyAcceptLog() {}
        public jobResizeNotifyAcceptLog(Pointer p) { super(p); read(); }


        /**
         * <  JobId
         */
        public int jobId;

        /**
         * <  Index
         */
        public int idx;

        /**
         * <  Notify Id
         */
        public int notifyId;

        /**
         * <  Resize Notify command pid
         */
        public int resizeNotifyCmdPid;

        /**
         * <  Resize Notify command pgid
         */
        public int resizeNotifyCmdPGid;

        /**
         * <  Status
         */
        public int status;
    }



    /**
     * \brief  job resize done notify log.
     */
    public static class jobResizeNotifyDoneLog extends Structure {
        public static class ByReference extends jobResizeNotifyDoneLog implements Structure.ByReference {}
        public static class ByValue extends jobResizeNotifyDoneLog implements Structure.ByValue {}
        public jobResizeNotifyDoneLog() {}
        public jobResizeNotifyDoneLog(Pointer p) { super(p); read(); }


        /**
         * <  JobId
         */
        public int jobId;

        /**
         * <  Index
         */
        public int idx;

        /**
         * <  Notify Id
         */
        public int notifyId;

        /**
         * <  Status
         */
        public int status;
    }



    /**
     * \brief  job resize release log.
     */
    public static class jobResizeReleaseLog extends Structure {
        public static class ByReference extends jobResizeReleaseLog implements Structure.ByReference {}
        public static class ByValue extends jobResizeReleaseLog implements Structure.ByValue {}
        public jobResizeReleaseLog() {}
        public jobResizeReleaseLog(Pointer p) { super(p); read(); }


        /**
         * <  JobId
         */
        public int jobId;

        /**
         * <  Index
         */
        public int idx;

        /**
         * <  Request Id
         */
        public int reqId;

        /**
         * <  Options
         */
        public int options;

        /**
         * <  User Id
         */
        public int userId;

        /**
         * <  User Name
         */
        public String userName;

        /**
         * <  Resize Notify command
         */
        public String resizeNotifyCmd;

        /**
         * <  Number of resized hosts
         */
        public int numResizeHosts;

        /**
         * <  Resized hosts
         */
        public Pointer resizeHosts;
    }



    /**
     * \brief  job resize cancel log.
     */
    public static class jobResizeCancelLog extends Structure {
        public static class ByReference extends jobResizeCancelLog implements Structure.ByReference {}
        public static class ByValue extends jobResizeCancelLog implements Structure.ByValue {}
        public jobResizeCancelLog() {}
        public jobResizeCancelLog(Pointer p) { super(p); read(); }


        /**
         * <  JobId
         */
        public int jobId;

        /**
         * <  Index
         */
        public int idx;

        /**
         * <  User Id
         */
        public int userId;

        /**
         * <  User name
         */
        public String userName;
    }



    /**
     * \brief log the running rusage of a job in the lsb.stream file
     */
    public static class jobRunRusageLog extends Structure {
        public static class ByReference extends jobRunRusageLog implements Structure.ByReference {}
        public static class ByValue extends jobRunRusageLog implements Structure.ByValue {}
        public jobRunRusageLog() {}
        public jobRunRusageLog(Pointer p) { super(p); read(); }


        /**
         * < The unique ID of the job
         */
        public int jobid;

        /**
         * < Job array index; must be 0 in JOB_NEW
         */
        public int idx;

        /**
         * < jrusage
         */
        public LibLsf.jRusage jrusage;
    }



    /**
     * \brief  SLA event log.
     */
    public static class slaLog extends Structure {
        public static class ByReference extends slaLog implements Structure.ByReference {}
        public static class ByValue extends slaLog implements Structure.ByValue {}
        public slaLog() {}
        public slaLog(Pointer p) { super(p); read(); }


        /**
         * < Service class name
         */
        public String name;

        /**
         * < Consumer name associated with the service class
         */
        public String consumer;

        /**
         * < Objectives
         */
        public int goaltype;

        /**
         * < The service class state (ontime, delayed)
         */
        public int state;

        /**
         * < Optimum number of job slots (or concurrently running jobs) needed for the  service class to meet its service-level goals
         */
        public int optimum;

        /**
         * < Job counters for the service class
         */
        public int[] counters = new int[NUM_JGRP_COUNTERS];
    }



    /**
     * \brief  a wrap of structure perfmonLog for performance metrics project
     */
    public static class perfmonLogInfo extends Structure {
        public static class ByReference extends perfmonLogInfo implements Structure.ByReference {}
        public static class ByValue extends perfmonLogInfo implements Structure.ByValue {}
        public perfmonLogInfo() {}
        public perfmonLogInfo(Pointer p) { super(p); read(); }


        /**
         * <  Sample period
         */
        public int samplePeriod;

        /**
         * <  Metrics
         */
        public IntByReference metrics;

        /**
         * <  Start time
         */
        public NativeLong startTime;

        /**
         * <  Log time
         */
        public NativeLong logTime;
    }



    /**
     * \brief performance metrics log in lsb.stream
     */
    public static class perfmonLog extends Structure {
        public static class ByReference extends perfmonLog implements Structure.ByReference {}
        public static class ByValue extends perfmonLog implements Structure.ByValue {}
        public perfmonLog() {}
        public perfmonLog(Pointer p) { super(p); read(); }


        /**
         * < Sample rate
         */
        public int samplePeriod;

        /**
         * < Number of Queries
         */
        public int totalQueries;

        /**
         * < Number of Job Query
         */
        public int jobQuries;

        /**
         * < Number of Queue Query
         */
        public int queueQuries;

        /**
         * < Number of Host Query
         */
        public int hostQuries;

        /**
         * < Number of Submission Requests
         */
        public int submissionRequest;

        /**
         * < Number of Jobs Submitted
         */
        public int jobSubmitted;

        /**
         * < Number of Dispatched Jobs
         */
        public int dispatchedjobs;

        /**
         * < Number of Job Completed
         */
        public int jobcompleted;

        /**
         * < Number of MultiCluster Jobs Sent
         */
        public int jobMCSend;

        /**
         * < Number of MultiCluster Jobs Received
         */
        public int jobMCReceive;

        /**
         * < Start Time
         */
        public NativeLong startTime;
    }



    /**
     * \brief task finish log.Task accounting record in ssched.acct
     */
    public static class taskFinishLog extends Structure {
        public static class ByReference extends taskFinishLog implements Structure.ByReference {}
        public static class ByValue extends taskFinishLog implements Structure.ByValue {}
        public taskFinishLog() {}
        public taskFinishLog(Pointer p) { super(p); read(); }


        /**
         * <  Job finish event
         */
        public jobFinishLog jobFinishLog;

        /**
         * < Task ID
         */
        public int taskId;

        /**
         * < Task index
         */
        public int taskIdx;

        /**
         * < Name of task
         */
        public String taskName;

        /**
         * < Bit mask of task options: \n TASK_IN_FILE (0x01)-specify input file \n TASK_OUT_FILE (0x02)-specify output file \n TASK_ERR_FILE (0x04)-specify error file \n TASK_PRE_EXEC (0x08)-specify pre-exec command \n TASK_POST_EXEC (0x10)-specify post-exec command \n TASK_NAME (0x20)-specify task name
         */
        public int taskOptions;

        /**
         * < Task Exit Reason \n TASK_EXIT_NORMAL = 0- normal exit \n TASK_EXIT_INIT = 1-generic task initialization failure \n TASK_EXIT_PATH = 2-failed to initialize path \n TASK_EXIT_NO_FILE = 3-failed to create task file \n TASK_EXIT_PRE_EXEC = 4- task pre-exec failed \n TASK_EXIT_NO_PROCESS = 5-fork failed \n TASK_EXIT_XDR = 6-xdr communication error \n TASK_EXIT_NOMEM = 7- no memory \n TASK_EXIT_SYS = 8-system call failed \n TASK_EXIT_TSCHILD_EXEC = 9-failed to run sschild \n TASK_EXIT_RUNLIMIT = 10-task reaches run limit \n TASK_EXIT_IO = 11-I/O failure \n TASK_EXIT_RSRC_LIMIT = 12-set task resource limit failed
         */
        public int taskExitReason;
    }



    /**
     * \brief End of stream event. The stream is moved to lsb.stream.0 and
     * a new lsb.stream is opened. Readers of lsb.stream when encounter
     * the event EVENT_END_OF_STREAM should close and reopen the
     * lsb.stream file.
     */
    public static class eventEOSLog extends Structure {
        public static class ByReference extends eventEOSLog implements Structure.ByReference {}
        public static class ByValue extends eventEOSLog implements Structure.ByValue {}
        public eventEOSLog() {}
        public eventEOSLog(Pointer p) { super(p); read(); }


        /**
         * < Event end of stream
         */
        public int eos;
    }



    /**
     * \brief job resize event: indicating a realized job allocation change
     */
    public static class jobResizeLog extends Structure {
        public static class ByReference extends jobResizeLog implements Structure.ByReference {}
        public static class ByValue extends jobResizeLog implements Structure.ByValue {}
        public jobResizeLog() {}
        public jobResizeLog(Pointer p) { super(p); read(); }


        /**
         * <  JobId
         */
        public int jobId;

        /**
         * <  Index
         */
        public int idx;

        /**
         * <  Start time
         */
        public NativeLong startTime;

        /**
         * <  User Id
         */
        public int userId;

        /**
         * <  User name
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < 0 grow, 1 shrink
         */
        public int resizeType;

        /**
         * < The start time of last allocation
         */
        public NativeLong lastResizeStartTime;

        /**
         * < The finish time of last allocation
         */
        public NativeLong lastResizeFinishTime;

        /**
         * < Allocation before the resize
         */
        public int numExecHosts;

        /**
         * <  Execute hosts
         */
        public Pointer execHosts;

        /**
         * < The delta of the allocation change
         */
        public int numResizeHosts;

        /**
         * <  Resize hosts
         */
        public Pointer resizeHosts;
    }



    /**
     * \brief  Log event types.
     */
    public static class eventLog extends Union {
        /**
         * <  Job new event
         */
        public jobNewLog jobNewLog;

        /**
         * <  Job start event
         */
        public jobStartLog jobStartLog;

        /**
         * <  Job status event
         */
        public jobStatusLog jobStatusLog;

        /**
         * <  sbatchd job status event
         */
        public sbdJobStatusLog sbdJobStatusLog;

        /**
         * <  Job switch event
         */
        public jobSwitchLog jobSwitchLog;

        /**
         * <  Job move event
         */
        public jobMoveLog jobMoveLog;

        /**
         * <  Queue control event
         */
        public queueCtrlLog queueCtrlLog;

/* New debug event*/
        public newDebugLog newDebugLog;

        /**
         * <  Host control event
         */
        public hostCtrlLog hostCtrlLog;

        /**
         * <  mbatchd start event
         */
        public mbdStartLog mbdStartLog;

        /**
         * <  mbatchd die event
         */
        public mbdDieLog mbdDieLog;

        /**
         * <  Unfulfill event
         */
        public unfulfillLog unfulfillLog;

        /**
         * <  Job finish event
         */
        public jobFinishLog jobFinishLog;

        /**
         * <  Load index event
         */
        public loadIndexLog loadIndexLog;

        /**
         * <  Migration initiated event
         */
        public migLog migLog;

        /**
         * <  Calendar event
         */
        public calendarLog calendarLog;

        /**
         * <  Job forward event
         */
        public jobForwardLog jobForwardLog;

        /**
         * <  Job accept event
         */
        public jobAcceptLog jobAcceptLog;

        /**
         * <  Job accepted from another  cluster event
         */
        public statusAckLog statusAckLog;

        /**
         * <  Job signal event
         */
        public signalLog signalLog;

        /**
         * <  Job execution event
         */
        public jobExecuteLog jobExecuteLog;

        /**
         * <  Job message event
         */
        public jobMsgLog jobMsgLog;

        /**
         * <  Job message ackknowledge event
         */
        public jobMsgAckLog jobMsgAckLog;

        /**
         * <  Job requeue event
         */
        public jobRequeueLog jobRequeueLog;

        /**
         * <  Checkpoint event
         */
        public chkpntLog chkpntLog;

        /**
         * <  Signal with action event
         */
        public sigactLog sigactLog;

        /**
         * <  Job occupy request event
         */
        public jobOccupyReqLog jobOccupyReqLog;

        /**
         * <  Job vacate event
         */
        public jobVacatedLog jobVacatedLog;

        /**
         * <  Job start accept event
         */
        public jobStartAcceptLog jobStartAcceptLog;

        /**
         * <  Job clean event
         */
        public jobCleanLog jobCleanLog;

        /**
         * <  Job exception event
         */
        public jobExceptionLog jobExceptionLog;

        /**
         * <  Job group new event
         */
        public jgrpNewLog jgrpNewLog;

        /**
         * <  Job group Ctrl event
         */
        public jgrpCtrlLog jgrpCtrlLog;

        /**
         * <  Job Force Request event
         */
        public jobForceRequestLog jobForceRequestLog;

        /**
         * <  Event switch event
         */
        public logSwitchLog logSwitchLog;

        /**
         * <  Job modify event
         */
        public jobModLog jobModLog;

        /**
         * <  Job group stratus event
         */
        public jgrpStatusLog jgrpStatusLog;

        /**
         * <  Job attribute setting event
         */
        public jobAttrSetLog jobAttrSetLog;

        /**
         * <  Job external message event
         */
        public jobExternalMsgLog jobExternalMsgLog;

        /**
         * <  Job chunk event
         */
        public jobChunkLog jobChunkLog;

        /**
         * < sbatchd  unreported status event
         */
        public sbdUnreportedStatusLog sbdUnreportedStatusLog;

        /**
         * <  Reservation finish event
         */
        public rsvFinishLog rsvFinishLog;

        /**
         * <  Host group control Log
         */
        public hgCtrlLog hgCtrlLog;

        /**
         * <  cpu profile event
         */
        public cpuProfileLog cpuProfileLog;

        /**
         * <  Data logging event
         */
        public dataLoggingLog dataLoggingLog;

        /**
         * <  Job run rusage event
         */
        public jobRunRusageLog jobRunRusageLog;

        /**
         * <  Event EOS event
         */
        public eventEOSLog eventEOSLog;

        /**
         * <  SLA event
         */
        public slaLog slaLog;

        /**
         * <  Performance event
         */
        public perfmonLog perfmonLog;

        /**
         * <  Task finish event
         */
        public taskFinishLog taskFinishLog;

        /**
         * <  Job resize notify start event
         */
        public jobResizeNotifyStartLog jobResizeNotifyStartLog;

        /**
         * <  Job resize notify accept event
         */
        public jobResizeNotifyAcceptLog jobResizeNotifyAcceptLog;

        /**
         * <  Job resize notify done event
         */
        public jobResizeNotifyDoneLog jobResizeNotifyDoneLog;

        /**
         * <  Job resize release event
         */
        public jobResizeReleaseLog jobResizeReleaseLog;

        /**
         * <  Job resize cancel event
         */
        public jobResizeCancelLog jobResizeCancelLog;

        /**
         * <  Job resize event
         */
        public jobResizeLog jobResizeLog;

/*#if defined(LSF_SIMULATOR)*/
/**< Job array element event */
        /*public jobArrayElementLog jobArrayElementLog;*/

        /**< LSF simulator status event */
        /*public mbdSimStatusLog   mbdSimStatusLog;*/
        /*#endif*/
    }




    /**
     * \brief  event records.
     */
    public static class eventRec extends Structure {
        public static class ByReference extends eventRec implements Structure.ByReference {}
        public static class ByValue extends eventRec implements Structure.ByValue {}
        public eventRec() {}
        public eventRec(Pointer p) { super(p); read(); }


        /**
         * < The mbatchd version number
         */
        public byte[] version = new byte[MAX_VERSION_LEN];

        /**
         * < Event type in \ref event_types
         */
        public int type;

        /**
         * < The time the event occurred
         */
        public NativeLong eventTime;

        /**
         * < The information for this type of event, contained in a structure  corresponding to type
         */
        public eventLog eventLog;
    }



    public static class eventLogFile extends Structure {
        public static class ByReference extends eventLogFile implements Structure.ByReference {}
        public static class ByValue extends eventLogFile implements Structure.ByValue {}
        public eventLogFile() {}
        public eventLogFile(Pointer p) { super(p); read(); }


/* event file directory */
        public byte[] eventDir = new byte[LibLsf.MAXFILENAMELEN];

/* start and end event time */
        public NativeLong beginTime, endTime;
    }



    public static class eventLogHandle extends Structure {
        public static class ByReference extends eventLogHandle implements Structure.ByReference {}
        public static class ByValue extends eventLogHandle implements Structure.ByValue {}
        public eventLogHandle() {}
        public eventLogHandle(Pointer p) { super(p); read(); }


/* open event file pointer */
        public Pointer fp;

/* current open events file name */
        public byte[] openEventFile = new byte[LibLsf.MAXFILENAMELEN];

/* current open event file number */
        public int curOpenFile;
        public int lastOpenFile;                   /* last open event file number, 0
                  means lsb.events */
    }




    public static final String LSF_JOBIDINDEX_FILENAME = "lsb.events.index";
    public static final String LSF_JOBIDINDEX_FILETAG = "#LSF_JOBID_INDEX_FILE";

/* structures used to handle jobId index file */

    public static class jobIdIndexS extends Structure {
        public static class ByReference extends jobIdIndexS implements Structure.ByReference {}
        public static class ByValue extends jobIdIndexS implements Structure.ByValue {}
        public jobIdIndexS() {}
        public jobIdIndexS(Pointer p) { super(p); read(); }


/* the index file name */
        public byte[] fileName = new byte[LibLsf.MAXFILENAMELEN];

/* open index file pointer */
        public Pointer fp;

/* version number for future use */
        public float version;

/* total number of rows(files) indices */
        public int totalRows;

/* last update time */
        public NativeLong lastUpdate;

/* current rows */
        public int curRow;
        /* the event file currently handled is */
        /* (totalRows - curRow + 1) */

/* time stamp of current row */
        public NativeLong timeStamp;

/* min jobId in that row */
        public long minJobId;

/* max jobId in that row */
        public long maxJobId;

/* total number of jobIds */
        public int totalJobIds;

/* jobId array of current row */
        public IntByReference jobIds;
    }



/* structures used to hold one element of sorted int list */

    public static class sortIntList extends Structure {
        public static class ByReference extends sortIntList implements Structure.ByReference {}
        public static class ByValue extends sortIntList implements Structure.ByValue {}
        public sortIntList() {}
        public sortIntList(Pointer p) { super(p); read(); }

        public int value;

/* points to next element */
        public sortIntList.ByReference forw;

/* points to prior element */
        public sortIntList.ByReference back;
    }



    public static class nqsStatusReq extends Structure {
        public static class ByReference extends nqsStatusReq implements Structure.ByReference {}
        public static class ByValue extends nqsStatusReq implements Structure.ByValue {}
        public nqsStatusReq() {}
        public nqsStatusReq(Pointer p) { super(p); read(); }

        public long jobId;
        public int opCode;
        public int reportCode;
        public String nqsQueue;
        public int fromUid;
        public String fromUserName;
        public String fromHostName;
        public int idx;
    }



    public static class nqsStatusReply extends Structure {
        public static class ByReference extends nqsStatusReply implements Structure.ByReference {}
        public static class ByValue extends nqsStatusReply implements Structure.ByValue {}
        public nqsStatusReply() {}
        public nqsStatusReply(Pointer p) { super(p); read(); }

        public String orgHost;
        public String orgUser;
        public NativeLong startTime;
        public String jobName;
        public String nqsQueue;
        public String lsbManager;
        public int options;
        public String outFile;
        public String errFile;
    }



/*
*  SBD uses the following data structure to communicate with
*  the resource manager.
*
 */
    public static final int LSB_MAX_SD_LENGTH = 128;

    public static class lsbMsgHdr extends Structure {
        public static class ByReference extends lsbMsgHdr implements Structure.ByReference {}
        public static class ByValue extends lsbMsgHdr implements Structure.ByValue {}
        public lsbMsgHdr() {}
        public lsbMsgHdr(Pointer p) { super(p); read(); }

        public int usrId;
        public long jobId;
        public int msgId;
        public int type;
        public String src;
        public String dest;
    }



    public static class lsbMsg extends Structure {
        public static class ByReference extends lsbMsg implements Structure.ByReference {}
        public static class ByValue extends lsbMsg implements Structure.ByValue {}
        public lsbMsg() {}
        public lsbMsg(Pointer p) { super(p); read(); }

        public lsbMsgHdr.ByReference header;
        public String msg;
    }



/* data structures related to API_CONF */

    public static final int CONF_NO_CHECK = 0x00;
    public static final int CONF_CHECK = 0x01;
    public static final int CONF_EXPAND = 0X02;
    public static final int CONF_RETURN_HOSTSPEC = 0X04;
    public static final int CONF_NO_EXPAND = 0X08;
    public static final int CONF_HAS_CU = 0X10;

    public static class paramConf extends Structure {
        public static class ByReference extends paramConf implements Structure.ByReference {}
        public static class ByValue extends paramConf implements Structure.ByValue {}
        public paramConf() {}
        public paramConf(Pointer p) { super(p); read(); }

        public parameterInfo.ByReference param;
    }



    public static class userConf extends Structure {
        public static class ByReference extends userConf implements Structure.ByReference {}
        public static class ByValue extends userConf implements Structure.ByValue {}
        public userConf() {}
        public userConf(Pointer p) { super(p); read(); }

        public int numUgroups;
        public Pointer /* groupInfoEnt.ByReference */ ugroups;
        public int numUsers;
        public Pointer /* userInfoEnt.ByReference */ users;
        public int numUserEquivalent;
        public Pointer /* userEquivalentInfoEnt.ByReference */ userEquivalent;
        public int numUserMapping;
        public Pointer /* userMappingInfoEnt.ByReference */ userMapping;
    }



    public static class hostConf extends Structure {
        public static class ByReference extends hostConf implements Structure.ByReference {}
        public static class ByValue extends hostConf implements Structure.ByValue {}
        public hostConf() {}
        public hostConf(Pointer p) { super(p); read(); }

        public int numHosts;
        public Pointer /* hostInfoEnt.ByReference */ hosts;
        public int numHparts;
        public Pointer /* hostPartInfoEnt.ByReference */ hparts;
        public int numHgroups;
        public Pointer /* groupInfoEnt.ByReference */ hgroups;
    }



    /**
     * \brief  lsb shared resource Instance.
     */
    public static class lsbSharedResourceInstance extends Structure {
        public static class ByReference extends lsbSharedResourceInstance implements Structure.ByReference {}
        public static class ByValue extends lsbSharedResourceInstance implements Structure.ByValue {}
        public lsbSharedResourceInstance() {}
        public lsbSharedResourceInstance(Pointer p) { super(p); read(); }


        /**
         * < Value used by mbatchd
         */
        public String totalValue;

        /**
         * < Reserved value
         */
        public String rsvValue;

        /**
         * < Number of Hosts associated with the resource.
         */
        public int nHosts;

        /**
         * < Hosts list
         */
        public Pointer hostList;
    }



    /**
     * \brief lsb shared resource information.
     */
    public static class lsbSharedResourceInfo extends Structure {
        public static class ByReference extends lsbSharedResourceInfo implements Structure.ByReference {}
        public static class ByValue extends lsbSharedResourceInfo implements Structure.ByValue {}
        public lsbSharedResourceInfo() {}
        public lsbSharedResourceInfo(Pointer p) { super(p); read(); }


        /**
         * < Resource name
         */
        public String resourceName;

        /**
         * < Number of instances
         */
        public int nInstances;

        /**
         * < List of instances
         */
        public Pointer /* lsbSharedResourceInstance.ByReference */ instances;
    }



    public static class queueConf extends Structure {
        public static class ByReference extends queueConf implements Structure.ByReference {}
        public static class ByValue extends queueConf implements Structure.ByValue {}
        public queueConf() {}
        public queueConf(Pointer p) { super(p); read(); }

        public int numQueues;
        public Pointer /* queueInfoEnt.ByReference */ queues;
    }



    /**
     * \brief  frame element information.
     */
    public static class frameElementInfo extends Structure {
        public static class ByReference extends frameElementInfo implements Structure.ByReference {}
        public static class ByValue extends frameElementInfo implements Structure.ByValue {}
        public frameElementInfo() {}
        public frameElementInfo(Pointer p) { super(p); read(); }


        /**
         * <  The job index in the frame job array.
         */
        public int jobindex;

        /**
         * <  The job status.
         */
        public int jobState;

        /**
         * <  The start frame of this frame job.
         */
        public int start;

        /**
         * <  The end frame of this frame job.
         */
        public int end;

        /**
         * <  The step of this frame job.
         */
        public int step;

        /**
         * <  The chunk size of this frame job.
         */
        public int chunk;
    }



    /**
     * \brief  frame job Infomation.
     */
    public static class frameJobInfo extends Structure {
        public static class ByReference extends frameJobInfo implements Structure.ByReference {}
        public static class ByValue extends frameJobInfo implements Structure.ByValue {}
        public frameJobInfo() {}
        public frameJobInfo(Pointer p) { super(p); read(); }


        /**
         * < The job ID that the LSF system assigned to the frame job array.
         */
        public long jobGid;

        /**
         * < The max job number in one frame job array.
         */
        public int maxJob;

        /**
         * < The user submitted the frame job array.
         */
        public byte[] userName = new byte[MAX_LSB_NAME_LEN];

        /**
         * < Full job name
         */
        public byte[] jobName = new byte[LibLsf.MAXLINELEN];

        /**
         * < The full job name of the frame job array.  frameElementPtr The pointer to frame ob array table.
         */
        public frameElementInfo.ByReference frameElementPtr;
    }



    public static class nqsRusageReq extends Structure {
        public static class ByReference extends nqsRusageReq implements Structure.ByReference {}
        public static class ByValue extends nqsRusageReq implements Structure.ByValue {}
        public nqsRusageReq() {}
        public nqsRusageReq(Pointer p) { super(p); read(); }

        public long jobId;
        public int mem;
        public float cpuTime;
    }



    public static class nqsRusageReply extends Structure {
        public static class ByReference extends nqsRusageReply implements Structure.ByReference {}
        public static class ByValue extends nqsRusageReply implements Structure.ByValue {}
        public nqsRusageReply() {}
        public nqsRusageReply(Pointer p) { super(p); read(); }

        public int status;
    }



/* end of data structures related to API_CONF */

/*
*  Structure used for the Advance Reservation API
*
*  MBD allows the LSF administration to make advance reservation on
*  behalf of a user, group or or for system maintenance purposes.
*  Clients can add a reservation, remove a reservation and show
*  reservation statuses.  The following data structures are used to
*  encapsulate these requests
*
*     addRsvRequest: to add a reservation
*     rmRsvRequest:  to remove a reservation
*     rsvInfoEnt:    to display reservation information
*
 */

    public static class _rsvEventInfo_prePost_t extends Structure {
        public static class ByReference extends _rsvEventInfo_prePost_t implements Structure.ByReference {}
        public static class ByValue extends _rsvEventInfo_prePost_t implements Structure.ByValue {}
        public _rsvEventInfo_prePost_t() {}
        public _rsvEventInfo_prePost_t(Pointer p) { super(p); read(); }

        public int shift;
    }



    public static final int RSV_EXECEVENTTYPE_PRE = 1;
    public static final int RSV_EXECEVENTTYPE_POST = 2;

    public static final String RSV_EXECEVENTNAME_PRE = "pre";
    public static final String RSV_EXECEVENTNAME_POST = "post";

    /**
     * \brief  reservation excution event
     */
    public static class _rsvExecEvent_t extends Structure {
        public static class ByReference extends _rsvExecEvent_t implements Structure.ByReference {}
        public static class ByValue extends _rsvExecEvent_t implements Structure.ByValue {}
        public _rsvExecEvent_t() {}
        public _rsvExecEvent_t(Pointer p) { super(p); read(); }


        /**
         * < Event type
         */
        public int type;

        /**
         * < Boolean: is there additional info?
         */
        public int infoAttached;

        /**
         * < Info pertaining to event, such as offset
         */
        public Pointer info;
    }



    /**
     * \brief  reservation excution command
     */
    public static class _rsvExecCmd_t extends Structure {
        public static class ByReference extends _rsvExecCmd_t implements Structure.ByReference {}
        public static class ByValue extends _rsvExecCmd_t implements Structure.ByValue {}
        public _rsvExecCmd_t() {}
        public _rsvExecCmd_t(Pointer p) { super(p); read(); }


        /**
         * < Full path to the command name
         */
        public String path;

        /**
         * < Size of events array
         */
        public int numEvents;

        /**
         * < Array of events that trigger -exec command
         */
        public Pointer /* _rsvExecEvent_t.ByReference */ events;
    }



    /**
     *  \addtogroup reservation_option reservation_option
     *    definitions of reservation options.
     */

    /**
     * <  User
     */
    public static final int RSV_OPTION_USER = 0x0001;

    /**
     * <  Group
     */
    public static final int RSV_OPTION_GROUP = 0x0002;

    /**
     * <  System
     */
    public static final int RSV_OPTION_SYSTEM = 0x0004;

    /**
     * <  Recur
     */
    public static final int RSV_OPTION_RECUR = 0x0008;

    /**
     * <  Resource requirement
     */
    public static final int RSV_OPTION_RESREQ = 0x0010;

    /**
     * <  Host
     */
    public static final int RSV_OPTION_HOST = 0x0020;

    /**
     * <  Open
     */
    public static final int RSV_OPTION_OPEN = 0x0040;

    /**
     * <  Delete
     */
    public static final int RSV_OPTION_DELETE = 0x0080;

    /**
     * <  Close
     */
    public static final int RSV_OPTION_CLOSED = 0x0100;

    /**
     * <  Execute
     */
    public static final int RSV_OPTION_EXEC = 0x0200;

    /**
     * <  Remote execute
     */
    public static final int RSV_OPTION_RMEXEC = 0x0400;

    /**
     * <  Next instance
     */
    public static final int RSV_OPTION_NEXTINSTANCE = 0x0800;

    /**
     * <  Disable
     */
    public static final int RSV_OPTION_DISABLE = 0x1000;

    /**
     * <  Add host
     */
    public static final int RSV_OPTION_ADDHOST = 0x2000;

    /**
     * <  Remote host
     */
    public static final int RSV_OPTION_RMHOST = 0x4000;

    /**
     * <  Description
     */
    public static final int RSV_OPTION_DESCRIPTION = 0x8000;

    /**
     * <  Timewindow mode
     */
    public static final int RSV_OPTION_TWMOD = 0x10000;

    /**
     * <  Switch open/close
     */
    public static final int RSV_OPTION_SWITCHOPENCLOSE = 0x20000;

    /**
     * <  User mode
     */
    public static final int RSV_OPTION_USERMOD = 0x40000;

    /**
     * <  Reservation name
     */
    public static final int RSV_OPTION_RSVNAME = 0x80000;

    /**
     * <  Expired
     */
    public static final int RSV_OPTION_EXPIRED = 0x100000;

    /**
     * \brief add reservation request.
     */
    public static class addRsvRequest extends Structure {
        public static class ByReference extends addRsvRequest implements Structure.ByReference {}
        public static class ByValue extends addRsvRequest implements Structure.ByValue {}
        public addRsvRequest() {}
        public addRsvRequest(Pointer p) { super(p); read(); }


        /**
         * <Reservation options \ref reservation_option
         */
        public int options;

        /**
         * < User or group for which the reservation is made
         */
        public String name;

        /**
         * < Minimum number of processors the required to run the job. See the -g option of brsvadd.
         */
        public int minNumProcs;

        /**
         * < Maximum number of processors the required to run the job.
         */
        public int maxNumProcs;

        /**< Range of number of processors */
        //struct procRange;

        /**
         * < The number of invoker specified hosts for the reservation. If numAskedHosts is 0, all qualified hosts will be considered.
         */
        public int numAskedHosts;

        /**
         * < The array of names of invoker specified hosts hosts for the reservation. The number of hosts is given by numAskedHosts. See the -m option of brsvadd.
         */
        public Pointer askedHosts;

        /**
         * < The resource requirements of the reservation. See the -R option of brsvadd.
         */
        public String resReq;

        /**
         * < Active time window for a recurring reservation. See the -t option of brsvadd.
         */
        public String timeWindow;

        /**
         * < Info for the -exec option.
         */
        public _rsvExecCmd_t.ByReference execCmd;

        /**
         * < Description for the reservation to be created. The description must be provided as a double quoted text string. The maximum length  is 512 chars.  Equivalent to the value of brsvadd -d.
         */
        public String desc;

        /**
         * < User-defined advance reservation name unique in an LSF cluster. The name is a string of letters, numeric chars, underscores, and dashes beginning with a letter. The maximum length of the name is 39 chars. Equivalent to the value of brsvadd -N.
         */
        public String rsvName;
    }



    /**
     * \brief  remove reservation request.
     */
    public static class rmRsvRequest extends Structure {
        public static class ByReference extends rmRsvRequest implements Structure.ByReference {}
        public static class ByValue extends rmRsvRequest implements Structure.ByValue {}
        public rmRsvRequest() {}
        public rmRsvRequest(Pointer p) { super(p); read(); }


        /**
         * < Reservation ID of the reservation that you wish to remove.
         */
        public String rsvId;
    }



    /**
     * \brief  modifiy reservation request
     */
    public static class modRsvRequest extends Structure {
        public static class ByReference extends modRsvRequest implements Structure.ByReference {}
        public static class ByValue extends modRsvRequest implements Structure.ByValue {}
        public modRsvRequest() {}
        public modRsvRequest(Pointer p) { super(p); read(); }


        /**
         * < Reservation ID of the reservation that you  wish to modify.
         */
        public String rsvId;

        /**
         * < LSF user name for the reservation.  See the -g option of brsvadd. .
         */
        public addRsvRequest fieldsFromAddReq;

        /**
         * < Disabled time duration
         */
        public String disabledDuration;
    }



    /**
     * \brief  host reservation infromation entry.
     */
    public static class hostRsvInfoEnt extends Structure {
        public static class ByReference extends hostRsvInfoEnt implements Structure.ByReference {}
        public static class ByValue extends hostRsvInfoEnt implements Structure.ByValue {}
        public hostRsvInfoEnt() {}
        public hostRsvInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < Host name.
         */
        public String host;

        /**
         * < Number of CPUs reserved on the host.
         */
        public int numCPUs;

        /**
         * < Number of job slots reserved on the host.
         */
        public int numSlots;

        /**
         * < Number of processors reserved on the host.
         */
        public int numRsvProcs;

        /**
         * < Count for used + suspended from reserved slots
         */
        public int numusedRsvProcs;

        /**
         * < Number of processors in use on the host.
         */
        public int numUsedProcs;
    }



    /**
     * \brief  reservation information entry.
     */
    public static class rsvInfoEnt extends Structure {
        public static class ByReference extends rsvInfoEnt implements Structure.ByReference {}
        public static class ByValue extends rsvInfoEnt implements Structure.ByValue {}
        public rsvInfoEnt() {}
        public rsvInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < Reservation options, see \ref reservation_option
         */
        public int options;

        /**
         * < Reservation ID returned from mbatchd. If the reservation fails, this is null. The memory for rsvid is allocated by the caller.
         */
        public String rsvId;

        /**
         * <  LSF user group name for the reservation. See the -g option of brsvadd.
         */
        public String name;

        /**
         * <  Number of hosts reserved
         */
        public int numRsvHosts;

        /**
         * <  Info about the reserved hosts
         */
        public Pointer /* hostRsvInfoEnt.ByReference */ rsvHosts;

        /**
         * < Active time window for a recurring reservation. See the -t option of  brsvadd.
         */
        public String timeWindow;

        /**
         * < Number of jobs running in the reservation.
         */
        public int numRsvJobs;

        /**
         * < Job IDs of jobs running in the reservation.
         */
        public LongByReference jobIds;

        /**
         * < Status of jobs running in the reservation.
         */
        public IntByReference jobStatus;

        /**
         * <  Description for the reservation to be created. The description must be provided as a double quoted text string. The maximum length is 512 chars. Equivalent to thevalue of brsvadd -d.
         */
        public String desc;

        /**
         * <  Null-terminated list of disabled durations
         */
        public Pointer disabledDurations;

        /**
         * <  The current state of the reservation - active or inactive.
         */
        public int state;

        /**
         * <  The time of the next instance of a recurring reservation.
         */
        public String nextInstance;

        /**
         * <  Creator of the reservation.
         */
        public String creator;
    }



/* backfill window related data structures and functions */

    public static class slotInfoRequest extends Structure {
        public static class ByReference extends slotInfoRequest implements Structure.ByReference {}
        public static class ByValue extends slotInfoRequest implements Structure.ByValue {}
        public slotInfoRequest() {}
        public slotInfoRequest(Pointer p) { super(p); read(); }

        /* options mask */

/* Option -R */
        public static int SLOT_OPTION_RESREQ = 0X001;

        public int options;

/* Resource requirement string */
        public String resReq;
    }



/*copy from SRInfo*/

    public static class SRInfoEnt extends Structure {
        public static class ByReference extends SRInfoEnt implements Structure.ByReference {}
        public static class ByValue extends SRInfoEnt implements Structure.ByValue {}
        public SRInfoEnt() {}
        public SRInfoEnt(Pointer p) { super(p); read(); }


/*number of reserved slots*/
        public int numReserved;

/* job's predicted start time */
        public NativeLong predictedStartTime;
    }



    public static class hostSRInfoEnt extends Structure {
        public static class ByReference extends hostSRInfoEnt implements Structure.ByReference {}
        public static class ByValue extends hostSRInfoEnt implements Structure.ByValue {}
        public hostSRInfoEnt() {}
        public hostSRInfoEnt(Pointer p) { super(p); read(); }

        public String host;
        public int hStatus;
        public int userJobLimit;
        public int maxJobs;
        public int numJobs;
        public int numRUN;
        public int numSSUSP;
        public int numUSUSP;
        public int numRESERVE;
        public int numSR;
        public Pointer /* SRInfoEnt.ByReference */ SRInfo;
    }



    public static class slotInfoReply extends Structure {
        public static class ByReference extends slotInfoReply implements Structure.ByReference {}
        public static class ByValue extends slotInfoReply implements Structure.ByValue {}
        public slotInfoReply() {}
        public slotInfoReply(Pointer p) { super(p); read(); }


/* to store the time of Master host */
        public NativeLong masterTime;
        public int numHosts;
        public Pointer /* hostSRInfoEnt.ByReference */ hostInfo;
        public int numAR;
        public Pointer /* rsvInfoEnt.ByReference */ ARInfo;
    }




/* the general limit related data structures and functions */


    public static final int LSB_RSRC_LIMIT_TYPE_SLOTS = 0;
    public static final int LSB_RSRC_LIMIT_TYPE_SLOT_PERPSR = 1;
    public static final int LSB_RSRC_LIMIT_TYPE_MEM = 2;
    public static final int LSB_RSRC_LIMIT_TYPE_MEM_PERCENT = 3;
    public static final int LSB_RSRC_LIMIT_TYPE_SWP = 4;
    public static final int LSB_RSRC_LIMIT_TYPE_SWP_PERCENT = 5;
    public static final int LSB_RSRC_LIMIT_TYPE_TMP = 6;
    public static final int LSB_RSRC_LIMIT_TYPE_TMP_PERCENT = 7;
    public static final int LSB_RSRC_LIMIT_TYPE_JOBS = 8;

/* all external resources */
    public static final int LSB_RSRC_LIMIT_TYPE_EXT_RSRC = 9;

    /**
     * \addtogroup _consumertype _consumertype
     * consumer types
     */
    public static interface consumerType {
        /**
         * < Queues
         */
        public static final int LIMIT_QUEUES = 1;

        /**
         * < Per-queue
         */
        public static final int LIMIT_PER_QUEUE = 2;

        /**
         * < Users
         */
        public static final int LIMIT_USERS = 3;

        /**
         * < Per-users
         */
        public static final int LIMIT_PER_USER = 4;

        /**
         * < Hosts
         */
        public static final int LIMIT_HOSTS = 5;

        /**
         * < Per-host
         */
        public static final int LIMIT_PER_HOST = 6;

        /**
         * < Projects
         */
        public static final int LIMIT_PROJECTS = 7;

        /**
         * < Per-project
         */
        public static final int LIMIT_PER_PROJECT = 8;
    }


    /**< Type definitions */

    /**
     * \brief  limit consumer
     */
    public static class _limitConsumer extends Structure {
        public static class ByReference extends _limitConsumer implements Structure.ByReference {}
        public static class ByValue extends _limitConsumer implements Structure.ByValue {}
        public _limitConsumer() {}
        public _limitConsumer(Pointer p) { super(p); read(); }


        /**
         * < Consumer type ( _consumertype ):  -  Queues per-queue -  Users and per-user -  Hosts and per-host -  Projects and per-project
         */
        public int type;

        /**
         * < Consumer name
         */
        public String name;
    }



    /**
     * \brief  limit resource.
     */
    public static class _limitResource extends Structure {
        public static class ByReference extends _limitResource implements Structure.ByReference {}
        public static class ByValue extends _limitResource implements Structure.ByValue {}
        public _limitResource() {}
        public _limitResource(Pointer p) { super(p); read(); }


        /**
         * < Resource name
         */
        public String name;

        /**
         * < Resource type
         */
        public int type;

        /**
         * < Resource val
         */
        public float val;
    }



    /**
     * \brief   limit information request
     */
    public static class _limitInfoReq extends Structure {
        public static class ByReference extends _limitInfoReq implements Structure.ByReference {}
        public static class ByValue extends _limitInfoReq implements Structure.ByValue {}
        public _limitInfoReq() {}
        public _limitInfoReq(Pointer p) { super(p); read(); }


        /**
         * < Limit policy name given by the user.
         */
        public String name;

        /**
         * < Number of consumers
         */
        public int consumerC;

        /**
         * < Consumer name, queue/host/user/project
         */
        public Pointer /* _limitConsumer.ByReference */ consumerV;
    }



    /**
     * \brief  limit item.
     */
    public static class _limitItem extends Structure {
        public static class ByReference extends _limitItem implements Structure.ByReference {}
        public static class ByValue extends _limitItem implements Structure.ByValue {}
        public _limitItem() {}
        public _limitItem(Pointer p) { super(p); read(); }


        /**
         * < Number of consumers
         */
        public int consumerC;

        /**
         * < Consumers, such as queue, host, user or project
         */
        public Pointer /* _limitConsumer.ByReference */ consumerV;

        /**
         * < Number of resources
         */
        public int resourceC;

        /**
         * < Resources list
         */
        public Pointer /* _limitResource.ByReference */ resourceV;
    }



    /**
     * \brief  limit information entry .
     */
    public static class _limitInfoEnt extends Structure {
        public static class ByReference extends _limitInfoEnt implements Structure.ByReference {}
        public static class ByValue extends _limitInfoEnt implements Structure.ByValue {}
        public _limitInfoEnt() {}
        public _limitInfoEnt(Pointer p) { super(p); read(); }


        /**
         * < Limit policy name given by the user
         */
        public String name;

        /**
         * < Limit configuration
         */
        public _limitItem confInfo;

        /**
         * < Size of limit dynamic usage info array
         */
        public int usageC;

        /**
         * < Limit dynamic usage info array
         */
        public Pointer /* _limitItem.ByReference */ usageInfo;

    }



/* action code for threshold based on type/model, is used for
*  predefinedThresholdTypeModel().
 */

    public static final int ADD_THRESHOLD = 1;
    public static final int GET_THRESHOLD = 2;
    public static final int DEL_THRESHOLD = 3;

/* Structure to hold thresholds defined based on host's type/model */

    public static class thresholdEntry extends Structure {
        public static class ByReference extends thresholdEntry implements Structure.ByReference {}
        public static class ByValue extends thresholdEntry implements Structure.ByValue {}
        public thresholdEntry() {}
        public thresholdEntry(Pointer p) { super(p); read(); }


/* Name of type or model */
        public String attr;

/* Pointer to hostInfo */
        public hostInfoEnt.ByReference hostEntryPtr;
    }



    /**
     * \page lsb_limitInfo lsb_limitInfo
     * \brief gets resource allocation limit configuration and dynamic usage
     * information.
     * <p/>
     * Displays current usage of resource allocation limits configured in Limit
     * sections in lsb.resources:
     * \li    Configured limit policy name
     * \li    Users
     * \li    Queues
     * \li    Hosts
     * \li    Project names
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_limitInfo( limitInfoReq.ByReference req,  limitInfoEnt.ByReference[] limitItemRef,
     * IntByReference size, lsInfo.ByReference lsInfo)</b>
     *
     * @return int:-1
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         blimits
     *         <p/>
     *         \b Files
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.queues \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.users \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.resources
     * @param req input, the user request for limit information
     * @param limitItemRef output, the limit information array
     * @param size output, the size of the limit information array
     * @param lsInfo Please refer to the \ref lsInfo structure.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * _limitInfoReq
     * \n _limitConsumer
     * \n _limitInfoEnt
     * \n _limitItem
     * \n _limitResource
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref _consumertype
     * #see \ref lsb_freeLimitInfoEnt
     */
    public static native int lsb_limitInfo(_limitInfoReq req, Pointer limitItemRef, IntByReference size, LibLsf.lsInfo lsInfo);

    /**
     * \page lsb_freeLimitInfoEnt lsb_freeLimitInfoEnt
     * \brief Frees the memory allocated by \ref lsb_limitInfo.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * void lsb_freeLimitInfoEnt(limitInfoEnt.ByReference  ent, int size)</b>
     *
     * @param size input, the size of the limit information array
     *             <p/>
     *             <b>Data Structures:</b>
     *             \par
     *             _limitInfoEnt
     *             \n _limitItem
     *             \n _limitConsumer
     *             \n _limitResource
     *             <p/>
     *             <b>Define Statements:</b>
     *             \par
     *             \ref _consumertype
     * return void
     *         \n There's no return value.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         blimits
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.queues \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.users \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.resources
     * @param ent input, the array of limit information
     * #see \ref lsb_limitInfo
     */

    public static native void lsb_freeLimitInfoEnt(_limitInfoEnt ent, int size);

    /**
     *  \addtogroup resizablejob_related resizablejob_related
     *  Resizable job related definitions.
     */

    /**
     * < Means release no slots
     */
    public static final int LSB_RESIZE_REL_NONE = 0x0;

    /**
     * < Means release all slots-In this case, nHosts, hosts and slots  indicate the slots that are not released
     */
    public static final int LSB_RESIZE_REL_ALL = 0x01;

    /**
     * < Means cancel any pending resize request
     */
    public static final int LSB_RESIZE_REL_CANCEL = 0x02;

    /**
     * < Means execute no resize notification command
     */
    public static final int LSB_RESIZE_REL_NO_NOTIFY = 0x04;

    /**
     * \brief  job resize release.
     */
    public static class job_resize_release extends Structure {
        public static class ByReference extends job_resize_release implements Structure.ByReference {}
        public static class ByValue extends job_resize_release implements Structure.ByValue {}
        public job_resize_release() {}
        public job_resize_release(Pointer p) { super(p); read(); }


        /**
         * < LSF job ID
         */
        public long jobId;

        /**
         * < Options is constructed from the bitwise inclusive OR of zero or more of the flags, as defined in \ref resizablejob_related .
         */
        public int options;

        /**
         * < Number of hosts in the hosts list, if no hosts are to be specified this should be zero
         */
        public int nHosts;

        /**
         * < Specified hosts list, nHosts number of elements
         */
        public Pointer hosts;

        /**
         * < Slots list, each element specifies the number of slots per corresponding host (0 implies all), nHosts number of elements
         */
        public IntByReference slots;

        /**
         * < Name and location of notification command
         */
        public String notifyCmd;
    }



    public static class job_resize_request extends Structure {
        public static class ByReference extends job_resize_request implements Structure.ByReference {}
        public static class ByValue extends job_resize_request implements Structure.ByValue {}
        public job_resize_request() {}
        public job_resize_request(Pointer p) { super(p); read(); }

        public long jobId;
        public int options;

/* array size */
        public int nHosts;

/* array of hosts */
        public Pointer hosts;

/* notifocation command */
        public String notifyCmd;
    }



/*
*  End of resizable job related definitions
 */

/* Job Dependency Display */


/* Job Dependency Display */
/* for options */
    /**
     *  \addtogroup query_depend query_depend
     *  Job Dependency Display for options
     */

    /**
     * <  Recursively
     */
    public static final int QUERY_DEPEND_RECURSIVELY = 0x1;

    /**
     * <  Detail
     */
    public static final int QUERY_DEPEND_DETAIL = 0x2;

    /**
     * <  Unsatisfied
     */
    public static final int QUERY_DEPEND_UNSATISFIED = 0x4;

    /**
     * <  Child
     */
    public static final int QUERY_DEPEND_CHILD = 0x8;

    /**
     * \brief  job dependent request.
     */

    public static class jobDepRequest extends Structure {
        public static class ByReference extends jobDepRequest implements Structure.ByReference {}
        public static class ByValue extends jobDepRequest implements Structure.ByValue {}
        public jobDepRequest() {}
        public jobDepRequest(Pointer p) { super(p); read(); }

        /**
         * < Job ID of the queried job or job array.
         */
        public long jobId;

        /**
         * < You can set the following bits into this field:\n QUERY_DEPEND_RECURSIVELY\n Query the dependency information recursively.\n QUERY_DEPEND_DETAIL\n Query the detailed dependency information.\n QUERY_DEPEND_UNSATISFIED\n Query the jobs that cause this job pend.\n QUERY_DEPEND_CHILD\n Query child jobs.
         */
        public int options;

        /**
         * < The level when you set QUERY_DEPEND_RECURSIVELY.
         */
        public int level;
    }




    /**
     * \brief  queried jobs.
     */
    public static class queriedJobs extends Structure {
        public static class ByReference extends queriedJobs implements Structure.ByReference {}
        public static class ByValue extends queriedJobs implements Structure.ByValue {}
        public queriedJobs() {}
        public queriedJobs(Pointer p) { super(p); read(); }


        /**
         * < Job ID of the queried job or job array.
         */
        public long jobId;

        /**
         * < The whole dependency condition of the job.
         */
        public String dependcondition;

        /**
         * < Whether the condition is satisfied.
         */
        public int satisfied;
    }



/* for hasDependency */
    /**
     *  \addtogroup job_has_depend job_has_depend
     *  options for hasDependency
     */

    /**
     * <  Job has dependency
     */
    public static final int JOB_HAS_DEPENDENCY = 0x1;

    /**
     * <  Job has individual  condition.
     */
    public static final int JOB_HAS_INDIVIDUAL_CONDITION = 0x2;

    /**
     * \brief  dependency jobs.
     */

    public static class dependJobs extends Structure {
        public static class ByReference extends dependJobs implements Structure.ByReference {}
        public static class ByValue extends dependJobs implements Structure.ByValue {}
        public dependJobs() {}
        public dependJobs(Pointer p) { super(p); read(); }

        /**
         * < Job ID. By default, it is the parent job of the queried job. Modify to child job by setting QUERY_DEPEND_CHILD in options of JobDepRequest.
         */
        public long jobId;

        /**
         * < The job name associated with the job ID.
         */
        public String jobname;

        /**
         * < The number of degrees of separation from the original job.
         */
        public int level;

        /**
         * < Job status of the job.
         */
        public int jobstatus;

        /**
         * < Whether the job ID has a dependency or not. When you set QUERY_DEPEND_RECURSIVELY in options of JobDepRequest, 0 indicates job ID does not have a dependency. Otherwise, one or more of the following bits displays:-  JOB_HAS_DEPENDENCY: Job has a dependency.-  JOB_HAS_INDIVIDUAL_CONDITION: Job has an individual dependency condition when it is an element of job array.
         */
        public int hasDependency;

        /**
         * < When you set "QUERY_DEPEND_DETAIL" into options, it is dependency condition of jobId. It is "" when you do not set "QUERY_DEPEND_DETAIL".
         */
        public String condition;

        /**
         * < Whether the condition is satisfied.
         */
        public int satisfied;

        /**
         * < Job ID. By default, it is the child job. Modify to parent job by setting QUERY_DEPEND_CHILD in options of JobDepRequest
         */
        public long depJobId;
    }



    /**
     * \brief  job dependent information.
     */

    public static class jobDependInfo extends Structure {
        public static class ByReference extends jobDependInfo implements Structure.ByReference {}
        public static class ByValue extends jobDependInfo implements Structure.ByValue {}
        public jobDependInfo() {}
        public jobDependInfo(Pointer p) { super(p); read(); }


        /**
         * < You can set the following bits into this field:\n QUERY_DEPEND_RECURSIVELY\n Query the dependency information recursively.\n QUERY_DEPEND_DETAIL\n Query the detailed dependency information.\n QUERY_DEPEND_UNSATISFIED\n Query the jobs that cause this job pend.\n QUERY_DEPEND_CHILD\n Query child jobs.
         */
        public int options;

        /**
         * < The number of jobs you queried. By default, the value is 1. However, when you set QUERY_DEPEND_DETAIL in the options and you query a job array where some elements have a dependency condition that has changed, the value is the number of the changed element + 1.
         */
        public int numQueriedJobs;

        /**
         * < The jobs you queried.
         */
        public Pointer /* queriedJobs.ByReference */ queriedJobs;

        /**
         * < The number of levels returned.
         */
        public int level;

        /**
         * < The number of jobs returned.
         */
        public int numJobs;

        /**
         * < The returned dependency jobs.
         */
        public Pointer /* dependJobs.ByReference */ depJobs;
    }




/*
*  Functional prototypes of the Advance Reservation API
 */


/* Macros */

    public static boolean IS_PEND(int s) {
        return (JNAUtils.toBoolean((s) & JOB_STAT_PEND) || JNAUtils.toBoolean((s) & JOB_STAT_PSUSP));
    }

/* Do not test JOB_STAT_UNKWN in IS_START() */

    public static boolean IS_START(int s) {
        return (JNAUtils.toBoolean((s) & JOB_STAT_RUN) || JNAUtils.toBoolean((s) & JOB_STAT_SSUSP) || JNAUtils.toBoolean((s) & JOB_STAT_USUSP));
    }

    public static boolean IS_FINISH(int s) {
        return (JNAUtils.toBoolean((s) & JOB_STAT_DONE) || JNAUtils.toBoolean((s) & JOB_STAT_EXIT));
    }

    public static boolean IS_SUSP(int s) {
        return (JNAUtils.toBoolean((s) & JOB_STAT_PSUSP) || JNAUtils.toBoolean((s) & JOB_STAT_SSUSP) || JNAUtils.toBoolean((s) & JOB_STAT_USUSP));
    }

/* Macro for checking post job process. (IO_SPOOL) */

    public static boolean IS_POST_DONE(int s) {
        return (((s) & JOB_STAT_PDONE) == JOB_STAT_PDONE);
    }

    public static boolean IS_POST_ERR(int s) {
        return (((s) & JOB_STAT_PERR) == JOB_STAT_PERR);
    }

    public static boolean IS_POST_FINISH(int s) {
        return (IS_POST_DONE(s) || IS_POST_ERR(s));
    }

/*On windows ,for dll library ,need to use _declspec(dllexport) to export
*a symbol .but if do so ,static library will can not work .so we are going
*to change lsberrno to a function.
*/

    public static int lsberrno() {
        return lsb_errno().getValue();
    }




/*
*  Version of the mbatchd that was last contacted.
*  -1 indicates the mbatchd has not been contacted.
 */
    //public int lsb_mbd_version;

/*
*  The data definition for host name list operations
 */
    public static final int PRINT_SHORT_NAMELIST = 0x01;
    public static final int PRINT_LONG_NAMELIST = 0x02;
    public static final int PRINT_MCPU_HOSTS = 0x04;

    public static class nameList extends Structure {
        public static class ByReference extends nameList implements Structure.ByReference {}
        public static class ByValue extends nameList implements Structure.ByValue {}
        public nameList() {}
        public nameList(Pointer p) { super(p); read(); }


/* number of names */
        public int listSize;

/* a group of names */
        public Pointer names;

/* the ocurrent of corresponding name */
        public IntByReference counter;
    }



    public static native nameList.ByReference lsb_parseShortStr(String string1, int int1);

    public static native nameList.ByReference lsb_parseLongStr(String string1);

    public static native String lsb_printNameList(nameList namelist1, int int1);

    public static native nameList.ByReference lsb_compressStrList(Pointer stringArray1, int int1);

    public static native String lsb_splitName(String string1, IntByReference int1);

    public static native IntByReference lsb_errno();


/* external routines related to API_CONF */

    public static native paramConf.ByReference lsb_readparam(LibLsf.lsConf lsConf1);

    public static native userConf.ByReference lsb_readuser(LibLsf.lsConf lsConf1, int int1, LibLsf.clusterConf clusterConf1);

    public static native userConf.ByReference lsb_readuser_ex(LibLsf.lsConf lsConf1, int int1, LibLsf.clusterConf clusterConf1, LibLsf.sharedConf sharedConf1);

    public static native hostConf.ByReference lsb_readhost(LibLsf.lsConf lsConf1, LibLsf.lsInfo lsInfo1, int int1, LibLsf.clusterConf clusterConf1);

    public static native queueConf.ByReference lsb_readqueue(LibLsf.lsConf lsConf1, LibLsf.lsInfo lsInfo1, int int1, LibLsf.sharedConf sharedConf1, LibLsf.clusterConf clusterConf1);

    public static native void updateClusterConf(LibLsf.clusterConf clusterConf1);

/* end of external routines related to API_CONF */

    /**
     * \page lsb_hostpartinfo lsb_hostpartinfo
     * Returns informaton about host partitions.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * hostPartInfoEnt.ByReference lsb_hostpartinfo (String[] hostParts,
     * IntByReference numHostParts)</b> @param hostParts An array of host partition names.
     *
     * @return null
     *         \n Function failed.
     *         <p/>
     *         <b>Errors:</b>
     *         \par
     *         If the function fails, lsberrno is set to indicate the error. If lsberrno is
     *         LSBE_BAD_HPART, (*hostParts)[*numHostParts] is not a host partition known
     *         to the LSF system. Otherwise, if.ByReference numHostParts is less than its original value,
     *         * numHostParts is the actual number of host partitions found.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts
     * @param numHostHosts The number of host partition names.
     * To get information on all host partitions, set hostParts to null;* numHostParts
     * will be the actual number of host partitions when this call returns.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * hostPartInfoEnt
     * \n hostPartUserInfo
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_usergrpinfo
     * #see \ref lsb_hostgrpinfo
     * @param stringArray1 stringArray1
     */
    public static native hostPartInfoEnt.ByReference lsb_hostpartinfo(Pointer stringArray1, IntByReference numHostHosts);

    /**
     * \page lsb_init lsb_init
     * \brief Initializes the LSF batch library (LSBLIB), and gets the
     * configuration environment.
     * <p/>
     * You must use \ref lsb_init before any other LSBLIB library routine in your
     * application.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_init(String appname)</b>
     *
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         <b>Errors:</b>
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param appName The name of your application.
     * If appName holds the name of your application, a logfile with the same
     * name as
     * your application receives LSBLIB transaction information.
     * If appName is null, the logfile $LSF_LOGDIR/bcmd receives LSBLIB
     * transaction information.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * see none
     */
    public static native int lsb_init(String appName);

    public static native int sch_lsb_init();

    /**
     * \page lsb_openjobinfo lsb_openjobinfo
     * \brief Returns the number of jobs in the master batch daemon.
     * <p/>
     * \ref lsb_openjobinfo accesses information about pending, running and
     * suspended jobs in the master batch daemon. Use \ref lsb_openjobinfo to
     * create a connection to the master batch daemon. Next, use \ref lsb_readjobinfo
     * to read job records.Close the connection using \ref lsb_closejobinfo.
     * <p/>
     * \ref lsb_openjobinfo opens a connection with mbatchd and returns the total
     * number of records in the connection on success.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_openjobinfo(long jobId, String jobName,
     * String userName, String queueName, String hostName,
     * int options)</b>
     *
     * @param jobId   Passes information about jobs with the given job ID.
     *                If jobId is 0, \ref lsb_openjobinfo looks to another parameter to return
     *                information about jobs.If a member of a job array is to be passed, use
     *                the array form jobID[ i ] where jobID is the job array name, and i is
     *                the index value.
     * @param options <lsf/lsbatch.h> defines the flags shown in
     *                \ref defs_lsb_openjobinfo constructed from bits. Use the bitwise OR to set more
     *                than one flag.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                none
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref defs_lsb_openjobinfo_a
     *                \n \ref defs_lsb_openjobinfo
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         bjobs
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param jobName Passes information about jobs with the given job name.
     * If jobName is null, \ref lsb_openjobinfo looks to another parameter to return
     * information about jobs.
     * @param userName Passes information about jobs submitted by the named user
     * or user group, or by all users if user is all. If user is null,
     * \ref lsb_openjobinfo assumes the user is invoking this call.
     * @param queueName Passes information about jobs belonging to the named
     * queue. If queue is null,jobs in all the queues of the batch system are counted.
     * @param hostName Passes information about jobs on the named host, host
     * group or cluster name. If host is null, jobs on all hosts of the batch
     * system will be considered.
     * #see \ref               lsb_openjobinfo_a
     * #see \ref               lsb_openjobinfo_a_ext
     * #see \ref               lsb_openjobinfo_req
     * #see \ref               lsb_closejobinfo
     * #see \ref               lsb_readjobinfo
     * #see \ref               lsb_readframejob
     */
    public static native int lsb_openjobinfo(long jobId, String jobName, String userName, String queueName, String hostName, int options);

    /**
     * \page lsb_openjobinfo_a lsb_openjobinfo_a
     * \brief Provides the name and number of jobs and hosts in the master batch
     * daemon.
     * <p/>
     * \ref lsb_openjobinfo_a provides more information on pending, running and
     * suspended jobs than \ref lsb_openjobinfo. Use \ref lsb_openjobinfo_a to create a
     * connection to the master batch daemon. Next, use \ref lsb_readjobinfo to read
     * job records. Close the connection using \ref lsb_closejobinfo.
     * <p/>
     * \ref lsb_openjobinfo_a passes information about jobs based on the value of
     * jobId,jobName, userName, queueName, or hostName. Only one parameter can be
     * chosen. The other parameters must be null or 0.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * jobInfoHead.ByReference lsb_openjobinfo_a(long jobId,
     * String jobName,
     * String userName,
     * String queueName,
     * String hostName,
     * int options)</b>
     *
     * @param jobId   Passes information about jobs with the given job ID. If jobId
     *                is 0, \ref lsb_openjobinfo looks to another parameter to return information
     *                about jobs.
     *                If information about a member of a job array is to be passed, use the array
     *                form jobID[ i ] where jobID is the job array name, and i is the index value.
     * @param options <lsf/lsbatch.h> defines the flags shown in def_lsb_openjobinfo_a
     *                constructed from bits. Use the bitwise OR to set more than one flag.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                jobInfoHead
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref defs_lsb_openjobinfo_a
     * @return null \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         bjobs
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param jobName Passes information about jobs with the given job name. If
     * jobName is null, \ref lsb_openjobinfo looks to another parameter to return
     * information about jobs.
     * @param userName Passes information about jobs submitted by the named user
     * or user group, or by all users if userName is all. If userName is null,
     * \ref lsb_openjobinfo_a assumes the user is invoking this call.
     * @param queueName Passes information about jobs belonging to the named queue.
     * If queueName is null, jobs in all queues of the batch system will be
     * considered.
     * @param hostName Passes information about jobs on the named host, host group
     * or cluster name. If hostName is null, jobs on all hosts of the batch system
     * will be considered.
     * #see \ref lsb_openjobinfo
     * #see \ref lsb_closejobinfo
     * #see \ref lsb_readjobinfo
     * #see \ref lsb_readframejob
     */
    public static native jobInfoHead.ByReference lsb_openjobinfo_a(long jobId, String jobName, String userName, String queueName, String hostName, int options);

    /**
     * \page lsb_openjobinfo_a_ext lsb_openjobinfo_a_ext
     * \brief  Returns the name and number of jobs and hosts in the master batch
     * daemon with additional host group information.
     * <p/>
     * \ref lsb_openjobinfo_a_ext is run from \ref lsb_openjobinfo_a using the same
     * parameters and provides the same information as \ref lsb_openjobinfo_a, but with
     * additional host group information.
     * <p/>
     * \ref lsb_openjobinfo_a_ext passes information about jobs based on the value of
     * jobId, jobName, userName, queueName, or hostName. Only one parameter can be
     * chosen. The other parameters must be null or 0.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * jobInfoHeadExt.ByReference
     * lsb_openjobinfo_a_ext (long jobId, String jobName,
     * String userName, String queueName,
     * String hostName, int options)</b>
     *
     * @param jobId   Passes information about jobs with the given job ID. If jobId
     *                is 0, \ref lsb_openjobinfo_a_ext looks to another parameter to return information
     *                about jobs. If information about a member of a job array is to be passed, use
     *                the array form jobID[ i ] where jobID is the job array name, and i is the
     *                index value.
     * @param options <lsf/lsbatch.h> defines the flags shown in
     *                def_lsb_openjobinfo_a constructed from bits. Use the bitwise OR to set more
     *                than one flag.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                jobInfoHeadExt
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref defs_lsb_openjobinfo_a
     * @return null \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         bjobs
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param jobName Passes information about jobs with the given job name. If
     * jobName is null, \ref lsb_openjobinfo_a_ext looks to another parameter to return
     * information about jobs.
     * @param userName Passes information about jobs submitted by the named user
     * or user group, or by all users if userName is all. If userName is null,
     * \ref lsb_openjobinfo_a_ext assumes the user is invoking this call.
     * @param queueName Passes information about jobs belonging to the named queue.
     * If queueName is null, jobs in all queues of the batch system will be considered.
     * @param hostName Passes information about jobs on the named host, host group
     * or cluster name. If hostName is null, jobs on all hosts of the batch system
     * will be considered.
     * #see \ref lsb_openjobinfo
     * #see \ref lsb_closejobinfo
     * #see \ref lsb_readjobinfo
     * #see \ref lsb_readframejob
     */
    public static native jobInfoHeadExt.ByReference lsb_openjobinfo_a_ext(long jobId, String jobName, String userName, String queueName, String hostName, int options);

    /**
     * \page lsb_openjobinfo_req lsb_openjobinfo_req
     * \brief  Extensible API.
     * <p/>
     * Instead of submitting individual requests this API defines
     * all job info requests as objects, and can easily be enhanced to include
     * additinal requests.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * jobInfoHeadExt.ByReference lsb_openjobinfo_req (jobInfoReq.ByReference req)</b>
     *
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         none
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param req  job information request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * jobInfoReq
     * \n \ref jobInfoHeadExt
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref defs_lsb_openjobinfo_a
     * \n \ref defs_lsb_openjobinfo
     * #see \ref               lsb_openjobinfo_a
     * #see \ref               lsb_openjobinfo_a_ext
     * #see \ref               lsb_closejobinfo
     * #see \ref               lsb_readjobinfo
     * #see \ref               lsb_readframejob
     */
    public static native jobInfoHeadExt.ByReference lsb_openjobinfo_req(jobInfoReq req);

    public static native int lsb_queryjobinfo(int int1, NativeLongByReference long1, String string1);

    public static native jobInfoEnt.ByReference lsb_fetchjobinfo(IntByReference int1, int int2, NativeLongByReference long1, String string1);

    public static native jobInfoEnt.ByReference lsb_fetchjobinfo_ext(IntByReference int1, int int2, NativeLongByReference long1, String string1, jobInfoHeadExt jobInfoHeadExt);

    /**
     * \page lsb_readjobinfo lsb_readjobinfo
     * \brief Returns the next job information record in mbatchd.
     * <p/>
     * \ref lsb_readjobinfo reads the number of records defined by the more parameter.
     * The more parameter receives its value from either \ref lsb_openjobinfo or
     * \ref lsb_openjobinfo_a. Each time \ref lsb_readjobinfo is called, it returns one
     * record from mbatchd. Use \ref lsb_readjobinfo in a loop and use more to
     * determine how many times to repeat the loop to retrieve job information records.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * \n \#include <time.h>
     * \n \#include <lsf/lsf.h>
     * <p/>
     * jobInfoEnt.ByReference lsb_readjobinfo(IntByReference more)</b>
     *
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If there are no more records, then lsberrno is set to LSBE_EOF.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.queues
     * @param more Number of job records in the master batch daemon.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * jobInfoEnt
     * \n jobExternalMsgReply
     * \n jRusage
     * \n pidInfo
     * \n reserveItem
     * \n submit
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref job_states
     * \n \ref jobgroup_counterIndex
     * \n \ref group_nodetypes
     * #see \ref lsb_openjobinfo
     * #see \ref lsb_openjobinfo_a
     * #see \ref lsb_closejobinfo
     * #see \ref lsb_hostinfo
     * #see \ref lsb_pendreason
     * #see \ref lsb_queueinfo
     * #see \ref lsb_suspreason
     */
    public static native jobInfoEnt.ByReference lsb_readjobinfo(IntByReference more);

    /**
     * \page  lsb_submit lsb_submit
     * Submits or restarts a job in the batch system.
     * <p/>
     * \ref lsb_submit submits or restarts a job in the batch system according to the
     * jobSubReq specification.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * long lsb_submit (submit.ByReference jobSubReq,
     * submitReply.ByReference jobSubReply)</b>
     *
     * @return long:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         If the environment variable BSUB_CHK_RESREQ is set, the value of lsberrno is
     *         either LSBE_RESREQ_OK or LSBE_RESREQ_ERR, depending on the result of
     *         resource requirement string checking. The badJobName field in the submitReply
     *         structure contains the detailed error message.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bsub
     *         \n brestart
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param jobSubReq
     * Describes the requirements for job submission to the batch system.
     * A job that does not meet these requirements is not submitted to the
     * batch system and an error is returned.
     * @param jobSubReply
     * Describes the results of the job submission to the batch system.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * submit
     * \n submitReply
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref lsb_submit_options
     * \n \ref lsb_submit_options2
     * \n \ref lsb_submit_options3
     * #see \ref lsb_modify
     * #see \ref ls_info
     * #see \ref lsb_queueinfo
     */
    public static native long lsb_submit(submit jobSubReq, submitReply jobSubReply);

    /**
     * \page lsb_readjobinfo_cond lsb_readjobinfo_cond
     * \brief Returns the next job information record for condensed host groups
     * in mbatchd.
     * <p/>
     * \ref lsb_readjobinfo_cond reads the number of records defined by the more
     * parameter. The more parameter receives its value from either \ref lsb_openjobinfo
     * or \ref lsb_openjobinfo_a. Each time \ref lsb_readjobinfo_cond is called, it
     * returns one record from mbatchd. Use \ref lsb_readjobinfo_cond in a loop and use
     * more to determine how many times to repeat the loop to retrieve job information
     * records.
     * <p/>
     * \ref lsb_readjobinfo_cond differs from \ref lsb_readjobinfo in that if jInfoHExt
     * is not null, \ref lsb_readjobinfo_cond substitutes hostGroup (if it is a condensed
     * host group) for job execution hosts.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * \n \#include <time.h>
     * \n \#include <lsf/lsf.h>
     * <p/>
     * jobInfoEnt.ByReference lsb_readjobinfo_cond(IntByReference more,
     * jobInfoHeadExt.ByReference jInfoHExt);</b>
     *
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If there are no more records, then lsberrno is set to LSBE_EOF.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.queues
     * @param more Number of job records in the master batch daemon.
     * @param jInfoHExt Job information header info for the condensed host group.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * jobInfoEnt
     * \n jobExternalMsgReply
     * \n jRusage
     * \n pidInfo
     * \n reserveItem
     * \n submit
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref external_msg_processing
     * \n \ref group_nodetypes
     * #see \ref lsb_openjobinfo
     * #see \ref lsb_openjobinfo_a
     * #see \ref lsb_closejobinfo
     * #see \ref lsb_hostinfo
     * #see \ref lsb_pendreason
     * #see \ref lsb_queueinfo
     * #see \ref lsb_readjobinfo
     * #see \ref lsb_suspreason
     */
    public static native jobInfoEnt.ByReference lsb_readjobinfo_cond(IntByReference more, jobInfoHeadExt jInfoHExt);

    /**
     * \page lsb_readframejob lsb_readframejob
     * \brief Returns all frame jobs information which matchs the specified
     * parameters and fills related information into the frame job information table.
     * <p/>
     * \ref lsb_readframejob gets all frame jobs information that matches the specified
     * parameters and fills related information into the frame job information table.
     * \ref lsb_readframejob is a wrapper of \ref lsb_openjobinfo, \ref lsb_readjobinfo, and
     * \ref lsb_closejobinfo. Memory allocated in frameJobInfoTbl will be freed by
     * user.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_readframejob(long jobId, String frameName,
     * String user, String queue, String host, int options,
     * frameJobInfo.ByReference[] frameJobInfoTbl)</b>
     *
     * @param jobId   Get information about the frame jobs with the given job ID.
     *                If jobID is 0, get information about frame jobs which satisfy the other
     *                specifications. If a job in a job array is to be modified, use the array
     *                form jobID[i] where jobID is the job array name, and i is the index value.
     * @param options <lsf/lsbatch.h> defines the following flags \ref defs_lsb_openjobinfo_a
     *                constructed from bits. Use the bitwise OR to set more than one flag.
     * @return int:-1
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param frameName Get information about frame jobs with the given frame name.
     * @param user Get information about frame jobs submitted by the named user
     * or user group, or by all users if user is all. If user is null, the user
     * invoking this routine is assumed.
     * @param queue Get information about frame jobs belonging to the named queue.
     * If queue is null,jobs in all queues of the batch system will be considered.
     * @param host Get information about frame jobs on the named host, host
     * group or cluster name.If host is null, jobs on all hosts of the batch
     * system will be considered.
     * @param frameJobInfoTbl The result of all frame jobs information.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * \n frameJobInfo
     * \n frameElementInfo
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_openjobinfo
     * #see \ref lsb_readjobinfo
     * #see \ref lsb_closejobinfo
     */

    public static native int lsb_readframejob(long jobId, String frameName, String user, String queue, String host, int options, Pointer frameJobInfoTbl);

    /**
     * \page lsb_closejobinfo lsb_closejobinfo
     * \brief Closes job information connection with the master batch daemon.
     * <p/>
     * Use \ref lsb_closejobinfo to close the connection to the master batch daemon
     * after opening a job information connection with \ref lsb_openjobinfo and reading
     * job records with \ref lsb_readjobinfo.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * void lsb_closejobinfo()</b>
     *
     * param void \n
     *             <p/>
     *             <b>Data Structures:</b>
     *             \par
     *             none
     *             <p/>
     *             <b>Define Statements:</b>
     *             \par
     *             none
     * return void
     *         \n There's no returns value.
     *         <p/>
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * #see \ref lsb_openjobinfo
     * #see \ref lsb_openjobinfo_a
     * #see \ref lsb_readjobinfo
     */

    public static native void lsb_closejobinfo();

    /**
     * \page  lsb_hostcontrol lsb_hostcontrol
     * Opens or closes a host, or restarts or shuts down its slave batch daemon.
     * <p/>
     * \ref lsb_hostcontrol opens or closes a host, or restarts or shuts down its
     * slave batch daemon. Any program using this API must be setuid to root if
     * LSF_AUTH is not defined in the lsf.conf file.
     * <p/>
     * To restart the master batch daemon, mbatchd, in order to use updated
     * batch LSF configuration files, use \ref lsb_reconfig.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_hostcontrol (hostCtrlReq.ByReference req)</b>
     *
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         <b>Errors:</b>
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param req The host control request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * hostCtrlReq
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref host_ctrl_option
     * #see \ref lsb_reconfig
     */
    public static native int lsb_hostcontrol(hostCtrlReq req);

    public static native int lsb_hghostcontrol(hgCtrlReq hostCtrlReq1, hgCtrlReply reply);

    /**
     * \page lsb_queueinfo lsb_queueinfo
     * \brief Returns information about batch queues.
     * <p/>
     * \ref lsb_queueinfo gets information about batch queues. See lsb.queues for more
     * information about queue parameters.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * queueInfoEnt.ByReference lsb_queueinfo(String[] queues,
     * IntByReference numQueues, String hosts, String users,
     * int options)</b>
     *
     * @param options Reserved for future use; supply 0.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                queueInfoEnt
     *                \n shareAcctInfoEnt
     *                \n apsFactorInfo
     *                \n apsFactorMap
     *                \n apsLongNameMap
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref queue_status
     *                \n \ref queue_attribute
     * @return null
     *         \n Function Failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         If lsberrno is LSBE_BAD_QUEUE, (*queues)[*numQueues] is not a queue known
     *         to the LSF system. Otherwise, if.ByReference numQueues is less than its original value,
     *         * numQueues is the actual number of queues found.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bqueues
     *         <p/>
     *         \b Files:
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.queues
     * @param queues An array of names of queues of interest.
     * @param numQueues The number of queue names. To get information on all queues,
     * set.ByReference numQueues to 0;* numQueues will be updated to the actual number of
     * queues when this call returns.If.ByReference numQueues is 1 and queues is null,
     * information on the system default queue is returned.
     * @param hosts The host or cluster names. If hosts is not null, then only
     * the queues that are enabled for the hosts are of interest.
     * @param user The name of user. If user is not null, then only the queues
     * that are enabled for the user are of interest.
     * #see \ref lsb_hostinfo
     * #see \ref lsb_userinfo
     * #see \ref lsb_usergrpinfo
     */
    public static native queueInfoEnt.ByReference lsb_queueinfo(Pointer queues, IntByReference numQueues, String hosts, String user, int options);

    /**
     * \page lsb_reconfig lsb_reconfig
     * \brief Dynamically reconfigures an LSF batch system.
     * <p/>
     * \ref lsb_reconfig dynamically reconfigures an LSF batch system to pick up new
     * configuration parameters and changes to the job queue setup since system
     * startup or the last reconfiguration (see lsb.queues).
     * <p/>
     * To restart a slave batch daemon, use \ref lsb_hostcontrol. This call is
     * successfully invoked only by root or by the LSF administrator.
     * <p/>
     * Any program using this API must be setuid to root if LSF_AUTH is not
     * defined in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_reconfig (mbdCtrlReq.ByReference req)</b>
     *
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         badmin reconfig
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param req mbatchd control request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * mbdCtrlReq
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref mbd_operation
     * #see \ref lsb_openjobinfo
     */
    public static native int lsb_reconfig(mbdCtrlReq req);

    /**
     * \page lsb_signaljob lsb_signaljob
     * \brief Sends a signal to a job.
     * <p/>
     * Use \ref lsb_signaljob when migrating a job from one host to another. Use
     * \ref lsb_signaljob to stop or kill a job on a host before using \ref lsb_mig to
     * migrate the job. Next, use \ref lsb_signaljob to continue the stopped job at
     * the specified host.
     * <p/>
     * Generally, use \ref lsb_signaljob to apply any UNIX signal to a job or process.
     * <p/>
     * Any program using this API must be setuid to root if LSF_AUTH is not defined
     * in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_signaljob (long jobId, int sigValue)</b>
     *
     * @param jobId    The job to be signaled. If a job in a job array is to be
     *                 signaled, use the array form jobID[ i ] where jobID is the job array name,
     *                 and i is the index value.
     * @param sigValue SIGSTOP, SIGCONT, SIGKILL or some other UNIX signal.
     *                 <p/>
     *                 <b>Data Structures:</b>
     *                 \par
     *                 none
     *                 <p/>
     *                 <b>Define Statements:</b>
     *                 \par
     *                 none
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bkill \n
     *         bstop \n
     *         bresume
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * #see \ref lsb_chkpntjob
     * #see \ref lsb_forcekilljob
     * #see \ref lsb_mig
     */

    public static native int lsb_signaljob(long jobId, int sigValue);

    /**
     * \page lsb_killbulkjobs lsb_killbulkjobs
     * \brief Kills bulk jobs as soon as possible.
     * <p/>
     * Use \ref lsb_killbulkjobs to kill bulk jobs on a local host immediately, or
     * to kill other jobs as soon as possible. If mbatchd rejects the request, it
     * issues null as the reservation ID.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_killbulkjobs(signalBulkJobs.ByReference s)</b>
     *
     * @return int:-1 \n
     *         The bulk jobs were not killed.
     *         <p/>
     *         \b Error:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         bkill -b
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param s The signal to a group of jobs.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * signalBulkJobs
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * see none
     */

    public static native int lsb_killbulkjobs(signalBulkJobs s);

    public static native int lsb_msgjob(long long1, String s);

    /**
     * \page lsb_chkpntjob lsb_chkpntjob
     * \brief Checkpoints a job.
     * <p/>
     * Checkpoints a job.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_chkpntjob(long jobId, int period, int options)</b>
     *
     * @param jobId   The job to be checkpointed.
     * @param period  The checkpoint period in seconds. The value 0
     *                disables periodic checkpointing.
     * @param options The bitwise inclusive OR of some of the following:
     *                \n LSB_CHKPNT_KILL
     *                Checkpoint and kill the job as an atomic action.
     *                \n LSB_CHKPNT_FORCE
     *                Checkpoint the job even if non-checkpointable conditions exist.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                none
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref chkpnt_job_option
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \note Any program using this API must be setuid to root if LSF_AUTH
     *         is not defined in the lsf.conf file.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bchkpnt
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * see none
     */
    public static native int lsb_chkpntjob(long jobId, NativeLong period, int options);

    /**
     * \page lsb_deletejob lsb_deletejob
     * \brief Kills a job in a queue
     * <p/>
     * Use \ref lsb_deletejob to send a signal to kill a running, user-suspended,
     * or system-suspended job. The job can be requeued or deleted from the batch
     * system.If the job is requeued, it retains its submit time but it is dispatched
     * according to its requeue time. When the job is requeued, it is assigned the
     * PEND status and re-run.If the job is deleted from the batch system, it is
     * no longer available to be requeued.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_deletejob (long jobId, int times, int options)</b>
     *
     * @param jobId   The job to be killed. If an element of a job array is to be
     *                killed, use the array form jobID[i] where jobID is the job array name,
     *                and i is the index value.
     * @param times   Original job submit time.
     * @param options If the preprocessor macro LSB_KILL_REQUEUE in lsbatch.h is
     *                compared with options and found true, then requeue the job using the same job ID.
     *                If the preprocessor macro LSB_KILL_REQUEUE in lsbatch.h is compared with
     *                options and found false, then the job is deleted from the batch system.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                none
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref kill_requeue
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \note Any program using this API must be setuid to root if LSF_AUTH is not defined in the
     *         \n lsf.conf file.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bkill
     *         \n brequeue -J
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * #see \ref lsb_signaljob
     * #see \ref lsb_chkpntjob
     */
    public static native int lsb_deletejob(long jobId, int times, int options);

    /**
     * \page lsb_forcekilljob lsb_forcekilljob
     * \brief This function is used to send special force kill signal.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_forcekilljob(long jobId)</b>
     *
     * @param jobId which job is to be killed.
     *              <p/>
     *              <b>Data Structures:</b>
     *              \par
     *              none
     *              <p/>
     *              <b>Define Statements:</b>
     *              \par
     *              none
     * @return int:-1
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * #see \ref lsb_signaljob
     */
    public static native int lsb_forcekilljob(long jobId);

    /**
     * \page lsb_submitframe lsb_submitframe
     * \brief Submits a frame job to the batch system.
     * <p/>
     * \ref lsb_submitframe submits a frame job to the batch system according to the
     * jobSubReq specification and frameExp.
     * <p/>
     * Any program using this API must be setuid to root if LSF_AUTH is not defined
     * in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_submitframe (submit.ByReference jobSubReq, String frameExp,
     * submitReply.ByReference jobSubReply)</b>
     *
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error and jobSubReply gives
     *         additional information about the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param jobSubReq Describes the requirements for job submission to the
     * batch system. A job that does not meet these requirements is not submitted
     * to the batch system and an error is returned. \n
     * See \ref lsb_submit for descriptions of the submit structure fields.
     * @param frameExp The syntax of frameExp is: \n
     * <b>frame_name[indexlist]</b> \n
     * frame_name is any name consisting of alphanumerics, periods, forward slashes,
     * dashes or underscores. indexlist is a list of one or more frame indexes,
     * separated by commas. These indexes can each be either a single integer or
     * a range, specified in the following format: \n
     * <b>start-end[xstep[:chunk]]</b> \n
     * start, end, step, and chunk are integers, but chunk must be positive.
     * If step and
     * chunk are ommitted, the default value is 1.\n
     * An example of a valid expression for frameExp is:\n
     * <b>Frame_job_1[5,10-15,20-30x2:3]</b>
     * @param jobSubReply Describes the results of the job submission to the
     * batch system. \n
     * See \ref lsb_submit for descriptions of the submitReply structure
     * fields.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * submit
     * \n submitReply
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref lsb_submit_options
     * \n \ref lsb_submit_options2
     * \n \ref lsb_submit_options3
     * see none
     */
    public static native int lsb_submitframe(submit jobSubReq, String frameExp, submitReply jobSubReply);

    /**
     * \page lsb_requeuejob lsb_requeuejob
     * \brief Requeues job arrays, jobs in job arrays, and individual jobs.
     * <p/>
     * Use \ref lsb_requeuejob to requeue job arrays, jobs in job arrays, and individual
     * jobs that are running, pending, done, or exited. In a job array, you can
     * requeue all the jobs or requeue individual jobs of the array.
     * <p/>
     * \ref lsb_requeuejob requeues jobs as if the jobs were in an array. A job not in an
     * array is considered to be a job array composed of one job.
     * <p/>
     * Jobs in a job array can be requeued independently of each other regardless of
     * any job's status (running, pending, exited, done). A requeued job is requeued
     * to the same queue it was originally submitted from or switched to. The job
     * submission time does not change so a requeued job is placed at the top of the
     * queue. Use \ref lsb_movejob to place a job at the bottom or any other position
     * in a queue.
     * <p/>
     * If a clean period is reached before \ref lsb_requeuejob is called, the cleaned
     * jobs cannot be requeued. Set the variable CLEAN_PERIOD in your lsb.params file
     * to determine the amount of time that job records are kept in MBD core memory
     * after jobs have finished or terminated.
     * <p/>
     * To requeue a job assign values to the data members of the jobrequeue data
     * structure, process command line options in case the user has specified a
     * different job, and call \ref lsb_requeuejob to requeue the job array.
     * <p/>
     * Assign values to the jobID, status, and options data members of the jobrequeue
     * data structure. Assign the job identification number to jobID. Assign
     * JOB_STAT_PEND or JOB_STAT_PSUSP to status. Assign REQUEUE_DONE, REQUEUE_EXIT,
     * and or REQUEUE_RUN to requeue running jobs.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_requeuejob(jobrequeue.ByReference  reqPtr)</b>
     *
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         brequeue -d
     *         \n brequeue -e
     *         \n brequeue -a
     *         \n brequeue -r
     *         \n brequeue -H
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     *         \n $LSB_SHAREDIR
     * @param reqPtr This structure contains the information required to requeue a job.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * jobrequeue
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref requeuejob_options
     * #see \ref lsb_movejob
     * #see \ref lsb_pendreason
     */
    public static native int lsb_requeuejob(jobrequeue reqPtr);

    /**
     * \page lsb_sysmsg lsb_sysmsg
     * \brief Returns a pointer to static data.
     * <p/>
     * \ref lsb_sysmsg returns a pointer to static data which stores the batch error
     * message corresponding to lsberrno. The global variable lsberrno maintained
     * by LSBLIB holds the error number from the most recent LSBLIB call that caused
     * an error. If lsberrno == LSBE_SYS_CALL, then the system error message defined
     * by errno is also returned. If lsberrno == LSBE_LSLIB, then the error message
     * returned by \ref ls_sysmsg is returned.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * String lsb_sysmsg ()</b>
     *
     * param void \n
     *             <p/>
     *             <b>Data Structures:</b>
     *             \par
     *             none
     *             <p/>
     *             <b>Define Statements:</b>
     *             \par
     *             none
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * #see \ref ls_perror
     * #see \ref ls_sysmsg
     */
    public static native String lsb_sysmsg();

    /**
     * \page lsb_perror lsb_perror
     * \brief Prints a batch LSF error message on stderr.
     * <p/>
     * \ref lsb_perror prints a batch LSF error message on stderr. The usrMsg is
     * printed out first, followed by a ":" and the batch error message corresponding
     * to lsberrno.
     * <p/>
     * \ref lsb_perror - Print LSBATCH error message on stderr. In addition
     * to the error message defined by lsberrno, user supplied message usrMsg1
     * is printed out first and a ':' is added to separate.ByReference  usrMsg1 and LSBATCH
     * error message.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * void lsb_perror (String usrMsg)</b>
     *
     * return void \n
     *         Prints out the user supplied error message.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         none
     *         <p/>
     *         \b Files:
     *         \par
     *         none
     * @param usrMsg A user supplied error message.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * see none
     */
    public static native void lsb_perror(String usrMsg);

    public static native void lsb_errorByCmd(String string1, String string2, int int1);

    public static native String lsb_sperror(String string1);

    /**
     * \page lsb_peekjob lsb_peekjob
     * \brief Returns the base name of the file related to the job ID
     * <p/>
     * \ref lsb_peekjob retrieves the name of a job file.
     * <p/>
     * Only the submitter can peek at job output.
     * <p/>
     * The storage for the file name will be reused by the next call.
     * <p/>
     * Any program using this API must be setuid to root if LSF_AUTH
     * is not defined in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * String  lsb_peekjob (long jobId)</b>
     *
     * @param jobId The job ID that the LSF system assigned to the job. If a job
     *              in a job array is to be returned, use the array form jobID[i] where jobID
     *              is the job array name, and i is the index value.
     *              <p/>
     *              <b>Data Structures:</b>
     *              \par
     *              none
     *              <p/>
     *              <b>Define Statements:</b>
     *              \par
     *              none
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         bpeek
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * see none
     */
    public static native String lsb_peekjob(long jobId);

    /**
     * \page lsb_mig lsb_mig
     * \brief Migrates a job from one host to another.
     * <p/>
     * \ref lsb_mig migrates a job from one host to another. Any program using
     * this API must be setuid to root if LSF_AUTH is not defined
     * in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_mig(submig.ByReference mig, IntByReference badHostIdx)</b>
     *
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error and badHostIdx indicates
     *         which askedHost is not acceptable.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         none
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param mig The job to be migrated.
     * @param badHostIdx If the call fails, (**askedHosts)[*badHostIdx] is not a
     * host known to the LSF system.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * submig
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_submit
     */
    public static native int lsb_mig(submig mig, IntByReference badHostIdx);

    public static native clusterInfoEnt.ByReference lsb_clusterinfo(IntByReference int1, Pointer stringArray1, int int2);

    public static native clusterInfoEntEx.ByReference lsb_clusterinfoEx(IntByReference int1, Pointer stringArray1, int int2);

    /**
     * \page lsb_hostinfo lsb_hostinfo
     * Returns information about job server hosts.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * hostInfoEnt.ByReference lsb_hostinfo(String[] hosts, IntByReference numHosts)</b>
     *
     * @return hostInfoEnt.ByReference :null
     *         \n Function failed.
     *         <p/>
     *         <b>Errors:</b>
     *         \par
     *         If the function fails, lsberrno is set to indicate the error. If lsberrno is
     *         LSBE_BAD_HOST, (*hosts)[*numHosts] is not a host known to the batch system.
     *         Otherwise, if.ByReference numHosts is less than its original value,* numHosts is the actual
     *         number of hosts found.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bhosts
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts
     * @param hosts
     * An array of host or cluster names.
     * @param numHosts
     * The number of host names.
     * To get information on all hosts, set.ByReference numHosts to 0;* numHosts will be set to the
     * actual number of hostInfoEnt structures when this call returns.
     * If.ByReference numHosts is 1 and hosts is null, information on the local host is returned.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * hostInfoEnt
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref host_status
     * \n \ref host_load_BusyReason
     * \n \ref host_attributes
     * #see \ref lsb_hostinfo_ex
     * #see \ref ls_info
     * #see \ref ls_loadofhosts
     * #see \ref lsb_queueinfo
     * #see \ref lsb_userinfo
     */
    public static native hostInfoEnt.ByReference lsb_hostinfo(Pointer hosts, IntByReference numHosts);

    /**
     * \page lsb_hostinfo_ex lsb_hostinfo_ex
     * Returns informaton about job server hosts that satisfy specified resource
     * requirements. \ref lsb_hostinfo_ex returns information about job server hosts
     * that satisfy the specified resource requirements.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * hostInfoEnt.ByReference lsb_hostinfo_ex(String[] hosts,
     * IntByReference numHosts, String resReq, int options)</b> @param hosts An array of host or cluster names.
     *
     * @param options Options is reserved for the future use.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                hostInfoEnt
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref host_status
     *                \n \ref host_load_BusyReason
     *                \n \ref host_attributes
     * @return null
     *         \n Function failed.
     *         <p/>
     *         <b>Errors:</b>
     *         \par
     *         If the function fails, lsberrno is set to indicate the error. If lsberrno is
     *         LSBE_BAD_HOST, (*hosts)[*numHosts] is not a host known to the batch system.
     *         Otherwise, if.ByReference numHosts is less than its original value,* numHosts is the actual
     *         number of hosts found.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts
     * @param numHosts The number of host names.
     * To get information on all hosts, set.ByReference numHosts to 0;* numHosts will be set
     * to the actual number of hostInfoEnt structures when this call returns.
     * If.ByReference numHosts is 1 and hosts is null, information on the local host is returned.
     * @param resReq Resource requirements.
     * If this option is specified, then only host information for those hosts
     * that satisfy the resource requirements is returned. Returned hosts are
     * sorted according to the load on the resource() given in resReq, or by
     * default according to CPU and paging load.
     * #see \ref ls_info
     * #see \ref ls_loadofhosts
     * #see \ref lsb_hostinfo
     * #see \ref lsb_queueinfo
     * #see \ref lsb_userinfo
     * @param string1 string1
     */

    public static native hostInfoEnt.ByReference lsb_hostinfo_ex(Pointer resReq, IntByReference numHosts, String string1, int options);

    /**
     * \page lsb_hostinfo_cond lsb_hostinfo_cond
     * Returns condensed information about job server hosts.
     * <p/>
     * \ref lsb_hostinfo_cond returns condensed information about job server hosts.
     * While \ref lsb_hostinfo returns specific information about individual hosts,
     * \ref lsb_hostinfo_cond returns the number of jobs in each state within the
     * entire host group. The condHostInfoEnt structure contains counters that
     * indicate how many hosts are in the ok, busy, closed, full, unreach, and
     * unavail states and an array of hostInfoEnt structures that indicate the
     * status of each host in the host
     * group.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * condHostInfoEnt.ByReference  lsb_hostinfo_cond
     * (String[] hosts, IntByReference numHosts,
     * String resReq, int options)</b>
     *
     * @param options Any options called with the function.
     *                <p/>
     *                <b>Data Structures</b>
     *                \par
     *                condHostInfoEnt
     *                \n hostInfoEnt
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                none
     * @return null
     *         \n Function failed.
     *         <p/>
     *         <b Errors:</b>
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param hosts An array of host names belonging to the host group.
     * @param numHosts The number of host names in the host group.
     * To get information on all hosts in the host group, set.ByReference numHosts to 0;
     * * numHosts will be set to the actual number of hostInfoEnt structures in
     * the host group when this call returns.
     * @param resReq Any resource requirements called with the function.
     * #see \ref lsb_hostinfo
     */
    public static native condHostInfoEnt.ByReference lsb_hostinfo_cond(Pointer hosts, IntByReference numHosts, String resReq, int options);

    /**
     * \page lsb_movejob lsb_movejob
     * \brief Changes the position of a pending job in a queue.
     * <p/>
     * Use \ref lsb_movejob to move a pending job to a new position that you specify
     * in a queue. Position the job in a queue by first specifying the job ID.
     * Next, count, beginning at 1, from either the top or the bottom of the queue,
     * to the position you want to place the job.
     * <p/>
     * To position a job at the top of a queue, choose the top of a queue parameter
     * and a postion of 1.To position a job at the bottom of a queue, choose the
     * bottom of the queue parameter and a position of 1.
     * <p/>
     * By default, LSF dispatches
     * jobs in a queue in order of their arrival (such as first-come-first-served),
     * subject to the availability of suitable server hosts.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_movejob (long jobId, IntByReference position, int opCode)</b>
     *
     * @param jobId  The job ID that the LSF system assigns to the job. If a job
     *               in a job array is to be moved, use the array form jobID[ i ] where jobID is
     *               the job array name, and i is the index value.
     * @param opCode The top or bottom position of a queue.
     *               \n \b TO_TOP
     *               \n The top position of a queue.
     *               \n \b TO_BOTTOM
     *               \n The bottom position of a queue.
     *               \n If an opCode is not specified for the top or bottom position, the
     *               function fails.
     *               <p/>
     *               <b>Data Structures:</b>
     *               \par
     *               none
     *               <p/>
     *               <b>Define Statements:</b>
     *               \par
     *               \ref movejob_options
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         btop
     *         \n bbot
     *         \n bjobs -q
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param position The new position of the job in a queue. position must be
     * a value of 1 or more.
     * #see \ref lsb_pendreason
     */

    public static native int lsb_movejob(long jobId, IntByReference opCode, int position);

    /**
     * \page lsb_switchjob lsb_switchjob
     * \brief Switches an unfinished job to another queue.
     * <p/>
     * \ref lsb_switchjob switches an unfinished job to another queue. Effectively,
     * the job is removed from its current queue and re-queued in the new queue.
     * <p/>
     * The switch operation can be performed only when the job is acceptable to
     * the new queue. If the switch operation is unsuccessful, the job will stay
     * where it is.A user can only switch his/her own unfinished jobs, but root
     * and the LSF administrator can switch any unfinished job.
     * <p/>
     * Any program using this API must be setuid to root if LSF_AUTH is not defined
     * in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_switchjob (long jobId, String queue)</b>
     *
     * @param jobId The job to be switched. If an element of a job array is to
     *              be switched, use the array form jobID[i] where jobID is the job array name,
     *              and i is the index value.
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bswitch
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param queue The new queue for the job.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * see none
     */
    public static native int lsb_switchjob(long jobId, String queue);

    /**
     * \page lsb_queuecontrol lsb_queuecontrol
     * \brief Changes the status of a queue.
     * <p/>
     * \ref lsb_queuecontrol changes the status of a queue.
     * <p/>
     * Any program using this API must be setuid to root if LSF_AUTH is not defined
     * in the lsf.conf file.
     * <p/>
     * If a queue is inactivated by its dispatch window (see lsb.queues), then it
     * cannot be re-activated by this call.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_queuecontrol (queueCtrlReq.ByReference req)</b>
     *
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param req queue control request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * queueCtrlReq
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref queue_ctrl_option
     * #see \ref lsb_queueinfo
     */
    public static native int lsb_queuecontrol(queueCtrlReq req);

    /**
     * \page lsb_userinfo lsb_userinfo
     * \brief Returns the maximum number of job slots that a user can use
     * simultaneously on any host and in the whole local LSF cluster.
     * <p/>
     * \ref lsb_userinfo gets the maximum number of job slots that a user can use
     * simultaneously on any host and in the whole local LSF cluster, as well as
     * the current number of job slots used by running and suspended jobs or
     * reserved for pending jobs. The maximum numbers of job slots are defined
     * in the LSF configuration file lsb.users (see lsb.users). The reserved
     * user name default, defined in the lsb.users configuration file, matches
     * users not listed in the lsb.users file who have no jobs started in the
     * system.
     * <p/>
     * The returned array will be overwritten by the next call.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * userInfoEnt.ByReference lsb_userinfo(String[] users, IntByReference numUsers)</b>
     *
     * @return null \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error. If lsberrno is
     *         LSBE_BAD_USER, (*users)[*numUsers] is not a user known to the LSF system.
     *         Otherwise, if.ByReference numUsers is less than its original value,* numUsers is the actual
     *         number of users found.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         busers
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.users
     * @param users An array of user names.
     * @param numUsers The number of user names.
     * To get information about all users, set.ByReference numUsers = 0;* numUsers will
     * be updated to the actual number of users when this call returns. To get
     * information on the invoker, set users = null,* numUsers = 1.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * userInfoEnt
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_hostinfo
     * #see \ref lsb_queueinfo
     */
    public static native userInfoEnt.ByReference lsb_userinfo(Pointer users, IntByReference numUsers);

    /**
     * \page lsb_hostgrpinfo lsb_hostgrpinfo
     * Returns LSF host group membership.
     * <p/>
     * \ref lsb_hostgrpinfo gets LSF host group membership.
     * <p/>
     * LSF host group is defined in the configuration file lsb.hosts.
     * <p/>
     * The storage for the array of groupInfoEnt structures will be reused by
     * the next call.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * groupInfoEnt.ByReference lsb_hostgrpinfo (String[] groups,IntByReference numGroups,
     * int options)</b>
     *
     * @param options The bitwise inclusive OR of some of the following flags:
     *                \n GRP_RECURSIVE
     *                \n Expand the group membership recursively. That is, if a member of a
     *                group is itself a group, give the names of its members recursively, rather
     *                than its name, which is the default.
     *                \n GRP_ALL
     *                \n Get membership of all groups.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                groupInfoEnt
     *                \n userShares
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref group_membership_option
     *                \n \ref group_define
     * @return null \n
     *         Function failed.
     *         <p/>
     *         <b>Errors:</b>
     *         \par
     *         On failure, returns null and sets lsberrno to indicate the error. If there
     *         are invalid groups specified, the function returns the groups up to the
     *         invalid ones and then sets lsberrno to LSBE_BAD_GROUP, which means that
     *         the specified (*groups)[*numGroups] is not a group known to the LSF system.
     *         If the first group specified is invalid, the function returns null.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts \n
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.users
     * @param groups An array of group names.
     * @param numGroups The number of group names.* numGroups will be updated
     * to the actual number of groups when this call returns.
     * #see \ref lsb_usergrpinfo
     */
    public static native groupInfoEnt.ByReference lsb_hostgrpinfo(Pointer groups, IntByReference numGroups, int options);

    /**
     * \page lsb_usergrpinfo lsb_usergrpinfo
     * \brief Returns LSF user group membership.
     * <p/>
     * \ref lsb_usergrpinfo gets LSF user group membership.
     * LSF user group is defined in the configuration file lsb.users.
     * The storage for the array of groupInfoEnt structures will be reused by
     * the next call.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * groupInfoEnt.ByReference lsb_usergrpinfo (String[] groups,
     * IntByReference numGroups, int options)</b>
     *
     * @param options The bitwise inclusive OR of some of flags in \ref group_membership_option
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                groupInfoEnt
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                \ref group_membership_option
     *                \n \ref group_define
     * @return null \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, returns null and sets lsberrno to indicate the error. If there
     *         are invalid groups specified, the function returns the groups up to the
     *         invalid ones. It then set lsberrno to LSBE_BAD_GROUP, that is the specified
     *         (*groups)[*numGroups] is not a group known to the LSF system. If the first
     *         group is invalid, the function returns null.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.hosts
     *         \n $LSB_CONFDIR/cluster_name/configdir/lsb.users
     * @param groups An array of group names.
     * @param numGroups The number of group names.* numGroups will be updated
     * to the actual number of groups when this call returns.
     * #see \ref lsb_hostgrpinfo
     */
    public static native groupInfoEnt.ByReference lsb_usergrpinfo(Pointer groups, IntByReference numGroups, int options);

    /**
     * \page lsb_parameterinfo lsb_parameterinfo
     * \brief Returns information about the LSF cluster.
     * <p/>
     * \ref lsb_parameterinfo gets information about the LSF cluster.
     * <p/>
     * The static storage for the parameterInfo structure is reused on the next call.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * parameterInfo.ByReference lsb_parameterinfo(String[] names,
     * IntByReference numUsers, int options)</b>
     *
     * @param options Reserved but not used; supply 0.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                \ref parameterInfo
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                none
     * @return null \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         none
     *         <p/>
     *         \b Files:
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param names Reserved but not used; supply null.
     * @param numUsers Reserved but not used; supply null.
     * see none
     */
    public static native parameterInfo.ByReference lsb_parameterinfo(Pointer names, IntByReference numUsers, int options);

    /**
     * \page lsb_modify lsb_modify
     * \brief  Modifies a submitted job's parameters.
     * <p/>
     * lsb_modify() allows for the modification of a submitted job's parameters.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * long lsb_modify (submit.ByReference jobsubReq,
     * submitReply.ByReference jobSubReply,
     * long jobId)</b>
     *
     * @param jobId The job to be modified. If an element of a job array is to
     *              be modified, use the array form jobID[i] where jobID is the job array name,
     *              and i is the index value.
     *              <p/>
     *              <b>Data Structures:</b>
     *              \par
     *              \ref submit
     *              \n \ref submitReply
     *              <p/>
     *              <b>Define Statements:</b>
     *              \par
     *              none
     * @return long:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command :</b>
     *         \par
     *         bmod
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param jobSubReq Describes the requirements for job modification to the
     * batch system. A job that does not meet these requirements is not submitted
     * to the batch system and an error is returned.
     * @param jobSubReply Describes the results of the job modification to the
     * batch system.
     * #see \ref lsb_submit
     * #see \ref ls_info
     * #see \ref ls_rtask
     * #see \ref lsb_queueinfo
     */
    public static native long lsb_modify(submit jobSubReq, submitReply jobSubReply, long jobId);

    public static native FloatByReference getCpuFactor(String string1, int int1);

    /**
     * \page lsb_suspreason lsb_suspreason
     * \brief Explains why a job was suspended.
     * <p/>
     * Using the SBD, \ref lsb_suspreason explains why system-suspended and
     * user-suspended jobs were suspended.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * String lsb_suspreason(int reasons, int subreasons,
     * loadIndexLog.ByReference ld)</b>
     *
     * @param reasons    Reasons a job suspends.
     * @param subreasons If reasons is SUSP_LOAD_REASON, subreasons indicates
     *                   the load indices that are out of bounds. The integer values for the load
     *                   indices are found in lsf.h.If reasons is SUSP_RES_LIMIT, subreasons
     *                   indicates the job's requirements for resource reservation are not satisfied.
     *                   The integer values for the job's requirements for resource reservation are
     *                   found in lsbatch.h.
     *                   \n Subreasons a job suspends if reasons is SUSP_LOAD_REASON:
     *                   - \b  R15S
     *                   \n 15 second CPU run queue length
     *                   - \b  R1M
     *                   \n 1 minute CPU run queue length
     *                   - \b  R15M
     *                   \n 15 minute CPU run queue length
     *                   - \b  UT
     *                   \n 1 minute CPU utilization
     *                   - \b  PG
     *                   \n Paging rate
     *                   - \b  IO
     *                   \n Disk IO rate
     *                   - \b LS
     *                   \n Number of log in sessions
     *                   - \b IT
     *                   \n Idle time
     *                   - \b TMP
     *                   \n Available temporary space
     *                   - \b SWP
     *                   \n Available swap space
     *                   - \b MEM
     *                   \n Available memory
     *                   - \b USR1
     *                   \n USR1 is used to describe unavailable or out of bounds user defined load
     *                   information of an external dynamic load indice on execution hosts.
     *                   - \b USR2
     *                   \n USR2 is used to describe unavailable or out of bounds user defined load
     *                   information of an external dynamic load indice on execution hosts.
     * @return null \n
     *         The function failed. The reason code is bad.
     *         <p/>
     *         \b Errors:
     *         \par
     *         No error handling
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bjobs -s
     *         <p/>
     *         <b>Environment Variable:</b>
     *         \par
     *         LSB_SUSP_REASONS
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.queues \n
     *         $LSB_SHAREDIR/cluster_name/logdir/lsb.events
     * @param ld When reasons is SUSP_LOAD_REASON, ld is used to determine the
     * name of any external load indices. ld uses the most recent load index log
     * in the lsb.events file.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * loadIndexLog
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref suspending_reasons \n
     * \ref suspending_subreasons
     * #see \ref lsb_pendreason
     */
    public static native String lsb_suspreason(int reasons, int subreasons, loadIndexLog ld);

    /**
     * \page lsb_pendreason  lsb_pendreason
     * \brief Explains why a job is pending.
     * <p/>
     * Use \ref lsb_pendreason to determine why a job is pending. Each pending reason is
     * associated with one or more hosts.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * String lsb_pendreason (int numReasons, IntByReference rsTb,
     * jobInfoHead.ByReference jInfoH,
     * loadIndexLog.ByReference ld, int clusterId)</b>
     *
     * @param numReasons The number of reasons in the rsTb reason table.
     * @param clusterId  MultiCluster cluster ID. If clusterId is greater than or
     *                   equal to 0, the job is a pending remote job, and \ref lsb_pendreason checks for
     *                   host_name\@cluster_name. If host name is needed, it should be found in
     *                   jInfoH->remoteHosts. If the remote host name is not available, the constant
     *                   string remoteHost is used.
     *                   <p/>
     *                   <b>Data Structures:</b>
     *                   \par
     *                   \ref jobInfoHead
     *                   \n \ref loadIndexLog
     *                   <p/>
     *                   <b>Define Statements:</b>
     *                   \par
     *                   \ref pending_reasons
     *                   \n \ref suspending_reasons
     *                   \n \ref suspending_subreasons
     * @return null \n
     *         The function fails. The reason code is bad.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If no PEND reason is found, the function fails and lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         bjobs -p
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * @param rsTb The reason table. Each entry in the table contains one of \ref pending_reasons
     * @param jInfoH jInfoH contains job information.
     * @param ld From \ref lsb_suspreason, when reasons is SUSP_LOAD_REASON, ld is used to
     * determine the name of any external load indices. ld uses the most recent load
     * index log in the lsb.events file.
     * #see \ref lsb_geteventrec
     */
    public static native String lsb_pendreason(int numReasons, IntByReference rsTb, jobInfoHead jInfoH, loadIndexLog ld, int clusterId);

    /**
     * \page lsb_calendarinfo lsb_calendarinfo
     * \brief Gets information about calendars defined in the batch system.
     * <p/>
     * \ref lsb_calendarinfo gets information about calendars defined in the batch system.
     * <p/>
     * On success, this routine returns a pointer to an array of calendarInfoEnt
     * structures which stores the information about the returned calendars and
     * numCalendars gives number of calendars returned. On failure null is returned
     * and lsberrno is set to indicate the error.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * calendarInfoEnt.ByReference lsb_calendarinfo(String[] calendars,
     * IntByReference numCalendars, String user)</b>
     *
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param calendars calendars is a pointer to an array of calendar names.
     * @param numCalendars numCalendars gives the number of calendar names. If
     * * numCalendars is 0, then information about all calendars is returned.
     * By default, only the invokers calendars are considered.
     * @param user Setting the user parameter will cause the given users calendars
     * to be considered.Use the reserved user name all to get calendars of all users.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * calendarInfoEnt
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_calendarop
     */
    public static native calendarInfoEnt.ByReference lsb_calendarinfo(Pointer calendars, IntByReference numCalendars, String user);

    public static native int lsb_calExprOccs(String string1, int int1, int int2, String string2, PointerByReference int3);

    /**
     * \page lsb_calendarop lsb_calendarop
     * \brief Adds, modifies or deletes a calendar.
     * <p/>
     * \ref lsb_calendarop is used to add, modify or delete a calendar. The oper
     * parameter is one of CALADD, CALMOD, or CALDEL. When the operation CALADD
     * is specified, the first element of the names array is used as the name of
     * the calendar to add. The desc and calExpr parameters should point to the
     * description string and the time expression list, respectively. See bcaladd()
     * for a description of time expressions.
     * <p/>
     * CALMOD permits the modification of the
     * description or time expression list associated with an existing calendar. The
     * first name in the names array indicates the calendar to be modified. The desc
     * and calExpr parameters can be set to the updated value or to null to
     * indicate that the existing value should be maintained.
     * <p/>
     * If the operation is
     * CALDEL then the names parameter points to an array of calendar names to be
     * deleted. numNames gives the number of names in the array. options is
     * reserved for the future use.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * IntByReference lsb_calendarop(int oper, int numNames, String[] names, byte
     * * desc, String calExpr, int options, String[] badStr)</b>
     *
     * @param oper     One of CALADD, CALMOD, or CALDEL. Depending on which one is
     *                 chosen, adds, modifies, or deletes a calendar.Defined in \ref calendar_command.
     * @param numNames The number of names in the array.
     * @param options  Currently unused.
     * @return int:-1
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error. If error
     *         is related to bad calendar name or time expression, the routine returns
     *         the name or expression in badStr.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param names Depending on oper, it defines the name of the calendar is going
     * to be added, modified or deleted.
     * @param desc The calendar's description list.
     * @param calExpr A calendar expression.
     * @param badStr Return from mbatchd indicating bad name or event time of calendar.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref calendar_command
     * #see \ref lsb_calendarinfo
     */
    public static native int lsb_calendarop(int oper, int numNames, Pointer names, String desc, String calExpr, int options, String badStr);

    /**
     * \page lsb_puteventrec lsb_puteventrec
     * \brief Puts information of an eventRec structure pointed to by logPtr
     * into a log file.
     * <p/>
     * \ref lsb_puteventrec puts information of an eventRec structure pointed to by
     * logPtr into a log file. log_fp is a pointer pointing to the log file name
     * that could be either event a log file or job log file.
     * <p/>
     * See \ref lsb_geteventrec for detailed information about parameters.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_puteventrec(Pointer log_fp, eventRec.ByReference logPtr)</b>
     *
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_SHAREDIR/cluster_name/logdir/lsb.events
     * @param logPtr The eventRec structure pointed to by logPtr into a log file.
     * @param log_fp A pointer pointing to the log file name that could be either
     * event a log file or job log file.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * eventRec
     * \n eventLog
     * \n xFile
     * \n jobAttrSetLog
     * \n logSwitchLog
     * \n dataLoggingLog
     * \n jgrpNewLog
     * \n jgrpCtrlLog
     * \n jgrpStatusLog
     * \n jobNewLog
     * \n jobModLog
     * \n jobStartLog
     * \n jobStartAcceptLog
     * \n jobExecuteLog
     * \n jobStatusLog
     * \n sbdJobStatusLog
     * \n sbdUnreportedStatusLog
     * \n jobSwitchLog
     * \n jobMoveLog
     * \n chkpntLog
     * \n jobRequeueLog
     * \n jobCleanLog
     * \n jobExceptionLog
     * \n sigactLog
     * \n migLog
     * \n signalLog
     * \n queueCtrlLog
     * \n hostCtrlLog
     * \n hgCtrlLog
     * \n mbdStartLog
     * \n mbdDieLog
     * \n unfulfillLog
     * \n jobFinishLog
     * \n loadIndexLog
     * \n calendarLog
     * \n jobForwardLog
     * \n jobAcceptLog
     * \n statusAckLog
     * \n jobMsgLog
     * \n jobMsgAckLog
     * \n jobOccupyReqLog
     * \n jobVacatedLog
     * \n jobForceRequestLog
     * \n jobChunkLog
     * \n jobExternalMsgLog
     * \n rsvRes
     * \n rsvFinishLog
     * \n cpuProfileLog
     * \n jobRunRusageLog
     * \n slaLog
     * \n perfmonLogInfo
     * \n perfmonLog
     * \n taskFinishLog
     * \n eventEOSLog
     * \n jobResizeNotifyStartLog
     * \n jobResizeNotifyAcceptLog
     * \n jobResizeNotifyDoneLog
     * \n jobResizeReleaseLog
     * \n jobResizeCancelLog
     * \n jobResizeLog
     * \n jRusage
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref event_types
     * \n \ref defs_lsb_XF_OP
     * \n \ref jobgroup_controltypes
     * \n \ref signal_action
     * #see \ref lsb_geteventrec
     */
    public static native int lsb_puteventrec(Pointer logPtr, eventRec log_fp);

    public static native int lsb_puteventrecRaw(Pointer pointer1, eventRec eventRec1, String string1);

    /**
     * \page lsb_geteventrec lsb_geteventrec
     * \brief Get an event record from a log file
     * <p/>
     * \ref lsb_geteventrec returns an eventRec from a log file.
     * <p/>
     * The storage for the eventRec structure returned by \ref lsb_geteventrec will be
     * reused by the next call.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * eventRec.ByReference lsb_geteventrec(Pointer  log_fp,IntByReference  lineNum)</b>
     *
     * @return null \n
     *         Function failed.If there are no more records, returns null and sets
     *         lsberrno to LSBE_EOF.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_SHAREDIR/cluster_name/logdir/lsb.acct
     *         \n $LSB_SHAREDIR/cluster_name/logdir/lsb.events
     *         \n $LSB_SHAREDIR/cluster_name/logdir/lsb.rsv.ids
     *         \n $LSB_SHAREDIR/cluster_name/logdir/lsb.rsv.state
     * @param log_fp Either an event log file or a job log file.
     * @param lineNum The number of the event record.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * eventRec
     * \n eventLog
     * \n xFile
     * \n jobAttrSetLog
     * \n logSwitchLog
     * \n dataLoggingLog
     * \n jgrpNewLog
     * \n jgrpCtrlLog
     * \n jgrpStatusLog
     * \n jobNewLog
     * \n jobModLog
     * \n jobStartLog
     * \n jobStartAcceptLog
     * \n jobExecuteLog
     * \n jobStatusLog
     * \n sbdJobStatusLog
     * \n sbdUnreportedStatusLog
     * \n jobSwitchLog
     * \n jobMoveLog
     * \n chkpntLog
     * \n jobRequeueLog
     * \n jobCleanLog
     * \n jobExceptionLog
     * \n sigactLog
     * \n migLog
     * \n signalLog
     * \n queueCtrlLog
     * \n hostCtrlLog
     * \n hgCtrlLog
     * \n mbdStartLog
     * \n mbdDieLog
     * \n unfulfillLog
     * \n jobFinishLog
     * \n loadIndexLog
     * \n calendarLog
     * \n jobForwardLog
     * \n jobAcceptLog
     * \n statusAckLog
     * \n jobMsgLog
     * \n jobMsgAckLog
     * \n jobOccupyReqLog
     * \n jobVacatedLog
     * \n jobForceRequestLog
     * \n jobChunkLog
     * \n jobExternalMsgLog
     * \n rsvRes
     * \n rsvFinishLog
     * \n cpuProfileLog
     * \n jobRunRusageLog
     * \n slaLog
     * \n perfmonLogInfo
     * \n perfmonLog
     * \n taskFinishLog
     * \n eventEOSLog
     * \n jobResizeNotifyStartLog
     * \n jobResizeNotifyAcceptLog
     * \n jobResizeNotifyDoneLog
     * \n jobResizeReleaseLog
     * \n jobResizeCancelLog
     * \n jobResizeLog
     * \n jRusage
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref event_types
     * \n \ref defs_lsb_XF_OP
     * \n \ref jobgroup_controltypes
     * \n \ref signal_action
     * #see \ref lsb_hostcontrol
     * #see \ref lsb_movejob
     * #see \ref lsb_pendreason
     * #see \ref lsb_puteventrec
     * #see \ref lsb_queuecontrol
     * #see \ref lsb_readjobinfo
     * #see \ref lsb_submit
     * #see \ref lsb_suspreason
     */
    public static native eventRec.ByReference lsb_geteventrec(Pointer log_fp, IntByReference lineNum);

    public static native eventRec.ByReference lsb_geteventrec_decrypt(Pointer pointer1, IntByReference int1);

    public static native eventRec.ByReference lsb_geteventrecord(Pointer pointer1, IntByReference int1);

    public static native eventRec.ByReference lsb_geteventrecordEx(Pointer pointer1, IntByReference int1, Pointer stringArray1);

    public static native eventRec.ByReference lsb_getnewjob_from_string(String string1);

    public static native eventInfoEnt.ByReference lsb_eventinfo(Pointer stringArray1, IntByReference int1, String string1);

    /**
     * \page lsb_sharedresourceinfo lsb_sharedresourceinfo
     * \brief Returns the requested shared resource information in dynamic values.
     * <p/>
     * \ref lsb_sharedresourceinfo returns the requested shared resource information in
     * dynamic values. The result of this call is a chained data structure as
     * defined in <lsf/lsbatch.h>, which contains requested information.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * LSB_SHARED_RESOURCE_INFO_T.ByReference lsb_sharedresourceinfo(
     * String[] resources,
     * IntByReference numResources,
     * String hostName, int options)</b>
     *
     * @param options options is reserved for future use. Currently, it should be set to 0.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                lsbSharedResourceInfo
     *                \n lsbSharedResourceInstance
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                none
     * @return null \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSF_CONFDIR/lsf.shared
     *         \n $LSF_CONFDIR/lsf.cluster.cluster_name
     * @param resources resources is an null terminated string array storing
     * requesting resource names.Setting resources to point to null returns all
     * shared resources.
     * @param numResources numResources is an input/output parameter. On input
     * it indicates how many resources are requested. Value 0 means requesting
     * all shared resources. On return it contains qualified resource number.
     * @param hostName hostName is a string containing a host name. Only shared resource
     * available on the specified host will be returned. If hostName is a null,
     * shared resource available on all hosts will be returned.
     * #see \ref ls_sharedresourceinfo
     */
    public static native Pointer lsb_sharedresourceinfo(Pointer resources, IntByReference numResources, String hostName, int options);

    /**
     * \page lsb_geteventrecbyline lsb_geteventrecbyline
     * Parse an event line and put the result in an event record structure.
     * The \ref lsb_geteventrecbyline function parses an event line and puts the result
     * in an event record structure.
     * <p/>
     * If the line to be parsed is a comment line, \ref lsb_geteventrecbyline sets errno to
     * bad event format and logs an error.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_geteventrecbyline(String line, eventRec.ByReference logRec)</b>
     *
     * @return int:-1
     *         \n Function failed and lserrno was set.
     *         <p/>
     *         \b Errors:
     *         \par
     *         none
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param line
     * Buffer containing a line of event text string
     * @param logRec
     * Pointer to an eventRec structure
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * eventRec
     * \n eventLog
     * \n xFile
     * \n jobAttrSetLog
     * \n logSwitchLog
     * \n dataLoggingLog
     * \n jgrpNewLog
     * \n jgrpCtrlLog
     * \n jgrpStatusLog
     * \n jobNewLog
     * \n jobModLog
     * \n jobStartLog
     * \n jobStartAcceptLog
     * \n jobExecuteLog
     * \n jobStatusLog
     * \n sbdJobStatusLog
     * \n sbdUnreportedStatusLog
     * \n jobSwitchLog
     * \n jobMoveLog
     * \n chkpntLog
     * \n jobRequeueLog
     * \n jobCleanLog
     * \n jobExceptionLog
     * \n sigactLog
     * \n migLog
     * \n signalLog
     * \n queueCtrlLog
     * \n hostCtrlLog
     * \n hgCtrlLog
     * \n mbdStartLog
     * \n mbdDieLog
     * \n unfulfillLog
     * \n jobFinishLog
     * \n loadIndexLog
     * \n calendarLog
     * \n jobForwardLog
     * \n jobAcceptLog
     * \n statusAckLog
     * \n jobMsgLog
     * \n jobMsgAckLog
     * \n jobOccupyReqLog
     * \n jobVacatedLog
     * \n jobForceRequestLog
     * \n jobChunkLog
     * \n jobExternalMsgLog
     * \n rsvRes
     * \n rsvFinishLog
     * \n cpuProfileLog
     * \n jobRunRusageLog
     * \n slaLog
     * \n perfmonLogInfo
     * \n perfmonLog
     * \n taskFinishLog
     * \n eventEOSLog
     * \n jobResizeNotifyStartLog
     * \n jobResizeNotifyAcceptLog
     * \n jobResizeNotifyDoneLog
     * \n jobResizeReleaseLog
     * \n jobResizeCancelLog
     * \n jobResizeLog
     * \n jRusage
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * <p/>
     * <b>Pre-Conditions:</b>
     * \par
     * The event record structure must have been initialized outside the
     * \ref lsb_geteventrecbyline function.
     * see none
     */

    public static native int lsb_geteventrecbyline(String line, eventRec logRec);
/* Retain lsb_connect for now */

    public static int lsb_connect(int a) {
        return lsb_rcvconnect();
    }

    public static native int lsb_rcvconnect();

    public static native int lsb_sndmsg(lsbMsgHdr lsbMsgHdr1, String string1, int int1);

    public static native int lsb_rcvmsg(lsbMsgHdr lsbMsgHdr1, Pointer stringArray1, int int1);

    /**
     * \page  lsb_runjob lsb_runjob
     * Starts a batch job immediately on a set of specified host().
     * \ref lsb_runjob starts a batch job immediately on a set of specified host().
     * The job must have been submitted and is in PEND or FINISHED status. Only
     * the LSF administrator or the owner of the job can start the job. If the
     * options is set to RUNJOB_OPT_NOSTOP, then the job will not be suspended by
     * the queue's RUNWINDOW,loadStop and STOP_COND and the hosts' RUNWINDOW and
     * loadStop conditions. By default, these conditions apply to the job as do
     * to other normal jobs.
     * <p/>
     * Any program using this API must be setuid to root
     * if LSF_AUTH is not defined in the lsf.conf file.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_runjob(runJobRequest.ByReference runJobRequest)</b>
     *
     * @param runJobRequest The job-starting request.
     *                      <p/>
     *                      <b>Data Structures:</b>
     *                      \par
     *                      runJobRequest
     *                      <p/>
     *                      <b>Define Statements:</b>
     *                      \par
     *                      \ref runjob_option
     * @return int:-1 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         brun
     *         <p/>
     *         \b Files:
     *         \par
     *         ${LSF_ENVDIR:-/etc}/lsf.conf
     * see none
     */
    public static native int lsb_runjob(runJobRequest runJobRequest);

/* API for job group */

    public static native int lsb_addjgrp(jgrpAdd jgrpAdd1, Pointer jgrpReply1);

    public static native int lsb_modjgrp(jgrpMod jgrpMod1, Pointer jgrpReply1);

    public static native int lsb_holdjgrp(String string1, int int1, Pointer jgrpReply1);

    public static native int lsb_reljgrp(String string1, int int1, Pointer jgrpReply1);

    public static native int lsb_deljgrp(String string1, int int1, Pointer jgrpReply1);

    public static native int lsb_deljgrp_ext(jgrpCtrl jgrpCtrl1, Pointer jgrpReply1);

    public static native jgrp.ByReference lsb_listjgrp(IntByReference int1);

    public static native serviceClass.ByReference lsb_serviceClassInfo(IntByReference int1);

/* API for Application Encapsulation */

    public static native appInfoEnt.ByReference lsb_appInfo(IntByReference int1);

    public static native void lsb_freeAppInfoEnts(int int1, appInfoEnt appInfoEnt1);

/* routine to convert the job id to string */

    public static native String lsb_jobid2str(long long1);

    public static native String lsb_jobid2str_r(long long1, byte[] byte1);

    public static native String lsb_jobidinstr(long long1);
/* routine to compose and decompose 64bit jobId */

    public static native void jobId32To64(LongByReference long1, int int1, int int2);

    public static native void jobId64To32(long long1, IntByReference int1, IntByReference int2);
/* API for job attribute operations */

    public static native int lsb_setjobattr(int int1, jobAttrInfoEnt jobAttrInfoEnt1);

/* API for remote task execution */

    public static native long lsb_rexecv(int int1, Pointer stringArray1, Pointer stringArray2, IntByReference int2, int int3);


    public static interface lsb_catchCallback extends Callback {
        int invoke(Pointer pointer);
    }

    public static native int lsb_catch(String string1, lsb_catchCallback callback);

    public static native void lsb_throw(String string1, Pointer pointer1);

/* API for job external message */

    /**
     *  \page lsb_postjobmsg lsb_postjobmsg
     *  \brief Sends messages and data posted to a job.
     *
     *  Use \ref lsb_postjobmsg to post a message and data to a job, open a TCP
     *  connection, and transfer attached message and data from the mbatchd. Use
     *  \ref lsb_readjobmsg to display messages and copy data files posted by
     *  \ref lsb_postjobmsg.
     *
     *  While you can post multiple messages and attached data files to a job,
     *  you must call \ref lsb_postjobmsg for each message and attached data file
     *  you want to post. By default, \ref lsb_postjobmsg posts a message to position
     *  0 of the message index (msgId) (see PARAMETERS) of the specified job.
     *  To post additional messages to a job, call \ref lsb_postjobmsg and increment
     *  the message index.
     *
     *  \ref lsb_readjobmsg reads posted job messages by their
     *  position in the message index.
     *
     *  If a data file is attached to a message and the flag EXT_ATTA_POST is set,
     *  use the JOB_ATTA_DIR parameter in lsb.params(5) to specify the directory
     *  where attachment data fies are saved. The directory must have at least 1MB
     *  of free space.The mbatchd checks for available space in the job attachment
     *  directory before transferring the file.
     *
     *  Use the MAX_JOB_ATTA_SIZE parameter in lsb.params(5) to set a maximum size
     *  for job message attachments.
     *
     *  Users can only send messages and data from their own jobs. Root and LSF
     *  administrators can also send messages of jobs submtted by other users, but
     *  they cannot attach data files to jobs owned by other users.
     *
     *  You can post messages and data to a job until it is cleaned from the system.
     *  You cannot send messages and data to finished or exited jobs.
     *
     *  <b>\#include <lsf/lsbatch.h> \n
     *     \#include <time.h>
     *
     *  int lsb_postjobmsg(jobExternalMsgReq.ByReference jobExternalMsg,
     *                    String filename)</b>
     *
     *  @param jobExternalMsg This structure contains the information required to
     *  define an external message of a job.
     *  @param filename Name of attached data file. If no file is attached, use null.
     *
     *  <b>Data Structures:</b>
     *  \par
     *  \ref jobExternalMsgReq
     *
     *  <b>Define Statements:</b>
     *  \par
     *  \ref external_msg_post
     *
     *  @return int:value \n
     *  The successful function returns a socket number.
     * return int:0 \n
     *  The EXT_ATTA_POST bit of options is not set or there is no attached data.
     *  return int:-1 \n
     *  The function failed.
     *
     *  \b Errors:
     *  \par
     *  If the function fails, lserrno is set to indicate the error.
     *
     *  <b>Equivalent line command:</b>
     *  \par
     *  bpost
     *
     *  \b Files:
     *  \par
     *  $LSB_CONFDIR/cluster_name/configdir/lsb.params
     *  \n $JOB_ATTA_DIR
     *  \n $LSB_SHAREDIR/info
     *
     * #see \ref lsb_readjobmsg
     *
     */

    public static native int lsb_postjobmsg(jobExternalMsgReq jobExternalMsg, String filename);
    /**
     *  \page lsb_readjobmsg lsb_readjobmsg
     *  \brief Reads messages and data posted to a job.
     *
     *  Use \ref lsb_readjobmsg to open a TCP connection, receive attached messages and
     *  data from the mbatchd, and display the messages posted by \ref lsb_postjobmsg.
     *
     *  By default, \ref lsb_readjobmsg displays the message "no description" or the
     *  message at index position 0 of the specified job. To read other messages,
     *  choose another index position. The index is populated by \ref lsb_postjobmsg.
     *
     *  If a data file is attached to a message and the flag EXT_ATTA_READ is set,
     *  \ref lsb_readjobmsg gets the message and copies its data file to the default
     *  directory JOB_ATTA_DIR, overwriting the specified file if it already exists.
     *  If there is no file attached, the system reports an error.
     *
     *  Users can only read messages and data from their own jobs. Root and LSF
     *  administrators can also read messages of jobs submtted by other users,
     *  but they cannot read data files attached to jobs owned by other users.
     *
     *  You can read messages and data from a job until it is cleaned from the
     *  system. You cannot read messages and data from done or exited jobs.
     *
     *  <b>\#include <lsf/lsbatch.h> \n
     *  \#include <time.h> \n
     *  int lsb_readjobmsg(jobExternalMsgReq.ByReference jobExternalMsg,
     *          jobExternalMsgReply.ByReference jobExternalMsgReply)</b>
     *
     *  @param jobExternalMsg the information required to define an external
     *  message of a job.
     *  @param jobExternalMsgReply the information required to define an
     *  external message reply.
     *
     *  <b>Data Structures:</b>
     *  \par
     *  jobExternalMsgReq
     *  \n jobExternalMsgReply
     *
     *  <b>Define Statements:</b>
     *  \par
     *  \ref external_msg_processing
     *  \n \ref ext_data_status
     *
     *  @return int:value \n
     *  The successful function returns a socket number.
     *  return int:0 \n
     *  The EXT_ATTA_READ bit of options is not set or there is no
     *  attached data.
     *  return int:-1 \n
     *  The function failed.
     *
     *  \b Errors:
     *  \par
     *  If the function fails, lserrno is set to indicate the error.
     *
     *  <b>Equivalent line commands:</b>
     *  \par
     *  bread
     *
     *  <b>Files:</b>
     *  \par
     *  $LSB_CONFDIR/cluster_name/configdir/lsb.params
     *  \n $JOB_ATTA_DIR
     *  \n $LSB_SHAREDIR/info
     * #see \ref lsb_postjobmsg
     */

    public static native int lsb_readjobmsg(jobExternalMsgReq jobExternalMsg, jobExternalMsgReply jobExternalMsgReply);

/* API for symphony job information update in bulk mode */

    public static native int lsb_bulkJobInfoUpdate(symJobStatusUpdateReqArray symJobStatusUpdateReqArray1, symJobStatusUpdateReplyArray symJobStatusUpdateReplyArray1);

/* API for advance reservation */

    /**
     * \page lsb_addreservation lsb_addreservation
     * \brief Makes an advance reservation.
     * <p/>
     * Use \ref lsb_addreservation to send a reservation request to mbatchd. If
     * mbatchd grants the reservation, it issues the reservation ID. If mbatchd
     * rejects the request, it issues null as the reservation ID.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_addreservation (addRsvRequest.ByReference request, String rsvId)</b>
     *
     * @return int:-1 \n
     *         The reservation failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         brsvadd
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param request The reservation request
     * @param rsvId Reservation ID returned from mbatchd. If the reservation
     * fails, this is null. The
     * memory for rsvid is allocated by the caller.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * addRsvRequest
     * \n _rsvExecCmd_t
     * \n _rsvExecEvent_t
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref reservation_option
     * #see \ref lsb_removereservation
     * #see \ref lsb_modreservation
     * #see \ref lsb_reservationinfo
     */
    public static native int lsb_addreservation(addRsvRequest request, String rsvId);

    /**
     * \page lsb_removereservation lsb_removereservation
     * \brief Removes a reservation.
     * <p/>
     * Use \ref lsb_removereservation to remove a reservation. mbatchd removes the
     * reservation with the specified reservation ID.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_removereservation(String rsvId)</b>
     *
     * @return int:-1 \n
     *         The reservation removal failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         brsvdel
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param rsvId Reservation ID of the reservation that you wish to remove.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_addreservation
     * #see \ref lsb_modreservation
     * #see \ref lsb_reservationinfo
     */
    public static native int lsb_removereservation(String rsvId);

    /**
     * \page lsb_reservationinfo lsb_reservationinfo
     * \brief Retrieve reservation information to display active advance reservations.
     * <p/>
     * Use \ref lsb_reservationinfo to retrieve reservation information from mbatchd.
     * This function allocates memory that the caller should free.
     * <p/>
     * If the \ref lsb_reservationinfo function succeeds, it returns the reservation
     * records pertaining to a particular reservation ID (rsvId) as an array of
     * rsvInfoEnt structs.
     * <p/>
     * If rsvId is null, all reservation information will be returned. If a
     * particular rsvId  is specified:
     * \li If found, the reservation record pertaining to a particular rsvId is
     * returned
     * \li If not found, the number of reservation records is set to zero and
     * the lsberrno  is set appropiately
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * rsvInfoEnt.ByReference lsb_reservationinfo(String rsvId, IntByReference numEnts,
     * int options)</b>
     *
     * @param options The parameter options is currently ignored.
     *                <p/>
     *                <b>Data Structures:</b>
     *                \par
     *                rsvInfoEnt
     *                \n hostRsvInfoEnt
     *                <p/>
     *                <b>Define Statements:</b>
     *                \par
     *                none
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         brsvs
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param rsvId Reservation ID of the requested reservation.
     * @param numEnts Number of reservation entries that mbatchd returns.
     * #see \ref lsb_addreservation
     * #see \ref lsb_modreservation
     * #see \ref lsb_removereservation
     */

    public static native rsvInfoEnt.ByReference lsb_reservationinfo(String rsvId, IntByReference numEnts, int options);

    public static native int lsb_freeRsvExecCmd(Pointer _rsvExecCmd_tArray1);

    public static native _rsvExecCmd_t.ByReference lsb_dupRsvExecCmd(_rsvExecCmd_t _rsvExecCmd_t1);

    public static native int lsb_parseRsvExecOption(String string1, Pointer _rsvExecCmd_tArray1);

    /**
     * \page lsb_modreservation lsb_modreservation
     * \brief Modifies an advance reservation.
     * <p/>
     * Use \ref lsb_modreservation to modify an advance reservation. mbatchd receives
     * the modification request and modifies the reservation with the specified
     * reservation ID.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_modreservation(modRsvRequest.ByReference request)</b>
     *
     * @return int:-1 \n
     *         The reservation modification failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         brsvmod
     *         <p/>
     *         \b Files:
     *         \par
     *         none
     * @param request modify reservation request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * modRsvRequest
     * \n addRsvRequest
     * \n _rsvExecCmd_t
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_addreservation
     * #see \ref lsb_removereservation
     * #see \ref lsb_reservationinfo
     */

    public static native int lsb_modreservation(modRsvRequest request);

/* routines for sorted integer list */
    /*
    sortIntList.ByReference  initSortIntList(int);
    int insertSortIntList(sortIntList.ByReference , int);
    sortIntList.ByReference  getNextSortIntList(sortIntList.ByReference , sortIntList.ByReference , IntByReference );
    void freeSortIntList(sortIntList.ByReference );
    int getMinSortIntList(sortIntList.ByReference , IntByReference );
    int getMaxSortIntList(sortIntList.ByReference , IntByReference );
    int getTotalSortIntList(sortIntList.ByReference );

    int updateJobIdIndexFile (String string1, String string1, int);
    */

/* Structures and routine for obtaining subset of info about jobs
*  This is being used by Maui integration.
 */

    public static class jobExtschInfoReq extends Structure {
        public static class ByReference extends jobExtschInfoReq implements Structure.ByReference {}
        public static class ByValue extends jobExtschInfoReq implements Structure.ByValue {}
        public jobExtschInfoReq() {}
        public jobExtschInfoReq(Pointer p) { super(p); read(); }

        public int qCnt;
        public Pointer queues;
    }



    public static class jobExtschInfo extends Structure {
        public static class ByReference extends jobExtschInfo implements Structure.ByReference {}
        public static class ByValue extends jobExtschInfo implements Structure.ByValue {}
        public jobExtschInfo() {}
        public jobExtschInfo(Pointer p) { super(p); read(); }

        public long jobId;
        public int status;
        public NativeLong jRusageUpdateTime;
        public LibLsf.jRusage runRusage;
    }



    public static class jobExtschInfoReply extends Structure {
        public static class ByReference extends jobExtschInfoReply implements Structure.ByReference {}
        public static class ByValue extends jobExtschInfoReply implements Structure.ByValue {}
        public jobExtschInfoReply() {}
        public jobExtschInfoReply(Pointer p) { super(p); read(); }

        public int jobCnt;
        public PointerByReference jobs;
    }



    public static native int getjobinfo4queues(jobExtschInfoReq jobExtschInfoReq1, jobExtschInfoReply jobExtschInfoReply1);

    public static native void free_jobExtschInfoReply(jobExtschInfoReply jobExtschInfoReply1);

    public static native void free_jobExtschInfoReq(jobExtschInfoReq jobExtschInfoReq1);

/* For RFC 725 */

    public static native String longer_strcpy(String dest, String src);

/* Structures and API for job diagnostics.  These are applicable only if
*  CONDENSE_PENDING_REASONS is enabled in lsb.params.
 */

    public static class diagnoseJobReq extends Structure {
        public static class ByReference extends diagnoseJobReq implements Structure.ByReference {}
        public static class ByValue extends diagnoseJobReq implements Structure.ByValue {}
        public diagnoseJobReq() {}
        public diagnoseJobReq(Pointer p) { super(p); read(); }

        public int jobCnt;
        public LongByReference jobId;
    }



    public static native int lsb_diagnosejob(diagnoseJobReq diagnoseJobReq1);

    public static final int SIM_STATUS_RUN = 0x01;
    public static final int SIM_STATUS_SUSPEND = 0x02;

/* simulator status reply
 */

    public static class simStatusReply extends Structure {
        public static class ByReference extends simStatusReply implements Structure.ByReference {}
        public static class ByValue extends simStatusReply implements Structure.ByValue {}
        public simStatusReply() {}
        public simStatusReply(Pointer p) { super(p); read(); }

        public int simStatus;
        public NativeLong curTime;
    }



    public static native simStatusReply.ByReference lsb_simstatus();

    public static native void free_simStatusReply(simStatusReply simStatusReply1);

/* batch command options flag for lease */
    public static final int LSB_HOST_OPTION_EXPORT = 0x1;
/* bhosts -x option */
    public static final int LSB_HOST_OPTION_EXCEPT = 0x2;
/* retrieve hosts that belong to batch partition */
    public static final int LSB_HOST_OPTION_BATCH = 0x4;


/* Display condensed host output */
    public static final int LSB_HOST_OPTION_CONDENSED = 0x08;

/* error codes, structures and routines for syntax check of RMS external scheduler options */

/*  non-rms option shown up in RMS[] */
    public static final int RMS_NON_RMS_OPTIONS_ERR = (-1);

/*  nodes and ptile co-exist */
    public static final int RMS_NODE_PTILE_ERR = (-2);

/*  rails and railmask co-exist */
    public static final int RMS_RAIL_RAILMASK_ERR = (-3);

/*  nodes is out of range 1..LSB_RMS_MAXNUMNODES */
    public static final int RMS_NODES_OUT_BOUND_ERR = (-4);

/*  ptile is out of range 1..LSB_RMS_MAXPTILE */
    public static final int RMS_PTILE_OUT_BOUND_ERR = (-5);

/*  rails is out of range 1..LSB_RMS_MAXNUMRAILS */
    public static final int RMS_RAIL_OUT_BOUND_ERR = (-6);

/*  railmask syntax error */
    public static final int RMS_RAILMASK_OUT_BOUND_ERR = (-7);

/*  nodes syntax error */
    public static final int RMS_NODES_SYNTAX_ERR = (-8);

/*  ptile syntax error */
    public static final int RMS_PTILE_SYNTAX_ERR = (-9);

/*  rails syntax error */
    public static final int RMS_RAIL_SYNTAX_ERR = (-10);

/*  railmask syntax error */
    public static final int RMS_RAILMASK_SYNTAX_ERR = (-11);

/*  base syntax error */
    public static final int RMS_BASE_SYNTAX_ERR = (-12);

/*  base string too NativeLong*/
    public static final int RMS_BASE_TOO_LONG = (-13);

/*  >=1 allocation types are specified */
    public static final int RMS_TOO_MANY_ALLOCTYPE_ERR = (-14);

/*  =1 allocation types are specified */
    public static final int RMS_NO_LSF_EXTSCHED_Y_ERR = (-15);

/*  error reading env from lsf.conf inside syntax check */
    public static final int RMS_READ_ENV_ERR = (-20);

/*  memory allocation problems inside syntax check function */
    public static final int RMS_MEM_ALLOC_ERR = (-21);

/*  [] mis-matched in RMS[] */
    public static final int RMS_BRACKETS_MISMATCH_ERR = (-22);

    public static interface rmsAllocType_t {
          public static final int RMS_ALLOC_TYPE_UNKNOWN = 0;
          public static final int RMS_ALLOC_TYPE_SLOAD = 1;
          public static final int RMS_ALLOC_TYPE_SNODE = 2;
          public static final int RMS_ALLOC_TYPE_MCONT = 3;
    }



    public static interface rmsTopology_t {
          public static final int RMS_TOPOLOGY_UNKNOWN = 0;
          public static final int RMS_TOPOLOGY_PTILE = 1;
          public static final int RMS_TOPOLOGY_NODES = 2;
    }



    public static interface rmsFlags_t {
          public static final int RMS_FLAGS_UNKNOWN = 0;
          public static final int RMS_FLAGS_RAILS = 1;
          public static final int RMS_FLAGS_RAILMASK = 2;
    }



    public static class rmsextschedoption extends Structure {
        public static class ByReference extends rmsextschedoption implements Structure.ByReference {}
        public static class ByValue extends rmsextschedoption implements Structure.ByValue {}
        public rmsextschedoption() {}
        public rmsextschedoption(Pointer p) { super(p); read(); }

        public /*rmsAllocType_t*/ int alloc_type;
        public /*rmsTopology_t*/ int topology;
        public int topology_value;
        public int set_base;
        public byte[] base = new byte[LibLsf.MAXHOSTNAMELEN];
        public /*rmsFlags_t*/ int flags;
        public int flags_value;
    }



    public static native int parseRmsOptions(String string1, rmsextschedoption rmsextschedoption1, LibLsf.config_param config_param1);

/* Stream interface.
*  By default the stream lsb.stream is located in a subdirectory
*  stream of the cluster working directory i.e.:
*  work/<clustername>/logdir/stream and the size of
*  lsb.stream is 1024MB
 */
    public static final int MBD_DEF_STREAM_SIZE = (1024 * 1024 * 1024);

/* default maximum number of backup stream.utc file */
    public static final int DEF_MAX_STREAM_FILE_NUMBER = 10;

    /**
     * \brief  Stream interface.
     */
    public static class lsbStream extends Structure {
        public static class ByReference extends lsbStream implements Structure.ByReference {}
        public static class ByValue extends lsbStream implements Structure.ByValue {}
        public lsbStream() {}
        public lsbStream(Pointer p) { super(p); read(); }

        public static interface trsFunc extends Callback {
            int invoke(String string1);
        }

        /**
         * < Pointer to full path to the stream file
         */
        public String streamFile;

        /**
         * < Max size of the stream file
         */
        public int maxStreamSize;

        /**
         * < Max number of backup stream files
         */
        public int maxStreamFileNum;

        /**
         * < Set to 1 to enable trace of the stream
         */
        public int trace;

        /**
         * < Pointer to a function that the library invokes, passing a trace buffer.
         */
        public trsFunc trs;
    }



     /**//*
     * \page lsb_openstream  lsb_openstream
     * \brief Open and create an lsb_stream file.
     * <p/>
     * \ref lsb_openstream opens the streamFile .
     * <p/>
     * This API function is inside liblsbstream.so.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_openstream(lsbStream.ByReference params)</b>
     *
     * @return int:-1 or null \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param params Parameters.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * \ref lsbStream
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_closestream
     * #see \ref lsb_readstreamline
     * #see \ref lsb_writestream
     * #see \ref lsb_readstream
     * #see \ref lsb_streamversion
     */
    // NOTE: Not in libbat
    //public static native int lsb_openstream(lsbStream params);

     /**//*
     * \page lsb_closestream lsb_closestream
     * \brief Close an lsb_stream file.
     * <p/>
     * \ref lsb_closestream closes the streamFile.
     * <p/>
     * This API function is inside liblsbstream.so.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_closestream(String config)</b>
     *
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param config Pointer to the handle of the stream file.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_openstream
     * #see \ref lsb_readstreamline
     * #see \ref lsb_writestream
     * #see \ref lsb_readstream
     * #see \ref lsb_streamversion
     */
    // NOTE: Not in libbat
    //public static native int lsb_closestream(String config);

     /**//*
     * \page lsb_streamversion lsb_streamversion
     * \brief Version of the current event type supported by mbatchd.
     * <p/>
     * \ref lsb_streamversion returns the event version number of mbatchd, which is the
     * version of the events to be written to the stream file. This API function
     * is inside liblsbstream.so.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * String  lsb_streamversion()</b>
     *
     * param void \n
     *             <p/>
     *             <b>Data Structures:</b>
     *             \par
     *             none
     *             <p/>
     *             <b>Define Statements:</b>
     *             \par
     *             none
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * #see \ref lsb_closestream
     * #see \ref lsb_geteventrec
     * #see \ref lsb_openstream
     * #see \ref lsb_puteventrec
     * #see \ref lsb_readstreamline
     * #see \ref lsb_writestream
     * #see \ref lsb_readstream
     */
    // NOTE: Not in libbat
    //public static native String lsb_streamversion();

     /**//*
     * \page lsb_writestream lsb_writestream
     * \brief Writes a current version eventRec structure into the lsb_stream file.
     * <p/>
     * \ref lsb_writestream writes an eventrRec to the open streamFile.
     * This API function is inside liblsbstream.so.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_writestream(eventRec.ByReference logPtr)</b>
     *
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param logPtr Pointer to the eventRec structure.
     * \n see \ref lsb_geteventrec for details on the eventRec structure.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * eventRec
     * \n eventLog
     * \n xFile
     * \n jobAttrSetLog
     * \n logSwitchLog
     * \n dataLoggingLog
     * \n jgrpNewLog
     * \n jgrpCtrlLog
     * \n jgrpStatusLog
     * \n jobNewLog
     * \n jobModLog
     * \n jobStartLog
     * \n jobStartAcceptLog
     * \n jobExecuteLog
     * \n jobStatusLog
     * \n sbdJobStatusLog
     * \n sbdUnreportedStatusLog
     * \n jobSwitchLog
     * \n jobMoveLog
     * \n chkpntLog
     * \n jobRequeueLog
     * \n jobCleanLog
     * \n jobExceptionLog
     * \n sigactLog
     * \n migLog
     * \n signalLog
     * \n queueCtrlLog
     * \n hostCtrlLog
     * \n hgCtrlLog
     * \n mbdStartLog
     * \n mbdDieLog
     * \n unfulfillLog
     * \n jobFinishLog
     * \n loadIndexLog
     * \n calendarLog
     * \n jobForwardLog
     * \n jobAcceptLog
     * \n statusAckLog
     * \n jobMsgLog
     * \n jobMsgAckLog
     * \n jobOccupyReqLog
     * \n jobVacatedLog
     * \n jobForceRequestLog
     * \n jobChunkLog
     * \n jobExternalMsgLog
     * \n rsvRes
     * \n rsvFinishLog
     * \n cpuProfileLog
     * \n jobRunRusageLog
     * \n slaLog
     * \n perfmonLogInfo
     * \n perfmonLog
     * \n taskFinishLog
     * \n eventEOSLog
     * \n jobResizeNotifyStartLog
     * \n jobResizeNotifyAcceptLog
     * \n jobResizeNotifyDoneLog
     * \n jobResizeReleaseLog
     * \n jobResizeCancelLog
     * \n jobResizeLog
     * \n jRusage
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref event_types
     * \n \ref defs_lsb_XF_OP
     * \n \ref jobgroup_controltypes
     * \n \ref signal_action
     * #see \ref lsb_closestream
     * #see \ref lsb_geteventrec
     * #see \ref lsb_openstream
     * #see \ref lsb_puteventrec
     * #see \ref lsb_readstreamline
     * #see \ref lsb_streamversion
     * #see \ref lsb_readstream
     */
    // NOTE: Not in libbat
    //public static native int lsb_writestream(eventRec logPtr);

     /**//*
     * \page lsb_readstream lsb_readstream
     * \brief Reads a current version eventRec structure from the lsb_stream file.
     * <p/>
     * \ref lsb_readstream reads an eventrRec from the open streamFile.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * eventRec lsb_readstream(IntByReference nline)</b>
     *
     * @return int:-1 \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param nline Line number in the stream file to be read.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * eventRec
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_closestream
     * #see \ref lsb_geteventrec
     * #see \ref lsb_openstream
     * #see \ref lsb_puteventrec
     * #see \ref lsb_readstreamline
     * #see \ref lsb_streamversion
     * #see \ref lsb_writestream
     */
    // NOTE: Not in libbat
    //public static native eventRec.ByReference lsb_readstream(IntByReference nline);

     /**//*
     * \page lsb_readstreamline lsb_readstreamline
     * \brief Reads a current version eventRec structure from the lsb_stream file.
     * <p/>
     * \ref lsb_readstreamline reads an eventrRec from the open streamFile
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * eventRec.ByReference lsb_readstreamline(String line)</b>
     *
     * @return null \n
     *         The function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         On failure, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         $LSB_CONFDIR/cluster_name/configdir/lsb.params
     * @param line Line number in the stream file to be read.
     * See \ref lsb_puteventrec and \ref lsb_geteventrec for details on the eventRec structure.
     * Additionally, there are three additional event types supported.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * eventRec
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * #see \ref lsb_closestream
     * #see \ref lsb_geteventrec
     * #see \ref lsb_openstream
     * #see \ref lsb_puteventrec
     * #see \ref lsb_readstream
     * #see \ref lsb_streamversion
     * #see \ref lsb_writestream
     */
    // NOTE: Not in libbat
    //public static native eventRec.ByReference lsb_readstreamline(String line);

    public static final int NUM_EXITRATE_TYPES = 4;

/* options for exit rate type */


/* all exited jobs */
    public static final int JOB_EXIT = 0x01;

/* jobs failed to start due to initialization problem on execution host*/
    public static final int JOB_INIT = 0x02;

/* jobs failed to start due to HPC specific initialization problem on execution host*/
    public static final int HPC_INIT = 0x04;

/* jobs exited not related to reasons set by LSF */
    public static final int JOB_EXIT_NONLSF = 0x08;

    /**
     * \brief  APS factor information
     */
    public static class apsFactorInfo extends Structure {
        public static class ByReference extends apsFactorInfo implements Structure.ByReference {}
        public static class ByValue extends apsFactorInfo implements Structure.ByValue {}
        public apsFactorInfo() {}
        public apsFactorInfo(Pointer p) { super(p); read(); }


        /**
         * <  Name
         */
        public String name;

        /**
         * <  Weight
         */
        public float weight;

        /**
         * <  Limit
         */
        public float limit;

        /**
         * <  Grace period
         */
        public int gracePeriod;
    }



/* options for job group delete */

/* delete the specified user's all empty job groups*/
    public static final int JGRP_DEL_USER_GROUPS = 0x01;

/* delete one job group's all empty children groups including itself*/
    public static final int JGRP_DEL_CHILD_GROUPS = 0x02;

/* delete all empty job groups */
    public static final int JGRP_DEL_ALL = 0x04;

    /**
     * ------------------------------------------------------------------------
     * lsb_getallocFromHhostfile
     * <p/>
     * Read the specified hostfile and return the host list. If path is null
     * then read the hostfile specified by LSB_DJOB_HOSTFILE. The hostfile
     * is assumed to be in simple format of one host per line. A host
     * can be repeated.
     * <p/>
     * This function will allocate the memory for hostlist.
     * It is the responsibility of the caller to free the lists when no longer
     * needed. On success hostlist will be a list of strings.
     * Before freeing hostlist the individual
     * elements must be freed.
     * <p/>
     * Parameters:
     * @param hostlist  [OUT]
     * @param path      [IN]    path to hostfile, if null check LSB_DJOB_HOSTFILE
     * <p/>
     * @return Value:
     * >0    success, length of hostlist not including the null last element
     * -1    failure, lsberrno is set
     * -------------------------------------------------------------------------
     */
    public static native int lsb_getallocFromHostfile(Pointer hostlist, String path);


    /**
     *  \addtogroup defs_lsb_launch defs_lsb_launch
     *  lsb_launch() Valid options are:
     */

    /**
     * < Disable standard input and redirect input from the special  device /dev/null. This is equivalent to blaunch -n.
     */
    public static final int LSF_DJOB_DISABLE_STDIN = 0x01;

    /**
     * < Replace existing enviornment variable values with envp.
     */
    public static final int LSF_DJOB_REPLACE_ENV = 0x02;

    /**
     * < Non-blocking mode; the parallel job does not wait once all tasks start.  This forces \ref lsb_launch not to wait for its tasks to finish.
     */
    public static final int LSF_DJOB_NOWAIT = 0x04;

    /**
     * < Display standard error messages with a corresponding host name where the message was generated.Cannot be specified with LSF_DJOB_NOWAIT.
     */
    public static final int LSF_DJOB_STDERR_WITH_HOSTNAME = 0x08;

    /**
     * < Display standard output messages with a corresponding host name  where the message was generated. Cannot be specified with LSF_DJOB_NOWAIT.
     */
    public static final int LSF_DJOB_STDOUT_WITH_HOSTNAME = 0x10;

    /**
     * < Use user's login shell to  launch tasks
     */
    public static final int LSF_DJOB_USE_LOGIN_SHELL = 0x20;

    /**
     * < Use /bin/sh to launch tasks
     */
    public static final int LSF_DJOB_USE_BOURNE_SHELL = 0x40;

    /**
     * < Separate stderr from stdout
     */
    public static final int LSF_DJOB_STDERR = 0x80;

/*
* -------------------------------------------------------------------------
*  lsb_launch (where, argv, options, envp)
*
*  DESCRIPTION:
*
*    The specified command (i.e., argv) will be launched on the remote
*    nodes in parallel
*
*  ARGUMENTS:
*    where [IN]:
*        A null terminated list of hosts.
*        If this parameter is null then the environment variable
*        LSB_MCPU_HOSTS will be used.
*        A task will be launched for each slot.
*    options [IN]:
*        options value obtained by ORing
*    Envp [IN]:
*        A Null terminated list of environment variables (in 'variable=value'
*        format).
*        The environment to set for each task.
*        If this parameter is null then the same environment used to start
*        the first task will be used.
*        If non-null, it is appended to the environment used for the
*        first task.
*        If LSF_DJOB_REPLACE_ENV is specified, Envp entries will overwrite
*        existing values except those LSF needs.
*
*  RETURN:
*    < 0 on failure
*    > 0 upon success (i.e., number of tasks issued)
*
 */

    /**
     * \page lsb_launch lsb_launch
     * \brief  Launch commands on remote hosts in parallel.
     * <p/>
     * \ref lsb_launch is a synchronous API call to allow source level integration with
     * vendor MPI implementations. This API will launch the specified command (argv)
     * on the remote nodes in parallel.
     * \n LSF must be installed before integrating your MPI implementation with
     * \ref lsb_launch. The \ref lsb_launch API requires the full set of liblsf.so,
     * libbat.so (or liblsf.a, libbat.a).
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_launch (String[] where, String[] argv, int userOptions, String[] envp)</b>
     *
     * @param userOptions [IN] Options to modify the behavior of \ref lsb_launch
     *                    Multiple option values can be specified. For example option values can be
     *                    separated by OR (|):
     *                    \n \ref lsb_launch (where, argv, LSF_DJOB_REPLACE_ENV | LSF_DJOB_DISABLE_STDIN, envp);
     * @return < 0 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line command:</b>
     *         \par
     *         blaunch
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param where [IN] A null-terminated list of hosts. A task will be launched
     * for each slot.If this parameter is null then the environment variable
     * LSB_MCPU_HOSTS will be used.
     * @param argv [IN] The command to be executed
     * @param envp [IN] A null-terminated list of environment variables specifying
     * the environment to set for each task.If envp is null, \ref lsb_launch uses the
     * same environment used to start the first task on the first execution host.
     * If non-null, envp values are appended to the environment used for the first
     * task.If the LSF_DJOB_REPLACE_ENV option is specified, envp entries will
     * overwrite all existing environment values except those needed by LSF.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref defs_lsb_launch
     * see none
     */
    public static native int lsb_launch(Pointer where, Pointer argv, int userOptions, Pointer envp);

/*
* -------------------------------------------------------------------------
*  lsb_getalloc
*
*  This function will allocate the memory for hostlist.
*
*  It is the responsibility of the caller to free the lists when no longer
*  needed. On success hostlist will be a list of strings.
*  Before freeing hostlist the individual
*  elements must be freed.
*
*  Parameters:
*     hostlist     [OUT]     null terminated list of hosts
*
*  Returns:
*    >0    success, length of hostlist not including the null last element
*    -1    failure, lsberrno is set
*
* -------------------------------------------------------------------------
 */

    /**
     * \page lsb_getalloc lsb_getalloc
     * \brief Allocates memory for a host list to be used for launching parallel
     * tasks through blaunch and the \ref lsb_launch API.
     * <p/>
     * It is the responsibility of the caller to free the host list when it is
     * no longer needed.On success, the host list will be a list of strings.
     * Before freeing host list, the individual elements must be freed.
     * <p/>
     * An application using the \ref lsb_getalloc API is assumed to be part of an
     * LSF job, and that LSB_MCPU_HOSTS is set in the environment.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_getalloc(String[][] hostlist)</b>
     *
     * @return < 0 \n
     *         Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         If the function fails, lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param hostlist [OUT] A null-terminated list of host names
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * none
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * none
     * see none
     */
    public static native int lsb_getalloc(Pointer hostlist);

    /**
     * \page lsb_resize_cancel lsb_resize_cancel
     * \brief Cancels a pending job resize allocation request.
     * <p/>
     * Use \ref lsb_resize_cancel to cancel a pending allocation request for a
     * resizable job. A running job can only have one pending request at any
     * particular time. If one request is still pending, additional requests
     * are rejected with a proper error code.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_resize_cancel(long jobId);</b>
     *
     * @param jobId LSF job ID
     *              <p/>
     *              <b>Data Structures:</b>
     *              \par
     *              none
     *              <p/>
     *              <b>Define Statements:</b>
     *              \par
     *              none
     * @return int:-1 \n
     *         On failure, returns -1.
     *         <p/>
     *         \b Errors:
     *         \par
     *         lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         bresize cancel job_ID
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * #see \ref lsb_resize_release
     */

    public static native int lsb_resize_cancel(long jobId);

    /**
     * \page lsb_resize_release lsb_resize_release
     * \brief Releases part of the allocation of a running resizable job.
     * <p/>
     * Use \ref lsb_resize_release to release part of a running job allocation.
     * A running job can only have one pending request at any particular time.
     * If one request is still pending, additional requests are rejected with
     * a proper error code.
     * <p/>
     * If a notification command is defined through job submission, application
     * profile,or the \ref lsb_resize_release API, the notification command is invoked
     * on the first execution host of the job allocation once allocation resize
     * requests have been satisfied.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * int lsb_resize_release(job_resize_release.ByReference req);</b>
     *
     * @return int:-1 \n
     *         On failure, returns -1.
     *         <p/>
     *         \b Errors:
     *         \par
     *         lsberrno is set to indicate the error.
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         release [-c] [-rnc resize_notification_cmd | -rncn] released_host_specification job_ID
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param req job resize release request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * job_resize_release
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref resizablejob_related
     * #see \ref lsb_resize_cancel
     */
    public static native int lsb_resize_release(job_resize_release req);

    public static native int lsb_resize_request(job_resize_request job_resize_request1);

    /**
     * \page  lsb_getjobdepinfo lsb_getjobdepinfo
     * Returns the job dependency information.
     * <p/>
     * \ref lsb_getjobdepinfo returns information about jobs (including job arrays) when
     * a job has one or more dependencies on it.
     * <p/>
     * <b>\#include <lsf/lsbatch.h>
     * <p/>
     * jobDependInfo.ByReference
     * lsb_getjobdepinfo(jobDepRequest.ByReference jobdepReq)</b>
     *
     * @return null
     *         \n Function failed.
     *         <p/>
     *         \b Errors:
     *         \par
     *         none
     *         <p/>
     *         <b>Equivalent line commands:</b>
     *         \par
     *         none
     *         <p/>
     *         <b>Files:</b>
     *         \par
     *         none
     * @param jobdepReq Job dependent Request.
     * <p/>
     * <b>Data Structures:</b>
     * \par
     * dependJobs
     * \n queriedJobs
     * \n jobDependInfo
     * \n jobDepRequest
     * <p/>
     * <b>Define Statements:</b>
     * \par
     * \ref job_has_depend
     * \n \ref query_depend
     */
    public static native jobDependInfo.ByReference lsb_getjobdepinfo(jobDepRequest jobdepReq);


    /**
     *  \page lsb_jsdl2submit lsb_jsdl2submit
     *  \brief  Accepts a JSDL job submission file as input and converts the file
     *   for use with LSF.
     *
     *  \ref lsb_jsdl2submit converts parameters specified in the JSDL file and merges
     *  them with the other command line and job script options. The merged submit
     *  request is then sent to mbatchd for processing.
     *
     *  Code must link to LSF_LIBDIR/libbat.jsdl.lib
     *
     *  <b>\#include <lsf/lsbatch.h>
     *
     *  int lsb_jsdl2submit(submit.ByReference req, String template);</b>
     *
     *  @param req Reads the specified JSDL options and maps them to the
     *  submitReq structure. Code must specify either jsdl or jsdl_strict.
     *  @param template The default template, which contains all of the bsub
     *  submission options.
     *
     *  <b>Data Structures:</b>
     *  \par
     *  submit
     *
     *  <b>Define Statements:</b>
     *  \par
     *  none
     *
     *  @return int:0 \n
     *  Function completed successfully.
     *  @return int:-1 \n
     *  Function failed.
     *
     *  <b>Errors:</b>
     *  \par
     *  On failure, sets lsberrno to indicate the error.
     *
     *  <b>Equivalent line command:</b>
     *  \par
     *   bsub with options
     *
     *  <b>Files:</b>
     *  \par
     *  $LSF_LIBDIR/jsdl.xsd
     *  \n $LSF_LIBDIR/jsdl-posix.xsd
     *  \n $LSF_LIBDIR/jsdl-lsf.xsd
     *
     *  @see \ref lsb_submit
     *  @see \ref lsb_modify
     */

    /**
     *  \page lsblib lsblib
     *  \brief Application Programming Interface (API) library functions for batch jobs
     *
     *  LSBLIB functions allow application programs to get information about the hosts,
     *  queues, users, jobs and configuration of the batch system. Application programs
     *  can also submit jobs and control hosts, queues and jobs. Finally, application
     *  programs can read batch log files and write batch error messages.
     *
     *  \note
     *  \par
     *  All LSBLIB APIs require that the batch header file <lsf/lsbatch.h> be included.
     *  \par
     *  Many LSBLIB APIs return a pointer to an array or structure. These data structures
     *  are in static storage or on the heap. The next time the API is called, the storage
     *  is overwritten or freed.
     *  \par
     *  Any program using LSBLIB APIs that change the state of the batch system (that
     *  is, except for APIs that just get information about the system) must be setuid
     *  to root if LSF_AUTH is not defined in the lsf.conf file.
     *  \par
     *  On systems which have both System V and BSD programming interfaces, LSBLIB
     *  typically requires the BSD programming interface. On System V-based versions of
     *  UNIX, for example SGI IRIX, it is normally necessary to link applications using
     *  LSBLIB with the BSD compatibility library.
     *  \par
     *  On AFS systems, the following needs to be added to the end of your linkage
     *  specifications when linking with LSBLIB (assuming your AFS library path is
     *  /usr/afsws):
     *  \par
     *  For HP-UX and Solaris,
     *  \par
     *  -lc -L/usr/afsws/lib -L/usr/afsws/lib/afs -lsys -lrx -llwp /usr/afsws/lib/afs/util.a
     *  \par
     *  For other platforms,
     *  \par
     *  -lc -L/usr/afsws/lib -L/usr/afsws/lib/afs -lsys -lrx -llwp
     *
     *  \b Files:
     *  \par
     *  ${LSF_ENVDIR:-/etc}/lsf.conf
     *  \n $LSF_CONFDIR/lsf.shared
     *  \n $LSF_CONFDIR/lsf.cluster.cluster_name
     *  \n $LSF_CONFDIR/lsf.task
     *  \n $LSF_CONFDIR/lsf.task.cluster_name
     *  \n $LSB_CONFDIR/cluster_name/configdir/lsb.hosts
     *  \n $$LSB_CONFDIR/cluster_name/configdir/lsb.params
     *  \n $LSB_CONFDIR/cluster_name/configdir/lsb.queues
     *  \n $LSB_CONFDIR/cluster_name/configdir/lsb.users
     *
     *  @see lsblibapis
     */
}
