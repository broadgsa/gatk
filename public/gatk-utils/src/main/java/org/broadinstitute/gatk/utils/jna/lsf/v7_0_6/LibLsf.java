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
import com.sun.jna.ptr.FloatByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import org.broadinstitute.gatk.utils.jna.clibrary.JNAUtils;
import org.broadinstitute.gatk.utils.jna.clibrary.LibC.timeval;

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
 * JNA wrappers for LSF's lsf.h and -llsf
 *
 * $Id: base.h,v 1.25.6.12.2.5.2.11.2.15 2009/08/17 07:25:05 qlnie Exp $
 ****************************************************************************
 *
 * Load Sharing Facility
 *
 * Header file for all components of load sharing facility.
 *
 ****************************************************************************/
@SuppressWarnings("unused")
public class LibLsf {

    static {
        /*
        LSF 7.0.6 on the mac is missing the unsatisfied exported symbol for environ which was removed on MacOS X 10.5+.
        nm $LSF_LIBDIR/liblsf.dylib | grep environ
        See "man environ" for more info, along with http://lists.apple.com/archives/java-dev/2007/Dec/msg00096.html
        For now, we export environ ourselves using libenvironhack.dylib available in c/libenvironhack.
        */
        if (Platform.isMac())
            NativeLibrary.getInstance("environhack");
        String lsfLibDir = System.getenv("LSF_LIBDIR");
        if (lsfLibDir != null) {
            NativeLibrary.addSearchPath("lsf", lsfLibDir);
        }
        Native.register("lsf");
    }

    public static final String PASSWD_FILE_LS = "passwd.lsfuser";
    public static final int PASSWORD_LEN = 64;
    public static final int MAXHOSTNAMELEN = JNAUtils.MAXHOSTNAMELEN;
    public static final int MAXPATHLEN = JNAUtils.MAXPATHLEN;


    public static final int LOG_EMERG = 0;
    public static final int LOG_ALERT = 1;
    public static final int LOG_CRIT = 2;
    public static final int LOG_ERR = 3;
    public static final int LOG_WARNING = 4;
    public static final int LOG_NOTICE = 5;
    public static final int LOG_INFO = 6;
    public static final int LOG_DEBUG = 7;

    public static final int INVALID_SOCKET = -1;

    public static boolean SOCK_INVALID(int c) {
        return ((c) == INVALID_SOCKET);
    }

    public static class rlimit extends Structure {
        public static class ByReference extends rlimit implements Structure.ByReference {}
        public static class ByValue extends rlimit implements Structure.ByValue {}
        public rlimit() {}
        public rlimit(Pointer p) { super(p); read(); }

        public NativeLong rlim_cur;
        public NativeLong rlim_max;
    }



    public static class rusage extends Structure {
        public static class ByReference extends rusage implements Structure.ByReference {}
        public static class ByValue extends rusage implements Structure.ByValue {}
        public rusage() {}
        public rusage(Pointer p) { super(p); read(); }

        public timeval ru_utime;
        public timeval ru_stime;


        public NativeLong ru_maxrss;
        //public static final int ru_first = ru_ixrss;
        public NativeLong ru_ixrss;
        public NativeLong ru_idrss;
        public NativeLong ru_isrss;
        public NativeLong ru_minflt;
        public NativeLong ru_majflt;
        public NativeLong ru_nswap;
        public NativeLong ru_inblock;
        public NativeLong ru_oublock;
        public NativeLong ru_msgsnd;
        public NativeLong ru_msgrcv;
        public NativeLong ru_nsignals;
        public NativeLong ru_nvcsw;
        public NativeLong ru_nivcsw;
        //public static final int ru_last = ru_nivcsw;
        // Listed in lsf.h but not present in structure.
        //public NativeLong ru_ioch;
    }




    public static final String _VERSION_STR_ = "Platform LSF 7.0";
    public static final String _WORKGROUP_STR_ = "";
    public static final String _MINOR_STR_ = "";
    public static final String _BUILD_STR_ = "";
    public static final String _NOTE_STR_ = "";
    public static final String _HOTFIX_STR_ = "";
    public static final String _OS_STR_ = "";

    public static final String _DATE_STR_ = "";
    public static final String _BUILD_INFO_ = _MINOR_STR_ + "" + _BUILD_STR_ + "" + _WORKGROUP_STR_ + ", " + _DATE_STR_ + "\nCopyright 1992-2009 Platform Computing Corporation\n\n" + _OS_STR_ + _NOTE_STR_ + _HOTFIX_STR_;
    public static final String _LS_VERSION_ = (_VERSION_STR_ + "" + _BUILD_INFO_);

    //public static int XDR_SETPOS (int xdrs, int pos)  { (*(xdrs)->x_ops->x_setpostn)(xdrs, 0); return (*(xdrs)->x_ops->x_setpostn)(xdrs, pos); }
    //public static int xdr_setpos (int xdrs, int pos)  { (*(xdrs)->x_ops->x_setpostn)(xdrs, 0); return (*(xdrs)->x_ops->x_setpostn)(xdrs, pos); }


    public static final int LSF_XDR_VERSION2_0 = 1;
    public static final int LSF_XDR_VERSION2_1 = 2;
    public static final int LSF_XDR_VERSION2_2 = 3;
    public static final int LSF_XDR_VERSION3_0 = 4;
    public static final int LSF_XDR_VERSION3_1 = 5;
    public static final int LSF_XDR_VERSION3_2 = 6;
    public static final int LSF_XDR_VERSION3_2_2 = 7;
    public static final int LSF_XDR_VERSION4_0 = 8;
    public static final int LSF_XDR_VERSION4_1 = 9;
    public static final int LSF_XDR_VERSION4_2 = 10;
    public static final int LSF_XDR_VERSION5_0 = 11;
    public static final int LSF_XDR_VERSION5_1 = 12;
    public static final int LSF_XDR_VERSION6_0 = 13;
    public static final int LSF_XDR_VERSION6_1 = 14;
    public static final int LSF_XDR_VERSION6_2 = 15;
    public static final int EGO_XDR_VERSION_1_1 = 16;
    public static final int LSF_XDR_VERSION7_0 = 17;
    public static final int EGO_XDR_VERSION_1_2 = LSF_XDR_VERSION7_0;
    public static final int LSF_XDR_VERSION7_0_EP1 = 18;
    public static final int LSF_XDR_VERSION7_0_EP2 = 19;
    public static final int LSF_XDR_VERSION7_0_EP3 = 20;
    public static final int LSF_XDR_VERSION7_0_EP4 = 21;
    public static final int LSF_XDR_VERSION7_0_EP5 = 22;
    public static final int LSF_XDR_VERSION7_0_EP6 = 23;
    public static final int EGO_XDR_VERSION_1_2_2 = LSF_XDR_VERSION7_0_EP1;
    public static final int EGO_XDR_VERSION_1_2_3 = LSF_XDR_VERSION7_0_EP2;

    public static final int EGO_XDR_VERSION = LSF_XDR_VERSION7_0_EP2;

    //public String LOG_VERSION;

    public static final int LSF_DEFAULT_SOCKS = 15;
    public static final int MAXLINELEN = 512;
    public static final int MAXLSFNAMELEN = 40;
    public static final int MAXLSFNAMELEN_70_EP1 = 128;

    public static final int MAXSRES = 32;
    public static final int MAXRESDESLEN = 256;
    public static final int NBUILTINDEX = 11;
    public static final int MAXTYPES = 128;
    public static final int MAXMODELS = 1024 + 2;
    public static final int MAXMODELS_70 = 128;
    public static final int MAXTYPES_31 = 25;
    public static final int MAXMODELS_31 = 30;
    public static final int MAXFILENAMELEN = 256;
    public static final int MAXEVARS = 30;

    public static final int GENMALLOCPACE = 1024;


    public static final int FIRST_RES_SOCK = 20;


    public static final int R15S = 0;
    public static final int R1M = 1;
    public static final int R15M = 2;
    public static final int UT = 3;
    public static final int PG = 4;
    public static final int IO = 5;
    public static final int LS = 6;
    public static final int IT = 7;
    public static final int TMP = 8;
    public static final int SWP = 9;
    public static final int MEM = 10;
    public static final int USR1 = 11;
    public static final int USR2 = 12;


    public static final float INFINIT_LOAD = (float) (0x7fffffff);
    public static final float INFINIT_FLOAT = (float) (0x7fffffff);

    public static final int INFINIT_INT = 0x7fffffff;
    public static final long INFINIT_LONG_INT = 0x7fffffffffffffffL;
    public static final short INFINIT_SHORT = 0x7fff;

    public static final int DEFAULT_RLIMIT = -1;

    public static final int LSF_RLIMIT_CPU = 0;
    public static final int LSF_RLIMIT_FSIZE = 1;
    public static final int LSF_RLIMIT_DATA = 2;
    public static final int LSF_RLIMIT_STACK = 3;
    public static final int LSF_RLIMIT_CORE = 4;
    public static final int LSF_RLIMIT_RSS = 5;
    public static final int LSF_RLIMIT_NOFILE = 6;
    public static final int LSF_RLIMIT_OPEN_MAX = 7;
    public static final int LSF_RLIMIT_VMEM = 8;
    public static final int LSF_RLIMIT_SWAP = LSF_RLIMIT_VMEM;
    public static final int LSF_RLIMIT_RUN = 9;
    public static final int LSF_RLIMIT_PROCESS = 10;
    public static final int LSF_RLIMIT_THREAD = 11;
    public static final int LSF_RLIM_NLIMITS = 12;

    public static final int LSF_RLIM_NLIMITS5_1 = 11;

    //public static int seteuid (int x) { return setresuid(-1,x,-1); }
    //public static int setegid (int x) { return setresgid(-1,x,-1); }

    public static final int LSF_NULL_MODE = 0;
    public static final int LSF_LOCAL_MODE = 1;
    public static final int LSF_REMOTE_MODE = 2;


    public static final int RF_MAXHOSTS = 5;


    public static final int RF_CMD_MAXHOSTS = 0;


    public static final int RF_CMD_RXFLAGS = 2;


    public static final int STATUS_TIMEOUT = 125;
    public static final int STATUS_IOERR = 124;
    public static final int STATUS_EXCESS = 123;
    public static final int STATUS_REX_NOMEM = 122;
    public static final int STATUS_REX_FATAL = 121;
    public static final int STATUS_REX_CWD = 120;
    public static final int STATUS_REX_PTY = 119;
    public static final int STATUS_REX_SP = 118;
    public static final int STATUS_REX_FORK = 117;
    public static final int STATUS_REX_AFS = 116;
    public static final int STATUS_REX_UNKNOWN = 115;
    public static final int STATUS_REX_NOVCL = 114;
    public static final int STATUS_REX_NOSYM = 113;
    public static final int STATUS_REX_VCL_INIT = 112;
    public static final int STATUS_REX_VCL_SPAWN = 111;
    public static final int STATUS_REX_EXEC = 110;
    public static final int STATUS_REX_MLS_INVAL = 109;
    public static final int STATUS_REX_MLS_CLEAR = 108;
    public static final int STATUS_REX_MLS_RHOST = 107;
    public static final int STATUS_REX_MLS_DOMIN = 106;
    public static final int STATUS_DENIED = 105;


    public static boolean REX_FATAL_ERROR(int s) {
        return (((s) == STATUS_REX_NOVCL) || ((s) == STATUS_REX_NOSYM) || ((s) == STATUS_REX_NOMEM) || ((s) == STATUS_REX_FATAL) || ((s) == STATUS_REX_CWD) || ((s) == STATUS_REX_PTY) || ((s) == STATUS_REX_VCL_INIT) || ((s) == STATUS_REX_VCL_SPAWN) || ((s) == STATUS_REX_MLS_INVAL) || ((s) == STATUS_REX_MLS_CLEAR) || ((s) == STATUS_REX_MLS_RHOST) || ((s) == STATUS_REX_MLS_DOMIN));
    }


    public static final int REXF_USEPTY = 0x00000001;
    public static final int REXF_CLNTDIR = 0x00000002;
    public static final int REXF_TASKPORT = 0x00000004;
    public static final int REXF_SHMODE = 0x00000008;
    public static final int REXF_TASKINFO = 0x00000010;
    public static final int REXF_REQVCL = 0x00000020;
    public static final int REXF_SYNCNIOS = 0x00000040;
    public static final int REXF_TTYASYNC = 0x00000080;
    public static final int REXF_STDERR = 0x00000100;


    public static final int EXACT = 0x01;
    public static final int OK_ONLY = 0x02;
    public static final int NORMALIZE = 0x04;
    public static final int LOCALITY = 0x08;
    public static final int IGNORE_RES = 0x10;
    public static final int LOCAL_ONLY = 0x20;
    public static final int DFT_FROMTYPE = 0x40;
    public static final int ALL_CLUSTERS = 0x80;
    public static final int EFFECTIVE = 0x100;


    public static final int RECV_FROM_CLUSTERS = 0x200;
    public static final int NEED_MY_CLUSTER_NAME = 0x400;


    public static final int SEND_TO_CLUSTERS = 0x400;


    public static final int NO_SORT = 0x800;


    public static final int EXCLUSIVE_RESOURCE = 0x1000;

    public static final int DT_CLUSTER_LOAD = 0x2000;


    public static final int FROM_MASTER = 0x01;


    public static final int KEEPUID = 0x01;


    public static final int RES_CMD_REBOOT = 1;

    public static final int RES_CMD_SHUTDOWN = 2;

    public static final int RES_CMD_LOGON = 3;

    public static final int RES_CMD_LOGOFF = 4;


    public static final int LIM_CMD_REBOOT = 1;
    public static final int LIM_CMD_SHUTDOWN = 2;
    public static final int LIM_CMD_REMOVEHOST = 3;
    public static final int LIM_CMD_ACTIVATE = 4;
    public static final int LIM_CMD_DEACTIVATE = 5;
    public static final int LIM_CMD_ELIM_ENV = 6;


    public static class connectEnt extends Structure {
        public static class ByReference extends connectEnt implements Structure.ByReference {}
        public static class ByValue extends connectEnt implements Structure.ByValue {}
        public connectEnt() {}
        public connectEnt(Pointer p) { super(p); read(); }

        public String hostname;
        public int[] csock = new int[2];
    }



    public static final int INTEGER_BITS = 32;

    public static int GET_INTNUM(int i) {
        return ((i) / INTEGER_BITS + 1);
    }


    public static final int LIM_UNAVAIL = 0x00010000;
    public static final int LIM_LOCKEDU = 0x00020000;
    public static final int LIM_LOCKEDW = 0x00040000;
    public static final int LIM_BUSY = 0x00080000;
    public static final int LIM_RESDOWN = 0x00100000;
    public static final int LIM_UNLICENSED = 0x00200000;
    public static final int LIM_SBDDOWN = 0x00400000;
    public static final int LIM_LOCKEDM = 0x00800000;

    public static final int LIM_OK_MASK = 0x00bf0000;
    public static final int LIM_PEMDOWN = 0x01000000;
    public static final int LIM_LOCKEDU_RMS = 0x80000000;


    public static boolean LS_ISUNAVAIL(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_UNAVAIL) != 0));
    }


    public static boolean LS_ISBUSYON(int[] status, int index) {
        return (((status) != null) && (((status[1 + (index) / INTEGER_BITS]) & (1 << (index) % INTEGER_BITS)) != 0));
    }

    public static boolean LS_ISBUSY(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_BUSY) != 0));
    }


    public static boolean LS_ISRMSLOCK(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_LOCKEDU_RMS) != 0));
    }


    public static boolean LS_ISLOCKEDU(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_LOCKEDU) != 0));
    }


    public static boolean LS_ISLOCKEDW(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_LOCKEDW) != 0));
    }


    public static boolean LS_ISLOCKEDM(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_LOCKEDM) != 0));
    }


    public static boolean LS_ISLOCKED(int[] status) {
        return (((status) != null) && (((status[0]) & (LIM_LOCKEDU | LIM_LOCKEDW | LIM_LOCKEDM)) != 0));
    }


    public static boolean LS_ISRESDOWN(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_RESDOWN) != 0));
    }


    public static boolean LS_ISSBDDOWN(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_SBDDOWN) != 0));
    }

    public static boolean LS_ISPEMDOWN(int[] status) {
        return (((status[0]) & LIM_PEMDOWN) != 0);
    }


    public static boolean LS_ISUNLICENSED(int[] status) {
        return (((status) != null) && (((status[0]) & LIM_UNLICENSED) != 0));
    }


    public static boolean LS_ISOK(int[] status) {
        return (((status) != null) && ((status[0] & LIM_OK_MASK) == 0));
    }


    public static boolean LS_ISOKNRES(int[] status) {
        return (((status) != null) && (((status[0] & ~(LIM_LOCKEDU_RMS)) & ~(LIM_RESDOWN | LIM_SBDDOWN | LIM_PEMDOWN)) == 0));
    }


    public static class placeInfo extends Structure {
        public static class ByReference extends placeInfo implements Structure.ByReference {}
        public static class ByValue extends placeInfo implements Structure.ByValue {}
        public placeInfo() {}
        public placeInfo(Pointer p) { super(p); read(); }

        public byte[] hostName = new byte[MAXHOSTNAMELEN];
        public int numtask;
    }




    public static class hostLoad extends Structure {
        public static class ByReference extends hostLoad implements Structure.ByReference {}
        public static class ByValue extends hostLoad implements Structure.ByValue {}
        public hostLoad() {}
        public hostLoad(Pointer p) { super(p); read(); }

        public byte[] hostName = new byte[MAXHOSTNAMELEN];
        public IntByReference status;
        public FloatByReference li;
    }




    public static interface valueType {
          public static final int LS_BOOLEAN = 0;
          public static final int LS_NUMERIC = 1;
          public static final int LS_STRING = 2;
          public static final int LS_EXTERNAL = 3;
    }



    public static interface orderType {
          public static final int INCR = 0;
          public static final int DECR = 1;
          public static final int NA = 2;
    }




    public static final int RESF_BUILTIN = 0x01;
    public static final int RESF_DYNAMIC = 0x02;
    public static final int RESF_GLOBAL = 0x04;
    public static final int RESF_SHARED = 0x08;
    public static final int RESF_LIC = 0x10;
    public static final int RESF_EXTERNAL = 0x20;
    public static final int RESF_RELEASE = 0x40;
    public static final int RESF_DEFINED_IN_RESOURCEMAP = 0x80;

    public static final int RESF_NON_CONSUMABLE = 0x100;
    public static final int RESF_REDEFINABLE = 0x200;
    public static final int RESF_ESRES = 0x400;


    public static class resItem extends Structure {
        public static class ByReference extends resItem implements Structure.ByReference {}
        public static class ByValue extends resItem implements Structure.ByValue {}
        public resItem() {}
        public resItem(Pointer p) { super(p); read(); }

        public byte[] name = new byte[MAXLSFNAMELEN];
        public byte[] des = new byte[MAXRESDESLEN];
        public /*valueType*/ int valueType;
        public /*orderType*/ int orderType;
        public int flags;
        public int interval;
    }




    public static class lsInfo extends Structure {
        public static class ByReference extends lsInfo implements Structure.ByReference {}
        public static class ByValue extends lsInfo implements Structure.ByValue {}
        public lsInfo() {}
        public lsInfo(Pointer p) { super(p); read(); }

        // The current version of JNA's Structure.getNativeAlignment passes a "null" to
        // Native.getNativeSize() when accessing the contents of a 2D array.
        // Although the method is marked as protected, there are also multiple "TO DO"
        // comments so when we upgrade don't want to have specialized code floating around.

        public int nRes;
        public Pointer /* resItem.ByReference */ resTable;
        public int nTypes;
        public byte[] hostTypes = new byte[MAXTYPES * MAXLSFNAMELEN];
        public int nModels;
        public byte[] hostModels = new byte[MAXMODELS * MAXLSFNAMELEN];
        public byte[] hostArchs = new byte[MAXMODELS * MAXLSFNAMELEN_70_EP1];
        public int[] modelRefs = new int[MAXMODELS];
        public float[] cpuFactor = new float[MAXMODELS];
        public int numIndx;
        public int numUsrIndx;
    }




    public static final int CLUST_STAT_OK = 0x01;
    public static final int CLUST_STAT_UNAVAIL = 0x02;
    public static final int CLUST_STAT_RECV_FROM = 0x04;
    public static final int CLUST_STAT_SEND_TO = 0x08;


    public static boolean IS_DEFAULT_AUTH(byte[] auth) {
        return (auth == null || auth[0] == '\0');
    }


    public static class clusterInfo extends Structure {
        public static class ByReference extends clusterInfo implements Structure.ByReference {}
        public static class ByValue extends clusterInfo implements Structure.ByValue {}
        public clusterInfo() {}
        public clusterInfo(Pointer p) { super(p); read(); }

        public byte[] clusterName = new byte[MAXLSFNAMELEN];
        public int status;
        public byte[] masterName = new byte[MAXHOSTNAMELEN];
        public byte[] managerName = new byte[MAXLSFNAMELEN];
        public int managerId;
        public int numServers;
        public int numClients;
        public int nRes;
        public Pointer resources;
        public int nTypes;
        public Pointer hostTypes;
        public int nModels;
        public Pointer hostModels;
        public int nAdmins;
        public IntByReference adminIds;
        public Pointer admins;
        public int analyzerLicFlag;
        public int jsLicFlag;
        public byte[] afterHoursWindow = new byte[MAXLINELEN];
        public byte[] preferAuthName = new byte[MAXLSFNAMELEN];
        public byte[] inUseAuthName = new byte[MAXLSFNAMELEN];
    }


    public static class hostInfo extends Structure {
        public static class ByReference extends hostInfo implements Structure.ByReference {}
        public static class ByValue extends hostInfo implements Structure.ByValue {}
        public hostInfo() {}
        public hostInfo(Pointer p) { super(p); read(); }

        public byte[] hostName = new byte[MAXHOSTNAMELEN];
        public String hostType;
        public String hostModel;
        public float cpuFactor;
        public int maxCpus;
        public int maxMem;
        public int maxSwap;
        public int maxTmp;
        public int nDisks;
        public int nRes;
        public Pointer resources;
        public int nDRes;
        public Pointer DResources;
        public String windows;
        public int numIndx;
        public FloatByReference busyThreshold;
        public byte isServer;
        public byte licensed;
        public int rexPriority;
        public int licFeaturesNeeded;


        public static final int LSF_BASE_LIC = 0;
        public static final int LSF_BATCH_LIC_OBSOLETE = 1;
        public static final int LSF_JS_SCHEDULER_LIC = 2;
        public static final int LSF_JS_LIC = 3;
        public static final int LSF_CLIENT_LIC = 4;
        public static final int LSF_MC_LIC = 5;
        public static final int LSF_ANALYZER_SERVER_LIC = 6;
        public static final int LSF_MAKE_LIC = 7;

        public static final int LSF_PARALLEL_LIC = 8;
        public static final int LSF_FLOAT_CLIENT_LIC = 9;
        public static final int LSF_FTA_LIC = 10;
        public static final int LSF_AFTER_HOURS_LIC = 11;
        public static final int LSF_RESOURCE_PREEMPT_LIC = 12;
        public static final int LSF_BACCT_LIC = 13;
        public static final int LSF_SCHED_FAIRSHARE_LIC = 14;
        public static final int LSF_SCHED_RESERVE_LIC = 15;
        public static final int LSF_SCHED_PREEMPTION_LIC = 16;
        public static final int LSF_SCHED_PARALLEL_LIC = 17;
        public static final int LSF_SCHED_ADVRSV_LIC = 18;
        public static final int LSF_API_CLIENT_LIC = 19;

        public static final int CLUSTERWARE_MANAGER_LIC = 20;
        public static final int LSF_MANAGER_LIC = 21;
        public static final int LSF_PCC_HPC_LIC = 22;
        public static final int sCLUSTERWARE_LIC = 23;
        public static final int OTTAWA_MANAGER_LIC = 24;

        public static final int SYMPHONY_MANAGER_ONLINE_LIC = 25;
        public static final int SYMPHONY_MANAGER_BATCH_LIC = 26;
        public static final int SYMPHONY_SCHED_JOB_PRIORITY_LIC = 27;
        public static final int LSF_DUALCORE_X86_LIC = 28;
        public static final int LSF_TSCHED_LIC = 29;
        public static final int LSF_WORKGROUP_LIC = 30;
        public static final int LSF_NUM_LIC_TYPE = 31;
        public static final int LSF_WG_NUM_LIC_TYPE = 2;
        public static final int LSF_NO_NEED_LIC = 32;

        public int licClass;
        public int cores;
        public static final int INET6_ADDRSTRLEN = 46;
        public byte[] hostAddr = new byte[INET6_ADDRSTRLEN];
        public int pprocs;

        public int cores_per_proc;
        public int threads_per_core;
    }

    public static boolean HAS_BATCH_LICENSES(int featureEnabled) {
        return (JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.CLUSTERWARE_MANAGER_LIC)) || JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.LSF_MANAGER_LIC)) || JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.LSF_WORKGROUP_LIC)) || JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.SYMPHONY_MANAGER_ONLINE_LIC)) || JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.SYMPHONY_MANAGER_BATCH_LIC)));
    }

    public static boolean HAS_SYMPHONY_LICENSES(int featureEnabled) {
        return (JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.SYMPHONY_MANAGER_ONLINE_LIC)) || JNAUtils.toBoolean(featureEnabled & (1 << hostInfo.SYMPHONY_MANAGER_BATCH_LIC)));
    }


    public static class config_param extends Structure {
        public static class ByReference extends config_param implements Structure.ByReference {}
        public static class ByValue extends config_param implements Structure.ByValue {}
        public config_param() {}
        public config_param(Pointer p) { super(p); read(); }

        public String paramName;
        public String paramValue;
    }



    public static class lsfRusage extends Structure {
        public static class ByReference extends lsfRusage implements Structure.ByReference {}
        public static class ByValue extends lsfRusage implements Structure.ByValue {}
        public lsfRusage() {}
        public lsfRusage(Pointer p) { super(p); read(); }

        public double ru_utime;
        public double ru_stime;
        public double ru_maxrss;
        public double ru_ixrss;
        public double ru_ismrss;
        public double ru_idrss;
        public double ru_isrss;
        public double ru_minflt;
        public double ru_majflt;
        public double ru_nswap;
        public double ru_inblock;
        public double ru_oublock;
        public double ru_ioch;
        public double ru_msgsnd;
        public double ru_msgrcv;
        public double ru_nsignals;
        public double ru_nvcsw;
        public double ru_nivcsw;
        public double ru_exutime;
    }




    public static class lsfAcctRec extends Structure {
        public static class ByReference extends lsfAcctRec implements Structure.ByReference {}
        public static class ByValue extends lsfAcctRec implements Structure.ByValue {}
        public lsfAcctRec() {}
        public lsfAcctRec(Pointer p) { super(p); read(); }

        public int pid;
        public String username;
        public int exitStatus;
        public NativeLong dispTime;
        public NativeLong termTime;
        public String fromHost;
        public String execHost;
        public String cwd;
        public String cmdln;
        public lsfRusage lsfRu;
    }




    public static class confNode extends Structure {
        public static class ByReference extends confNode implements Structure.ByReference {}
        public static class ByValue extends confNode implements Structure.ByValue {}
        public confNode() {}
        public confNode(Pointer p) { super(p); read(); }

        public confNode.ByReference leftPtr;
        public confNode.ByReference rightPtr;
        public confNode.ByReference fwPtr;
        public String cond;
        public int beginLineNum;
        public int numLines;
        public Pointer lines;
        public byte tag;
    }



    public static class pStack extends Structure {
        public static class ByReference extends pStack implements Structure.ByReference {}
        public static class ByValue extends pStack implements Structure.ByValue {}
        public pStack() {}
        public pStack(Pointer p) { super(p); read(); }

        public int top;
        public int size;
        public PointerByReference nodes;
    }



    public static class confHandle extends Structure {
        public static class ByReference extends confHandle implements Structure.ByReference {}
        public static class ByValue extends confHandle implements Structure.ByValue {}
        public confHandle() {}
        public confHandle(Pointer p) { super(p); read(); }

        public confNode.ByReference rootNode;
        public String fname;
        public confNode.ByReference curNode;
        public int lineCount;
        public pStack.ByReference ptrStack;
    }



    public static class lsConf extends Structure {
        public static class ByReference extends lsConf implements Structure.ByReference {}
        public static class ByValue extends lsConf implements Structure.ByValue {}
        public lsConf() {}
        public lsConf(Pointer p) { super(p); read(); }

        public confHandle.ByReference confhandle;
        public int numConds;
        public Pointer conds;
        public IntByReference values;
    }



    public static class sharedConf extends Structure {
        public static class ByReference extends sharedConf implements Structure.ByReference {}
        public static class ByValue extends sharedConf implements Structure.ByValue {}
        public sharedConf() {}
        public sharedConf(Pointer p) { super(p); read(); }

        public lsInfo.ByReference lsinfo;
        public int numCls;
        public Pointer clusterNames;
        public Pointer servers;
    }




    public static class lsSharedResourceInstance extends Structure {
        public static class ByReference extends lsSharedResourceInstance implements Structure.ByReference {}
        public static class ByValue extends lsSharedResourceInstance implements Structure.ByValue {}
        public lsSharedResourceInstance() {}
        public lsSharedResourceInstance(Pointer p) { super(p); read(); }

        public String value;
        public int nHosts;
        public Pointer hostList;

    }




    public static class lsSharedResourceInfo extends Structure {
        public static class ByReference extends lsSharedResourceInfo implements Structure.ByReference {}
        public static class ByValue extends lsSharedResourceInfo implements Structure.ByValue {}
        public lsSharedResourceInfo() {}
        public lsSharedResourceInfo(Pointer p) { super(p); read(); }

        public String resourceName;
        public int nInstances;
        public Pointer /* lsSharedResourceInstance.ByReference */ instances;
    }



    public static class clusterConf extends Structure {
        public static class ByReference extends clusterConf implements Structure.ByReference {}
        public static class ByValue extends clusterConf implements Structure.ByValue {}
        public clusterConf() {}
        public clusterConf(Pointer p) { super(p); read(); }

        public clusterInfo.ByReference clinfo;
        public int numHosts;
        public Pointer /* hostInfo.ByReference */ hosts;
        public int defaultFeatures;
        public int numShareRes;
        public Pointer /* lsSharedResourceInfo.ByReference */ shareRes;
    }




    public static class pidInfo extends Structure {
        public static class ByReference extends pidInfo implements Structure.ByReference {}
        public static class ByValue extends pidInfo implements Structure.ByValue {}
        public pidInfo() {}
        public pidInfo(Pointer p) { super(p); read(); }

        public int pid;
        public int ppid;
        public int pgid;
        public int jobid;
    }




    public static class jRusage extends Structure {
        public static class ByReference extends jRusage implements Structure.ByReference {}
        public static class ByValue extends jRusage implements Structure.ByValue {}
        public jRusage() {}
        public jRusage(Pointer p) { super(p); read(); }

        public int mem;
        public int swap;
        public int utime;
        public int stime;
        public int npids;
        public Pointer /* pidInfo.ByReference */ pidInfo;

        public int npgids;
        public IntByReference pgid;
        public int nthreads;
    }




    public static final int NUM_SUBS = 2;
    public static final int LEN_SUBS = 64;
    public static final int NUM_CLASS_TYPE = 3;

    public static class licUsage extends Structure {
        public static class ByReference extends licUsage implements Structure.ByReference {}
        public static class ByValue extends licUsage implements Structure.ByValue {}
        public licUsage() {}
        public licUsage(Pointer p) { super(p); read(); }

        public int licDisplayMask;
        public int usingDemoLicense;
        public float[] total = new float[hostInfo.LSF_NUM_LIC_TYPE];
        public float[] inUse = new float[hostInfo.LSF_NUM_LIC_TYPE];
    }



    public static class hostClassInfo extends Structure {
        public static class ByReference extends hostClassInfo implements Structure.ByReference {}
        public static class ByValue extends hostClassInfo implements Structure.ByValue {}
        public hostClassInfo() {}
        public hostClassInfo(Pointer p) { super(p); read(); }

        public int numHosts;
        public int numCpus;
        public int numCores;
    }



    public static class lsfLicUsage extends Structure {
        public static class ByReference extends lsfLicUsage implements Structure.ByReference {}
        public static class ByValue extends lsfLicUsage implements Structure.ByValue {}
        public lsfLicUsage() {}
        public lsfLicUsage(Pointer p) { super(p); read(); }

        public licUsage licUsage;
        public hostClassInfo[] hostInfo = new hostClassInfo[NUM_CLASS_TYPE];
        // The current version of JNA's Structure.getNativeAlignment passes a "null" to
        // Native.getNativeSize() when accessing the contents of a 2D array.
        // Although the method is marked as protected, there are also multiple "TO DO"
        // comments so when we upgrade don't want to have specialized code floating around.
        public byte[] substitution = new byte[NUM_SUBS * LEN_SUBS];
        public byte[] cluster = new byte[MAXFILENAMELEN];
    }


    public static class param_entry extends Structure {
        public static class ByReference extends param_entry implements Structure.ByReference {}
        public static class ByValue extends param_entry implements Structure.ByValue {}
        public param_entry() {}
        public param_entry(Pointer p) { super(p); read(); }

        public static int HAS_PARAM_VALUE = 0x001;
        public static final int HAS_PARAM_DEFAULT = 0x002;

        public int flags;
        public String key;
        public String value;
        public String default_value;
    }



    public static class params_key_value_pair extends Structure {
        public static class ByReference extends params_key_value_pair implements Structure.ByReference {}
        public static class ByValue extends params_key_value_pair implements Structure.ByValue {}
        public params_key_value_pair() {}
        public params_key_value_pair(Pointer p) { super(p); read(); }

        public int num_params;
        public String daemon_time;
        public Pointer /* param_entry.ByReference */ param;
    }




    public static final int LSE_NO_ERR = 0;
    public static final int LSE_BAD_XDR = 1;
    public static final int LSE_MSG_SYS = 2;
    public static final int LSE_BAD_ARGS = 3;
    public static final int LSE_MASTR_UNKNW = 4;
    public static final int LSE_LIM_DOWN = 5;
    public static final int LSE_PROTOC_LIM = 6;
    public static final int LSE_SOCK_SYS = 7;
    public static final int LSE_ACCEPT_SYS = 8;
    public static final int LSE_BAD_TASKF = 9;
    public static final int LSE_NO_HOST = 10;
    public static final int LSE_NO_ELHOST = 11;
    public static final int LSE_TIME_OUT = 12;
    public static final int LSE_NIOS_DOWN = 13;
    public static final int LSE_LIM_DENIED = 14;
    public static final int LSE_LIM_IGNORE = 15;
    public static final int LSE_LIM_BADHOST = 16;
    public static final int LSE_LIM_ALOCKED = 17;
    public static final int LSE_LIM_NLOCKED = 18;
    public static final int LSE_LIM_BADMOD = 19;
    public static final int LSE_SIG_SYS = 20;
    public static final int LSE_BAD_EXP = 21;
    public static final int LSE_NORCHILD = 22;
    public static final int LSE_MALLOC = 23;
    public static final int LSE_LSFCONF = 24;
    public static final int LSE_BAD_ENV = 25;
    public static final int LSE_LIM_NREG = 26;
    public static final int LSE_RES_NREG = 27;
    public static final int LSE_RES_NOMORECONN = 28;
    public static final int LSE_BADUSER = 29;
    public static final int LSE_RES_ROOTSECURE = 30;
    public static final int LSE_RES_DENIED = 31;
    public static final int LSE_BAD_OPCODE = 32;
    public static final int LSE_PROTOC_RES = 33;
    public static final int LSE_RES_CALLBACK = 34;
    public static final int LSE_RES_NOMEM = 35;
    public static final int LSE_RES_FATAL = 36;
    public static final int LSE_RES_PTY = 37;
    public static final int LSE_RES_SOCK = 38;
    public static final int LSE_RES_FORK = 39;
    public static final int LSE_NOMORE_SOCK = 40;
    public static final int LSE_WDIR = 41;
    public static final int LSE_LOSTCON = 42;
    public static final int LSE_RES_INVCHILD = 43;
    public static final int LSE_RES_KILL = 44;
    public static final int LSE_PTYMODE = 45;
    public static final int LSE_BAD_HOST = 46;
    public static final int LSE_PROTOC_NIOS = 47;
    public static final int LSE_WAIT_SYS = 48;
    public static final int LSE_SETPARAM = 49;
    public static final int LSE_RPIDLISTLEN = 50;
    public static final int LSE_BAD_CLUSTER = 51;
    public static final int LSE_RES_VERSION = 52;
    public static final int LSE_EXECV_SYS = 53;
    public static final int LSE_RES_DIR = 54;
    public static final int LSE_RES_DIRW = 55;
    public static final int LSE_BAD_SERVID = 56;
    public static final int LSE_NLSF_HOST = 57;
    public static final int LSE_UNKWN_RESNAME = 58;
    public static final int LSE_UNKWN_RESVALUE = 59;
    public static final int LSE_TASKEXIST = 60;
    public static final int LSE_BAD_TID = 61;
    public static final int LSE_TOOMANYTASK = 62;
    public static final int LSE_LIMIT_SYS = 63;
    public static final int LSE_BAD_NAMELIST = 64;
    public static final int LSE_NO_LICENSE = 65;
    public static final int LSE_LIM_NOMEM = 66;
    public static final int LSE_NIO_INIT = 67;
    public static final int LSE_CONF_SYNTAX = 68;
    public static final int LSE_FILE_SYS = 69;
    public static final int LSE_CONN_SYS = 70;
    public static final int LSE_SELECT_SYS = 71;
    public static final int LSE_EOF = 72;
    public static final int LSE_ACCT_FORMAT = 73;
    public static final int LSE_BAD_TIME = 74;
    public static final int LSE_FORK = 75;
    public static final int LSE_PIPE = 76;
    public static final int LSE_ESUB = 77;
    public static final int LSE_DCE_EXEC = 78;
    public static final int LSE_EAUTH = 79;
    public static final int LSE_NO_FILE = 80;
    public static final int LSE_NO_CHAN = 81;
    public static final int LSE_BAD_CHAN = 82;
    public static final int LSE_INTERNAL = 83;
    public static final int LSE_PROTOCOL = 84;
    public static final int LSE_THRD_SYS = 85;
    public static final int LSE_MISC_SYS = 86;
    public static final int LSE_LOGON_FAIL = 87;
    public static final int LSE_RES_RUSAGE = 88;
    public static final int LSE_NO_RESOURCE = 89;
    public static final int LSE_BAD_RESOURCE = 90;
    public static final int LSE_RES_PARENT = 91;
    public static final int LSE_NO_PASSWD = 92;
    public static final int LSE_SUDOERS_CONF = 93;
    public static final int LSE_SUDOERS_ROOT = 94;
    public static final int LSE_I18N_SETLC = 95;
    public static final int LSE_I18N_CATOPEN = 96;
    public static final int LSE_I18N_NOMEM = 97;
    public static final int LSE_NO_MEM = 98;
    public static final int LSE_REGISTRY_SYS = 99;
    public static final int LSE_FILE_CLOSE = 100;
    public static final int LSE_LIMCONF_NOTREADY = 101;
    public static final int LSE_MASTER_LIM_DOWN = 102;
    public static final int LSE_MLS_INVALID = 103;
    public static final int LSE_MLS_CLEARANCE = 104;
    public static final int LSE_MLS_RHOST = 105;
    public static final int LSE_MLS_DOMINATE = 106;
    public static final int LSE_NO_CAL = 107;
    public static final int LSE_NO_NETWORK = 108;
    public static final int LSE_GETCONF_FAILED = 109;
    public static final int LSE_TSSINIT = 110;
    public static final int LSE_DYNM_DENIED = 111;
    public static final int LSE_LIC_OVERUSE = 112;
    public static final int LSE_EGOCONF = 113;
    public static final int LSE_BAD_EGO_ENV = 114;
    public static final int LSE_EGO_CONF_SYNTAX = 115;
    public static final int LSE_EGO_GETCONF_FAILED = 116;
    public static final int LSE_NS_LOOKUP = 117;
    public static final int LSE_BAD_PASSWD = 118;

    public static final int LSE_UNKWN_USER = 119;
    public static final int LSE_NOT_WINHOST = 120;
    public static final int LSE_NOT_MASTERCAND = 121;
    public static final int LSE_HOST_UNAUTH = 122;
    public static final int LSE_UNRESOLVALBE_HOST = 123;
    public static final int LSE_RESOURCE_NOT_CONSUMABLE = 124;
    public static final int LSE_SHUTDOWN = 125;
    public static final int LSE_BAD_SYNTAX = 126;
    public static final int LSE_NERR = 127;


    public static boolean LSE_ISBAD_RESREQ(int s) {
        return (((s) == LSE_BAD_EXP) || ((s) == LSE_UNKWN_RESNAME) || ((s) == LSE_UNKWN_RESVALUE));
    }

    public static boolean LSE_SYSCALL(int s) {
        return (((s) == LSE_SELECT_SYS) || ((s) == LSE_CONN_SYS) || ((s) == LSE_FILE_SYS) || ((s) == LSE_MSG_SYS) || ((s) == LSE_SOCK_SYS) || ((s) == LSE_ACCEPT_SYS) || ((s) == LSE_SIG_SYS) || ((s) == LSE_WAIT_SYS) || ((s) == LSE_EXECV_SYS) || ((s) == LSE_LIMIT_SYS) || ((s) == LSE_PIPE) || ((s) == LSE_ESUB) || ((s) == LSE_REGISTRY_SYS) || ((s) == LSE_MISC_SYS));
    }


    /*
    public static void TIMEVAL (int level, int func, int val)  {
        if (timinglevel > level) {
            timeval before, after;
            timezone tz;
            gettimeofday(&before, &tz);
            func;
            gettimeofday(&after, &tz);
            val = (int)((after.tv_sec - before.tv_sec)*1000 +  (after.tv_usec-before.tv_usec)/1000);
        } else {
            func;
            val = 0;
        }
    }
    */

    public static class ls_timeval extends Structure {
        public static class ByReference extends ls_timeval implements Structure.ByReference {}
        public static class ByValue extends ls_timeval implements Structure.ByValue {}
        public ls_timeval() {}
        public ls_timeval(Pointer p) { super(p); read(); }

        public float rtime;
        public float utime;
        public float stime;
    }



    /*
    public static void LS_TIMEVAL_ZERO(ls_timeval tv) {                            tv.rtime = 0.0;          tv.utime = 0.0;          tv.stime = 0.0;      }

    public static int LS_TIMEVAL_INC (ls_timeval tv, int newtv) {                                  tv.rtime += newtv.rtime;       tv.utime += newtv.utime;       tv.stime += newtv.stime;      }

    public static void LOG_TIME_MSG(int level, String name, ls_timeval tv, int count, String msg) { if (timinglevel > level) {  ls_syslog(LOG_INFO, "L%d %s rtime %.2f ms, utime %.2f ms, stime %.2f ms, count %d %s",  level, name, tv.rtime, tv.utime, tv.stime, count, msg);  } }; }

    public static void TIMEIT (int level, String func, String name) {
        if  (timinglevel > level && clockticks > 0) {
            timeval _before, _after;
            timezone _tz;
            tms _buf, _buf2;
            gettimeofday(&_before, &_tz);
            times(&_buf);
            func;
            gettimeofday(&_after, &_tz);
            times(&_buf2);
            ls_syslog(LOG_INFO,"L%d %s rtime %.2f ms, utime %.2f ms, stime %.2f ms",  level,  name,  (_after.tv_sec - _before.tv_sec)*1000.0 +  (_after.tv_usec - _before.tv_usec)/1000.0,  1000.0*((_buf2.tms_utime - _buf.tms_utime)/clockticks),  1000.0*((_buf2.tms_stime - _buf.tms_stime)/clockticks));
        } else {
            func;
        }
    }

    public static int TIMEVAL2 (int level, String func, ls_timeval tv) {
        if (timinglevel > level && clockticks > 0) {
            timeval _before, _after;
            timezone _tz;
            tms _buf, _buf2;
            gettimeofday(&_before, &_tz);
            times(&_buf);
            func;
            gettimeofday(&_after, &_tz);
            times(&_buf2);
            tv.rtime = (_after.tv_sec - _before.tv_sec)*1000.0 +  (_after.tv_usec - _before.tv_usec)/1000.0;
            tv.utime = 1000.0*((_buf2.tms_utime - _buf.tms_utime)/clockticks);
            tv.stime = 1000.0*((_buf2.tms_stime - _buf.tms_stime)/clockticks);
        } else {
            func;
            tv.rtime = 0.0;
            tv.utime = 0.0;
            tv.stime = 0.0;
        }
    }

    public static int TIMEIT_START_BLOCK (int level) {
        tms _buf, _buf2;
        timeval _before, _after;
        timezone _tz;
        if  (timinglevel > level) {
            gettimeofday(&_before, &_tz);
            times(&_buf);
        }
    }

    public static int TIMEIT_END_BLOCK (int level, String name)  {
        if  (timinglevel > level) {
            float rt, ut, st;
            gettimeofday(&_after, &_tz);
            times(&_buf2);
            rt = (_after.tv_sec - _before.tv_sec)*1000.0 +  (_after.tv_usec - _before.tv_usec)/1000.0;
            ut = 1000.0*((_buf2.tms_utime - _buf.tms_utime)/clockticks);
            st = 1000.0*((_buf2.tms_stime - _buf.tms_stime)/clockticks);
            ls_syslog(LOG_INFO,"L%d %s rtime %.2f ms, utime %.2f ms, stime %.2f ms",  level, name, rt, ut, st);
        }
    }
    */

    public static final int LC_SCHED = 0x00000001;
    public static final int LC_EXEC = 0x00000002;
    public static final int LC_TRACE = 0x00000004;
    public static final int LC_COMM = 0x00000008;
    public static final int LC_XDR = 0x00000010;
    public static final int LC_CHKPNT = 0x00000020;
    public static final int LC_LICENCE = 0x00000040;
    public static final int LC_LICENSE = 0x00000040;
    public static final int LC_FILE = 0x00000080;
    public static final int LC_AFS = 0x00000100;
    public static final int LC_AUTH = 0x00000200;
    public static final int LC_HANG = 0x00000400;
    public static final int LC_MULTI = 0x00000800;
    public static final int LC_SIGNAL = 0x00001000;
    public static final int LC_DCE = 0x00002000;
    public static final int LC_PIM = 0x00004000;
    public static final int LC_MEMORY = 0x00004000;
    public static final int LC_SYS = 0x00008000;
    public static final int LC_JLIMIT = 0x00010000;
    public static final int LC_FAIR = 0x00020000;
    public static final int LC_PREEMPT = 0x00040000;
    public static final int LC_PEND = 0x00080000;
    public static final int LC_EEVENTD = 0x00100000;
    public static final int LC_LOADINDX = 0x00200000;
    public static final int LC_RESOURCE = 0x00200000;

    public static final int LC_JGRP = 0x00400000;
    public static final int LC_JARRAY = 0x00800000;
    public static final int LC_MPI = 0x01000000;
    public static final int LC_ELIM = 0x02000000;
    public static final int LC_M_LOG = 0x04000000;
    public static final int LC_PERFM = 0x08000000;
    public static final int LC_DLOG = 0x10000000;
    public static final int LC_HPC = 0x20000000;
    public static final int LC_LICSCHED = 0x40000000;

    public static final int LC_XDRVERSION = 0x80000000;
    public static final int LC_FLEX = 0x80000000;

    public static final int LC_ADVRSV = LC_DLOG;
    public static final int LC_RESREQ = LC_M_LOG;


    public static final int LOG_DEBUG1 = LOG_DEBUG + 1;
    public static final int LOG_DEBUG2 = LOG_DEBUG + 2;
    public static final int LOG_DEBUG3 = LOG_DEBUG + 3;


    public static final int LSF_EVENT_LIM_DOWN = 1;
    public static final int LSF_EVENT_RES_DOWN = 2;
    public static final int LSF_EVENT_SBD_DOWN = 3;
    public static final int LSF_EVENT_HOST_UNLIC = 4;
    public static final int LSF_EVENT_MASTER_ELECT = 5;
    public static final int LSF_EVENT_MASTER_RESIGN = 6;
    public static final int LSF_EVENT_MBD_UP = 7;
    public static final int LSF_EVENT_MBD_DOWN = 8;
    public static final int LSF_EVENT_MBD_RECONFIG = 9;
    public static final int LSF_EVENT_WORKDIR_FULL = 10;
    public static final int LSF_EVENT_HOST_OPENED = 11;
    public static final int LSF_EVENT_HOST_CLOSED = 12;
    public static final int LSF_EVENT_QUEUE_OPENED = 13;
    public static final int LSF_EVENT_QUEUE_CLOSED = 14;
    public static final int LSF_EVENT_SCH_DOWN = 15;
    public static final int LSF_EVENT_LIC_OVERUSE = 16;

    public static final int LSF_NIOS_REQUEUE = 127;


    /*
    public int lserrno;
    public int masterLimDown;
    public int ls_nerr;
    public String[] ls_errmsg;
    public int logclass;
    public int timinglevel;
    public int clockticks;


    public int lsf_lim_version;
    */


    public static native int ls_readconfenv(config_param config_param1, String string);


    public static native Pointer ls_placereq(String resreq, IntByReference numhosts, int options, String fromhost);


    public static native Pointer ls_placeofhosts(String resreq, IntByReference numhosts, int options, String fromhost, Pointer hostlist, int listsize);

    // NOTE: Not in liblsf
    //public static native Pointer ls_placeoftype(String resreq, IntByReference numhosts, int options, String fromhost, String hosttype);


    public static native hostLoad.ByReference ls_load(String resreq, IntByReference numhosts, int options, String fromhost);


    public static native hostLoad.ByReference ls_loadofhosts(String resreq, IntByReference numhosts, int options, String fromhost, Pointer hostlist, int listsize);

    // NOTE: Not in liblsf
    //public static native hostLoad.ByReference ls_loadoftype(String resreq, IntByReference numhosts, int options, String fromhost, String hosttype);


    public static native hostLoad.ByReference ls_loadinfo(String resreq, IntByReference numhosts, int options, String fromhost, Pointer hostlist, int listsize, Pointer indxnamelist);


    public static native int ls_loadadj(String resreq, placeInfo hostlist, int listsize);


    public static native int ls_eligible(String task, String resreqstr, byte mode);


    public static native String ls_resreq(String task);


    public static native int ls_insertrtask(String task);


    public static native int ls_insertltask(String task);


    public static native int ls_deletertask(String task);


    public static native int ls_deleteltask(String task);


    public static native int ls_listrtask(Pointer taskList, int sortflag);


    public static native int ls_listltask(Pointer taskList, int sortflag);


    public static native Pointer ls_findmyconnections();


    public static native int ls_isconnected(String hostName);

    // NOTE: Not in liblsf
    //public static native int ls_lostconnection();


    public static native String ls_getclustername();


    public static native clusterInfo.ByReference ls_clusterinfo(String string1, IntByReference int1, Pointer stringArray1, int int2, int int3);


    public static native lsSharedResourceInfo.ByReference ls_sharedresourceinfo(Pointer stringArray1, IntByReference int1, String string1, int int2);


    public static native String ls_getmastername();


    public static native String ls_getmyhostname();


    public static native String ls_getmyhostname2();


    public static native hostInfo.ByReference ls_gethostinfo(String string1, IntByReference int1, Pointer stringArray1, int int2, int int3);

    public static native String ls_getISVmode();

    public static native int ls_isshutdown();

    public static native int ls_isPartialLicensingEnabled();

    /* NOTE: ls_getLicenseUsage() is not supported by LSF v8.x
    *  Wei Xing, ICR
    */
//    public static native lsfLicUsage.ByReference ls_getLicenseUsage();

    public static native lsInfo.ByReference ls_info();

    public static native Pointer ls_indexnames(lsInfo lsInfo1);

    public static native int ls_isclustername(String string);


    public static native String ls_gethosttype(String hostname);


    public static native FloatByReference ls_getmodelfactor(String modelname);


    public static native FloatByReference ls_gethostfactor(String hostname);


    public static native String ls_gethostmodel(String hostname);

    // NOTE: Not in liblsf
    //public static native IntByReference ls_gethostrespriority(String hostname);


    public static native int ls_lockhost(NativeLong duration);


    public static native int ls_unlockhost();


    public static native int ls_limcontrol(String hostname, int opCode);

    public static native void ls_remtty(int ind, int enableIntSus);

    public static native void ls_loctty(int ind);


    public static native String ls_sysmsg();


    public static native void ls_perror(String usrMsg);


    public static native lsConf.ByReference ls_getconf(String string);

    public static native void ls_freeconf(lsConf lsConf1);

    public static native sharedConf.ByReference ls_readshared(String string1);

    public static native clusterConf.ByReference ls_readcluster(String string1, lsInfo lsInfo1);

    public static native clusterConf.ByReference ls_readcluster_ex(String string1, lsInfo lsInfo1, int int1);


    public static native int _ls_initdebug(String appName);

    public static native void ls_syslog(int level, String fmt, Pointer args);

    public static native void ls_errlog(Pointer fp, String fmt, Pointer args);

    // NOTE: va_list is too compiler specific.  Skipping this function.
    //public static native void  ls_verrlog (Pointer fp, String fmt, va_list ap);

    public static native int ls_fdbusy(int fd);


    public static native String ls_getmnthost(String fn);

    public static native int ls_servavail(int int1, int int2);

    public static native int ls_getpriority(IntByReference priority);

    public static native int ls_setpriority(int newPriority);

    public static native void ls_ruunix2lsf(rusage rusage, lsfRusage lsfRusage);

    public static native void ls_rulsf2unix(lsfRusage lsfRusage, rusage rusage);

    public static native void cleanLsfRusage(lsfRusage lsfRusage1);

    public static native void cleanRusage(rusage rusage1);


    // NOTE: Not in liblsf
    //public static native int getBEtime(String string1, byte byte1, NativeLongByReference long1);


    public static native int ls_postevent(int int1, String string1, Pointer stringArray1, int int2);

    public static native int ls_postmultievent(int int1, String string1, Pointer stringArray1, int int2, int int3);

    public static class extResInfo extends Structure {
        public static class ByReference extends extResInfo implements Structure.ByReference {}
        public static class ByValue extends extResInfo implements Structure.ByValue {}
        public extResInfo() {}
        public extResInfo(Pointer p) { super(p); read(); }

        public String name;
        public String type;
        public String interval;
        public String increasing;
        public String des;
    }




    // NOTE: Not in liblsf
    //public static native int lim_vcl_get_eres_version();

    // NOTE: Not in liblsf
    //public static native extResInfo.ByReference lim_vcl_get_eres_def(String string1);

    // NOTE: Not in liblsf
    //public static native String lim_vcl_get_eres_loc(String string1);

    // NOTE: Not in liblsf
    //public static native String lim_vcl_get_eres_val(String string1);


    public static int isspace(byte c) {
        return ((c == 0x20 || c == 0x09 || c == 0x0a || c == 0x0b || c == 0x0c || c == 0x0d) ? 8 : 0);
    }

    public static final int LSF_VERSION = LSF_XDR_VERSION7_0_EP6;
    public static final String LSF_CURRENT_VERSION = "7.06";


    public static final String LSF_PRODUCT_COPYRIGHT_STR = "Copyright 1992-2009 Platform Computing Corp.";


    public static final String LSF_NAME_STR = "Platform LSF";
    public static final String LSF_IDENTIFIER_STR = "";
    public static final String LSF_PRODUCT_NAME_STR = LSF_NAME_STR + LSF_IDENTIFIER_STR;


    public static final String LSF_PRODUCT_COMMENT_STR = "";


    public static final String LSF_PRODUCT_BUILD_STR = "";


    public static final String LSF_PRODUCT_BUILD_DATE_STR = "";


    public static final int LSF_PRODUCT_MAJOR_VERSION = 7;
    public static final int LSF_PRODUCT_MINOR_VERSION = 0;
    public static final int LSF_PRODUCT_MAINTAIN_VERSION = 6;

    public static final String LSF_PRODUCT_MAJOR_VERSION_STR = "7";
    public static final String LSF_PRODUCT_MINOR_VERSION_STR = "0";
    public static final String LSF_PRODUCT_MAINTAIN_VERSION_STR = "6";

    public static final String LSF_PRODUCT_VERSION_STR = LSF_PRODUCT_MAJOR_VERSION_STR + "." + LSF_PRODUCT_MINOR_VERSION_STR + "." + LSF_PRODUCT_MAINTAIN_VERSION_STR;
    public static final String LSF_FILE_VERSION_STR = LSF_PRODUCT_MAJOR_VERSION_STR + "." + LSF_PRODUCT_MINOR_VERSION_STR + "." + LSF_PRODUCT_MAINTAIN_VERSION_STR;


    public static final String _VERSION_STR_LSID_ = "Platform LSF HPC 7";
    public static final String _LSID_VERSION_ = (_VERSION_STR_LSID_ + " Update " + _MINOR_STR_ + ", " + _DATE_STR_ + "\nCopyright 1992-2009 Platform Computing Corporation\n");


    /* Removing since the ls_nio functions which use fd_set, etc. are not in liblsf.

    public static final int NIO_STDIN_ON = 0x01;
    public static final int NIO_STDIN_OFF = 0x02;
    public static final int NIO_TAGSTDOUT_ON = 0x03;
    public static final int NIO_TAGSTDOUT_OFF = 0x04;

    public static final int NIO_TASK_STDINON = 0x01;
    public static final int NIO_TASK_STDINOFF = 0x02;
    public static final int NIO_TASK_ALL = 0x03;
    public static final int NIO_TASK_CONNECTED = 0x04;

    public static interface nioType {
          public static final int NIO_STATUS = 0;
          public static final int NIO_STDOUT = 1;
          public static final int NIO_EOF = 2;
          public static final int NIO_IOERR = 3;
          public static final int NIO_REQUEUE = 4;
          public static final int NIO_STDERR = 5;
    }



    public static class nioEvent extends Structure {
        public static class ByReference extends nioEvent implements Structure.ByReference {}
        public static class ByValue extends nioEvent implements Structure.ByValue {}
        public nioEvent() {}
        public nioEvent(Pointer p) { super(p); read(); }

        public int tid;
        public *//*nioType*//* int type;
        public int status;
    }



    public static class nioInfo extends Structure {
        public static class ByReference extends nioInfo implements Structure.ByReference {}
        public static class ByValue extends nioInfo implements Structure.ByValue {}
        public nioInfo() {}
        public nioInfo(Pointer p) { super(p); read(); }

        public int num;
        public Pointer / * nioEvent.ByReference * / ioTask;
    }


    public static final int FD_SETSIZE = 64;

    public static class fd_set extends Structure {
        public static class ByReference extends fd_set implements Structure.ByReference {}
        public static class ByValue extends fd_set implements Structure.ByValue {}
        public fd_set() {}
        public fd_set(Pointer p) { super(p); read(); }

        public int count;
        public int[] fd = new int[FD_SETSIZE];
    }
    */

    public static native int ls_initdebug(String appName);

    // NOTE: Not in liblsf
    //public static native int ls_nioinit(int sock);

    // NOTE: Not in liblsf
    //public static native int ls_nioselect(int int1, fd_set fd_set1, fd_set fd_set2, fd_set fd_set3, Pointer nioInfoArray1, timeval timeval1);

    // NOTE: Not in liblsf
    //public static native int ls_nioctl(int int1, int int2);

    // NOTE: Not in liblsf
    //public static native int ls_nionewtask(int int1, int int2);

    // NOTE: Not in liblsf
    //public static native int ls_nioremovetask(int int1);

    // NOTE: Not in liblsf
    //public static native int ls_niowrite(String string1, int int1);

    // NOTE: Not in liblsf
    //public static native int ls_nioclose();

    // NOTE: Not in liblsf
    //public static native int ls_nioread(int int1, String string1, int int2);

    // NOTE: Not in liblsf
    //public static native int ls_niotasks(int int1, IntByReference int2, int int3);

    // NOTE: Not in liblsf
    //public static native int ls_niostatus(int int1, IntByReference int2, rusage rusage1);

    // NOTE: Not in liblsf
    //public static native int ls_niokill(int int1);

    // NOTE: Not in liblsf
    //public static native int ls_niosetdebug(int int2);

    // NOTE: Not in liblsf
    //public static native int ls_niodump(int int1, int int2, int int3, String string1);


    public int lsf_res_version;


    public static native int ls_initrex(int a, int b);

    public static int ls_init(int a, int b) {
        return ls_initrex(a, b);
    }


    public static native int ls_donerex();

    public static native int ls_niossync(int int1);


    public static native int ls_setstdin(int on, IntByReference rpidlist, int len);


    public static native int ls_getstdin(int on, IntByReference rpidlist, int maxlen);

    public static native int ls_setstdout(int on, String format);


    public static native int ls_stdinmode(int onoff);


    public static native int ls_stoprex();


    public static native int ls_chdir(String string1, String string2);


    public static native int ls_connect(String string1);


    public static native int ls_rkill(int int1, int int2);


    public static native int ls_rsetenv(String host, Pointer env);

    public static native int ls_rsetenv_async(String host, Pointer env);


    public static native int ls_rescontrol(String host, int opcode, int options);


    public static native lsfAcctRec.ByReference ls_getacctrec(Pointer pointer1, IntByReference int1);

    public static native int ls_putacctrec(Pointer pointer1, lsfAcctRec lsfAcctRec1);


    // NOTE: No idea what resLogRecord is.
    //public static native resLogRecord.ByReference ls_readrexlog (Pointer );


    public static native int ls_rexecv(String string1, Pointer string2, int int1);


    public static native int ls_rexecve(String string1, Pointer stringArray1, int int1, Pointer stringArray2);

    public static native int ls_rexecv2(String string1, Pointer stringArray1, int int1);

    public static native int ls_startserver(String string1, Pointer stringArray1, int int1);


    public static native int ls_rtask(String string1, Pointer stringArray1, int int1);


    public static native int ls_rtaske(String string1, Pointer stringArray1, int int1, Pointer stringArray2);

    public static native int ls_rtask2(String string1, Pointer stringArray1, int int1, Pointer stringArray2);


    public static native int ls_rwait(IntByReference int1, int int2, rusage rusage1);


    public static native int ls_rwaittid(int int1, IntByReference int2, int int3, rusage rusage1);


    public static native int ls_conntaskport(int tid);


    public static native int ls_ropen(String host, String fn, int flags, int mode);


    public static native int ls_rclose(int rfd);


    public static native int ls_rwrite(int rfd, String buf, int len);


    public static native int ls_rread(int rfd, String buf, int len);


    public static native NativeLong ls_rlseek(int rfd, NativeLong offset, int whence);


    public static native int ls_runlink(String host, String fn);

    public static native int ls_rfstat(int rfd, Pointer buf);

    public static native int ls_rstat(String host, String fn, Pointer buf);


    public static native String ls_rgetmnthost(String host, String fn);


    public static native int ls_rfcontrol(int command, int arg);


    public static native int ls_rfterminate(String host);
}
