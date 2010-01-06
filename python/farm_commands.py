#!/usr/bin/env python

import os
import sys
import subprocess
import re

#justPrintCommands = False

def cmd(cmd_str_from_user, farm_queue=False, output_head=None, just_print_commands=False, outputFile = None, waitID = None, jobName = None, die_on_fail = False):
    """if farm_queue != False, submits to queue, other
die_on_fail_msg: if != None, die on command failure (non-zero return) and show die_on_fail_msg"""

    if farm_queue:
        if outputFile <> None:
            farm_stdout = outputFile
        elif output_head <> None:
            farm_stdout = output_head+".stdout"
        else:
            #farm_stdout = None
            farm_stdout = "%J.lsf.output"
            
        cmd_str = "bsub -q "+farm_queue 
        if farm_stdout <> None:
            cmd_str += " -o " + farm_stdout

        if waitID <> None:
            cmd_str += " -w \"ended(%s)\"" % (str(waitID))

        if jobName <> None:
            cmd_str += " -J %s" % (jobName)

        cmd_str += " '"+cmd_str_from_user + "'"
        
        print ">>> Farming via "+cmd_str
    else:
        cmd_str = cmd_str_from_user
        print ">>> Executing "+cmd_str

    if just_print_commands or (globals().has_key("justPrintCommands") and globals().justPrintCommands):
        return -1
    elif farm_queue:
        result = subprocess.Popen([cmd_str, ""], shell=True, stdout=subprocess.PIPE).communicate()[0]
        p = re.compile('Job <(\d+)> is submitted to queue')
        jobid = p.match(result).group(1)
        return jobid        
    else:
        # Actually execute the command if we're not just in debugging output mode
        status = os.system(cmd_str)
        if not farm_queue:
            print "<<< Exit code:", status,"\n"
        if die_on_fail != None and status != 0:
            print "### Failed with exit code "+str(status)+" while executing command "+cmd_str_from_user
            sys.exit()
        return int(status)
