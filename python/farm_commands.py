#!/usr/bin/env python

import os
import subprocess
import re

#justPrintCommands = False

def cmd(cmd_str_from_user, farm_queue=False, output_head=None, just_print_commands=False, outputFile = None, waitID = None):
    # if farm_queue is non-False, submits to queue, other

    if farm_queue:
        if outputFile <> None:
            farm_stdout = outputFile
        elif output_head <> None:
            farm_stdout = output_head+".stdout"
        else:
            farm_stdout = None
            
        cmd_str = "bsub -q "+farm_queue 
        if farm_stdout <> None:
            cmd_str += " -o " + farm_stdout

        if waitID <> None:
            cmd_str += " -w \"ended(%s)\"" % (str(waitID))

        cmd_str += " \""+cmd_str_from_user + "\""
        
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
        return int(status)
        
