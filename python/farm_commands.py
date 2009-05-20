#!/usr/bin/env python

import os
#justPrintCommands = False

def cmd(cmd_str, farm_queue=False, output_head=None, just_print_commands=False, outputFile = None):
    # if farm_queue is non-False, submits to queue, other

    if farm_queue:
        if outputFile <> None:
            farm_stdout = outputFile
        elif output_head <> None:
            farm_stdout = output_head+".stdout"
        else:
            farm_stdout = None
            
        cmd_str = "bsub -q "+farm_queue+" "+cmd_str 
        if farm_stdout <> None:
            cmd_str += " -o " + farm_stdout
        
        print ">>> Farming via "+cmd_str
    else:
        print ">>> Executing "+cmd_str

    if just_print_commands or (globals().has_key("justPrintCommands") and globals().justPrintCommands):
        return -1
    else:
        # Actually execute the command if we're not just in debugging output mode
        status = os.system(cmd_str)
        if not farm_queue:
            print "<<< Exit code:", status,"\n"
        return status
        
