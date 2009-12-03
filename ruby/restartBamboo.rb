#!/util/bin/ruby
# #########################################################
# 
# Dec 3rd, 2009 
# Aaron
#
# This script will restart Bamboo, if it's running; or if 
# it's left a pid file without stopping the server
# 
# It has two arguments:
# - the log file to store any results to (required)
# - DRYRUN, do not actually execute commands (optional)
# #########################################################

# require file utils, lets us change directory
require 'fileutils'

# get the time and date, and take off any endlines
td = `date`.chomp!

# what machine we have to run this on?
machine = "gsa2"

# check the command line arguments
if (ARGV.size < 1 || ARGV.size > 2)
  log_file.puts "#{td}: restartBamboo.rb <logFileLocation> <DRYRUN>\n"
  log_file.puts "#{td}: logFileLocation: where to put the output"
  log_file.puts "#{td}: DRYRUN (optional): do not execute any Bamboo commands, only echo them\nexiting..."
  log_file.close()
  exit(1)
end

# bamboo location, shell script, and PID file
bamboo_location = "/humgen/gsa-scr1/gsa-engineering/bamboo/Bamboo"
bamboo_script   = "bamboo.sh"
pid_file_loc    = "bamboo.pid"

# where we want to log the results, taken from the command line
log_file = File.open(ARGV[0],"a")

if (`uname -n`.chomp! != machine)
  log_file.puts "#{td}: You can only run this script on machine #{machine}...exiting!"
  log_file.close()
  exit(1)
end

echo = ""
# are we a dryrun?
if (ARGV[1] == "DRYRUN")
  echo = "echo "
elsif (ARGV.size == 2)
  log_file.puts "#{td}: The second parameter must be DRYRUN or nothing!"
  log_file.close()
  exit(1)
end

# first check if the bamboo dir exists
if (!File.exists?(bamboo_location)) 
  log_file.puts "#{td}: Bamboo location: #{bamboo_location} does not exist"
  endlog_file.close()
  exit(1)
end

# we have to cd to the directory, since bamboo uses some relative paths (stinks!)
FileUtils.cd(bamboo_location)

# output our current location
currentDir = FileUtils.pwd()
log_file.puts "current dir = #{currentDir}"

# check to see if the bamboo location has a pid file
if (File.exists?(bamboo_location + File::SEPARATOR + "bamboo.pid")) 
  log_file.puts "#{td}: PID file exists! Is bamboo running with that PID..."
  
  # get the pid
  pid = `cat #{pid_file_loc}`.to_i
  log_file.puts "#{td}: got a PID value of #{pid}"

  # check for that pid
  if (`ps --no-headers #{pid}` == "")
    log_file.puts "#{td}: unable to find process #{pid}"
    log_file.puts "#{td}: trying to remove pid file: #{pid_file_loc}"
    FileUtils.rm(pid_file_loc)
  else
    log_file.puts "#{td}: found process #{pid}, attempting to shut down..."
    shutdown_result = `#{echo}./#{bamboo_script} stop`
    log_file.puts "#{td}: shutdown result:\n#{shutdown_result}"
  end
end

# regardless of what happened above, restart the server
startup_result = `#{echo}./#{bamboo_script} start`
log_file.puts "#{td}: startup result:\n#{startup_result}"

# close the log file 
log_file.close()
