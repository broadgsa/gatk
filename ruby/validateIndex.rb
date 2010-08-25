# this ruby files validates a linear index
# set the include path to include the current directory
$LOAD_PATH << File.dirname(__FILE__)

# require a couple of files
require "index/Index.rb"
require "index/LinearIndex.rb"
require "index/IntervalIndex.rb"
require "optparse"
require "yaml"

# This hash will hold all of the options
# parsed from the command-line by
# OptionParser.
options = {}

optparse = OptionParser.new do|opts|
  # Set a banner, displayed at the top
  # of the help screen.
  opts.banner = "Usage: ruby validateIndex.rb [options] file1 file2 ..."

  # Define the options, and what they do
  options[:index] = [] if options[:index] == nil
  opts.on( '-i', '--index INDEX (REQUIRED)', 'Specify the index.  Multiple are allowed' ) do |file| options[:index].push(file) end

  options[:verbose] = false
  opts.on( '-v', '--verbose', 'Output more information' ) do options[:verbose] = true end

  options[:validate] = false
  opts.on( '-c', '--check', 'Check (Validate) the index(es) passed in as parameters' ) do options[:check] = true end
  
  options[:diff] = false
  opts.on( '-d', '--diff', 'Diff two indexes' ) do options[:diff] = true end
  
  options[:print] = false
  opts.on( '-p', '--print', 'Print all of the information about the file' ) do options[:print] = true end
	
  # This displays the help screen, all programs are
  # assumed to have this option.
  opts.on_tail( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
end
# parse the command line
optparse.parse!

#Now raise an exception if we have not found a host option
if options[:index].size == 0
	puts "you must at least specify an index file!"
	puts optparse
end

# a function to load an index
def loadIndex(file)
	indexTry = Index.new(file)
	indexTry.close()
	if (indexTry.type == 1)
		puts "Linear index..."
		index = LinearIndex.new(file)
	else
		puts "Interval index..."
		index = IntervalIndex.new(file)
	end
	index
end

#################### Control Block ####################

# load all of the indexes
indexes = []
options[:index].each {|indexFile|
	indexes.push(loadIndex(indexFile))
}

# switch on the flags supplied
if (options[:diff])
	if (options[:index].size != 2)
		print "Unable to diff indexes if you don't supply two and only two indexes\n";
		exit(1)
	else
		indexes[0].diff(indexes[1])
	end
elsif (options[:validate])
	indexes.each {|index| index.validate() }
elsif (options[:print])
	indexes.each {|index| puts YAML::dump( index ) }
end




# if they specified validate
if (options[:check])
	options[:index].each {|index|
		idx = Index.new(index).validate()
	}
end




exit
