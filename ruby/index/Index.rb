require './util/BinaryFileReader'
# this base class for all index types (at least linear and tree)
class Index
	attr_reader :type, :headerVersion, :fileName, :fileSize, :t5, :md5, :flags 

	# construct, given a file name
	def initialize(fileName)
		@inputFile = fileName 
		@file = BinaryFileReader.new(fileName)
		magic = @file.readBytes(4)
		if (magic != "TIDX")
			print "#{@inputFile}: !! Magic number is not what we expected, TIDX, instead we saw #{magic} !!\n"
			exit(1)
		end
		@type = @file.readInt()
		@headerVersion = @file.readInt()
		@fileName = @file.readString()
		@fileSize = @file.readLong()
		@ts = @file.readLong()
		@md5 = @file.readString()
		@flags = @file.readUInt()
		@seqDict = readSeqDictionary() if (@flags == 32768)
		@propCount = readPropertyDictionary() if (@headerVersion >= 3) 
	end
	
	def validate()
		 f = Proc.new{ print "#{@inputFile}:\t\terror: invalid type, we saw #{@type} but expected [1-2]\n"; return false} if @type < 1 or @type > 2
		 f = Proc.new{ print "#{@inputFile}:\t\terror: invalid header version, we saw #{@headerVersion} but expected [1-3]\n"; return false} if @headerVersion < 1 or @headerVersion > 3
		 f = Proc.new{ print "#{@inputFile}:\t\twarning: on fileName, we saw '#{@fileName}' but expected actual text\n"; return false} if @fileName == ""
		 f = Proc.new{ print "#{@inputFile}:\t\twarning: on TS, we saw '#{@ts}' but expected actual text\n"; return false} if @ts == ""
		 f = Proc.new{ print "#{@inputFile}:\t\twarning: on md5, we saw '#{@md5}' but expected actual text\n"; return false} if @md5 == ""
		 f.call if f != nil
		 return true
	end 
	
	# diff two headers
	def diffHeader(otherIndex)
		self.instance_variables.each { |var|
			next if "#{var}" == "@file" or  "#{var}" == "@sequences"
			puts "Other header doesn't define #{var}" if !(otherIndex.instance_variable_defined?(var))
			one = (self.instance_variable_get(var)).to_s
			two = (otherIndex.instance_variable_get(var)).to_s
			puts "type #{var} not equal, #{one} != #{two}" if one != two
		}
	end
	
	# read the sequence dictionary, assuming we have one
	def readSeqDictionary()
		sequences = []
		count = @file.readInt()
		count.times {|index|
			sequences.add(@file.readString())
			@file.readInt() # drop the sizes for now
		}
		sequences # return sequences
	end
	
	# read the sequence dictionary, assuming we have one
	def readPropertyDictionary()
		sequences = {}
		count = @file.readInt()
		count.times {|index|
			sequences.put(@file.readString(),@file.readString()) }
		sequences # return sequences
	end
	
	def close()
		@file.close()
	end
end