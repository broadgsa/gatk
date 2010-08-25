# the linear index class implementation
$LOAD_PATH << File.dirname(__FILE__)
require "Index.rb"

class LinearIndex < Index
	attr_accessor :nSeq, :sequences
	def initialize(file)
		super(file)
		@nSeq = @file.readInt()
		@sequences = Array.new()
		@nSeq.times {|index|
			@sequences.push(LISeqEntry.new(@file))
		}
	end
	
	def diff(otherIndex)
		diffHeader(otherIndex)
		if (otherIndex.type != @type)
			print "Indexes are not the same type (#{otherIndex.type} != #{@type})\n"
			return false
		end
		ret = false
		notInOther = @sequences.reject {|item|
								return true if !otherIndex.sequences.include?(item) 
						 		item.diff(otherIndex.sequences[otherIndex.sequences.index(item)])
						 	}
		notInOther.pretty_print
	end	
		
			
end

class LISeqEntry
	def initialize(file)
		@contig = file.readString()
		@binWidth = file.readInt()
 		@binCount = file.readInt()
 		@longestFeature = file.readInt()
 		@maxBin = file.readInt()
 		@totalBin = file.readInt()
 		@startPositions = Array.new()
 		@binCount.times { |index|
 			@startPositions.push(file.readLong())
 		}
 		@finalPos = file.readLong()
 	end
 	
 	# print a summary of the index characteristics
	def diff(otherLISeqEntry)
		self.instance_variables.each { |var|
			next if "#{var}" == "@file" or  "#{var}" == "@sequences"
			puts "Other LISeqEntry doesn't define #{var}" if !(otherLISeqEntry.instance_variable_defined?(var))
			one = (self.instance_variable_get(var)).to_s
			two = (otherLISeqEntry.instance_variable_get(var)).to_s
			puts "otherLISeqEntry: type #{var} not equal, #{one} != #{two}" if one != two
		}
	end
end

