# the implementation of the interval index class
$LOAD_PATH << File.dirname(__FILE__)
require "Index.rb"

class IntervalIndex < Index
	def initialize(file)
		super(file)
		@nSeq = @file.readInt()
		@sequences = Array.new()
		@nSeq.times {|index|
			@sequences.push(TISeqEntry.new(@file))
		}
	end
	
	def diff(otherIndex)
		diffHeader(otherIndex)
		if (otherIndex.type != @type)
			print "Indexes are not the same type (#{otherIndex.type} != #{@type})\n"
			return false
		end
		ret = false
	end	
end

class TISeqEntry
	def initialize(file)
		@contig = file.readString()
		@binCount = file.readInt()
 		@startPositions = Array.new()
 		@endPositions = Array.new()
 		@positions = Array.new()
 		@sizes = Array.new()
 		@binCount.times { |index|
 			@startPositions.push(file.readInt())
 			@endPositions.push(file.readInt())
 			@positions.push(file.readLong())
 			@sizes.push(file.readInt())
 		} 		
 	end
end

