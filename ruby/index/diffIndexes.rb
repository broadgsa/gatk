# this ruby files takes two indexes (of the same type) and diff's them.  If they're different types,
# it'll stop after the header

# a function to exit, printing a message
def exitWithError(message) 
	puts ###########################################
	puts message
	puts ###########################################
	exit(1)
end


if (ARGV.size != 2)
	exitWithError("We take two files as input, try again!")	
end

# open the indexes
index1 = File.new(ARGV[0])
index2 = File.new(ARGV[1])

# a helper function for comparing values
def compValues(file1, file2, byteCount, type)
	index1Type = file1.sysread(byteCount)
	index2Type = file2.sysread(byteCount)
	if (index1Type != index2Type)
		print "#{type}, index1 (#{index1Type}) != index2 (#{index2Type})\n"
	else 
		print "#{type}, index1 (#{index1Type}) == index2 (#{index2Type})\n"
	end
end

# a helper function for comparing values
def compInts(file1, file2, byteCount, type)
	if (byteCount > 4) 
		upack = "q"
	else 
		upack = "i"
	end
	index1Type = file1.sysread(byteCount).unpack(upack)
	index2Type = file2.sysread(byteCount).unpack(upack)
	if (index1Type != index2Type)
		print "#{type}, index1 (#{index1Type}) != index2 (#{index2Type})\n"
	else 
		print "#{type}, index1 (#{index1Type}) == index2 (#{index2Type})\n"
	end
end


# a helper function for reading strings
def readString(index)
	buffer = []
	ch = index.sysread(1)
	while (ch != "\0")
		buffer.push(ch)
		ch = index.sysread(1)
	end
	buffer.to_s
end

# validate the magic number from both
exitWithError("Magic number not valid for index 1") if index1.sysread(4) != "TIDX"
exitWithError("Magic number not valid for index 2") if index2.sysread(4) != "TIDX"

# validate the types
compInts(index1,index2,4,"types")

# validate the versions
v1 = index1.sysread(4).unpack("i")
v2 = index2.sysread(4).unpack("i")
if (v1 != v2)
	print "version, index1 (#{v1}) != index2 (#{v2})\n"
else 
	print "version, index1 (#{v1}) == index2 (#{v2})\n"
end

# validate the filenames
fl1 = readString(index1)
fl2 = readString(index2)
if (fl1 != fl2)
	print "filename, index1 (#{fl1}) != index2 (#{fl2})\n"
else 
	print "filename, index1 (#{fl1}) == index2 (#{fl2})\n"
end
# validate the sizes
compInts(index1,index2,8,"sizes")


# validate the T5?
compValues(index1,index2,8,"T5")

# validate the MD5 - just a byte, we don't write the MD5 sums in yet
# validate the filenames
fl1 = readString(index1)
fl2 = readString(index2)
if (fl1 != fl2)
	print "md5, index1 (#{fl1}) != index2 (#{fl2})\n"
else 
	print "md5, index1 (#{fl1}) == index2 (#{fl2})\n"
end

# validate the flags
index1Flags = (index1.sysread(4)).unpack("L")
index2Flags = (index2.sysread(4)).unpack("L")
if (index1Flags != index2Flags)
		print "Flags are different, index1 = #{index1Flags[0]}, index2 = #{index2Flags[0]}\n"
end
		
def readSeqDictionary(file)
	puts "reading seq dict"
	sequences = []
	count = (file.sysread(4)).unpack("i")	
	puts count
	count[0].times {|index|
		sequences.add(readString(file).to_s)
		file.sysread(4) # drop the sizes for now
	}
	sequences # return sequences
end

if (index1Flags[0] == 32768)
	puts readSeqDictionary(index1)
elsif (index2Flags[0] == 32768)
	puts readSeqDictionary(index2)
end

# bump off the prop dictionary
index1.sysread(4) if (v1[0] == 3) 
index2.sysread(4) if (v2[0] == 3) 

def readSeqEntry(i1)
    puts "--------------------------------------------" 
	puts "Contig --> #{readString(i1).to_s}"
	print "bin width         = #{(i1.sysread(4)).unpack("i")}\n"
	binCount = (i1.sysread(4)).unpack("i")[0]
	print "number of bins    = #{binCount}\n"		
	print "longest feature   = #{(i1.sysread(4)).unpack("i")}\n"
	print "max bin size      = #{(i1.sysread(4)).unpack("i")}\n"
	print "total bin size    = #{(i1.sysread(4)).unpack("i")}\n"
	lastStartPos = -1
	binCount.times { |index|
		startPos = (i1.sysread(8)).unpack("q")[0]	
		#if (startPos < lastStartPos)
			puts "bin at index #{index}, this start = #{startPos}, last start = #{lastStartPos}"
		#end
		lastStartPos = startPos
	}
		
end

def compContigLinear(i1, i2)
	seqC1 = (i1.sysread(4)).unpack("L")
	seqC2 = (i2.sysread(4)).unpack("L")	
	print "seq count 1 = #{seqC1}, count 2 = #{seqC2}\n"
	puts "\nentries for index 1"
	seqC1[0].times { |index|
		readSeqEntry(i1)			
	} 
	puts "\nentries for index 2"	
	seqC2[0].times { |index|
		readSeqEntry(i2)	
	} 
end

compContigLinear(index1,index2)		
	

print "Done!\n"
# close the files
index1.close()
index2.close()