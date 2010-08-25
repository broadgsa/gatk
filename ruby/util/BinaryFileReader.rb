# a ruby class for reading in binary files; really this just adds some conv. methods like readInt(), readLong(), etc.
class BinaryFileReader
	# constructor
	def initialize(fileName)
		@file = File.open(fileName,"r")
	end

	# read and return an int (4 byte, signed, machine based endian)
	def readInt()
		(@file.sysread(4)).unpack("i")[0]
	end
	
	# read and return an int (4 byte, unsigned, machine based endian)
	def readUInt()
		(@file.sysread(4)).unpack("L")[0]
	end
	
	# read and return an long (8 byte, signed, machine based endian)
	def readLong()
		(@file.sysread(8)).unpack("q")[0]
	end
	
	# read and return a set number of bytes as a string
	def readBytes(count)
		(@file.sysread(count)).to_s
	end
	
	# read and return a null terminated string
	def readString()
		buffer = []
		ch = @file.sysread(1)
		while (ch != "\0")
			buffer.push(ch)
			ch = @file.sysread(1)
		end
		buffer.to_s
	end

	# close the file
	def close()
		@file.close()
	end
end

