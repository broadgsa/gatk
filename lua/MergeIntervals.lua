-- Merges all intervals into a sorted list of unoverlapping intervals
-- 
-- ps: this utility is similar to -im on GATK, but for queue scripts we need the intervals
--     file fixed before we can run Unified Genotyper because of scatter-gather.
--
--	Usage:
--		lua MergeIntervals.lua raw.intervals > merged.intervals
--

assert(table.getn(arg) > 0, "\n\nMissing input file\n\nUsage:\n\tlua MergeIntervals.lua raw.intervals > merged.intervals\n\n")


local function comp(a, b)
   if tonumber(a) and tonumber(b) then return tonumber(a) < tonumber(b) end
   return a<b
end


-- read file and create chromosome tables with raw intervals marked as "x"
local intervals = {}
local intervalKeys = {}
local sortedChrs = {}
for l in io.lines(arg[1]) do
	local chr, a, b = l:match("([%a%d]+):(%d+)-(%d+)")
	a,b = tonumber(a), tonumber(b)
	if not intervals[chr] then 
		intervals[chr] = {}
		intervals[chr].min = a
		intervals[chr].max = b
		intervalKeys[chr] = {}
		table.insert(sortedChrs, chr)
	end
	for i=a,b do intervals[chr][i] = "x" end
	intervals[chr].min = math.min(intervals[chr].min, a)
	intervals[chr].max = math.max(intervals[chr].max, b)
	table.insert(intervalKeys[chr], a)
end

local intervalSize = {}
local intervalsIndex = {}
for c, chrTable in pairs(intervals) do
	intervalSize[c] = {}
	intervalsIndex[c] = {}
	local intervalStart = -1
	local lastPositionRead = -1
	table.sort(intervalKeys[c])
	for _,startKey in ipairs(intervalKeys[c]) do
		local i = startKey
		while chrTable[i] and i > lastPositionRead do
			if (intervalStart < 0 or i ~= (intervalStart + intervalSize[c][intervalStart] + 1)) then 
				intervalStart = i
				intervalSize[c][intervalStart] = 0 
				table.insert(intervalsIndex[c], intervalStart)
			else 
				intervalSize[c][intervalStart] = intervalSize[c][intervalStart] + 1
			end
			lastPositionRead = i
			i = i + 1
		end
	end
end

table.sort(sortedChrs, comp)
for _,c in pairs(sortedChrs) do
	table.sort(intervalsIndex[c])
	for _,intervalStart in ipairs(intervalsIndex[c]) do
		if intervalSize[c][intervalStart] > 0 then 
			print(c..":"..intervalStart.."-"..intervalStart + intervalSize[c][intervalStart])
		end
	end
end
