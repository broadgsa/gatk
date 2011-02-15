local ampFile = arg[1]
local samFile = arg[2]
local headerFile = io.open(arg[2]:match("(.*).sam")..".header.sam", "w")
local bodyFile = io.open(arg[2]:match("(.*).sam")..".body.sam", "w")


-- These sizes are hardcoded for hg19, but future versions of this script should optionally take a .fai file to build this table.
chrlength = {}
chrlength["1"] =249250621
chrlength["2"] =243199373
chrlength["3"] =198022430
chrlength["4"] =191154276
chrlength["5"] =180915260
chrlength["6"] =171115067
chrlength["7"] =159138663
chrlength["8"] =146364022
chrlength["9"] =141213431
chrlength["10"]=135534747
chrlength["11"]=135006516
chrlength["12"]=133851895
chrlength["13"]=115169878
chrlength["14"]=107349540
chrlength["15"]=102531392
chrlength["16"]=90354753
chrlength["17"]=81195210
chrlength["18"]=78077248
chrlength["19"]=59128983
chrlength["20"]=63025520
chrlength["21"]=48129895
chrlength["22"]=51304566
chrlength["X"] =155270560
chrlength["Y"] =59373566
chrlength["MT"]=16569

local header = {}
local reads = {}

-- reads the amplicon file (global var) and returns an amplicon hash table indexed by amplicon name
local function readAmpliconFile()
  local amplicons = {}
  for l in io.lines(ampFile) do
    local chr, startPos, endPos, ampName = l:match("([%w]+)%s+(%d+)%s+(%d+)%s+([%w%p_]+)")
    amplicons[ampName] = {chr=chr, startPos=tonumber(startPos), endPos=tonumber(endPos)}
  end
  return amplicons
end

-- updates the global header table with the entries
local function processSamHeaderLine(l, amplicons)
  if l:sub(1,3) == "@HD" then header.hd = l
  elseif l:sub(1,3) == "@CO" then header.co = l
  else
    if not header.sq then header.sq = {} end
    local ampName = l:match("@SQ%s+SN:ps%d+_([%w%p_]+).*")
    local chr = amplicons[ampName].chr
    header.sq[chr] = chrlength[chr]
  end
end

local function printHeader()
  if header.hd then headerFile:write(header.hd.."\n") else headerFile:write("@HD VN:PB_v0.1") end
  for chr, len in pairs(header.sq) do
    headerFile:write("@SQ\tSN:"..chr.."\tLN:"..len.."\n")
  end
  if header.co then headerFile:write(header.co.."\n") end
end

local function printReads()
  for _, v in ipairs(reads) do bodyFile:write(v.."\n") end
end


local amplicons = readAmpliconFile()  -- amplicons indexed by name

--for i,v in pairs(amplicons) do print("'"..i.."'") end

for l in io.lines(samFile) do
  if l:sub(1,1) == "@" then processSamHeaderLine(l, amplicons)
  else
    local before, amp, mapStart, after = l:match("(%d+%s+%d+%s+)ps%d+_([%w%p_]+)%s+(%d+)(.*)")
    table.insert(reads, before..amplicons[amp].chr.."\t"..amplicons[amp].startPos + mapStart - 1 ..after)
  end
end

printHeader()
printReads()
