------------------------------------------------------------------------------------------------------------------------
-- Creates a new interval table
--
-- Return values:
-- 1: Interval table
------------------------------------------------------------------------------------------------------------------------

function newInterval(chr, start, finish, strand, info)
  return {chr= chr, start=tonumber(start), finish=tonumber(finish), strand=strand, info=info}
end


------------------------------------------------------------------------------------------------------------------------
-- Parses a line from an interval list file (not a header line!)
--
-- Return values:
-- 1: chromosome
-- 2: interval start
-- 3: interval end
-- 4: strand (+/-)
-- 5: info field
------------------------------------------------------------------------------------------------------------------------
local function parseIntervalLine(l)
  return l:match("(%w+)%s+(%d+)%s+(%d+)%s+([%+%-])%s+(.*)")
end

------------------------------------------------------------------------------------------------------------------------
-- Reads an interval list file into a table
--
-- Return values:
-- 1: table of intervals
-- 2: header string
------------------------------------------------------------------------------------------------------------------------
local function readIntervalList(filename)
  t = {}
  header = ""
  for l in io.lines(filename) do
    if l:sub(1,1) == "@" then header = header .. l .."\n"
    else
        local chr, start, finish, strand, info = parseIntervalLine(l)
        table.insert(t, newInterval(chr, start, finish, strand, info))
    end
  end
  return t, header
end

------------------------------------------------------------------------------------------------------------------------
-- Checks if two intervals have the same chromosome, start and end.
--
-- Return values:
-- 1: true/false
------------------------------------------------------------------------------------------------------------------------
local function isSameInterval (i1, i2)
  return i1.chr == i2.chr and i1.start == i2.start and i1.finish == i2.finish
end

------------------------------------------------------------------------------------------------------------------------
-- Checks if the line from an interval list file is a header line
--
-- Return values:
-- 1: true/false
------------------------------------------------------------------------------------------------------------------------
local function isIntervalHeaderLine(l)
  return l:sub(1,1) == "@"
end


------------------------------------------------------------------------------------------------------------------------
-- Compares the start position of two intervals
--
-- Return values:
-- 1: -1, 0 or 1 (respectively for a < b, a == b, a > b)
------------------------------------------------------------------------------------------------------------------------
local function compIntervals(a, b)

  local function c(x,y)
    if x < y then return -1
    elseif x > y then return 1
    else return 0 end
  end
  -- same chromosomes
  if a.chr == b.chr then return c(a.start, b.start)
  else
    x = tonumber(a.chr)
    y = tonumber(b.chr)
    if x and y then return c(x,y)
    else return c(a.chr, b.chr) end
  end
end

------------------------------------------------------------------------------------------------------------------------
-- Compare function to sort a list of intervals (use with table.sort)
--
-- Return values:
-- 1: true if a < b, false otherwise.
------------------------------------------------------------------------------------------------------------------------
local function sortCompInterval(a, b)
  if a.chr == b.chr then return a.start < b.start end
  local x = tonumber(a.chr)
  local y = tonumber(b.chr)
  if x and y then
    return x < y end
  return a.chr < b.chr
end

------------------------------------------------------------------------------------------------------------------------
-- Checks if the interval is a valid human genome interval
--
-- Return values:
-- 1: true/false
------------------------------------------------------------------------------------------------------------------------
local function isValidInterval(interval)
  local x
  if interval.chr == "X" then x = 23
  elseif interval.chr == "Y" then x = 24
  elseif interval.chr == "MT" then x = 25
  else x = tonumber(interval.chr) end
  return x >= 1 and x <= 25 and interval.start < interval.finish and chr_limits[x] > interval.finish
end


------------------------------------------------------------------------------------------------------------------------
-- Checks if the intervals are overlapping. Intervals are said to overlap if one of the following is true:
-- i1i2: i1 starts before i2, but ends inside i2.
-- i2i1: i2 starts before i1, but ends inside i1.
-- i1_inside: i1 is fully contained inside i2.
-- i2_inside: i2 is fully contained inside i1.
--
-- Return values:
-- 1: true/false
-- 2: if true, returns "i1i2", "i2i1", "i1_inside", "i2_inside"
------------------------------------------------------------------------------------------------------------------------
local function isOverlappingInterval(i1, i2)
  if i1.chr ~= i2.chr then return false
  elseif i1.start < i2.start and i1.finish < i2.finish then return true, "i1i2"
  elseif i2.start < i1.start and i2.finish < i1.finish then return true, "i2i1"
  elseif i1.start > i2.start and i1.finish < i2.finish then return true, "i1_inside"
  elseif i2.start > i1.start and i2.finish < i1.finish then return true, "i2_inside"
  else return false end
end