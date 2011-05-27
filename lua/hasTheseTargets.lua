-- This script takes two interval lists (a and b) and creates another interval list (c) with the intervals in a that are
-- inside the targets in b. For example, if a is a list of exon targets and b a list of genes, then c would be all exons
-- covered by the genes in file b.
--
-- Author: carneiro
-- Date: 5/26/2011

local targetList = arg[1]  -- list of targets to check if they are part of the set
local targetSet  = arg[2]  -- set of targets


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
-- 1: index pairs table with each interval (each interval is a table of the form {c, s, e, strand, info})
-- 2: header string
------------------------------------------------------------------------------------------------------------------------
local function readIntervals(filename)
  t = {}
  header = ""
  for l in io.lines(filename) do
    if l:sub(1,1) == "@" then header = header .. l .."\n"
    else
        local c, s, e, strand, info = parseIntervalLine(l)
        table.insert(t, {c=c, s=tonumber(s), e=tonumber(e), strand=strand, info=info})
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
local function isSameInterval(i, c, s, e)
  return i.c == c and s == i.s and e == i.e
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
  if a.c == b.c then return c(a.s, b.s)
  else
    x = tonumber(a.c)
    y = tonumber(b.c)
    if x and y then return c(x,y)
    else return c(a.c, b.c) end
  end
end

------------------------------------------------------------------------------------------------------------------------
-- Creates a new interval
--
-- Return values:
-- 1: interval table {c, s, e, st, info}
------------------------------------------------------------------------------------------------------------------------
local function newInterval(c,s,e,st,info)
  return {c=c, s=tonumber(s), e=tonumber(e), st=st, info=info }
end

local function printInterval(i)
  print(i.c, i.s, i.e, i.st, i.info)
end

local function intervalContainsInterval(a, b)
  return a.c == b.c and a.s <= b.s and a.e >= b.e
end

local function isInterceptingInterval(a, b)
  return a.s < b.s and a.e < b.e and a.e > b.s
end

local function findInterval(i, intervals)
  local start = 1
  local finish = table.getn(intervals)
  local current = math.floor((start + finish) / 2)

  while start < finish and not intervalContainsInterval(i, intervals[current]) and not isInterceptingInterval(i, intervals[current]) do
    if compIntervals(i, intervals[current]) < 0 then
      finish = current - 1
    else
      start = current + 1
    end
    current = math.floor((start + finish) / 2)
  end
  return intervalContainsInterval(i, intervals[current]), current
end




local a, header = readIntervals(targetList)
io.write(header)

for l in io.lines(targetSet) do
  if not isIntervalHeaderLine(l) then
    local c, s, e, st, info = parseIntervalLine(l)
    local intA = newInterval(c,s,e,st,info)
    local intervalExists, i = findInterval(intA, a)
    if intervalExists then
      print(a[i].c, a[i].s, a[i].e, st, info)
    end
  end
end
