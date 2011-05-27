-- This script takes the interval summary output from GATK's Depth of Coverage and the
-- interval list for 'genes of interest' and generates an intervals file sorted by
-- coverage.
--
-- Author: carneiro
-- Date: 5/26/2011
-------------------------------------------------------------------------------
-- Global script variables
-------------------------------------------------------------------------------
local coverage = arg[1]
local genesOfInterest= arg[2]

local genes = {}
local header = ""

-------------------------------------------------------------------------------
-- Helper functions
-------------------------------------------------------------------------------

local function isSameGene(i, c, s, e)
  return genes[i].c == c and s == genes[i].s and e == genes[i].e
end

local function isMergedGene(i, c, s, e)
  return genes[i].c == c and s <= genes[i].s and e >= genes[i].e
end

local function comp(a,b)
  return a.avgCoverage < b.avgCoverage
end

for l in io.lines(genesOfInterest) do
  if l:sub(1,1) == "@" then header = header .. l .."\n"
  else
    local c, s, e, info = l:match("(%w+)%s+(%d+)%s+(%d+)%s+%+%s+(.*)")
    table.insert(genes, {c=c, s=tonumber(s), e=tonumber(e), info=info})
  end
end

local i = 1
for l in io.lines(coverage) do
  local geneOk = false
  if l:match("(%w+)") ~= "Target" then    -- skip the first line (header)
    local c, s, e, totalCoverage, avgCoverage = l:match("(%w+):(%d+)%-(%d+)%s+(%d+)%s+([%d%.]+)")
    s = tonumber(s)
    e = tonumber(e)
    while genes[i] ~= nil and (isSameGene(i, c, s, e) or isMergedGene(i, c, s, e)) do
      genes[i].totalCoverage = tonumber(totalCoverage)
      genes[i].avgCoverage = tonumber(avgCoverage)
      geneOk = true
      i = i + 1
    end
    if not geneOk then
      error("Warning: Gene mismatch! Crazy!", c,s,e,"--",genes[i].c, genes[i].s, genes[i].e)
      break
    end
  end
end

table.sort(genes, comp)
io.write(header)
for _, g in ipairs(genes) do
  print(g.c,g.s,g.e,"+",g.avgCoverage..":::"..g.totalCoverage..":::"..g.info)
end
