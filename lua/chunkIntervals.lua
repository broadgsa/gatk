local infile = arg[1] or io.stdin
local outfile = (arg[2] and io.open(arg[2], "w")) or io.stdout


local function cutcommas(x)
  local n = ""
  for v in x:gmatch("(%d+)") do
    n = n .. v
  end
  if n == "" then return tonumber(x) else return tonumber(n) end
end

for l in io.lines(infile) do
  local chr, startPos, endPos = l:match("(.*):(%d+)%-([%d,]+)")
  startPos = cutcommas(startPos)
  endPos = cutcommas(endPos)
  for i=startPos,endPos,1000000 do
    if endPos > i+999999 then outfile:write(chr..":"..i.."-"..i+999999 .."\n")
    else outfile:write(chr..":"..i.."-"..endPos.."\n") end
  end
end