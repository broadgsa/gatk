--[[
-- Updates a list of BAM files to the latest version in the picard aggregation path
-- Usage:
--
-- lua updateSampleList.lua samples.list > updated_samples.list
 ]]
function latestVersion(sample)
  local version = tonumber(sample:match("/v(%d+)/"))
  f = io.open(sample)
  while (f == nil) do
    version = version + 1
    sample = sample:gsub("/v(%d+)/", "/v"..version.."/")
    f = io.open(sample)
  end
  return(sample)
end

for sample in io.lines(arg[1]) do
  print(latestVersion(sample))
end
