count_zeros = function(list) {
  zeros = 0
  for (x in list) {
    if (x == 0.0) {
      zeros = zeros + 1
    }
  }
  zeros
}
  

load = function(max_rows) {
  files = list.files(path=".", pattern="304NA.*")
  #max_rows = -1
  
  #FREESTANDING as a filter
  #HIT_TWICE for ZEROS...
  
  print ("Parsing file 1")
  t = read.table(files[1],header=T, nrows = max_rows)
  f = data.frame(loc=t$location,gc=t$gc,freestanding=t$freestanding)
  ht = data.frame(1:nrow(f))

  for (file in files) {
    print (file)
    t = read.table(file, header=T, nrows = max_rows)
    norm_cov = t$normalized_coverage
    #names(norm_cov) = c("norm_cov.1")
    f=cbind (f, norm_cov)
    ht=cbind (ht, t$hit_twice)
  }


  wgs = read.table("/seq/dirseq/analysis/agilent/rt-pcr/perfdata//OV-0751-WGS.baits.coverage.txt", header=T, nrows = max_rows)
  f=cbind (f, wgs_norm_cov = wgs$normalized_coverage)
  
  f=cbind(f,ht)
  
  # Compute normalized variance
  print("Calculating variance")
  var = apply(f[4:10], 1, var)
  print("Calculating std. dev.")
  sd = apply(f[4:10], 1, sd)
  print("Calculating mean")
  mean = apply(f[4:10], 1, mean)
  print("Binding normalized variance")
  f=cbind (f, normvar=var/mean/mean)
  print("Binding normalized std. dev.")
  f=cbind (f, normsd=sd/mean)
  print("Binding mean")
  f=cbind (f, mean=mean)
  print("Binding std. dev.")
  f=cbind (f, sd=sd)
  print("Binding variance")
  f=cbind (f, var=var)
  
  print("Calculating and binding number of zeros")
  count_zeros = apply(f[4:10], 1, count_zeros)
  num_not_hit_twice = apply(f[12:18], 1, count_zeros)
  f=cbind(f, count_zeros, num_not_hit_twice)
  
  print ("Parsing sequences file")
  seqs = read.table("whole_exome_agilent_designed_120.design.1line.sorted2",header=T,nrows=max_rows)
  f=cbind (f, seqs)
  
  #of = f[order(f$normvar),]
}

write_splits = function(f) {
  set.seed(0987123409)

  # Low variance
  nz = f[f$count_zeros < 1 & f$freestanding==1,] # Take reads with no zeros
  d =         write_split(nz, "Low_GC_Norm_Coverage",  0.0, 0.35, 0.8, 1.2, 0.0, 0.3, 0.0)  
  d = rbind(d,write_split(nz, "Mid_GC_Norm_Coverage",  0.45, 0.55, 0.8, 1.2, 0.0, 0.1, 0.0))
  d = rbind(d,write_split(nz, "High_GC_Norm_Coverage", 0.63, 1.0, 0.8, 1.2, 0.0, 0.3, 0.0))
  d = rbind(d,write_split(nz, "Low_GC_Undercovered",  0.0, 0.35, 0.2, 0.3, 0.0, 0.3, 0.0))
  d = rbind(d,write_split(nz, "Mid_GC_Undercovered",  0.45, 0.55, 0.2, 0.3, 0.0, 0.3, 0.0))
  d = rbind(d,write_split(nz, "High_GC_Undercovored", 0.63, 1.0, 0.2, 0.3, 0.0, 0.3, 0.0))
  az = f[f$count_zeros == 7 & f$freestanding==1,] # Take reads with all zeros
  d = rbind(d,write_split(az, "Low_GC_No_Coverage",  0.0, 0.35, 0.0, 0.1, -1.0, -1.0, 0.1))
  d = rbind(d,write_split(az, "Mid_GC_No_Coverage",  0.45, 0.55, 0.0, 0.1, -1.0, -1.0, 0.1))
  d = rbind(d,write_split(az, "High_GC_No_Coverage", 0.63, 1.0, 0.0, 0.1, -1.0, -1.0, 0.01))

  # High variance
  d = rbind(d,write_split(nz, "Mid_GC_Norm_Coverage_High_Variation",  0.45, 0.55, 0.8, 1.2, 0.355, 1000.0))
  d
}

write_split = function(data, label, gc_low, gc_high, cov_low, cov_high, normsd_low, normsd_high, wgs_cov_low) {
  if (normsd_high < 0.0) {
    # We have no coverage samples
    s = data[data$gc >= gc_low & data$gc <= gc_high  & data$mean >= cov_low & data$mean <= cov_high & data$wgs_norm_cov >= wgs_cov_low,]
    #s = s[order(runif(nrow(s))),] # Randomize rows
    s = s[order(s$wgs_norm_cov, decreasing = T),] # order according to norm SD
  }else{
    # We have low or normal coverage samples, so take those with tightest norm SDs
    s = data[data$gc >= gc_low & data$gc <= gc_high  & data$mean >= cov_low & data$mean <= cov_high & data$normsd >= normsd_low & data$normsd <= normsd_high ,]
    s = s[order(s$normsd),] # order according to norm SD
  }
  # & data$mean < 1.1 & data$mean > 0.9,]
  # & data$mean >= cov_low & data$mean <= cov_high 
  #print(s)
  print(nrow(s))
  s = s[1:50, ] #-c(3,11,12:18,19,23:25)]
  s = cbind(class=rep(label,50), s)
  s
}

#f=load()
#nz=f[f$count_zeros < 1,]
#print(summary(nz))

create_500 = function(f) {

  f = load(-1)
  s = write_splits(f)
  write.csv(s, "500_exome_baits_for_nanostring.csv")

}
