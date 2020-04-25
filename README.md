# hybridStructures

# Download data from IDS

# pre-filter to obtain machine-readable word lists

# use jAligns GotohBigram version to obtain alignments

# Obtain thresholds using R

library(mixtools)
files = list.files(getwd(), pattern="sm")

for(i in files){
 t = read.table(i)
 mix = normalmixEM(unlist(t))
 write(paste(mix$mu[1], " ", mix$mu[2]), file="mu", append=TRUE)
 write(paste(mix$lambda[1], " ", mix$lambda[2]), file="lambda", append=TRUE)
 write(paste(mix$sigma[1], " ", mix$sigma[2]), file="sigma", append=TRUE)
 pdf(paste(i, ".pdf", sep=""))
 plot(mix, which=2)
 dev.off()
}


m = read.table("mu")
s = read.table("sigma")
l = read.table("lambda")

intersection_plus <- vector(, dim(s)[1]) 
intersection_minus <- vector(, dim(s)[1])
intersection <- vector(, dim(s)[1])  
for(i in 1:dim(s)[1]){
  a = (1/(2*s[i,1]^2)) - 1/(2*s[i,2]^2)
  b = (m[i,2]/s[i,2]^2) - (m[i,1]/s[i,1]^2)
  c= (m[i,1]^2)/(2*s[i,1]^2) - (m[i,2]^2)/(2*s[i,2]^2) - log((l[i,1]*s[i,2])/(l[i,2]*s[i,1]))

  delta = b^2 - 4*a*c
 
  if((delta>0) && (a!=0)){
    intersection_plus[i] = (-b + sqrt(delta))/(2*a)
    intersection_minus[i] = (-b - sqrt(delta))/(2*a)
  } else {
    intersection_plus[i] = -20
    intersection_minus[i] = -20
  }
  if((intersection_plus[i] < m[i,1] && intersection_plus[i] > m[i,2]) || (intersection_plus[i] > m[i,1] && intersection_plus[i] < m[i,2])){
    intersection[i] = intersection_plus[i]
  } else if ((intersection_minus[i] < m[i,1] && intersection_minus[i] > m[i,2]) || (intersection_minus[i] > m[i,1] && intersection_minus[i] < m[i,2])){
    intersection[i] = intersection_minus[i]
  } else {
    intersection[i] = -20
  }
}

for(i in 1:length(intersection)){
  write(paste(unlist(strsplit(files[i], ".sm.nscores")), " ", intersection[i]), file="threshold", append=TRUE)
}


#Filtering on command line

for f in *.gz 
do
t=$(grep $f threshold | awk '{print $2}')                      
zgrep "IDS" $f | LC_ALL=C awk '$7 >= theta' theta=$t  > $f.aboveT
done

#obtain distances using R
