# hybridStructures

## 1. Download data from IDS

## 2. pre-filter to obtain machine-readable word lists

## 3. use jAligns GotohBigram version to obtain alignments

## 4. Obtain thresholds using R

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


## 5. Filtering on command line

for f in *.gz 
do
t=$(grep $f threshold | awk '{print $2}')                      
zgrep "IDS" $f | LC_ALL=C awk 'if($7 >= theta) print $7' theta=$t  > $f.aboveT
done

## 6. obtain distances using R

files = list.files(getwd(), pattern="aboveT")
sims = numeric(length(files))
for(i in 1:length(files)){
 t = read.table(files[i])
 sims[i] = sum(unlist(t))/length(unlist(t))
}

sim = matrix(0,14,14)
sim[upper.tri(sim,diag=TRUE)] <- sims

for (i in 2:14){
 for (j in 1:i-1){
  sim[i,j] = sim[j,i]}}

dist = matrix(0,14,14)

for (i in 1:14){
 for (j in 1:14){
  dist[i,j] = sim[i,i] + sim[j,j] - 2*sim[i,j]
 }
}

colnames(dist) <- c("Catalan", "Danish", "Dutch", "English", "French", "German", "Italian", "Latin", "OldEnglish", "OldHighGerman", "Portuguese", "Romanian", "Spanish", "Swedish")
rownames(dist) <- c("Catalan", "Danish", "Dutch", "English", "French", "German", "Italian", "Latin", "OldEnglish", "OldHighGerman", "Portuguese", "Romanian", "Spanish", "Swedish")

write.table(dist, file="distances", row.names=TRUE, col.names=TRUE)

