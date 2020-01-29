x = 5 
x %/% 4 ## division integer 
x %% 4 ## remainder

# "hello
# my name is Kaiya
# I am 26 years old
# I am a PhD student"
## there
## hello

path_to_file = "/Users/kprovost/Documents/Dissertation/CHAPTER3_TRAITS/Internships/Master_Spreadsheet_Morphology_round2 (Forcellati_Meghan) - KLP_Master_Spreadsheet_Morpholo.csv"
file = "Master_Spreadsheet_Morphology_round2 (Forcellati_Meghan) - KLP_Master_Spreadsheet_Morpholo.csv"

getwd()
setwd("/Users/kprovost/Documents/Dissertation/CHAPTER3_TRAITS/Internships/")
getwd()

Morphology_Data = read.csv(file)
names(Morphology_Data)
colnames(Morphology_Data)
rownames(Morphology_Data)
summary(Morphology_Data)

catalog = Morphology_Data$CATALOG.NUMBER
catalog_unique = unique(Morphology_Data$CATALOG.NUMBER)

Morphology_Data[,5]

Morphology_Data$COUNTRY == "MEXICO"
sum(Morphology_Data$COUNTRY == "MEXICO")
which(Morphology_Data$COUNTRY == "MEXICO")
which(Morphology_Data$STATE == "ARIZONA")

mexican_birds = Morphology_Data[Morphology_Data$COUNTRY == "MEXICO",]

chihuahuan_birds = Morphology_Data[Morphology_Data$WHICH.SIDE.OF.CFB == "CHIHUAHUAN",]
sonoran_birds = Morphology_Data[Morphology_Data$WHICH.SIDE.OF.CFB == "SONORAN",]

mexican_fusca = Morphology_Data[Morphology_Data$COUNTRY=="MEXICO" & Morphology_Data$SPP == "FUSCA",]

nrow(mexican_birds)
nrow(mexican_fusca)

ncol(mexican_birds)
ncol(mexican_fusca)
length(mexican_birds)
length(mexican_fusca)

bill_length = Morphology_Data$BILL.LENGTH
bill_width = Morphology_Data$BILL.WIDTH
bill_height = Morphology_Data$BILL.HEIGHT
spp = Morphology_Data$SPP

plot(bill_length)
plot(bill_length,bill_width)

boxplot(bill_length~spp,las=2)
boxplot(Morphology_Data$BILL.LENGTH~Morphology_Data$SPP)
boxplot(bill_length~spp,las=2,col=c("red","blue"))

hist(bill_length,breaks=50)

par(mfrow=c(1,2))
plot(bill_length,bill_width)
plot(bill_length,bill_height)

par(mfrow=c(1,1))

Morphology_Data[,c(11:13)]

bill_correlation = cor(Morphology_Data[,c(11:13)],use="pairwise.complete.obs")

library("corrplot")

corrplot(bill_correlation,method="number")

# c("circle", "square", "ellipse", "number", "shade","color", "pie")


pca = prcomp(na.omit(Morphology_Data[,c(11:13)]))
pca = prcomp(na.omit(Morphology_Data[,c(11:13)]),center=TRUE,scale.=TRUE)

results = pca$x
summary(pca)
plot(pca,type="l")

colnames(results)
PC1 = results[,1]
PC2 = results[,2]
PC3 = results[,3]

par(mfrow=c(1,2))
plot(PC1,PC2)
plot(bill_width,bill_length)

plot(PC1,PC3)

boxplot(PC1~spp[c(1:648)],las=2)

for(i in 1:10) {
  print(i*2)
}

list_of_species = c("bird","mammal","lizard")
for(species in list_of_species) {
  print(species)
}

animal = "bird"
if(animal=="lizard"){
  print("hooray")
} else if(animal=="bird") {
  print("woohoo")
} else {
  print("boo")
}

count = 0 
while(count < 5){
  print(count)
  count = count + 1 
}


## day 2

#kipps index calculated (primary_length-secondary_length/primary_length*100)

kipps = function(primary,secondary) {
  
  temp = ((primary-secondary)/primary)*100
  
  return(temp)
  
}

## paste strings of text together

"/Users/kprovost/Documents/Dissertation/CHAPTER3_TRAITS/Internships/morphology/Master_Spreadsheet_Morphology_round2_6feb2019.csv"

directory = "/Users/kprovost/Documents/Dissertation/CHAPTER3_TRAITS/Internships/morphology/"
file_prefix = "Master_Spreadsheet_Morphology_round2_"
the_date = "6feb2019"
file_suffix = ".csv"

our_file = paste(directory,file_prefix,the_date,file_suffix,sep="")

our_csv = read.csv(our_file)

our_csv$CATALOG.NUMBER
our_csv$STATE == "TEXAS"

names(our_csv)

agg = aggregate(cbind(BILL.LENGTH,BILL.WIDTH) ~ CATALOG.NUMBER,
          data = our_csv,
          FUN=mean,
          na.action=na.pass)


plot(our_csv$BILL.LENGTH,our_csv$BILL.WIDTH,
       col="red")
points(agg$BILL.LENGTH,agg$BILL.WIDTH,
       pch="+")





our_csv_good = our_csv[our_csv$CONDITION=="",]
nrow(our_csv)
nrow(our_csv_good)

fusca = our_csv_good[our_csv_good$SPP=="FUSCA",]

t.test(fusca$BILL.LENGTH[fusca$WHICH.SIDE.OF.CFB=="SONORAN"],
       fusca$BILL.LENGTH[fusca$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])

fusca_son = fusca[fusca$WHICH.SIDE.OF.CFB=="SONORAN",]
fusca_chi = fusca[fusca$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]

fusca_chi_agg = aggregate(cbind(BILL.LENGTH,BILL.WIDTH) ~ CATALOG.NUMBER,
                          data = fusca_chi,
                          FUN=mean,
                          na.action=na.pass)

fusca_son_agg = aggregate(cbind(BILL.LENGTH,BILL.WIDTH) ~ CATALOG.NUMBER,
                data = fusca_son,
                FUN=mean,
                na.action=na.pass)

t.test(fusca_son_agg$BILL.LENGTH,fusca_chi_agg$BILL.LENGTH)



library(lme4)

model1 = glm(BILL.LENGTH ~ SPP, data = our_csv_good) ## AIC Akaike Information Criterion: 3287

summary(model1)

plot(model1)


model2 = glm(BILL.LENGTH ~ TARSUS.LENGTH + SPP, data = our_csv_good) ## AIC Akaike Information Criterion: 3253
model2
summary(model2)
plot(model2)


model3 = glm(BILL.LENGTH ~ TARSUS.LENGTH, data = fusca) ## AIC Akaike Information Criterion: 5028
model3
summary(model3)
plot(model3)

 
res = model3$residuals
res = unname(res)
residual_data = merge(model3$data,res,by="row.names",all.x=TRUE)

plot(residual_data$BILL.LENGTH,residual_data$y)
plot(fusca$BILL.LENGTH,fusca$TARSUS.LENGTH)
plot(residual_data$TARSUS.LENGTH,residual_data$y)




model5 = glm(BILL.LENGTH ~ TARSUS.LENGTH + SPP + BILL.WIDTH + SEX + AGE + DATE, 
             data = our_csv_good) ## AIC Akaike Information Criterion: 2977
model5
summary(model5)
car::Anova(model5)
plot(model5)



model6 = glm(BILL.LENGTH ~ TARSUS.LENGTH + SPP + SEX + AGE + DATE, 
             data = our_csv_good) ## AIC Akaike Information Criterion: 2976
model6
summary(model6)
car::Anova(model6)
plot(model6)

model7 = glm(BILL.LENGTH ~ TARSUS.LENGTH + SPP + SEX*AGE + DATE, 
             data = our_csv_good) ## AIC Akaike Information Criterion: 2938
summary(model7)
car::Anova(model7)

