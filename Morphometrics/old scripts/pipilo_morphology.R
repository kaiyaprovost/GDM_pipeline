Pipilo <- read.csv("~/Documents/Classes/Stephanie/Pipilo_fuscus_measurements_27September2017.csv")

names(Pipilo)

hist(x=Pipilo$BILL.LENGTH..mm.)
hist(x=Pipilo$BILL.WIDTH..mm.)

summary(Pipilo)

plot(x=Pipilo$BILL.WIDTH..mm.,y=Pipilo$BILL.LENGTH..mm.)
trendline <- lm(Pipilo$BILL.LENGTH..mm. ~ Pipilo$BILL.WIDTH..mm.)
trendline
summary(trendline)
abline(trendline)

boxplot(Pipilo$BILL.LENGTH..mm. ~ Pipilo$CATALOG.NUMBER)

barplot(Pipilo$BILL.LENGTH..mm.)

class(Pipilo) = "data.frame"

Pipilo[,2] == "SONORAN"
Pipilo$WHICH.SIDE.OF.CFB == "SONORAN"
sum(Pipilo$WHICH.SIDE.OF.CFB == "SONORAN")
which(Pipilo$WHICH.SIDE.OF.CFB == "SONORAN")

SonoranBirds <- Pipilo[Pipilo$WHICH.SIDE.OF.CFB=="SONORAN",]
ChihuahuanBirds <- Pipilo[Pipilo$WHICH.SIDE.OF.CFB=="CHIHUAHUAN?",]
UnclearBirds <- Pipilo[Pipilo$WHICH.SIDE.OF.CFB=="UNCLEAR",]

mean(SonoranBirds$BILL.LENGTH..mm.)
sd(SonoranBirds$BILL.LENGTH..mm.)
median(SonoranBirds$BILL.LENGTH..mm.)
IQR(SonoranBirds$BILL.LENGTH..mm.)
quantile(SonoranBirds$BILL.LENGTH..mm. , 0.33)

t.test(x=SonoranBirds$BILL.LENGTH..mm.,y=ChihuahuanBirds$BILL.LENGTH..mm.)
t.test(x=SonoranBirds$BILL.WIDTH..mm.,
       y=ChihuahuanBirds$BILL.WIDTH..mm.)
t.test(x=SonoranBirds$BILL.HEIGHT..mm.,y=ChihuahuanBirds$BILL.HEIGHT..mm.)

boxplot(SonoranBirds$BILL.WIDTH..mm.,ChihuahuanBirds$BILL.WIDTH..mm.,
        names=c("SON","CHI"))

boxplot(SonoranBirds$BILL.WIDTH..mm.,ChihuahuanBirds$BILL.WIDTH..mm.,
        UnclearBirds$BILL.WIDTH..mm.,
        names=c("SON","CHI","UNC"))
