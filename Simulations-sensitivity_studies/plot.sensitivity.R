
## Plot results of the power study

setwd("/nfs/users/rg/dgarrido/sqtlseeker/simulations/Ha.dgarrido.extended2")
load("Ha.dgarrido.RData")

## Sensitivity vs angle 
## 1 factor (A), 1 level

## Between different Ho's

S <- subset(storage,scenario == 1 & ho == 0.2 & ha == "Ha2", select=c(t,angle,potA))
plot(as.numeric(as.character(S$angle)),as.numeric(as.character(S$potA)),col="blue",
    ylim=c(0,1),pch=20, main= "Comparison between different Ho's (1 factor, 1 level)",
    ylab="Potency of factor A", xlab="Angle (deg)")

S2 <- subset(storage,scenario == 1 & ho == 0.8 & ha == "Ha2", select=c(t,angle,potA))
points(as.numeric(as.character(S2$angle)),as.numeric(as.character(S2$potA)),col="red",pch=20)
points(S$t,S2$potA,col="red",pch=20)

leg.text=c("(0.2,0.114, ...)","(0.8,0.028, ...)") # legend
colors= c("blue","red")
legend(x='bottomleft',lty=1,legend=leg.text,col=colors,title="Ho")


## Between different Ha's

S <- subset(storage,scenario == 1 & ho == 0.2 & ha == "Ha1", select=c(t,angle,potA))
plot(S$t,S$potA,col="blue",pch=20,main= "Comparison between different Ha's (1 factor, 1 level)",
    ylab="Potency of factor A", xlab="t value")
points(as.numeric(as.character(S$angle)),as.numeric(as.character(S$potA)),col="blue",pch=20,
    main= "Comparison between different Ha's (1 factor, 1 level)",ylab="Potency of factor A", 
    xlab="Angle (deg)")

S2 <- subset(storage,scenario == 1 & ho == 0.2 & ha == "Ha2", select=c(t,angle,potA))
plot(as.numeric(as.character(S2$angle)),as.numeric(as.character(S2$potA)),col="red",pch=20,
    ylim=c(0,1),main= "Comparison between different Ha's (1 factor, 1 level)",ylab="Potency of factor A", 
    xlab="Angle (deg)")
points(S$t,S6$potA,col="red",pch=20)

leg.text=c("Ha1","Ha2") # legend
colors= c("blue","red")
legend(x='bottomright',lty=1,legend=leg.text,col=colors)


## Positive vs neative angles 

S <- subset(storage,scenario == 1 & ho == 0.8 & ha == "Ha2" & as.numeric(as.character(t))<=0, select=c(potA)) #-
S2 <- subset(storage,scenario == 1 & ho == 0.8 & ha == "Ha2" & as.numeric(as.character(t))>=0, select=c(potA)) #+

plot(as.numeric(as.character(S2$potA)),rev(as.numeric(as.character(S$potA))),col="blue",pch=20, ylim=c(0,1),xlim=c(0,1),,main="(+) VS (-) angles",ylab="Potency of A (-)",xlab="Potency of A (+)")
summary(lm(as.numeric(as.character(S2$potA))~rev(as.numeric(as.character(S$potA)))))
lines(c(0,1),c(0,1))

## Different scenarios in terms of one/more factors/levels

## One factor/one level
S <- subset(storage,scenario == 1 & ho == 0.8 & ha == "Ha1", select=c(angle,potA,potB,potAB))
plot(as.numeric(as.character(S$angle)),as.numeric(as.character(S$potA)),col="blue",ylim=c(0,1),pch=20,main= "Scenario: 1 factor (A), 1 level",ylab="Potency", xlab="Angle (deg)")
points(as.numeric(as.character(S$angle)),as.numeric(as.character(S$potB)),col="red",pch=20)
points(as.numeric(as.character(S$angle)),as.numeric(as.character(S$potAB)),col="green",pch=20)

leg.text=c("A","B","AB") # legend
colors= c("blue","red","green")
legend(x=0,y=0.5,lty=1,legend=leg.text,col=colors,title="Factors")

## One factor/two levels
S2 <- subset(storage,scenario == 2 & ho == 0.8 & ha =="Ha1", select=c(angle,potA,potB,potAB))
plot(as.numeric(as.character(S2$angle)),as.numeric(as.character(S2$potA)),col="blue",pch=20,ylim=c(0,1),main= "Scenario: 1 factor (A), 2 levels",ylab="Potency", xlab="Angle (deg)")
points(as.numeric(as.character(S2$angle)),as.numeric(as.character(S2$potB)),col="red",pch=20)
points(as.numeric(as.character(S2$angle)),as.numeric(as.character(S2$potAB)),col="green",pch=20)

leg.text=c("A","B","AB") # legend
colors= c("blue","red","green")
legend(x=-0.5,y=0.5,lty=1,legend=leg.text,col=colors,title="Factors")

## Two factors/one level each
S3 <- subset(storage,scenario == 3 & ho == 0.8 & ha == "Ha1", select=c(angle,potA,potB,potAB))
plot(as.numeric(as.character(S3$angle)),as.numeric(as.character(S3$potA)),col="blue",pch=20,ylim=c(0,1),main= "Scenario: 2 factors (A,B), 1 level",ylab="Potency", xlab="Angle (deg)")
points(as.numeric(as.character(S3$angle))*0.5,as.numeric(as.character(S3$potB)),col="red",pch=20)
points(as.numeric(as.character(S3$angle))*0.5,as.numeric(as.character(S3$potAB)),col="green",pch=20)

leg.text=c("A","B","AB") # legend
colors= c("blue","red","green")
legend(x="center",lty=1,legend=leg.text,col=colors,title="Factors")

