#Evaluación

#Cargamos las preguntas

library(readxl)
eval <- read_excel("Evaluación.xlsx")

x <- eval[sample(1:22, 10, replace = F),]
x <- rbind.data.frame(x[order(x$Número),],
       eval[23,]) 
ranE <- data.frame(Numero=1:1044, pregunta=1:1044)

x <- seq(11,1044, by=11)
y <- seq(1,1033, by=11)

for (i in 1:95){
ranE[y[i]:x[i],] <- rbind(eval[sample(1:22, 10, replace = F),], c("",""))
}

write.csv(ranE, "eval1.csv")
