library(raster)
library(spatstat)
library(vegan) 

########## sub-sampling plots in spp ###############
######################################
ventana <- function(pp,size){
  punto<-rpoint(1,win=pp$window)
  x<- punto$x
  y<- punto$y
  s2<-size/2
  return(owin(c(x-s2, x+s2), c(y-s2, y+s2)))
}
############################################


shinyApp(
  ui = fluidPage(
    
    ##Hacemos el programa ----
    
    mainPanel(
      plotOutput("plot")
    ) ),
  
  server = function(input, output) {
    
    inputPanel(
      selectInput("Dens_D", label = "Density of disturbed habitat:",
                  choices = seq(200, 1000, by=100) , selected = 20),
      
      sliderInput("Sample", label = "Number of samples:",
                  min = 5, max = 100, value = 10, step = 1)
    )
    
    output$plot <- renderPlot({
      
      # CONSTRUCT ANALYSIS WINDOW USING THE FOLLOWING:
      xrange=c(0, 5000)
      yrange=c(0, 5000)
      window<-owin(xrange, yrange)
      
      # Build maps from random points and interpole in same line
      set.seed(25)
      elev   <- density(rpoispp(lambda=0.05, win=window)) #
      elev <- elev*2000
      
      set.seed(3)
      hum <- density(rpoispp(lambda=0.05, win=window))
      
      
      
      sppE_d <- rpoint(input$Dens_D, elev) 
      sppH_d <- rpoint(input$Dens_D, hum)
      
      sppE_nd <- rpoint(input$Dens_D*2, elev) 
      sppH_nd <- rpoint(input$Dens_D*2, hum)
      
      
      ## Genero ppp de los dos  grupos de especies
      
      xE_d <- sppE_d$x
      yE_d <- sppE_d$y
      NsppE_d <- as.factor(c(sample(paste("sp", 1:20, sep=""), input$Dens_D*0.16, replace=T),
                             sample(paste("sp", 21:30, sep=""), input$Dens_D*0.84, replace=T)))
      
      ComE <- ppp(x=xE_d,y=yE_d, marks= NsppE_d, 
                  window=window) 
      
      xH_d <- sppH_d$x
      yH_d <- sppH_d$y
      NsppH_d <- as.factor(c(sample(paste("sp", 31:40, sep=""), input$n_breaks*0.08, replace=T), 
                             sample(paste("sp", 41:50, sep=""), input$n_breaks*0.92, replace=T)))
      
      
      ComH_d <- ppp(x=xH_d,y=yH_d, marks= NsppH_d, 
                    window=window) 
      ComEH_d<-superimpose(ComE_d,ComH_d)
      
      
      k=30
      
      unif.sample_d<- vector(mode = "list", length= k)
      
      for (i in 1:k){
        unif.sample_d[[i]]<-ComEH_d[ventana(ComEH_d, 300)]
      }
      
      
      plot(ComEH_d$x, ComEH_d$y, type="p", pch=21, bg=as.numeric(ComEH_d$marks))
      
      for (i in 1:k){
        plot(unif.sample_d[[i]]$window, add=T, lwd=2)
      }
      
    })
    
  }
  
)



library(vegan)

richness_d <- matrix(NA,k+1,specnumber(table(ComEH_d$marks)))
colnames(richness_d) <- levels(ComEH_d$marks)

for (i in 1:k) {
  richness_d[i,] <- table(unif.sample_d[[i]]$marks)
}

richness[nrow(richness_d),] <- table(ComEH_d$marks)

specnumber(colSums(richness_d[1:k,]))

#############