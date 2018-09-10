## app.R ##
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      box(plotOutput("distPlot")),
      
      box(
        title = "Controls",
        style = "padding-bottom: 20px;",
        column(4, sliderInput("k",
                              "Samples:",
                              min = 2,
                              max = 50,
                              value = 0)),
        column(4, selectInput("bins", label = "Density:",
                              choices =c( 30, 50, 100, 200, 500) , 
                              selected = 50)),
        column(4, checkboxInput(inputId = "Total",
                                label = strong("Show yotal richness"),
                                value = FALSE))
      ))
      )
    )



server <- function(input, output) {
 
  output$distPlot <- renderPlot({
    # generate bins based on as.numeric(input$bins) from ui.R
    library(spatstat)
    library(raster)
    library(vegan) 
    
    ########## sub-sampling plots in spp ###############
    ######################################
    ventana <- function(pp,size){
      punto<-rpoint(1,win=c((pp$window$xrange)-c(-15,15),
                            (pp$window$yrange)-c(-15,15)))
      x<- punto$x
      y<- punto$y
      s2<-size/2
      return(owin(c(x-s2, x+s2), c(y-s2, y+s2)))
    }
    ############################################
    
    xrange=c(0, 500)
    yrange=c(0, 500)
    window<-owin(xrange, yrange)
    
    # Build maps from random points ----
    set.seed(25)
    elev   <- density(rpoispp(lambda=0.05, win=window)) #
    elev <- elev*2000
    set.seed(3)
    hum <- density(rpoispp(lambda=0.05, win=window))
    
    sppE_d <- rpoint(as.numeric(input$bins), elev) 
    sppH_d <- rpoint(as.numeric(input$bins), hum)
    
    sppE_nd <- rpoint(as.numeric(input$bins)*2.5, elev) 
    sppH_nd <- rpoint(as.numeric(input$bins)*2.5, hum)
    
    
    ##Ensamble of communities----
    
    ###Disturbed ----
    xE_d <- sppE_d$x
    yE_d <- sppE_d$y
    NsppE_d <- as.factor(c(sample(paste("sp", 1:20, sep=""), as.numeric(input$bins)*0.16, replace=T),
                           sample(paste("sp", 21:30, sep=""), as.numeric(input$bins)*0.84, replace=T), 
                           sample(paste("sp", 21:30, sep=""),as.numeric(input$bins)-length(c(sample(paste("sp", 1:20, sep=""), as.numeric(input$bins)*0.16, replace=T),
                                                                                             sample(paste("sp", 21:30, sep=""), as.numeric(input$bins)*0.84, replace=T))))))
    
    ComE_d <- ppp(x=xE_d,y=yE_d, marks= NsppE_d, 
                  window=window) 
    
    xH_d <- sppH_d$x
    yH_d <- sppH_d$y
    NsppH_d <- as.factor(c(sample(paste("sp", 31:40, sep=""), as.numeric(input$bins)*0.08, replace=T), 
                           sample(paste("sp", 41:50, sep=""), as.numeric(input$bins)*0.92, replace=T), 
                           sample(paste("sp", 31:40, sep=""),as.numeric(input$bins)-length(c(sample(paste("sp", 31:40, sep=""),as.numeric(input$bins)*0.08, replace=T), 
                                                                                             sample(paste("sp", 41:50, sep=""), as.numeric(input$bins)*0.92, replace=T))))))
    ComH_d <- ppp(x=xH_d,y=yH_d, marks= NsppH_d, 
                  window=window) 
    ComEH_d<-superimpose(ComE_d,ComH_d)
    
    
    ###Undisturbed ----
    xE_nd <- sppE_nd$x
    yE_nd <- sppE_nd$y
    NsppE_nd <- as.factor(c(sample(paste("sp", 1:20, sep=""), (2.5*as.numeric(input$bins))*0.16, replace=T),
                            sample(paste("sp", 21:30, sep=""), (2.5*as.numeric(input$bins))*0.84, replace=T), 
                            sample(paste("sp", 21:30, sep=""),2.5*as.numeric(input$bins)-length(c(sample(paste("sp", 1:20, sep=""), (2.5*as.numeric(input$bins))*0.16, replace=T),
                                                                                                  sample(paste("sp", 21:30, sep=""), (2.5*as.numeric(input$bins))*0.84, replace=T))))))
    
    ComE_nd <- ppp(x=xE_nd,y=yE_nd, marks= NsppE_nd, 
                   window=window) 
    
    xH_nd <- sppH_nd$x
    yH_nd <- sppH_nd$y
    NsppH_nd <- as.factor(c(sample(paste("sp", 31:40, sep=""), as.numeric(input$bins)*0.08*2.5, replace=T), 
                            sample(paste("sp", 41:50, sep=""), as.numeric(input$bins)*0.92*2.5, replace=T), 
                            sample(paste("sp", 31:40, sep=""),2.5*as.numeric(input$bins)-length(c(sample(paste("sp", 31:40, sep=""),as.numeric(input$bins)*0.08*2.5, replace=T), 
                                                                                                  sample(paste("sp", 41:50, sep=""), as.numeric(input$bins)*0.92*2.5, replace=T))))))
    ComH_nd <- ppp(x=xH_nd,y=yH_nd, marks= NsppH_nd, 
                   window=window) 
    ComEH_nd<-superimpose(ComE_nd,ComH_nd)
    
    ##Hacemos el muestreo de las comunidades----
    
    k <- input$k
    
    unif.sample_d<- vector(mode = "list", length= k)
    
    for (i in 1:k){
      unif.sample_d[[i]]<-ComEH_d[ventana(ComEH_d, 30)]
    }
    
    unif.sample_nd<- vector(mode = "list", length= k)
    
    for (i in 1:k){
      unif.sample_nd[[i]]<-ComEH_nd[ventana(ComEH_nd, 30)]
    }
    
    ###Calculamos la riqueza----
    
    ##Riqueza disturbada
    richness_d <- matrix(NA,k+1,specnumber(table(ComEH_d$marks)))
    colnames(richness_d) <- levels(ComEH_d$marks)
    
    for (i in 1:k) {
      richness_d[i,] <- table(unif.sample_d[[i]]$marks)
    }
    richness_d[nrow(richness_d),] <- table(ComEH_d$marks)
    
    ##Riqueza no disturbada
    
    richness_nd <- matrix(NA,k+1,specnumber(table(ComEH_nd$marks)))
    colnames(richness_nd) <- levels(ComEH_nd$marks)
    
    for (i in 1:k) {
      richness_nd[i,] <- table(unif.sample_nd[[i]]$marks)
    }
    richness_nd[nrow(richness_nd),] <- table(ComEH_nd$marks)
    
    
    ###Graficamos las comunidades -----
    
    par(mfcol= c(1,3), mar=c(2,2,1,1))
    
    plot(ComEH_d$x, ComEH_d$y, pch=c(21:25), bg=as.numeric(ComEH_d$marks), cex=0.8, yaxt="n", xaxt="n")
    plot(elev, col=grey.colors(5), add=T)
    points(ComEH_d$x, ComEH_d$y, pch=c(21:25), 
           bg=as.numeric(ComEH_d$marks), col=as.numeric(ComEH_d$marks), cex=0.9)
    for (i in 1:k){
      plot(unif.sample_d[[i]]$window, add=T, lwd=1)
    }
    
    barplot(c(specnumber(colSums(richness_d[1:k,])),
              specnumber(colSums(richness_nd[1:k,]))), ylim=c(0,55), 
            las=2, col="grey20", border="white")
    
    plot(ComEH_nd$x, ComEH_nd$y, pch=c(21:25), bg=as.numeric(ComEH_nd$marks), cex=0.8, yaxt="n", xaxt="n")
    plot(elev, col=grey.colors(5), add=T)
    points(ComEH_nd, pch=c(21:25), 
           bg=as.numeric(ComEH_nd$marks),col=as.numeric(ComEH_nd$marks), cex=0.9)
    for (i in 1:k){
      plot(unif.sample_d[[i]]$window, add=T, lwd=1)
    }
    
    
    # draw the histogram with the specified number of bins
  })
  
}

shinyApp(ui, server)