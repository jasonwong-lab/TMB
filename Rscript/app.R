library(shinydashboard)
library(shiny)
library(colourpicker)
library(GenomicRanges)
library(VariantAnnotation)
source("plot_tmb.R")

#====================================================#
## Sidebar content ####
#====================================================#
sidebar <- dashboardSidebar(
  width = 150,
  sidebarMenu(
    menuItem("Dashboard", tabName = "main", icon = icon("dashboard")),
    menuItem("TCGA data", tabName = "tcga", icon = icon("th")),
    menuItem("My data", tabName = "my", icon = icon("bar-chart"))
  )
)

#====================================================#
## Dashboard Home ####
#====================================================#
bodyHome <- tabItem(tabName = "main", value="main_panel",
                    fluidRow(
                      box(
                        title = "Welcome to improve TMB evaluation!", width = 12, status = "primary",
                        HTML("This TMB (tumor mutational burden) prediction tool can optimize TMB that is derived from target region sequencing data. Linear regression is employed to modeling the relationship between whole coding regions derived TMB and panel derived mutations in this tool, and it is based on pan-cancer data of 10,179 samples across 33 cancer types from The Cancer Genome Atlas (TCGA). From the TCGA data module, you can get the landscape of TMB across different cancer types. By uploading mutation file and panel bed file, you can obtain the optimized TMB prediction. For more information to use the tool, please check the <a href='https://github.com/fanghu-hku/TMB/blob/main/TMB_shiny_app_tutorial.rst' target='_blank'>tutorial</a>."),
                      )
                    ),
                    fluidRow(
                      box(
                        title = "TCGA data explore", width = 6, status = "warning",
                        p("Mutation burden profile of each cancer type from TCGA can be checked in this module. You can get a landscape of TMB across various panels and indicate the TMB of interest among the patients."),
                        img(src='tmb_profile.png', align = "center", width="100%")
                      ),
                      box(
                        title = "My data explore", width = 6,status = "warning",
                        p("This module can optimize panel derived TMB by linear regression model. The correlation of simulated panel TMB and whole exome TMB that are based on TCGA data is provided."),
                        img(src='tmb_cor.png', align = "center", width="100%")
                      ),
                     ),
                    fluidRow(
                      box(
                        title = "References", width = 12, status = "success",
                        #h4("If you use intervene, please cite this paper:"),
                        HTML("<h5>Buchhalter, Ivo, et al. \"Size matters: Dissecting key parameters for panel‐based tumor mutational burden analysis.\" International Journal of Cancer 144.4 (2019): 848-858.</h5>"),
                        HTML("<h5>Fancello, Laura, et al. \"Tumor mutational burden quantification from targeted gene panels: major advancements and challenges.\" Journal for immunotherapy of cancer 7.1 (2019): 183.</h5>"),
                        HTML("Khan, Aziz, and Anthony Mathelier. \"Intervene: a tool for intersection and visualization of multiple gene or genomic region sets.\" BMC bioinformatics 18.1 (2017): 287.</h5>")
                      )
                    )
)
#====================================================#
## tcga data module ####
#====================================================#
bodytcga <- tabItem(tabName = "tcga",
                    h2("TCGA data"),
                    
                    fluidRow(
                      box( width = 4, status = "warning",height = "100%",
                           selectInput("ttype",label="Cancer type",choices = tty_list,selected = "Choose cancer type"),
                           checkboxGroupInput("size_type",label="Panel type",choices = list("WES" = 1, "F1CDX (simulated)" = 2,"MSK (simulated)" = 3),selected=1,inline = F),
                           sliderInput("range",label="Y-axis range",min = 0, max = 1000, value = c(0, 100),ticks=T,),numericInput("num","Input TMB of interest",NULL)),
                      box(
                        status = "warning", width = 8,height = "100%",
                        plotOutput("plot_bar")
                      )
                    )
)

#====================================================#
## my data module ####
#====================================================#
bodymy <- tabItem(tabName = "my", value="upset_plot",
                     h2("My data explore"),
                     fluidRow(
                       box(
                         status = "warning", width = 4,
                         selectInput("ttype.usr",label="Cancer type",choices = tty_list,selected = "Choose cancer type"),
                         fileInput("mut", "Input mutation vcf file"),
                         fileInput("panel", "Input panel bed file"),
                         actionButton("submit","Submit"),
                         br(),
                         HTML("<hr> <a href='COAD_test_sap.vcf'> <i class='fa fa-download'> </i> Mutation example data</a>"),
                         HTML("<br> <a href='msk_coding.bed'> <i class='fa fa-download'> </i> Panel example data</a>")
                       ),
                       box(
                         status = "warning", width = 8, height = "100%",
                         plotOutput("plot_cor"),
                         br(),
                         tableOutput("write_num"),
                         plotOutput("plot_pos")
                       
                       )
                     )
                     
)

## Body content ####
ui<-(dashboardPage(
  title="Optimize TMB evaluation",
  skin = "black",
  dashboardHeader(
    # Set height of dashboardHeader
    title = tags$a(tags$img(src='3189359-200.png', height = "100%"),tags$img(src='top_icon.png', height = "100%")),
    titleWidth = 150),
  sidebar,
  dashboardBody(
    tabItems(
      bodyHome,
      #tcga plots
      bodytcga,
      #my Plots
      bodymy
    )
  )
))


panel.type<-list(exome.tmb.list,f1cdx.tmb.list,msk.tmb.list)
# all coding region
exome.bed<-read.table("data/region_bed/exome.final.bed")
gr.exome<-GRanges(seqnames = Rle(exome.bed$V1),ranges = IRanges(exome.bed$V2,exome.bed$V3))
#==================================================
#      server
#==================================================
server<-function(input,output){

  # ===tcga data===
  ylim<-reactive(input$range)
  number<-reactive(input$num)
  output$plot_bar<-renderPlot({
  # order by wes data
    tmp1<-panel.type[[1]]
    y<-tmp1[[input$ttype]]
    y<-as.numeric(y)
    index.order<-order(y)
    
    # deal with check boxplot
    if (is.null(input$size_type) ){return(10)}
    p1 <- 1 %in% input$size_type
    p2 <- 2 %in% input$size_type
    p3 <- 3 %in% input$size_type
    
    tmp1<-panel.type[[1]]
    y1<-as.numeric(tmp1[[input$ttype]])
    y1<-y1[index.order]
    
    tmp1<-panel.type[[2]]
    y2<-as.numeric(tmp1[[input$ttype]])
    y2<-y2[index.order]
    
    tmp1<-panel.type[[3]]
    y3<-as.numeric(tmp1[[input$ttype]])
    y3<-y3[index.order]
    
    n<-length(y1)
    index<-1:n
    col<-c("#66c2a5","#fc8d62","#8da0cb")
    par(mar=c(4,4,0,1),mgp=c(2,.5,0))
    
    if (p1 & p2 & p3){
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col=col[1])
      points(y2,col=col[2])
      points(y3,col=col[3])
    } else if (p1 & p2){
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col=col[1])
      points(y2,col=col[2])
    } else if (p1 & p3){
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col=col[1])
      points(y3,col=col[3])
    } else if (p2 & p3){
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col="white")
      points(y2,col=col[2])
      points(y3,col=col[3])
    } else if (p1){
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col=col[1])
    } else if (p2) {
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col="white")
      points(y2,col=col[2])
    } else if (p3) {
      plot(y1,ylab = "Mutations/Mb",xlab="Sample index",ylim=ylim(),col="white")
      points(y3,col=col[3])
    }
    mtext(paste0(input$ttype," (n=",n,")"),3,cex=1,line=-1)
    tmp1<-abs(y1-number())
    number<-index[which(tmp1==min(tmp1),arr.ind=TRUE)][1]
    abline(v=number,lty=2,col="#bdbdbd")
    abline(h=10,lty=2,col="#bdbdbd")
  })

  # ===my data===
  tumor<-eventReactive(input$submit,{input$ttype.usr})
  p.panel<-observeEvent(input$submit,{
    tumor<-tumor()
    
    # upload panel region
    infile1 <- input$panel
    f1cdx.bed<-read.table(infile1$datapath)
    gr.f1cdx<-GRanges(seqnames = Rle(f1cdx.bed$V1),ranges = IRanges(f1cdx.bed$V2,f1cdx.bed$V3))
    
    # restrict panel within coding regions
    tmp<-findOverlapPairs(gr.f1cdx,gr.exome)
    gr.f1cdx<-pintersect(tmp)
    
    # upload mutations and calculate mutations in the given panel
    infile2 <- input$mut
    # mut.panel.bed<-read.table(infile2$datapath)
    # gr.panel.mut<-GRanges(seqnames = Rle(mut.panel.bed$V1),ranges = IRanges(mut.panel.bed$V2,mut.panel.bed$V2),sample=mut.panel.bed$V6)
    vcf <- readVcf(infile2$datapath, "hg19")
    gr.panel.mut<-rowRanges(vcf)
    seqlevelsStyle(gr.panel.mut) <- "UCSC"
    
    
    # calculate tcga muations in given panel (all mutations)
    mut.exome.bed<-read.table(paste0("data/exome_all_bed/",input$ttype.usr,".bed"))
    gr.exome.mut<-GRanges(seqnames = Rle(mut.exome.bed$V1),ranges = IRanges(mut.exome.bed$V2,mut.exome.bed$V2),sample=mut.exome.bed$V6)
    tmp<-subsetByOverlaps(gr.exome.mut,gr.f1cdx)
    panel.mut.cout<-as.data.frame(table(mcols(tmp)$sample))
    rownames(panel.mut.cout)<-panel.mut.cout$Var1
    
    # caluculate tcga muations in all coding region (non-syn mutations
    exome.mut.cout<-read.table(paste0("data/exome/",input$ttype.usr,".exome.tmb.txt"),row.names = 1)
    
    # prepare correlation data
    sname<-intersect(rownames(panel.mut.cout),rownames(exome.mut.cout))
    final.data<-data.frame(panel.mut.cout[sname,]$Freq,as.numeric(exome.mut.cout[sname,]))
    colnames(final.data)<-c("panel","wes")
    rownames(final.data)<-sname
    
    # model estimate
    region.exome=36.747178
    tmp<-data.frame(x=gr.f1cdx)
    region.panel=sum(tmp[,3]-tmp[,2])/1000000
    x<-final.data$panel/region.panel
    y<-final.data$wes/region.exome
    fit<-lm(y~x)
    max<-max(x,y)
    z<-predict(fit,newdata=data.frame(x=obs.panel))
    if(z<0){z<-"Too few mut to estimate"}
    
    # plot cor of panel and wes
    output$plot_cor<-renderPlot({
      par(mar=c(4,4,0,1),mgp=c(2,.5,0))
      plot(x,y,xlim=c(0,max),ylim=c(0,max),xlab="Panel (mut/Mb)",ylab="WES (mut/Mb)")
      fit<-lm(y~x)
      abline(fit,lty=2,col="#fc8d62")
      p<-signif(summary(fit)$coefficients[,4][2],5)
      tmp<-cor.test(x,y)
      r<-signif(tmp$estimate,3)
      mtext(paste0(tumor()," (R=",r,", P=",p,")"),3,cex=1,line=-1)
    }
    )
    
    # output results
    obs.panel<-sum(countOverlaps(gr.panel.mut, gr.f1cdx))/region.panel
    write.out<-data.frame(PANEL=obs.panel,Predicted_WES=z)
    colnames(write.out)<-c("Observed mutations","Predicted TMB")
    write.out[1,]<-paste(round(write.out[1,],3),"mut/Mb")
    output$write_num<-renderTable(write.out,width="100%",align="c")
    
    # output position
    if (is.character(z)){
      return(10)
    }else{
      tmp1<-panel.type[[1]]
      y1<-as.numeric(tmp1[[input$ttype.usr]])
      index.order<-order(y1)
      y1<-y1[index.order]

      output$plot_pos<-renderPlot({
        par(mar=c(3,1,1.5,3),mgp=c(2,.5,0))
        n<-length(y1)
        label<-paste("TCGA",input$ttype.usr,"index")
        plot(rep(1,n),yaxt="n",ylab="",xlab=label,cex=.7,col="#cbd5e8")
        index=1:n
        
        # where the point is
        tmp1<-abs(y1-z)
        number<-index[which(tmp1==min(tmp1),arr.ind=TRUE)][1]
        #abline(v=number,lty=2,col="#bdbdbd")
        points(number,1,pch=8,col="#e41a1c")
        mtext(paste0("Rank of the predicted TMB in TCGA ",input$ttype.usr," data: ",number,"/",n),3,adj=0,line=0.35)
      }, height = 85, width = 650)
      }
  })
}
shinyApp(ui = ui, server = server)