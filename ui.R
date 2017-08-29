library(shiny)
require(rCharts)
options(RCHART_WIDTH = 800)

## do your own stuff
source("/home/albrecht/Rfun/powerQuantitative.R")

## loads popNames
load("data/info.Rdata")

   
## Define UI for dataset viewer application
shinyUI(pageWithSidebar(
 
  ## Application title
  headerPanel("PBS calculator"),
  sidebarPanel( 
      
      textInput("pw", "Enter password", ""),
      
############### PBS plots
      conditionalPanel(condition="input.conditionedPanels==1", ##
                      
                       actionButton("runPlot", em("Run")),
                       selectInput("together", "Plot windows and single marker values together", choices = c("YES","NO"),selected="NO"),
                       selectInput("plotType", "Plot points (p) or lines (l)", choices = c("p","l"),selected="p"),
                       numericInput("quantileLine","Display x (0-100) % PBS (whole genome) quantile line (red line):", 99.9),
                       numericInput("posPBS", "Get PBS for a specific position", 0)

                       
                       ),
############## PBS table
      conditionalPanel(condition="input.conditionedPanels==2",
                       
                       actionButton("runTable", em("Run"))
                      

                       ),  
          ############## PBS genes table
conditionalPanel(condition="input.conditionedPanels==3",


        actionButton("runPBSgenes", em("Run"))
                       

    ),
        ############## PBS manhattan plot
conditionalPanel(condition="input.conditionedPanels==4",

                 actionButton("runManhattan", em("Run")),
                 numericInput("quantileLine","Display x (0-100) % PBS (whole genome) quantile line (red line):", 99.9)
                 
                 ),
################### whole genome PBS table 
conditionalPanel(condition="input.conditionedPanels==5",
                 
                 actionButton("runWG", em("Run"))
                 
),
############### download results
      conditionalPanel(condition="input.conditionedPanels==6", ##

                 actionButton("runDownload", em("Run"))

                 ),


############### display error
conditionalPanel(condition="input.conditionedPanels==7", ##

                 actionButton("runError", em("Run"))
                 ),



selectInput("pop1", "Choose first population:", 
            choices = popNames,selected=popNames[1]),
selectInput("pop2", "Choose second population:", 
            choices = popNames,selected=popNames[2]),
selectInput("pop3", "Choose outgroup population:", 
            choices = popNames,selected=popNames[3]),
numericInput("chr", "Choose chromosome:", 11),
numericInput("start", "Starting position (Mb):", 60),
numericInput("end", "Ending position (Mb):", 62),
numericInput("winSize", "Window size (Bp):", 50000),
selectInput("ifWindows", "With windows", choices = c("YES","NO"),selected="NO"),
selectInput("FstOnly", "Plot Fst values beeween first and second pop, instead of PBS values", choices = c("YES","NO"),selected="NO"),
textInput("clearData","Write 'clear' to clear shiny data", "NA"),
numericInput("minWin", "Exclude windows or genes (gene PBS) below number SNPs:", 10),
helpText(""),
downloadButton('downloadData', 'Download top 1000 in interval result table'), 
downloadButton('downloadgenePBS', 'Download PBS for genes result table (chr=-1 for whole genome)'), 
helpText("")

  
#################shared panel bottom
   
  ),
   
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
 

   
  mainPanel(

      tabsetPanel(
          tabPanel("PBS Plot",value=1,plotOutput("PBSsinglePlot", height = "1000px")),
          tabPanel("table", value=2,tableOutput("tableTop")),
          tabPanel("PBS genes", value=3,chartOutput("genePBS","datatables")),
          tabPanel("manhattan PBS Plot",value=4,plotOutput("manhattanPlot", height = "700px")),
          tabPanel("top WG table ", value=5,chartOutput("tableTopWg","datatables")),
          tabPanel("download whole genome",value=6,chartOutput("wholeGenome","datatables")),
          tabPanel("Where is the error", value=7,textOutput("mistake")),
          
          id = "conditionedPanels"
          
      )
    
 
  ) 
))












