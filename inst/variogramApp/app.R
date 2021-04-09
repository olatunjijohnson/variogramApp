#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(geoR)
library(ggplot2)
library(dplyr)
# source('~/Documents/Lancaster/Lancaster Job/projects/LOALOA/rcode/loaloa/myvariogramplot.R')
source('myvariogramplot.R')

##### set the names of the correlation function ##########
cor.names <- c("matern", "exponential", "gaussian", "spherical", "circular", 
               "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", 
               "gneiting", "gneiting.matern", "pure.nugget")
choices = data.frame(
    var = cor.names,
    num = 1:length(cor.names)
)
# List of choices for selectInput
mylist <- as.list(choices$num)
# Name it
names(mylist) <- choices$var
#####################################################


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Interactive Variogram Estimation"),

    # Sidebar with a slider input for number of bins 
    
    sidebarLayout(
        sidebarPanel(
                fileInput(inputId = "SDALGCPoutput", label = "Upload a csv file:"), 
                radioButtons("model", 'choose the model', 
                             choices=c("Continuous data" ='continuous', "Prevalence data" = 'prevalence', "Count data" = 'count'), 
                             selected = NULL),
                sliderInput(inputId = "dist",
                            label = "Distance:",
                            min = 1,
                            max = 8,
                            value = 5, step=1), actionButton("change", "Change slider max value"),
                selectInput(
                    inputId = "functions",
                    label = "Correlation functions",
                    choices = mylist,
                    selected = 2
                ),
                selectInput(
                    inputId = "xaxis",
                    label = "X axis",
                    choices = "",
                    selected = 1
                ),
                selectInput(
                    inputId = "yaxis",
                    label = "Y axis",
                    choices = "",
                    selected = 5
                ),
                conditionalPanel(condition = "input.model=='continuous'",
                                 selectInput(
                                     inputId = "y",
                                     label = "continuous variable",
                                     choices = ""
                                 )
                                 ),
                conditionalPanel(condition = "input.model=='prevalence'",
                                 selectInput(
                                     inputId = "p",
                                     label = "Number of postives",
                                     choices = ""
                                 ),
                                 selectInput(
                                     inputId = "m",
                                     label = "Number Examined",
                                     choices = ""
                                 )
                ),
                conditionalPanel(condition = "input.model=='count'",
                                 selectInput(
                                     inputId = "c",
                                     label = "count variable",
                                     choices = ""
                                 ),
                                 selectInput(
                                     inputId = "e",
                                     label = "offset",
                                     choices = ""
                                 )
                                 ),
                selectInput(
                    inputId = "D",
                    label = "Covariates",
                    choices = "",
                    multiple = TRUE
                ) 
                        
        ),

        # Show a plot of the generated distribution
        mainPanel(
            # tabsetPanel(type = "tabs",
            #             tabPanel("Continuous data", value=21, plotOutput("distPlot"), dataTableOutput('mytable')),
            #             tabPanel("Prevalence data", value=22, plotOutput("distPlot"), dataTableOutput('mytable')),
            #             tabPanel("Count data", value=23, plotOutput("distPlot"), dataTableOutput('mytable')),
            #             id= "tabselected", selected = 21
            # )
            plotOutput("distPlot")  , dataTableOutput('mytable')
        )
    )
)

# Define server logic required to draw a histogram
server <- function(session, input, output) {
    options(shiny.maxRequestSize=200*1024^2)    # to increase the size of the input
    data_all <- reactive({
        req(input$SDALGCPoutput)
        dff <- input$SDALGCPoutput
        if (is.null(dff))
            return(NULL)
        x <- readRDS(dff$datapath)
        x
    })
    
    observeEvent(input$change,{
        updateSliderInput(session, "dist", max = 200, step = round(200/7))
    })
    
    # observeEvent(input$input$xaxis, ignoreNULL=T, ignoreInit=T, {
    #     updateSliderInput(session, "change", min = 0, max = as.numeric(dist(cbind(range(selxaxis), range(selyaxis))))/2, 
    #                       step =  min(round(as.numeric(dist(cbind(range(selxaxis), range(selyaxis))))/40),1))
    # })
    
    
    observe({
        datachoices <- names(data_all())
        updateSelectInput(session, "xaxis", choices = datachoices)
        updateSelectInput(session, "yaxis", choices = datachoices)
        updateSelectInput(session, "y", choices = datachoices)
        updateSelectInput(session, "p", choices = datachoices)
        updateSelectInput(session, "m", choices = datachoices)
        updateSelectInput(session, "c", choices = datachoices)
        updateSelectInput(session, "e", choices = datachoices)
        updateSelectInput(session, "D", choices = datachoices)
    })
    
    # output$mytable = renderDataTable({
    #     data_all()[, c(input$xaxis,input$yaxis, input$p, input$m)]
    # })
    
    output$mytable = renderDataTable({
        data_all()[, input$D, drop=FALSE]
    })
    
    
    output$distPlot <- renderPlot({
        df <- data_all()
        if(input$model=='continuous'){
            if(is.null(input$D)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                temp.fit <- lm(df[, input$y] ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                                data = residd, max.dist = input$dist)
            } else{
                xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                temp.fit <- lm(df[, input$y] ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                                data = residd, max.dist = input$dist)
            }

        } else if(input$model=='prevalence'){
            if(is.null(input$D)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(logit ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                                data = residd, max.dist = input$dist)
            } else{
                xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(logit ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                                data = residd, max.dist = input$dist)
            }
        }else{
            if(is.null(input$D)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logc <- log((df[, input$c])/(df[, input$e]))
                temp.fit <- lm(logc ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                                data = residd, max.dist = input$dist)
            } else{
                xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                logc <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(logc ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                                data = residd, max.dist = input$dist)
            }
        }
        myvariogramplot(vario)
    })
    

}

# Run the application 
shinyApp(ui = ui, server = server)
