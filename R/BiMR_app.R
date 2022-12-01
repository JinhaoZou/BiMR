#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# 1. ADD run to run
# 2. Add the distribution together
# 3. add a place people can estimate their estimations together
# 4. write a R Package to understand the estimations


library(shiny)
library(ggplot2)
library(dplyr)
library(MendelianRandomization)
library(ivmodel)
#source("Data_produce04182022.R")


# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Bidirectional Mendelian Randomization Estimations"),

    # generate tabset Panel to seperate the app in two functions: Simulation and Analysis
    tabsetPanel( id = "tabset",
      ########### Simulation ############################
      tabPanel("Simulation",
               sidebarLayout(
                 sidebarPanel(

                   selectInput("scenario", h4("Under which scenario the data is generated?"),
                               choices = list("UMR" = 1, "BMR with infinite feedback cycles" = 2,
                                              "BMR with finite feedback cycles" = 3),
                               selected = 2),

                   selectInput("method", h4("Which method you would like to use for estimate the causal effects?"),
                               choices = list("Ratio" = 1, "BiRatio" = 2, "Liml" = 3, "BiLiml" = 4, "All methods" = 5),
                               selected = 5),


                   numericInput("n_sim",
                                "Number of simulations:",
                                value = 5),

                   sliderInput("bins",
                               "Number of bins:",
                               min = 0,
                               max = 100,
                               value = 20),

                   numericInput("nIVs",
                                "Number of strong IVs",
                                value = 1),

                   numericInput("nIVw",
                                "Number of weak IVs:",
                                value = 0,),


                   numericInput("g12",
                                "Gamma12:",
                                value = 0.9),

                   numericInput("g21",
                                "Gamma21:",
                                value = -0.9,),

                   sliderInput("epsilon",
                               "Error variance:",
                               min = 0,
                               max = 1,
                               value = 0.1, step = 0.1),

                   sliderInput("cy",
                               "Confounder Effects:",
                               min = 0,
                               max = 1,
                               value = 0.1, step = 0.1),

                   sliderInput("C_var",
                               "Confounder variance:",
                               min = 0,
                               max = 2,
                               value = 1, step = 0.2),

                   # Run the code
                   actionButton("run", "Run")

                 ),#end of sidebarPanel
                 mainPanel(
                   plotOutput("dist.g12"),
                   plotOutput("dist.g21"),
                   #textOutput("text1"),
                   htmlOutput("text2"),
                   #tableOutput("table1"),
                   #tableOutput("table2"),
                   tableOutput("table3")
                   #htmlOutput("text3")
                 )#end of mainPanel
               )#end of sidebarLayout
      ),#end of tabPanel; Simulation

      ########### Data Analysis #########################
      tabPanel("Data Analysis",
               sidebarLayout(
                 sidebarPanel(

                   # Input: Select files
                   fileInput("file_X1", "Choose csv file for X1", multiple = FALSE,
                             accept = c(".csv")),

                   fileInput("file_X2", "Choose csv file for X2", multiple = FALSE,
                             accept = c(".csv")),

                   fileInput("file_Y1", "Choose csv file for Y1", multiple = FALSE,
                             accept = c(".csv")),

                   fileInput("file_Y2", "Choose csv file for Y2", multiple = FALSE,
                             accept = c(".csv")),

                   fileInput("file_con", "Choose csv file for confounder", multiple = FALSE,
                             accept = c(".csv")),

                   # Horizontal line
                   tags$hr(),

                   # Run the code
                   actionButton("run2", "Run")

                 ),#end of sidebarPanel
                 mainPanel(

                   tableOutput("AlyResult")

                 )#end of mainPanel
               )#end of sidebarLayout
      )#end of tabPanel; Data Analysis

    )#end of tabsetPanel
)# end of ui


# Define server logic required
server <- function(input, output) {

    ############## Simulation Results ############################
    Method <- c("Ratio",  "BiRatio", "Liml", "BiLiml", "all")
    Causal <- c("uni", "bi_infi", "bi_fl")

    methodr <- reactive({input$method})
    scenarior <- reactive({input$scenario})

    sim_result <- eventReactive(input$run, {Sim_one(method = Method[as.numeric(methodr())], causal = Causal[as.numeric(scenarior())],
                                    g12 = input$g12, g21 = input$g21, nIVs = input$nIVs, nIVw = input$nIVw,
                                    con1_var = input$C_var, bcy1 = input$cy, bcy2 = input$cy,
                                    Ve1 = input$epsilon, Ve2 = input$epsilon,
                                    n_sim = input$n_sim)})


    # The distribution figure for the
    output$dist.g12 <- renderPlot({

        # Plot the hist based on the methods, if there is only one method for estimation
        if(as.numeric(methodr()) <= 4){
        # generate bins based on input$bins from ui.R
          g12 <- sim_result()$Est_slt[1,]
          bins.g12 <- seq(min(g12), max(g12), length.out = input$bins + 1)
          # draw the histogram with the specified number of bins
          hist(g12, breaks = bins.g12, col = 'darkgray', border = 'white', main = "Distribution of estimations of Gamma12")
        }else{
          data.G12 <- data.frame(method = rep(Method[-5], each = input$n_sim),
                                 value = c(sim_result()$Est_slt[1*8-7,], sim_result()$Est_slt[2*8-7,],
                                           sim_result()$Est_slt[3*8-7,], sim_result()$Est_slt[4*8-7,]))
          #print(data.G12)
          data.G12 %>%
          ggplot( aes(x=value, fill=method)) +
            geom_histogram( color="#e9ecef", alpha=1, position = 'dodge') +
            scale_fill_manual(values=c("Ratio" = "darkslategray1", "BiRatio" =  "deepskyblue3",
                                         "Liml"= "salmon1 ", "BiLiML" = "brown3")) +
            labs(fill="", title = "Distribution of estimations of Gamma21") +
            theme(plot.title = element_text(size = 18, hjust = 0.5),
                  panel.grid.major = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "Black"))
        }
    })


    output$dist.g21 <- renderPlot({
      if(as.numeric(methodr()) <= 4){
      # generate bins based on input$bins from ui.R
        g21 <- sim_result()$Est_slt[5,]
        bins.g21 <- seq(min(g21), max(g21), length.out = input$bins + 1)
         # draw the histogram with the specified number of bins
        hist(g21, breaks = bins.g21, col = 'darkgray', border = 'white', main = "Distribution of estimations of Gamma21")
      }else{
        data.G21 <- data.frame(method = rep(Method[-5], each = input$n_sim),
                               value = c(sim_result()$Est_slt[1*8-3,], sim_result()$Est_slt[2*8-3,],
                                         sim_result()$Est_slt[3*8-3,], sim_result()$Est_slt[4*8-3,]))

        data.G21 %>%
          ggplot( aes(x=value, fill=method)) +
          geom_histogram( color="#e9ecef", alpha=1, position = 'dodge') +
          scale_fill_manual(values=c("Ratio" = "darkslategray1", "BiRatio" =  "deepskyblue3",
                                     "Liml"= "salmon1 ", "BiLiML" = "brown3")) +
          labs(fill="", title = "Distribution of estimations of Gamma21") +
          theme(plot.title = element_text(size = 18, hjust = 0.5),
                panel.grid.major = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "Black"))

      }
    })



    output$text2 <- renderUI({
      Scenarios <- c("UMR", "BMR with infinite feedback cycles", "BMR with finite feedback cycles")
      HTML(paste("You have selected scenario as", Scenarios[as.numeric(scenarior())],"<br/>Please make sure gamma12 = 0 when senario is UMR"))
    })

    #output$table1 <- renderTable({sim_result()$result1}, rownames = TRUE)
    #output$table2 <- renderTable({sim_result()$result2}, rownames = TRUE)
    output$table3 <- renderTable({sim_result()$result3}, rownames = TRUE)

    # output$text3 <- renderUI({HTML(paste("Conclution about results: ",
    #                                 "<br/> - For strong SNPs, BiRatio/BiLIML works better than Ratio/LIML",
    #                                 "<br/> - For weak SNPs, BiRatio works better than Ratio, differen between LIML and BiLIML is not obvious",
    #                                 "<br/> - when the g12 and g21 has different direction and large, the difference between Ratio/Liml and BiRatio/BiLiml is large",
    #                                 "<br/>                            ",
    #                                 "<br/> Conclution about parameters:",
    #                                 "<br/> - Variance increase, the estimation bias increase, the difference between Ratio/Liml and BiRatio/BiLiml decreace",
    #                                 "<br/> - Confounder effect and variance increase, the estimation bias increase, the difference between LIML and BiLIML decrease"))
    # })

    ############## Data Analysis #####################

    output$AlyResult <- renderTable({

      req(input$file_X1)
      req(input$file_X2)
      req(input$file_Y1)
      req(input$file_Y2)
      req(input$file_con)
      #X1r <- reactive({input$file_X1})

      X1 <- read.csv2(input$file_X1$datapath)
      X2 <- read.csv2(input$file_X2$datapath)
      Y1 <- read.csv2(input$file_Y1$datapath)
      Y2 <- read.csv2(input$file_Y2$datapath)
      con <- read.csv2(input$file_con$datapath)

      example_data <- list(X1 = as.matrix(X1), Y1 = as.matrix(Y1), X2 = as.matrix(X2), Y2 = as.matrix(Y2), con = as.matrix(con))

      result <- eventReactive(input$run2, {Est_all(data = example_data, method = "all")})
      result2 <- result()
      rownames(result2) <- c("Estimated Gamma12", "S.D. of estimated Gamma12", "Lower bound of CI of Gamma12", "Upper bound of CI of Gamma12",
                             "Estimated Gamma21", "S.D. of estimated Gamma21", "Lower bound of CI of Gamma21", "Upper bound of CI of Gamma21")

      return(result2)


      },

      rownames = TRUE)

}


# Run the application
shinyApp(ui = ui, server = server)
