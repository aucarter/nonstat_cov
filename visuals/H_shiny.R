###
# Interactive visualization of the anisotropic correlation function
###
library(shiny)
# Define UI for application that draws an SIS model
ui <- fluidPage(
    # Application title
    titlePanel("Iso-correlation curve"),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("gamma",
                        "gamma",
                        min = 0.01,
                        max = 5,
                        value = 1),
            sliderInput("v_x",
                        "v_x",
                        min = -10,
                        max = 10,
                        value = 0.0),
            sliderInput("v_y",
                        "v_y",
                        min = -10,
                        max = 10,
                        value = 0.0),

            width = 3
        ),
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("HPlot")
        )
    )
)

# Define server logic required to draw plot
server <- function(input, output) {
    ### Define the anistotropy function
    makeH <- function(gamma, v_x, v_y) {
        H <- gamma*diag(2) + c(v_x, v_y) %*% t(c(v_x, v_y))
        return(H)
    }
    
    plotH <- function(input) {
        H <- makeH(input$gamma, input$v_x, input$v_y)
        eigenH <- eigen(H)
        eigVal <- eigenH$values
        eigVec <- eigenH$vectors
        eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
        ctr    <- c(0, 0)
        angles <- seq(0, 2*pi, length.out=200)
        xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
        yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
        ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
        ellRot  <- eigVec %*% t(ellBase)                                          # rotated ellipse
        plot((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, type="l", lwd=2, xlab = "x", ylab = "y")
        matlines(xMat, yMat, lty=1, lwd=2, col="green")
        points(ctr[1], ctr[2], pch=4, col="red", lwd=3)
    }

    output$HPlot <- renderPlot({
        plotH(input)
    }, width = 600, height = 500)
}

# Run the application
shinyApp(ui = ui, server = server)