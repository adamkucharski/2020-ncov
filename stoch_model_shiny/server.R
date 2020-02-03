

# Define server logic for slider examples
shinyServer(function(input, output,session) {
  
  numerical_solver <- function(r0, k){
    
    fun <- function (s) (1 + (r0/k)*(1 - s))^(-k) - s
    solutions <- rootSolve::multiroot(fun, c(0, 1))$root
    
    realistic_sol <- min(solutions)
    return(realistic_sol)
    
  }
  
  # observe({
  #   input$R0value
  # })
  
  # output_estimate <- reactive({
  #   
  #   ss <- seq(0.001,0.999,0.001)
  #   calculate_prob <- abs((1 + (input$R0value/input$k)*(1 - ss))^(-input$k) - ss)
  #   
  #   prob <- ss[min(calculate_prob)==calculate_prob]
  #   
  #   prob
  #   
  # })
  # 
  output$plot    <- renderPlot({ 
    
    if(input$k=="SARS-like"){kk=0.16}
    if(input$k=="MERS-like"){kk=0.25}
    if(input$k=="Random-mixing"){kk=1}
    
    # calculate probability of outbreak based on branching process
    ss <- seq(0.001,0.999,0.001)
    calculate_prob <- abs((1 + (input$R0value/kk)*(1 - ss))^(-kk) - ss)
    
    prob <- 1-ss[min(calculate_prob)==calculate_prob]
    
    par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))
    
    n_seq <- seq(0,10,1)
    prob_seq <- 1-(1-prob)^n_seq
    
    plot(n_seq,prob_seq,type="l",ylim=c(0,1),xlim=c(-0.5,10.5),xlab=c("number of independently introduced cases"),ylab="probability of large outbreak",col="white",xaxs="i")
    grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
    points(n_seq,prob_seq,col="blue",pch=19)

    })

  
})








