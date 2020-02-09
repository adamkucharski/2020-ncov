

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
    
    if(input$k=="SARS-like"){kk=c(0.16,0.11,0.64)}
    if(input$k=="MERS-like"){kk=c(0.25,0.09,0.91)} # TO CHECK
    if(input$k=="Random-mixing"){kk=c(1,1,1)}
    if(input$k=="nCoV-like (early estimate)"){kk=c(0.54,0.014,1.75)} # TO CHECK
    
    # calculate probability of outbreak based on branching process
    ss <- seq(0.001,0.999,0.001)
    calculate_prob1 <- abs((1 + (input$R0value/kk[1])*(1 - ss))^(-kk[1]) - ss)
    calculate_prob2 <- abs((1 + (input$R0value/kk[2])*(1 - ss))^(-kk[2]) - ss)
    calculate_prob3 <- abs((1 + (input$R0value/kk[3])*(1 - ss))^(-kk[3]) - ss)
    
    prob1 <- 1-ss[min(calculate_prob1)==calculate_prob1]
    prob2 <- 1-ss[min(calculate_prob2)==calculate_prob2]
    prob3 <- 1-ss[min(calculate_prob3)==calculate_prob3]
    
    par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))
    
    n_seq <- seq(0,10,1)
    prob_seq1 <- 1-(1-prob1)^n_seq
    prob_seq2 <- 1-(1-prob2)^n_seq
    prob_seq3 <- 1-(1-prob3)^n_seq
    
    plot(n_seq,prob_seq1,type="l",ylim=c(0,1),xlim=c(-0.5,10.5),xlab=c("number of independently introduced cases"),ylab="probability of large outbreak",col="white",xaxs="i")
    grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
    points(n_seq,prob_seq1,col="blue",pch=19)
    for(ii in 1:11){
      lines(c(n_seq[ii],n_seq[ii]),c(prob_seq2[ii],prob_seq3[ii]),col="blue")
    }


    })

  
})








