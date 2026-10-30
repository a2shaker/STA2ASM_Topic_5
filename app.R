# STA2ASM Week 5 — Hypothesis testing: One-sample inference
# Standalone Shiny app; requires shiny and base R.
library(shiny)

PAL <- list(teal="#236A6E", blue="#174F8A", coral="#C0504D", green="#548235",
            rose="#C65D7B", lightteal="#EEF5F5", lightrose="#F5D8E0")
math <- function(x) withMathJax(HTML(paste0("\\(",x,"\\)")))
fmt <- function(x, digits=3) formatC(x, format="f", digits=digits)
shade <- function(x, y, selection, colour, alpha=.4) {
  ids <- which(selection)
  if (!length(ids)) return(invisible(NULL))
  runs <- split(ids, cumsum(c(1, as.integer(diff(ids)!=1))))
  for (r in runs) polygon(c(x[min(r)],x[r],x[max(r)]),c(0,y[r],0),
                          col=adjustcolor(colour, alpha.f=alpha),border=NA)
}
get_t <- function(n, mean, sd, mu0, alpha, alternative) {
  df <- n-1
  stat <- (mean-mu0)/(sd/sqrt(n))
  p <- switch(alternative, two=2*pt(-abs(stat),df), less=pt(stat,df), greater=pt(stat,df,lower.tail=FALSE))
  crit <- switch(alternative, two=c(qt(alpha/2,df),qt(1-alpha/2,df)),
                 less=qt(alpha,df), greater=qt(1-alpha,df))
  ci <- switch(alternative,
               two=mean+c(-1,1)*qt(1-alpha/2,df)*sd/sqrt(n),
               less=c(-Inf,mean+qt(1-alpha,df)*sd/sqrt(n)),
               greater=c(mean-qt(1-alpha,df)*sd/sqrt(n),Inf))
  list(n=n,df=df,stat=stat,p=p,crit=crit,ci=ci,d=(mean-mu0)/sd,mean=mean,sd=sd,mu0=mu0,alpha=alpha,alternative=alternative)
}
# Exact t-test power under a normal population; power.t.test uses the noncentral t.
power_calc <- function(n, delta, sd, alpha, alternative) {
  if (!is.finite(n) || n<=1 || !is.finite(sd) || sd<=0 || !is.finite(delta) || !is.finite(alpha) || alpha<=0 || alpha>=1) return(NA_real_)
  tryCatch(power.t.test(n=n,delta=abs(delta),sd=sd,sig.level=alpha,
                        type="one.sample",alternative=alternative)$power,
           error=function(e) NA_real_)
}
required_n <- function(delta, sd, alpha, alternative, target) {
  if (!is.finite(delta) || abs(delta)<1e-8 || !is.finite(sd) || sd<=0 || !is.finite(target) || target<=alpha || target>=1) return(NA_real_)
  tryCatch(ceiling(power.t.test(power=target,delta=abs(delta),sd=sd,
                                sig.level=alpha,type="one.sample",alternative=alternative)$n),
           error=function(e) NA_real_)
}
ui <- fluidPage(
  tags$head(tags$style(HTML('
    body {font-family:"Helvetica Neue",Arial,sans-serif;color:#333;background:white;}
    h1,h2,h3,h4 {color:#236A6E;font-weight:500;}
    .well {background:#eef5f5;border:1px solid #cfdfe0;border-left:4px solid #4d8386;}
    .navbar-default {background:white;border:0;border-bottom:1px solid #dce6e6;border-radius:0;box-shadow:none;margin-bottom:16px;}
    .navbar-default .navbar-nav>li>a {color:#4a5b5c;padding:12px 14px 9px;border-bottom:3px solid transparent;}
    .navbar-default .navbar-nav>.active>a,.navbar-default .navbar-nav>.active>a:hover {color:#236A6E;font-weight:600;background:transparent;border-bottom:3px solid #236A6E;}
    .btn-primary {background:#236A6E;border-color:#236A6E;}
    .btn-primary:hover {background:#174F54;border-color:#174F54;}
    .result {background:#eef5f5;border-left:4px solid #236A6E;padding:12px 15px;margin:12px 0;}
    .note {background:#EEF8EF;border-left:4px solid #548235;padding:10px 14px;margin-top:12px;}
    .subtle {color:#58696a;font-size:13px;}
    .statgrid {display:flex;flex-wrap:wrap;gap:12px;margin:10px 0;}
    .statcard {flex:1 1 145px;border:1px solid #cfdfe0;border-radius:4px;padding:10px;background:white;}
    .statcard .value {font-size:20px;font-weight:600;color:#236A6E;}
    .statcard .key {font-size:12px;color:#58696a;}
    .well .control-label {font-weight:600;}
    .plotnote {font-size:13px;color:#58696a;margin-bottom:15px;}
    .shiny-input-container {max-width:100%;}
  '))),
  titlePanel("STA2ASM — Week 5: One-sample inference"),
  withMathJax(),
  p(class="subtle",withMathJax(HTML("Explore one-sample \\(t\\)-tests, confidence intervals, power and sample-size planning using the Week 5 readings."))),
  navbarPage(title=NULL,id="main_tab",collapsible=TRUE,
             tabPanel(HTML("One-sample <i>t</i>-test"),value="ttest",
                      fluidRow(column(12,h2("One-sample t-test"),
                                      p(withMathJax(HTML("Change the sample statistics to explore how the test statistic, critical values, rejection regions, \\(p\\)-value and confidence interval change."))))),
                      fluidRow(
                        column(3,wellPanel(
                          numericInput("tt_n",HTML("Sample size \\(n\\)"),30,min=2,max=10000,step=1),
                          numericInput("tt_mean",HTML("Sample mean change \\(\\bar{x}\\) (mmHg)"),-8,step=1),
                          numericInput("tt_sd",HTML("Sample standard deviation \\(s\\) (mmHg)"),15,min=.01,step=1),
                          numericInput("tt_mu0",HTML("Hypothesised mean change \\(\\mu_0\\) (mmHg)"),0,step=1),
                          selectInput("tt_alpha",HTML("Significance level \\(\\alpha\\)"),choices=c("0.01"="0.01","0.05"="0.05","0.10"="0.10"),selected="0.05"),
                          selectInput("tt_alt",HTML("Alternative hypothesis \\(H_1\\)"),choices=c("Two-sided: mean differs from null"="two","Lower-sided: mean less than null"="less","Upper-sided: mean greater than null"="greater")),
                          radioButtons("tt_shading","Show shaded area",choices=c("Rejection region"="rejection","\\(p\\)-value"="p","Both"="both"),selected="both"),
                          checkboxInput("tt_clinical","Show clinical threshold (−5 mmHg)",TRUE),
                          actionButton("tt_reset","Restore reading defaults"))),
                        column(9,
                               h3(withMathJax(HTML("\\(t\\)-distribution and confidence interval"))),
                               plotOutput("tt_combined",height="570px"),
                               uiOutput("tt_results")
                        )
                      )
             ),
             tabPanel("Power and sample size",value="power",
                      fluidRow(column(12,h2("Power and sample size"),
                                      p("Calculate power for a chosen sample size, or find the sample size required for a target power."))),
                      fluidRow(
                        column(3,wellPanel(
                          radioButtons("pw_goal","What would you like to calculate?",
                                       choices=c("Power (given sample size)"="power","Required sample size (given target power)"="n"),selected="power"),
                          conditionalPanel("input.pw_goal == 'power'",numericInput("pw_n",HTML("Sample size \\(n\\)"),12,min=2,max=10000,step=1)),
                          conditionalPanel("input.pw_goal == 'n'",numericInput("pw_target","Target power",.80,min=.06,max=.999,step=.05)),
                          numericInput("pw_delta",HTML("Assumed mean change \\(\\mu_1-\\mu_0\\) (mmHg)"),-8,step=1),
                          numericInput("pw_sd",HTML("Assumed population SD \\(\\sigma\\) (mmHg)"),15,min=.01,step=1),
                          selectInput("pw_alpha",HTML("Significance level \\(\\alpha\\)"),choices=c("0.01"="0.01","0.05"="0.05","0.10"="0.10"),selected="0.05"),
                          div(class="pw-null-hypothesis", style="color: #222; margin: 10px 0 8px;", withMathJax(HTML("Null hypothesis: \\(H_0: \\mu = 0\\) mmHg."))),
                          uiOutput("pw_alternative_ui"),
                          actionButton("pw_reset","Restore reading defaults"),
                          tags$hr(),h4("Simulation"),
                          selectInput("sim_m","Number of simulated samples",choices=c("500"=500,"1,000"=1000,"2,000"=2000),selected=1000)
                        )),
                        column(9,
                               uiOutput("pw_results"),
                               h3("Power curve"),plotOutput("pw_curve",height="320px"),
                               h3("Sampling distributions under the null and alternative"),
                               plotOutput("pw_distributions",height="330px"),
                               div(class="plotnote","This visual uses normal sampling distributions with known population SD. Numerical power calculations use the noncentral t-distribution, so the shaded normal-approximation power may differ, especially for small samples."),
                               h3("Simulation comparison"),uiOutput("sim_results")
                        )
                      )
             )
  )
)
server <- function(input,output,session) {
  observeEvent(input$tt_reset,{
    updateNumericInput(session,"tt_n",value=30)
    updateNumericInput(session,"tt_mean",value=-8)
    updateNumericInput(session,"tt_sd",value=15)
    updateNumericInput(session,"tt_mu0",value=0)
    updateSelectInput(session,"tt_alpha",selected="0.05")
    updateSelectInput(session,"tt_alt",selected="two")
    updateRadioButtons(session,"tt_shading",selected="both")
    updateCheckboxInput(session,"tt_clinical",value=TRUE)
  })
  tt <- reactive({
    req(is.finite(input$tt_n),input$tt_n>=2,is.finite(input$tt_sd),input$tt_sd>0,
        is.finite(input$tt_mean),is.finite(input$tt_mu0))
    get_t(as.integer(input$tt_n),input$tt_mean,input$tt_sd,input$tt_mu0,
          as.numeric(input$tt_alpha),input$tt_alt)
  })
  output$tt_combined <- renderPlot({
    d <- tt(); df <- d$df; ci <- d$ci
    # Separate horizontal scales: t statistic above, mean change below.
    old <- par(no.readonly=TRUE)
    on.exit({layout(1);par(old)}, add=TRUE)
    layout(matrix(c(1,2),ncol=1), heights=c(1.65,1))
    
    # Upper panel: null distribution, rejection region and p-value area.
    lim <- max(4.3,abs(d$stat)*1.17,abs(d$crit)*1.15)
    lim <- min(lim,80)
    x <- seq(-lim,lim,length.out=1800); y <- dt(x,df)
    par(mar=c(4.2,4.7,1.5,1.5))
    plot(x,y,type="n",xlab="t-statistic",ylab="Density",las=1,bty="l",
         ylim=c(0,max(y)*1.18))
    reject <- switch(d$alternative,two=abs(x)>=d$crit[2],less=x<=d$crit,greater=x>=d$crit)
    p_area <- switch(d$alternative,two=abs(x)>=abs(d$stat),less=x<=d$stat,greater=x>=d$stat)
    if(input$tt_shading %in% c("rejection","both")) shade(x,y,reject,PAL$rose,.40)
    if(input$tt_shading %in% c("p","both")) {
      ids <- which(p_area)
      if(length(ids)) {
        runs <- split(ids,cumsum(c(1,as.integer(diff(ids)!=1))))
        for(r in runs) polygon(c(x[min(r)],x[r],x[max(r)]),c(0,y[r],0),
                               density=18,angle=45,col=PAL$green,border=NA)
      }
    }
    lines(x,y,col=PAL$teal,lwd=2.8)
    abline(v=d$crit,col=PAL$rose,lty=2,lwd=1.8)
    abline(v=d$stat,col=PAL$blue,lty=3,lwd=2.3)
    legend("topright",
           legend=expression("Distribution under "*H[0], "Critical value(s)", "Observed "*t,
                             "Rejection region", italic(p)*"-value area"),
           col=c(PAL$teal,PAL$rose,PAL$blue,PAL$rose,PAL$green),
           lty=c(1,2,3,NA,NA),lwd=c(2.8,1.8,2.3,NA,NA),
           pch=c(NA,NA,NA,15,15),pt.cex=c(NA,NA,NA,1.4,1.4),bty="n",cex=.85)
    if(abs(d$stat)>lim) mtext("Observed t lies outside display window",side=3,line=0,cex=.75,col=PAL$blue)
    
    # Lower panel: confidence interval on the original mean-change scale.
    finite <- ci[is.finite(ci)]
    pts <- c(d$mu0,d$mean,finite,if(isTRUE(input$tt_clinical)) -5)
    # Keep the lower null value aligned with t = 0 in the upper panel.
    # Both panels use symmetric horizontal limits and identical left/right margins.
    half_range <- max(abs(pts-d$mu0),na.rm=TRUE)
    half_range <- max(half_range*1.25, 3)
    limits <- d$mu0 + c(-half_range,half_range)
    par(mar=c(5.0,4.7,2.0,1.5))
    plot(NA,xlim=limits,ylim=c(0,1),xlab="Mean change in SBP (mmHg)",
         ylab="",yaxt="n",las=1,bty="n")
    abline(h=.48,col="#CFDFE0",lwd=1)
    abline(v=d$mu0,col=PAL$rose,lty=2,lwd=2)
    if(isTRUE(input$tt_clinical)) abline(v=-5,col=PAL$green,lty=3,lwd=2)
    lo <- if(is.finite(ci[1])) ci[1] else limits[1]
    hi <- if(is.finite(ci[2])) ci[2] else limits[2]
    segments(lo,.48,hi,.48,col=PAL$blue,lwd=3)
    if(is.finite(ci[1])) segments(lo,.40,lo,.56,col=PAL$blue,lwd=2) else arrows(lo+.10*half_range,.48,lo,.48,length=.1,col=PAL$blue,lwd=2)
    if(is.finite(ci[2])) segments(hi,.40,hi,.56,col=PAL$blue,lwd=2) else arrows(hi-.10*half_range,.48,hi,.48,length=.1,col=PAL$blue,lwd=2)
    points(d$mean,.48,pch=21,bg=PAL$teal,cex=1.55)
    mtext("Confidence interval and clinical importance",side=3,adj=0,line=.6,
          cex=1.02,col=PAL$teal,font=2)
    
    # Labels sit beside their reference lines, not on top of them.
    # Positions are chosen in plot coordinates so they remain legible as inputs change.
    offset <- .035 * diff(limits)
    text(d$mean,.74,labels=expression(bold(bar(x))),col=PAL$teal,cex=1.35)
    text(d$mu0-offset,.18,labels=expression(bold(mu[0])),
         adj=c(1,.5),col=PAL$rose,cex=1.35)
    if(isTRUE(input$tt_clinical)) {
      # Put the clinical label on the side with more room, avoiding the reference line.
      right_room <- limits[2]-(-5)
      left_room <- (-5)-limits[1]
      if(right_room >= left_room) {
        text(-5+offset,.87,"Clinical threshold",adj=c(0,.5),
             col=PAL$green,cex=1.02,font=2)
      } else {
        text(-5-offset,.87,"Clinical threshold",adj=c(1,.5),
             col=PAL$green,cex=1.02,font=2)
      }
    }
  })
  output$tt_results <- renderUI({
    d <- tt(); ci <- d$ci
    ci_text <- if(d$alternative=="two") sprintf("[%.2f, %.2f]",ci[1],ci[2]) else if(d$alternative=="less") sprintf("(−∞, %.2f]",ci[2]) else sprintf("[%.2f, ∞)",ci[1])
    withMathJax(div(class="result",
                    div(class="statgrid",
                        div(class="statcard",div(class="key",HTML("\\(t\\)-statistic")),div(class="value",fmt(d$stat,3))),
                        div(class="statcard",div(class="key","Degrees of freedom"),div(class="value",d$df)),
                        div(class="statcard",div(class="key",HTML("\\(p\\)-value")),div(class="value",if(d$p<.0001) "< 0.0001" else fmt(d$p,4))),
                        div(class="statcard",div(class="key",sprintf("%.0f%% confidence interval",100*(1-d$alpha))),div(class="value",ci_text)),
                        div(class="statcard",div(class="key",HTML("Cohen's \\(d\\)")),div(class="value",fmt(d$d,3)))),
                    p(HTML(sprintf("At \\(\\alpha = %.2f\\), the result %s statistically significant.",
                                   d$alpha,if(d$p<d$alpha) "is" else "is not")))
    ))
  })
  # Alternative hypothesis follows the direction of the assumed mean change.
  # Keep the existing selection when the labels are refreshed.
  output$pw_alternative_ui <- renderUI({
    delta <- input$pw_delta
    direction <- if (is.null(delta) || !is.finite(delta) || delta == 0) {
      "\\(H_1: \\mu \\ne 0\\) (choose a nonzero change)"
    } else if (delta < 0) {
      "\\(H_1: \\mu < 0\\)"
    } else {
      "\\(H_1: \\mu > 0\\)"
    }
    selected <- isolate(input$pw_alt)
    if (is.null(selected) || !selected %in% c("two.sided", "one.sided")) selected <- "two.sided"
    withMathJax(radioButtons(
      "pw_alt", "Alternative hypothesis",
      choiceNames = list(
        HTML("Two-sided: \\(H_1: \\mu \\ne 0\\)"),
        HTML(paste0("One-sided: ", direction))
      ),
      choiceValues = c("two.sided", "one.sided"),
      selected = selected
    ))
  })
  
  observeEvent(input$pw_reset,{
    updateRadioButtons(session,"pw_goal",selected="power")
    updateNumericInput(session,"pw_n",value=12)
    updateNumericInput(session,"pw_target",value=.80)
    updateNumericInput(session,"pw_delta",value=-8)
    updateNumericInput(session,"pw_sd",value=15)
    updateSelectInput(session,"pw_alpha",selected="0.05")
    updateRadioButtons(session,"pw_alt",selected="two.sided")
    updateSelectInput(session,"sim_m",selected="1000")
  })
  pw <- reactive({
    req(is.finite(input$pw_delta),is.finite(input$pw_sd),input$pw_sd>0)
    alpha <- as.numeric(input$pw_alpha)
    if(input$pw_goal=="power") {
      req(is.finite(input$pw_n),input$pw_n>=2)
      n <- as.integer(input$pw_n)
      power <- power_calc(n,input$pw_delta,input$pw_sd,alpha,input$pw_alt)
      target <- .8
    } else {
      req(is.finite(input$pw_target),input$pw_target>alpha,input$pw_target<1)
      target <- input$pw_target
      n <- required_n(input$pw_delta,input$pw_sd,alpha,input$pw_alt,target)
      power <- if(is.finite(n)) power_calc(n,input$pw_delta,input$pw_sd,alpha,input$pw_alt) else NA_real_
    }
    list(n=n,power=power,target=target,delta=input$pw_delta,sd=input$pw_sd,alpha=alpha,alternative=input$pw_alt,goal=input$pw_goal)
  })
  output$pw_results <- renderUI({
    d <- pw()
    if(!is.finite(d$n) || !is.finite(d$power)) return(div(class="note","Choose a nonzero assumed mean change and a target power above the significance level to calculate a finite sample size."))
    div(class="result",div(class="statgrid",
                           div(class="statcard",div(class="key",if(d$goal=="power") "Calculated power" else "Required sample size"),div(class="value",if(d$goal=="power") fmt(d$power,3) else format(d$n,big.mark=","))),
                           div(class="statcard",div(class="key",if(d$goal=="power") "Sample size" else "Achieved power"),div(class="value",if(d$goal=="power") d$n else fmt(d$power,3))),
                           div(class="statcard",div(class="key",HTML("Population effect size \\(|\\delta|\\)")),div(class="value",fmt(abs(d$delta)/d$sd,3)))),
        p(class="subtle",withMathJax(HTML("Calculated with <code>power.t.test</code> using the noncentral \\(t\\)-distribution."))))
  })
  output$pw_curve <- renderPlot({
    d <- pw()
    if(!is.finite(d$n)) {plot.new();text(.5,.5,"Enter a nonzero effect and valid settings.");return()}
    upper <- min(3000,max(60,ceiling(d$n*1.8)))
    ns <- unique(round(seq(2,upper,length.out=100)))
    vals <- vapply(ns,function(k) power_calc(k,d$delta,d$sd,d$alpha,d$alternative),numeric(1))
    old <- par(mar=c(4.5,4.7,1.1,1.2));on.exit(par(old))
    plot(ns,vals,type="l",col=PAL$teal,lwd=3,ylim=c(0,1),xlab="Sample size (n)",ylab="Power",las=1,bty="l")
    abline(h=d$target,col=PAL$green,lty=2,lwd=2)
    abline(v=d$n,col=PAL$blue,lty=3,lwd=2)
    points(d$n,d$power,pch=21,bg=PAL$teal,cex=1.5)
    legend("bottomright",c("Power","Target power","Selected / required n"),col=c(PAL$teal,PAL$green,PAL$blue),lty=c(1,2,3),lwd=2,bty="n",cex=.85)
    if(d$n>upper) mtext("Selected n is beyond display range",side=3,line=.1,cex=.75)
  })
  output$pw_distributions <- renderPlot({
    d <- pw()
    if(!is.finite(d$n)) {plot.new();text(.5,.5,"Enter valid power settings.");return()}
    se <- d$sd/sqrt(d$n)
    mu0 <- 0; mu1 <- d$delta
    if(d$alternative=="one.sided") {
      cut <- if(mu1<mu0) qnorm(d$alpha,mu0,se) else qnorm(1-d$alpha,mu0,se)
      reject <- function(x) if(mu1<mu0) x<=cut else x>=cut
    } else {
      cut <- qnorm(c(d$alpha/2,1-d$alpha/2),mu0,se)
      reject <- function(x) x<=cut[1] | x>=cut[2]
    }
    xmin <- min(mu0,mu1)-4.1*se; xmax <- max(mu0,mu1)+4.1*se
    x <- seq(xmin,xmax,length.out=1700)
    null <- dnorm(x,mu0,se); alt <- dnorm(x,mu1,se)
    old <- par(mar=c(4.4,4.7,1.1,1.2));on.exit(par(old))
    plot(x,null,type="n",ylim=c(0,max(c(null,alt))*1.2),xlab="Sample mean change (mmHg)",ylab="Density",las=1,bty="l")
    shade(x,alt,reject(x),PAL$green,.30)
    shade(x,null,reject(x),PAL$rose,.50)
    lines(x,alt,col=PAL$teal,lwd=2.8)
    lines(x,null,col=PAL$rose,lwd=2.8)
    abline(v=cut,col=PAL$blue,lty=2,lwd=1.8)
    legend("topright",c("Null distribution","Alternative distribution","Rejection cutoff","Significance level","Power"),
           col=c(PAL$rose,PAL$teal,PAL$blue,PAL$rose,PAL$green),
           lty=c(1,1,2,NA,NA),lwd=c(2.8,2.8,1.8,NA,NA),
           pch=c(NA,NA,NA,15,15),pt.cex=c(NA,NA,NA,1.4,1.4),bty="n",cex=.75)
  })
  # Recompute the simulation whenever its parameters or repetition count change.
  # A local seed makes the display reproducible and avoids unrelated RNG changes.
  sim <- reactive({
    d <- pw()
    m <- as.integer(input$sim_m)
    req(is.finite(m), m > 0)
    if(!is.finite(d$n) || !is.finite(d$power) || d$n > 5000)
      return(list(error="Simulation requires a finite sample size of at most 5,000."))
    set.seed(42)
    pvals <- numeric(m)
    direction <- if(d$alternative == "two.sided") "two.sided" else if(d$delta < 0) "less" else "greater"
    for(i in seq_len(m)) {
      x <- rnorm(d$n, mean=d$delta, sd=d$sd)
      pvals[i] <- t.test(x, mu=0, alternative=direction)$p.value
    }
    list(d=d,m=m,simulated=mean(pvals < d$alpha))
  })
  output$sim_results <- renderUI({
    r <- sim()
    if(!is.null(r$error)) return(div(class="note",r$error))
    d <- r$d
    withMathJax(div(class="result",
                    div(class="statgrid",
                        div(class="statcard",div(class="key","Simulated power"),div(class="value",fmt(r$simulated,3))),
                        div(class="statcard",div(class="key","Theoretical power"),div(class="value",fmt(d$power,3)))),
                    p(HTML(sprintf("Based on %s simulated samples of size %s, with mean change %.2f mmHg, SD %.2f mmHg, \\(\\alpha = %.2f\\) and a %s test.",
                                   format(r$m,big.mark=","),format(d$n,big.mark=","),d$delta,d$sd,d$alpha,d$alternative))),
                    h4("R code:"),
                    tags$pre(tags$code(sprintf(paste0(
                      "# Settings selected above\n",
                      "n <- %d\n",
                      "mean.change <- %s\n",
                      "sigma <- %s\n",
                      "alpha <- %s\n",
                      "alternative <- \"%s\"\n\n",
                      "%s",
                      "# Theoretical power using the noncentral t-distribution\n",
                      "power.t.test(n = n, delta = abs(mean.change), sd = sigma,\n",
                      "             sig.level = alpha, type = \"one.sample\",\n",
                      "             alternative = alternative)$power\n\n",
                      "# Estimate power by simulation\n",
                      "set.seed(42)\n",
                      "m <- %d\n",
                      "p.values <- numeric(m)\n",
                      "test.direction <- \"%s\"\n",
                      "for (i in 1:m) {\n",
                      "  x <- rnorm(n, mean = mean.change, sd = sigma)\n",
                      "  p.values[i] <- t.test(x, mu = 0,\n",
                      "                         alternative = test.direction)$p.value\n",
                      "}\n",
                      "mean(p.values < alpha)"),
                      d$n,format(d$delta,trim=TRUE),format(d$sd,trim=TRUE),
                      format(d$alpha,trim=TRUE),d$alternative,
                      if(d$goal=="n") sprintf(paste0(
                        "# Required sample size for the chosen target power\n",
                        "target.power <- %s\n",
                        "ceiling(power.t.test(power = target.power, delta = abs(mean.change),\n",
                        "                     sd = sigma, sig.level = alpha,\n",
                        "                     type = \"one.sample\", alternative = alternative)$n)\n\n"),
                        format(d$target,trim=TRUE)) else "",
                      r$m,
                      if(d$alternative=="two.sided") "two.sided" else if(d$delta<0) "less" else "greater")))
    ))
  })
}
shinyApp(ui,server)
