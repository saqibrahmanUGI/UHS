
operations<-c('Oesophagectomy','Gastrectomy','Wide.Local.Excision','Mastectomy',
              'Liver resection','Pancreatectomy','Right Hemicolectomy',
              'Anterior.Resection','Appendicectomy','Emergency.laparotomy')

specialties<-c('UGI','Breast','HPB','Colorectal','General')
df<-as.data.frame(1:10000)
colnames(df)[1]<-'id'
df$Operation<-rep(operations,each=1000)
df$Readmission<-1
df$Return.to.theatre<-1
df$Escalation.of.Care<-1
df$Complication<-1
df$Mortality30d<-1
df$Mortality90d<-1

df$Age<-rnorm(n = 10000, mean = 55, sd = 15)
df$BMI<-rnorm(n = 10000, mean = 25, sd = 5)
ASA<-c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
df$ASA<-sample(ASA,nrow(df),replace=T)
df$Gender<-sample(c('Male','Female'),nrow(df),replace=TRUE)
smoking=c('Never','Never','Never','Never','Ex Smoker','ExSmoker','Current Smoker','Current Smoker')
df$Smoking<-sample(smoking,nrow(df),replace=T)
df$Specialty<-rep(specialties,each=2000)

n10<-c(0,0,0,0,0,0,0,0,0,1)
n20<-c(0,0,0,0,0,0,0,0,1,1)
n30<-c(0,0,0,0,0,0,0,1,1,1)
n40<-c(0,0,0,0,0,0,1,1,1,1)
n50<-c(0,0,0,0,0,1,1,1,1,1)
nlist<-list(n10,n20,n30,n40,n50)
sample(n10,nrow(df),replace=TRUE)

for (q in 3:8){
for (i in 1:10){
  df[((i*1000)-999):(i*1000),q]<-sample(sample(nlist,1)[[1]],1000,replace=T)
}
}

outcomes<-colnames(df)[3:8]




proc<-operations[2]
outcome<-outcomes[2]
riskadjustment<-TRUE
riskadjustment<-FALSE



{
data<-subset(df,df$Operation==proc)

benchmarkDF<-data[1:(nrow(data)-100),]
analysisDF<-data[(nrow(data)-99):nrow(data),]
covariates<-colnames(df)[c(9:13)]
benchmarkDF<-droplevels(benchmarkDF)
if (riskadjustment==TRUE){
model <- glm(as.formula(paste(outcome,'~',paste(covariates,collapse='+'))),data= benchmarkDF,family  = binomial)
preds<-predict(model,analysisDF,type='response')
}


###########Choose  variable for analysis - here i have chosen 
###########unplanned admission to critical care
benchoutcome<-benchmarkDF[,outcome]
analysisoutcome<-analysisDF[,outcome]

##benchmark previous performance
benchmark<-mean(benchoutcome)
}
{
  #specify exponential weighting factor - higher values will detect problems earlier, but
  #increased risk of false negative. Should not be higher than 0.1, somewhere between 0.01
  #and 0.1 is probably most appropriate
  lambda<-0.05
  #set control limit at above 95% confidence interval
  L<-1.95996
  
  ##establish control limits
  llim = benchmark - L*sqrt((benchmark)*(1-benchmark)*(lambda/(2-lambda)))
  ulim = benchmark +  L*sqrt((benchmark)*(1-benchmark)*(lambda/(2-lambda)))
  
  ##create trace
  x = benchmark
  EWMA = rep(benchmark,length(analysisoutcome))
  for (i in 1:length(analysisoutcome)) {
    EWMA[i] = (1-lambda)*ifelse(riskadjustment==TRUE,preds[i],x) + lambda*analysisoutcome[i]
    x = EWMA[i]
  }
  
  ###a vector of the sequential weighted proportions is created
  EWMA

  #plot trace (plotly package)
  library(plotly)
  EWMAplot<-plot_ly(x=1:nrow(analysisDF),y=EWMA*100, type="scatter", mode="lines", name="Weighted Average")%>% 
    layout(legend=list(orientation='h'),title="EWMA Chart",xaxis=list(title = ""),yaxis=list(title = "Weighted percentage",hoverformat='.1f',range=c(((min(EWMA)*0.9)*100),ifelse((ulim>max(EWMA)),((ulim*1.1)*100),(max(EWMA*1.25)*100)))))%>%
    add_trace(x=c(1:length(EWMA)), y=ulim*100, type="scatter", mode="lines", color=I("red"),name="Alert Line")%>%
    add_trace(x=c(1:length(EWMA)), y=benchmark*100, type="scatter", mode="lines", color=I("green"),name="Overall Average")%>%plotly::config(displayModeBar=FALSE)%>%
    layout(margin = list(b=100),annotations = 
             list(x = 1, y = -0.3, text = paste('Escalation of Care'), 
                  showarrow = F, xref='paper', yref='paper', 
                  xanchor='right', yanchor='auto', xshift=0, yshift=0.5
             ))
  
}
EWMAplot

{
  #Combinations of odds ratio of negative effect and tolerability set to 
  #calculated risk of false positive using Markov-chain monte carlo simulation
  
  ##set target
  OR<-2
  ##set limit
  limit<-2
  
  z = analysisoutcome
  
  s=(z*(log(OR)-log(1+ifelse(riskadjustment==TRUE,preds,benchmark)*(OR-1)))) -((1-z)*log(1+ifelse(riskadjustment==TRUE,preds,benchmark)*(OR-1)))
  
  CUSUM = rep(NA,length(z))
  x = 0
  
  for (i in 1:length(z)) {
    CUSUM[i] = ifelse(x + s[i]<0,0,x+s[i])
    CUSUM[i] = ifelse(x + s[i]>=limit,0,CUSUM[i])
    x = CUSUM[i]
  }
  

  #plot trace (plotly package)
  library(plotly)
  CUSUMplot<-plot_ly(x=1:nrow(analysisDF),y=CUSUM, type="scatter", mode="lines", name="Trace")%>% 
    layout(legend=list(orientation='h'),title="CUSUM Chart",xaxis=list(title = ""),yaxis=list(title = "Chart Statistic",hoverformat='.1f',range=c(0,ifelse((limit>max(CUSUM)),((limit*1.1)),(max(CUSUM*1.25))))))%>%
    add_trace(x=c(1:length(CUSUM)), y=limit, type="scatter", mode="lines", color=I("red"),name="Alert Line")%>%
    plotly::config(displayModeBar=FALSE)%>%
    layout(margin = list(b=100),annotations = 
             list(x = 1, y = -0.3, text = paste('Escalation of Care'), 
                  showarrow = F, xref='paper', yref='paper', 
                  xanchor='right', yanchor='auto', xshift=0, yshift=0.5
             ))
  
}

CUSUMplot








