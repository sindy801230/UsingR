# 本代码用来计算GM(1,1)灰度模型
# 输入：x 原始数据，adv为外推预测步长
# 输出：actual 原始数据，fit 拟合数据,degree 拟合精度，
#       C 后验差比值，P 小概率误差，predict 外推预测值

library(readxl)
dat <- read_excel('103-107指考PR.xlsx')

# 主函数
GM <- function(x,adv=0) {
  x0 <- x
  k = length(x0)
  # AGO
  x1 = cumsum(x0)
  # construct matrix B & Y
  B = cbind(-0.5*(x1[-k]+x1[-1]),rep(1,times=k-1))
  Y = x0[-1]
  # compute BtB...
  BtB = t(B)%*%B
  BtB.inv = solve(BtB)
  BtY = t(B)%*%Y
  
  # estimate
  alpha = BtB.inv%*%BtY
  
  # 建立预测模型
  
  predict <- function(k) {
    y = (x0[1] - alpha[2]/alpha[1])*exp(-alpha[1]*k)+alpha[2]/alpha[1]
    return(y)
  }
  pre <- sapply(X=0:(k-1),FUN=predict)
  prediff <- c(pre[1],diff(pre))
  # 计算残差
  error <- round(abs(prediff-x0),digits=6)
  emax <- max(error)
  emin <- min(error)
  # 模型评价
  incidence <- function(x) {
    return((emin + 0.5*emax)/(x+0.5*emax))
  }
  degree <- mean(sapply(error,incidence))
  
  s1 <- sqrt(sum((x0-mean(x0))^2)/5)
  s2 <- sqrt(sum((error-mean(error))^2)/5)
  
  C <- s2/s1
  
  e <- abs(error-mean(error))
  p <- length(e<0.6745)/length(e)
  
  result <- list(actual = x0,
                 fit = prediff,
                 degree = degree,
                 C = C,
                 P = p)
  
  # 外推预测第k+adv项
  if (adv > 0) {
    pre.adv <- predict(k+adv-1)-predict(k+adv-2)
    
    result$predict <- pre.adv
  }
  class(result) <- 'GM1.1'
  return(result)
}

# 增加对应gm1.1类的泛型函数 print & plot
print.GM1.1 <- function(mod){
  cat('the result of GM(1,1)\n')
  cat('Actual Data:', '\n',mod$actual ,'\n')
  cat('Fit Data:', '\n',round(mod$fit,2) ,'\n')
  cat('Degree:', round(mod$degree,3) ,'\n')
  cat('C:', round(mod$C,3) ,'\n')
  cat('P:', round(mod$P,3) ,'\n')
  if (!is.null(mod$predict)) {
    cat('Predict Data:', round(mod$predict,2), '\n')
  }
}

plot.GM1.1 <- function(mod,adv=5){
  prex <- numeric(adv)
  x <- mod$actual
  for (k in 1:adv){
    prex[k] <- GM(x,k)$predict    
  }
  
  value = c(x,prex)
  
  res <- data.frame(index = 1:length(value),
                    value = value,
                    type = factor(c(rep(1,length(x)),rep(2,length(prex)))))
  library(ggplot2)
  p <- ggplot(res,aes(x=index,y= value))
  p + geom_point(aes(color=type),size=3)+ 
    geom_path(linetype=2) +
    theme_bw()
}


# 原始数据
#x = c(26.7,31.5,32.8,34.1,35.8,37.5)
yea <- dat[grep('106',dat$學年度),]
x <- yea$pr[1:10]

# 预测第7项
res <- GM(x,length(x)+3)
print(res)
plot(res,3)

#方法二####
#编写应用于R软件的GM(1,1)模型
gm11<-function(x0,t){ #x0为输入学列，t為預測個數
  x1<-cumsum(x0) #一次累加生成序列1-AG0序列
  b<-numeric(length(x0)-1)
  n<-length(x0)-1
  for(i in 1:n){ #生成x1的紧邻均值生成序列
    b[i]<--(x1[i]+x1[i+1])/2 
    b} #得序列b，即为x1的紧邻均值生成序列
  D<-numeric(length(x0)-1)
  D[]<-1
  B<-cbind(b,D)
  BT<-t(B)#做逆矩阵
  M<-solve(BT%*%B)
  YN<-numeric(length(x0)-1)
  YN<-x0[2:length(x0)]
  alpha<-M%*%BT%*%YN  #模型的最小二乘估计参数列满足alpha尖
  alpha2<-matrix(alpha,ncol=1)
  a<-alpha2[1]
  u<-alpha2[2]
  cat("GM(1,1)参数估计值：",'\n',"发展系数-a=",-a,"  ","灰色作用量u=",u,'\n','\n') #利用最小二乘法求得参数估计值a,u
  y<-numeric(length(c(1:t)))
  y[1]<-x1[1]
  for(w in 1:(t-1)){  #将a,u的估计值代入时间响应序列函数计算x1拟合序列y
    y[w+1]<-(x1[1]-u/a)*exp(-a*w)+u/a 
  }
  cat("x(1)的模拟值：",'\n',y,'\n')
  xy<-numeric(length(y))
  xy[1]<-y[1]
  for(o in 2:t){ #运用后减运算还原得模型输入序列x0预测序列
    xy[o]<-y[o]-y[o-1] 
  } 
  cat("x(0)的模拟值：",'\n',xy,'\n','\n')                       
  #计算残差e
  e<-numeric(length(x0))
  for(l in 1:length(x0)){
    e[l]<-x0[l]-xy[l] #得残差
  }
  cat("残差：",'\n',e,'\n')
  #计算相对误差
  e2<-numeric(length(x0))
  for(s in 1:length(x0)){
    e2[s]<-(abs(e[s])/x0[s]) #得相对误差
  }
  cat("相对残差：",'\n',e2,'\n','\n')
  cat("残差平方和=",sum(e^2),'\n')
  cat("平均相对误差=",sum(e2)/(length(e2)-1)*100,"%",'\n')
  cat("相对精度=",(1-(sum(e2)/(length(e2)-1)))*100,"%",'\n','\n')
  #后验差比值检验
  avge<-mean(abs(e));esum<-sum((abs(e)-avge)^2);evar=esum/(length(e)-1);se=sqrt(evar)  #计算残差的方差se
  avgx0<-mean(x0);x0sum<-sum((x0-avgx0)^2);x0var=x0sum/(length(x0));sx=sqrt(x0var)  #计算原序列x0的方差sx
  cv<-se/sx  #得验差比值
  cat("后验差比值检验:",'\n',"C值=",cv,'\n')#对后验差比值进行检验，与一般标准进行比较判断预测结果好坏。
  if(cv < 0.35){     
    cat("C值<0.35, GM(1,1)预测精度等级为：好",'\n','\n')
  }else{
    if(cv<0.5){
      cat("C值属于[0.35,0.5), GM(1,1)模型预测精度等级为：合格",'\n','\n')
    }else{
      if(cv<0.65){
        cat("C值属于[0.5,0.65), GM(1,1)模型预测精度等级为：勉强合格",'\n','\n')
      }else{
        cat("C值>=0.65, GM(1,1)模型预测精度等级为：不合格",'\n','\n')
      }
    }
  }
  #画出输入序列x0的预测序列及x0的比较图像
  plot(xy,col='blue',type='b',pch=16,xlab='時間序列',ylab='值')
  #points(x0,col='red',type='b',pch=4)
  legend('topleft',c('預測值'),pch=c(16,4),lty=l,col=c('blue'))
  #legend('topleft',c('预测价格','原始价格'),pch=c(16,4),lty=l,col=c('blue','red'))
}

#a<-c(1.95,2.23,2.4,2.15,1.8,1.95)
yea <- dat[grep('106',dat$學年度),]
a <- yea$pr[1:10]

par(family = "STKaiti")
gm11(a,length(a)+6)


