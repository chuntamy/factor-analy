#################Lead action###############################
library(psych)
library(ppcor)
library(parallel)
library(nFactors)
library(GPArotation)

cor<-read.csv("/Users/chuntamy/Desktop/MVA/PCA/HW 5.9 TEST DATA.csv", row.name="cor")
Numbersubject=200 #設定樣本數，為得知取幾個因素
KMO(cor)  #KMO指數高則可選為因素 
ev<-eigen(cor) #特徵值、向量
ev #取3個變數(eigenvalue>1)

#平行步驟
ap<-parallel(subject=Numbersubject, var=ncol(cor),rep=100,cent=.05) #var=ncol(correlation)為指示變數個數
nS<-nScree(x=ev$values, aparallel=ap$eigen$qevpea) #x=eigenvalues
plotnScree(nS)  #show scree plot and parallel to determine no. of factor)
numberfactor=3 #設定因素個數結果

#------PCF method-------
evalue<-ev$values 
evalue

evaluematrixsquare<-diag(sqrt(evalue))  #eigenvalue對角矩陣(√)
evaluematrixsquare

evector<- -ev$vectors 
evectormatrix<-as.matrix(evector) #負號方便解釋
evectormatrix

evectormatrix_std<-t(evectormatrix%*%evaluematrixsquare) #標準化後轉置
evectormatrix_std
evectormatrix_sol<-t(t(evectormatrix%*%evaluematrixsquare))
pcfmloading<-evectormatrix_sol[,1:numberfactor] #取第一～所需因素個數的pattern loading
pcfmloading #「pattern loading」of PCF
#           [,1]        [,2]         [,3]
#[1,] 0.6982359  0.44430899  0.074717532
#[2,] 0.7211931  0.47820624 -0.018516710
#[3,] 0.6587904  0.54515498  0.007130233
#[4,] 0.5771257  0.55851545 -0.112381592
#[5,] 0.7374168 -0.23028612  0.403782375
#[6,] 0.6787844 -0.37030349  0.221647335
#[7,] 0.7548855 -0.34594313  0.204905627
#[8,] 0.7768674 -0.27011041  0.226959879
#[9,] 0.6489433 -0.21149481 -0.157063607
#[10,] 0.6158057 -0.09733703 -0.176401859
#[11,] 0.5779177 -0.24645072 -0.483258780
#[12,] 0.4037785 -0.25957409 -0.663783682


pcfmcom_f<-as.matrix(pcfmloading^2) #各因素communality
pcfmcom_f 

#[,1]        [,2]         [,3]
#[1,] 0.4875334 0.197410482 5.582710e-03
#[2,] 0.5201195 0.228681209 3.428686e-04
#[3,] 0.4340047 0.297193956 5.084023e-05
#[4,] 0.3330741 0.311939507 1.262962e-02
#[5,] 0.5437836 0.053031698 1.630402e-01
#[6,] 0.4607483 0.137124672 4.912754e-02
#[7,] 0.5698521 0.119676646 4.198632e-02
#[8,] 0.6035229 0.072959634 5.151079e-02
#[9,] 0.4211274 0.044730053 2.466898e-02
#[10,] 0.3792167 0.009474497 3.111762e-02
#[11,] 0.3339889 0.060737956 2.335390e-01
#[12,] 0.1630371 0.067378707 4.406088e-01

pcfmcom<-rowSums(pcfmcom_f) #各指示變數communality
pcfmcom  
#[1] 0.6905266 0.7491436 0.7312495 0.6576432 0.7598555 0.6470005 0.7315151 0.7279933
#[9] 0.4905264 0.4198088 0.6282659 0.6710246

sum(pcfmcom)  #total communality(將各因素communality加總)
# 7.904553

ncol(cor)-sum(pcfmcom)  #total unique variance 
#4.095447
1-pcfmcom #unique variance 
#[1] 0.3094734 0.2508564 0.2687505 0.3423568 0.2401445 0.3529995 0.2684849 0.2720067
#[9] 0.5094736 0.5801912 0.3717341 0.3289754

pcfmloading%*%t(pcfmloading)  #再生成相關矩陣
pcfmcomres<-as.matrix(cor-pcfmloading%*%t(pcfmloading))  #殘差矩陣(原相關矩陣－再生成矩陣)
pcfmcomres #殘差矩陣
sqrt((sum(pcfmcomres^2)-sum(diag(pcfmcomres)^2))/(ncol(cor)*(ncol(cor)-1))) #RMSR
# 0.06570678

#------PCF+varimax method-------
pcfv<-principal(cor, nfactors=numberfactor, rotate="varimax")
pcfv$loading
pcfvloading<-pcfv$loading[,1:3]
plot(pcfvloading)
pcfvcom_f<-matrix(pcfv$loading^2,nrow=ncol(cor), ncol=pcfv$factors)
pcfvcom_f   #各因素communality
#             [,1]        [,2]        [,3]
 #[1,] 0.088743225 0.595288531 0.006494814
 #[2,] 0.061564608 0.662396845 0.025182153
 #[3,] 0.034161623 0.689049153 0.008038766
 #[4,] 0.003534367 0.630464452 0.023644372
 #[5,] 0.697278764 0.061516816 0.001059896
 #[6,] 0.592697520 0.010700873 0.043602112
 #[7,] 0.642303519 0.028422696 0.060788860
 #[8,] 0.624185058 0.058684756 0.045123495
 #[9,] 0.222905596 0.047975607 0.219645196
#[10,] 0.145839013 0.084844997 0.189124777
#[11,] 0.073213630 0.023612506 0.531439722
#[12,] 0.003951109 0.001779381 0.665294102

pcfv$communality #各指示變數communality
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.6905266 0.7491436 0.7312495 0.6576432 0.7598555 0.6470005 0.7315151 0.7279933 
#t9       t10       t11       t12 
#0.4905264 0.4198088 0.6282659 0.6710246 

sum(pcfv$communality)   #show total communality
#7.904553

ncol(cor)-sum(pcfv$communality)    #total variance of unique
# 4.095447

pcfv$uniquenesses   # variance of unique
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.3094734 0.2508564 0.2687505 0.3423568 0.2401445 0.3529995 0.2684849 0.2720067 
#t9       t10       t11       t12 
#0.5094736 0.5801912 0.3717341 0.3289754 

pcfv$rms #RMSR
#0.06570678

#------PCF+quartimax method-------
pcfq<-principal(cor, nfactors=numberfactor, rotate="quartimax", scores=T, residuals=T)
pcfq$loading  #pattern loading of PCF
pcfqloading<-pcfq$loading[,1:3]
plot(pcfqloading)   #plot pattern loading of PCF
pcfqcom_f<-matrix(pcfq$loading^2,nrow=ncol(cor), ncol=pcfq$factors)
pcfqcom_f   #各因素communality
#             [,1]        [,2]         [,3]
#[1,] 0.105211853 0.585314546 1.715002e-07
#[2,] 0.081917008 0.660517690 6.708908e-03
#[3,] 0.046173779 0.684641567 4.341962e-04
#[4,] 0.009619855 0.637541578 1.048176e-02
#[5,] 0.702076261 0.049954593 7.824622e-03
#[6,] 0.628477757 0.008086201 1.043655e-02
#[7,] 0.689798070 0.024291736 1.742527e-02
#[8,] 0.666818354 0.052006083 9.168872e-03
#[9,] 0.285427975 0.051699884 1.533985e-01
#[10,] 0.195499787 0.090362617 1.339464e-01
#[11,] 0.134534204 0.032848544 4.608831e-01
#[12,] 0.028725724 0.006503956 6.357949e-01

pcfq$communality    #各指示變數communality
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.6905266 0.7491436 0.7312495 0.6576432 0.7598555 0.6470005 0.7315151 0.7279933 
#t9       t10       t11       t12 
#0.4905264 0.4198088 0.6282659 0.6710246 

sum(pcfq$communality)   #total communality
#7.904553

ncol(cor)-sum(pcfq$communality)    #total variance of unique
#4.095447

pcfq$uniquenesses   #show variance of unique
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.3094734 0.2508564 0.2687505 0.3423568 0.2401445 0.3529995 0.2684849 0.2720067 
#t9       t10       t11       t12 
#0.5094736 0.5801912 0.3717341 0.3289754 

pcfq$rms #show RMSR
#0.06570678

#------PAF method-------
paf<-fa(cor, nfactors=numberfactor, rotate="none", scores=T, residuals=T, nin.err=0.001, max.iter=30)
paf$loading   #show pattern loading of PAF
pafloading<-paf$loading[,1:3]
plot(pafloading)   #plot pattern loading of paf
pafcom_f<-matrix(paf$loading^2,nrow=ncol(cor), ncol=paf$factors)
pafcom_f   #各因素communality
#           [,1]        [,2]         [,3]
#[1,] 0.4552439 0.143119541 6.039897e-03
#[2,] 0.5048000 0.198747677 1.046702e-05
#[3,] 0.4117015 0.229915899 1.855311e-03
#[4,] 0.2962498 0.192182745 2.885860e-03
#[5,] 0.5274211 0.060723561 1.075071e-01
#[6,] 0.4128993 0.099097604 7.220797e-03
#[7,] 0.5470764 0.117036871 1.646885e-02
#[8,] 0.5792026 0.074085890 2.374089e-02
#[9,] 0.3588609 0.022692074 3.151518e-02
#[10,] 0.3143940 0.003923366 3.016708e-02
#[11,] 0.2943736 0.032379003 1.804845e-01
#[12,] 0.1297769 0.022346885 1.181126e-01

paf$communality    #各指示變數communality
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.6044033 0.7035582 0.6434727 0.4913185 0.6956518 0.5192177 0.6805821 0.6770294 
#t9       t10       t11       t12 
#0.4130682 0.3484844 0.5072371 0.2702363 

sum(paf$communality)   #show total communality
#6.55426

ncol(cor)-sum(paf$communality)    #total variance of unique
#5.44574

paf$uniquenesses   #show variance of unique
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.3955967 0.2964418 0.3565273 0.5086815 0.3043482 0.4807823 0.3194179 0.3229706 
#t9       t10       t11       t12 
#0.5869318 0.6515156 0.4927629 0.7297637 

paf$values         #communality of all factor
paf$Vaccounted     #communality of each factor
paf$residual   #show resdiual matrix 
paf$rms        #RMSR 
#0.02481592

#operation procedure of reproduced correlation matrix, residual matrix, RMSR
pafloadingmatrix<-matrix(paf$loading, nrow=ncol(cor), ncol=paf$factors)
pafloadingmatrix     #show pattern loading matrix
pafloadingmatrixtran<-t(pafloadingmatrix)
pafloadingmatrixtran  #show pattern loading matrix transpose
pafreproducecorr<-pafloadingmatrix%*%pafloadingmatrixtran  
pafreproducecorr    #show reproduced correlation matrix
correlationmatrix<-as.matrix(cor[,1:ncol(cor)])
correlationmatrix   #show correlation matrix
pafresidual<- correlationmatrix - pafreproducecorr
pafresidual    #show residual matrix
pafRMSR<-sqrt((sum(pafresidual^2)-sum(diag(pafresidual)^2))/(ncol(cor)*(ncol(cor)-1)))
pafRMSR
#0.02481592

#------PAF+varimax method-------
pafv<-fa(cor, nfactors=numberfactor, rotate="varimax", scores=T, residuals=T, nin.err=0.001, max.iter=30)
pafv$rot.mat   #orthgonal transformation matrix
pafv$loading   #show pattern loading of PAF+v
pafvloading<-pafv$loading[,1:3]
plot(pafvloading)   #plot pattern loading of paf
pafvcom_f<-matrix(pafv$loading^2,nrow=ncol(cor), ncol=pafv$factors)
pafvcom_f   #show communality of each correlation with factor

pafv$communality    #show communality of each indicator
#       t1        t2        t3        t4        t5        t6        t7        t8 
#0.6044033 0.7035582 0.6434727 0.4913185 0.6956518 0.5192177 0.6805821 0.6770294 
#t9       t10       t11       t12 
#0.4130682 0.3484844 0.5072371 0.2702363

sum(pafv$communality)   #show total communality
#6.55426
ncol(cor)-sum(pafv$communality)    #show total variance of unique
#5.44574
pafv$uniquenesses   #show variance of unique
pafv$values         #show communality of all factor
pafv$Vaccounted     #show communality of each factor
pafv$residual   #show resdiual matrix

pafv$rms        #show RMSR
#0.02481592

pafv$R2         # regression->R2
#MR1       MR2       MR3 
#0.7942377 0.8219524 0.5933992 

#------PAF+Quartimax method-------
pafq<-fa(cor, nfactors=numberfactor, rotate="quartimax", scores=T, residuals=T, max.iter=30)
pafq$loading   #show pattern loading of PAF+v
pafqloading<-pafq$loading[,1:3]
plot(pafqloading)   #plot pattern loading of paf
pafqcom_f<-matrix(pafq$loading^2,nrow=ncol(cor), ncol=pafq$factors)

pafqcom_f   #communality of each indicator with factor
#            [,1]        [,2]         [,3]
#[1,] 0.12096537 0.482980251 4.577256e-04
#[2,] 0.10375879 0.597626626 2.172754e-03
#[3,] 0.06621896 0.577173675 8.009954e-05
#[4,] 0.03250398 0.453136215 5.678257e-03
#[5,] 0.62782330 0.042503752 2.532472e-02
#[6,] 0.50169611 0.012053775 5.467767e-03
#[7,] 0.65769803 0.020153396 2.730714e-03
#[8,] 0.63227569 0.044357567 3.961342e-04
#[9,] 0.27345080 0.051607435 8.800992e-02
#[10,] 0.19737685 0.076679676 7.442787e-02
#[11,] 0.19180457 0.031865048 2.835675e-01
#[12,] 0.08613776 0.009203552 1.748950e-01

pafq$communality    #show communality of each indicator
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.6044033 0.7035582 0.6434727 0.4913185 0.6956518 0.5192177 0.6805821 0.6770294 
#t9       t10       t11       t12 
#0.4130682 0.3484844 0.5072371 0.2702363 

sum(pafq$communality)   #show total communality
#6.55426

ncol(cor)-sum(pafq$communality)    #total variance of unique
#5.44574

pafq$uniquenesses   #variance of unique
#t1        t2        t3        t4        t5        t6        t7        t8 
#0.3955967 0.2964418 0.3565273 0.5086815 0.3043482 0.4807823 0.3194179 0.3229706 
#t9       t10       t11       t12 
#0.5869318 0.6515156 0.4927629 0.7297637 

pafq$values         #communality of all factor
pafq$Vaccounted     #show communality of each factor
pafq$residual   #show resdiual matrix
pafq$rms        #show RMSR
#0.02481592

pafq$R2         # regression->R2
#      MR1       MR2       MR3 
#0.8658388 0.8184572 0.5252933 

