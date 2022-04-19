KM.survival = function(x,DATA,y,NAME,to.excl=NULL,...){
					library(survival)
					
					if(x == 1){
							TYPE = DATA[,which(colnames(DATA)==NAME)]
					}else{
						if(x == 2){
								TYPE = DATA$predict.label.weighted
						}else{
							if(x == 3){
									TYPE = DATA$predict.label.w.ambi.unweighted
							}else{
								if(x == 4){
										TYPE = DATA$predict.label.w.ambi.weighted
								}else{
									if(x == 5){	to.excl = to.excl
											TYPE.1 = DATA[,which(colnames(DATA)==NAME)]
											TYPE = TYPE.1[-which(TYPE.1 == to.excl)]
									}else{
										if(x == 6){
												to.excl = to.excl
												TYPE.1 = DATA[,which(colnames(DATA)==NAME)]
												TYPE = TYPE.1[-to.excl]
										}else{
											if(x == 7){
													TYPE = DATA$LAUREN.histo
											}else{
												TYPE = DATA$LAUREN.histo[-which(DATA$LAUREN.histo == 3)]
											}
										}	
									}
								}
							}
						}
					}
					
					if(x == 5){
							FIT = survfit(Surv(DATA$OS.months[-which(TYPE.1 == to.excl)],DATA$status.of.survival[-which(TYPE.1 == to.excl)]) ~ as.factor(TYPE))
							CHI.test = survdiff(Surv(DATA$OS.months[-which(TYPE.1 == to.excl)],DATA$status.of.survival[-which(TYPE.1 == to.excl)]) ~ as.factor(TYPE))
							COX.test = coxph(Surv(DATA$OS.months[-which(TYPE.1 == to.excl)],DATA$status.of.survival[-which(TYPE.1 == to.excl)]) ~ as.factor(TYPE))
					}else{
						if(x == 6){
								FIT = survfit(Surv(DATA$OS.months[-to.excl],DATA$status.of.survival[-to.excl]) ~ as.factor(TYPE))
								CHI.test = survdiff(Surv(DATA$OS.months[-to.excl],DATA$status.of.survival[-to.excl]) ~ as.factor(TYPE))
								COX.test = coxph(Surv(DATA$OS.months[-to.excl],DATA$status.of.survival[-to.excl]) ~ as.factor(TYPE))
						}else{
							if(x == 8){
									FIT = survfit(Surv(DATA$OS.months[-which(DATA$LAUREN.histo == 3)], DATA$status.of.survival[-which(DATA$LAUREN.histo == 3)]) ~ as.factor(TYPE))
									CHI.test = survdiff(Surv(DATA$OS.months[-which(DATA$LAUREN.histo == 3)], DATA$status.of.survival[-which(DATA$LAUREN.histo == 3)]) ~ as.factor(TYPE))
									COX.test = coxph(Surv(DATA$OS.months[-which(DATA$LAUREN.histo == 3)], DATA$status.of.survival[-which(DATA$LAUREN.histo == 3)]) ~ as.factor(TYPE))
							}else{		
								FIT = survfit(Surv(DATA$OS.months, DATA$status.of.survival) ~ as.factor(TYPE))
								CHI.test = survdiff(Surv(DATA$OS.months, DATA$status.of.survival) ~ as.factor(TYPE))
								COX.test = coxph(Surv(DATA$OS.months, DATA$status.of.survival) ~ as.factor(TYPE))
							}
						}
					}
					
					if(x != 3 & x != 4 & x != 7){
									pdf(paste("Kaplan_meier_",y,".pdf",sep=""))
									plot(FIT, main="Kaplan-Meier estimate", xlab="time (days)", ylab="survival function", col = c("red","blue","green","magenta"))
									legend('topright', levels(as.factor(TYPE)),  lty=1, col= c("red","blue","green","magenta"), bty='n', cex=.75)
                  dev.off()
					}else{
						pdf(paste("Kaplan_meier_basic_",y,".pdf",sep=""))
						plot(FIT, main="Kaplan-Meier estimate", xlab="time (days)", ylab="survival function", col = c("red","blue","green","magenta"))
						legend('topright', levels(as.factor(TYPE)),  lty=1, col= c("red","blue","green","magenta"), bty='n', cex=.75)
            dev.off()
					}
					
					
					sink(file=paste("summary_FIT_",y,".txt",sep=""))
					print(summary(FIT),times=c(0,50,100,150),extend=TRUE)
					sink()
					
					sink(file=paste("summary_COX_test_",y,".txt",sep=""))
					print(summary(COX.test))
					sink()
					
					sink(file=paste("summary_CHI_test_",y,".txt",sep=""))
					print(CHI.test)
					sink()
          
          tempopdir=getwd()
				
          
}
