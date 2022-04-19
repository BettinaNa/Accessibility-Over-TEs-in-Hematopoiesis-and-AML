COX.regression = function(Type = c("univariate","multivariate"),DATA,time.period,survival.status,VAR=1,excl=NULL,project_name, varname,model,special,...){
library(survival)

if(Type == "univariate"){
# UNIVARIATE
	TIME = DATA[,which(colnames(DATA) == time.period)]
	SURVIVAL = DATA[,which(colnames(DATA) == survival.status)]
	VARIABLE = DATA[,which(colnames(DATA) == varname)]

	if(is.null(excl)){
		COX.test = coxph(Surv(TIME, SURVIVAL)~VARIABLE,na.omit(DATA))
	}else{
		COX.test = coxph(Surv(TIME[-excl], SURVIVAL[-excl])~VARIABLE[-excl],na.omit(DATA))
	}

	cox_path = paste( project_name, '/', "summary_COX_test_",varname,"_",special,".txt", sep = '')
	sink(file= cox_path)
		print(colnames(VARIABLE))
		print(summary(COX.test))
	sink()

	ph_path = paste( project_name, '/', "PH_assumption_COX_test_",varname,"_",special,".txt",sep="")
    sink(file = ph_path)
	    print(colnames(VARIABLE))
	    print(cox.zph(COX.test))
    sink()
 # END UNIVARIATE
 }else{
 	# NOT UNIVARIATE
	DATA = na.omit(DATA)
	TIME = DATA[,which(colnames(DATA) == time.period)]
	SURVIVAL = DATA[,which(colnames(DATA) == survival.status)]
	VARIABLE = DATA[,which(match(colnames(DATA),varname,nomatch=0)!=0)]
	print(colnames(VARIABLE))


	if(model == "additive"){
		if(VAR == 2){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]),na.omit(DATA))
		}else{
		if(VAR == 3){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]),na.omit(DATA))
		}else{
		if(VAR == 4){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]),na.omit(DATA))
		}else{
		if(VAR == 5){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]) + (VARIABLE[,5]),na.omit(DATA))
		}else{
		if(VAR == 6){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]) + (VARIABLE[,5]) + (VARIABLE[,6]),na.omit(DATA))
		}else{
		if(VAR == 7){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]) + (VARIABLE[,5]) + (VARIABLE[,6]) + (VARIABLE[,7]),na.omit(DATA))
		}else{
		if(VAR == 8){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]) + (VARIABLE[,5]) + (VARIABLE[,6]) + (VARIABLE[,7]) + (VARIABLE[,8]),na.omit(DATA))
		}else{
		if(VAR == 9){
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]) + (VARIABLE[,5]) + (VARIABLE[,6]) + (VARIABLE[,7]) + (VARIABLE[,8]) + (VARIABLE[,9]),na.omit(DATA))
		}else{
		COX.test = coxph(Surv(TIME, SURVIVAL)~(VARIABLE[,1]) + (VARIABLE[,2]) + (VARIABLE[,3]) + (VARIABLE[,4]) + (VARIABLE[,5]) + (VARIABLE[,6]) + (VARIABLE[,7]) + (VARIABLE[,8]) + (VARIABLE[,9]) + (VARIABLE[,10]),na.omit(DATA))
		}
		}
		}
		}
		}
		}
		}
		} 

		cox_path = paste( project_name, '/', "summary_COX_test_",as.character(VAR),"_",special,".txt", sep = '')
		ph_path = paste( project_name, '/', "PH_assumption_COX_test_",as.character(VAR),"_",special,".txt" ,sep="")

		sink( file = cox_path )
			print(colnames(VARIABLE))
			print(summary(COX.test))
		sink()

	    sink( file = ph_path )
		    print(colnames(VARIABLE))
		    print(cox.zph(COX.test))
	    sink()
		#loglikelihoodScore = COX.test$loglik[2]
		loglikelihoodScore = COX.test$loglik
		log_path = paste(project_name, "/logLLScore_",as.character(VAR),"_",special,".txt",sep="")
		write.table(loglikelihoodScore, file= log_path,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	# END ADDITIVE
	}else{
	print("in progress")
	}
	}
}

