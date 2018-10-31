dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

MS_5000_0pcase = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_5000_beta_135_fast_migrate_0_MS_starting_cond_response_paper.out'), col.names=c('time', 'E&D'))
absRiskMat_5000_MS = read.csv(file=paste0(dir,'absolute_risk_matrix_N_5000_MS.csv'),header=FALSE)
absRiskMat_5000 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_5000.csv'),header=FALSE)
filter_5000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_5000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
absRiskMat_5000_MS_test = read.csv(file=paste0(dir,'absolute_risk_matrix_N_5000_MS_test.csv'),header=FALSE)

MS_20000_0pcase = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out'), col.names=c('time', 'E&D'))
absRiskMat_20000_MS = read.csv(file=paste0(dir,'absolute_risk_matrix_N_20000_MS.csv'),header=FALSE)
absRiskMat_20000 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_20000.csv'),header=FALSE)
filter_20000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

marker = list(color = brewer.pal(8,"Set2"))

png(paste0(dir,'ICA_vs_Non_ICA_starting_cond_N_5000.png'), width=1200, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(filter_5000_no_pcase, type='l',lty=1, col = marker$color[4],xlab = '',ylab = 'Probability of circulation',xlim=c(0,5),xaxt='n')
lines(MS_5000_0pcase,col=marker$color[5])
legend('topright',legend=c('EE (N=5000)','MS (N=5000)'),col=marker$color[4:5],lwd=2,bty='n')
plot(absRiskMat_5000[1:(nrow(absRiskMat_5000)-3),],type='l',col=marker$color[4],xlab="Time since last paralytic case (years)",ylab="Absolute risk of circulation",yaxt='n',ylim=c(-.03,.1))
lines(absRiskMat_5000_MS_test[1:(nrow(absRiskMat_5000_MS_test)-1),],col=marker$color[5])
abline(0,0)
ticks<-c(-.02,0,.02,.04,.06,.08)
axis(2,at=ticks,labels=ticks)
legend('topright',legend=c('EE (N=5000)','MS (N=5000)'),col=marker$color[4:5],lwd=2,bty='n')
dev.off()

png(paste0(dir,'ICA_vs_Non_ICA_starting_cond_N_20000.png'), width=1200, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(filter_20000_no_pcase, type='l',lty=1, col = marker$color[4],xlab = '',ylab = 'Probability of circulation',xlim=c(0,5),xaxt='n')
lines(MS_20000_0pcase,col=marker$color[5])
legend('topright',legend=c('EE (N=20000)','MS (N=20000)'),col=marker$color[4:5],lwd=2,bty='n')
plot(absRiskMat_20000[1:(nrow(absRiskMat_20000)-4),],type='l',col=marker$color[4],xlab="Time since last paralytic case (years)",ylab="Absolute risk of circulation",yaxt='n',ylim=c(-.03,.1))
lines(absRiskMat_20000_MS,col=marker$color[5])
abline(0,0)
ticks<-c(-.02,0,.02,.04,.06,.08)
axis(2,at=ticks,labels=ticks)
legend('topright',legend=c('EE (N=20000)','MS (N=20000)'),col=marker$color[4:5],lwd=2,bty='n')
dev.off()
# read_ed_data <- function(fn, dir) {
#   res <- read.table(paste0(dir, fn), col.names=c('time', 'E&D'))
#   N <- as.integer(gsub("^.+_N_(\\d+)_.+$","\\1",fn))
#   ICA <- gsub("^.+-([01])pcase_.+$","\\1",fn) != "0" # TRUE when not "0"
#   EE<-TRUE
#   if(grepl("mostly_sus_starting_cond.out",fn)){
#     EE<-FALSE #not EE if mostly sus starting cond
#   }
#   within(res,{
#     N = N
#     ICA = ICA
#     EE = EE
#   })
# }
# 
# fnlist <- c(
#   'e_and_d_values-0pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out',
#   'e_and_d_values-1pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out',
#   'e_and_d_value-0pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
#   'e_and_d_value-1pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
#   'e_and_d_value-0pcase_det_1_N_10000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
#   'e_and_d_value-1pcase_det_1_N_10000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
#   'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out',
#   'e_and_d_values-1pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out',
#   'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out',
#   'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out',
#   'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
#   'e_and_d_value-1pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out'
# )
# 
# thing <- lapply(fnlist, read_ed_data, dir=dir)
# thing2 <- Reduce(rbind, thing)
# require(ggplot2)
# 
# for(k in 0:2){
#   thing2 <- Reduce(rbind, thing[(4*k + 1):(4*k+4)])
#   ggsave(paste0("change_starting_cond_N_",thing2$N[1],".png"),ggplot(thing2) + theme_minimal() +
#   aes(x=time, y=E.D, linetype=ICA, color=EE, group=interaction(EE, ICA, N)) +
#   geom_line() +
#   coord_cartesian(xlim=c(0,5))+
#     ggtitle(paste0("Comparing starting conditions, N = ",thing2$N[1]))+
#     xlab("Years since paralytic case")+
#     ylab("Probability of circulation"))
# }

