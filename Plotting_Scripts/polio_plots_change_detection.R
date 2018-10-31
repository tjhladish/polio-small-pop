dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

library(RColorBrewer)

full_detect_3500 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_3500 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_0.5_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_3500 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_0.25_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

full_detect_3500_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_3500_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_0.5_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_3500_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_0.25_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

# full_detect_25000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
# half_detect_25000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_50_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
# quarter_detect_25000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_25_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
# 
# full_detect_25000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
# half_detect_25000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_50_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
# quarter_detect_25000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_25_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))

full_detect_20000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_20000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_0.5_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_20000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_0.25_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

full_detect_20000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_20000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_0.5_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_20000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_0.25_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

full_detect_10000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_10000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_0.5_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_10000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_0.25_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

full_detect_10000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_10000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_0.5_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_10000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_0.25_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))


marker = list(color = brewer.pal(8,"GnBu"))

png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_20000_response_paper.png'), width=1400, height=800, res=150)
plot(full_detect_20000, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Endemic potential statistic (probability of circulation)',xlim=c(0,10),col=marker$color[5])
lines(full_detect_20000_filter,type='l', col=marker$color[5],lty=4)
lines(half_detect_20000,type='l', col=marker$color[6])
lines(half_detect_20000_filter,type='l', col=marker$color[6],lty=4)
lines(quarter_detect_20000,type='l', col=marker$color[7])
lines(quarter_detect_20000_filter,type='l', col=marker$color[7],lty=4)
legend('topright', legend = c('ICA, 100% detection (N=20,000)','no ICA, 100% detection (N=20,000)', 'ICA, 50% detection (N=20,000)','no ICA, 50% detection (N=20,000)',
                              'ICA, 25% detection (N=20,000)', 'no ICA, 25% detection (N=20,000)'),col=c(marker$color[5],marker$color[5],marker$color[6],marker$color[6],marker$color[7],marker$color[7]),lty=c(1,4,1,4,1,4),lwd=2, bty='n' )

dev.off()

png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_3500_response_paper.png'), width=1400, height=800, res=150)
plot(full_detect_3500, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Endemic potential statistic (probability of circulation)',xlim=c(0,10),col=marker$color[5])
lines(full_detect_3500_filter, col=marker$color[5],lty=4)
lines(half_detect_3500, col=marker$color[6])
lines(half_detect_3500_filter, col=marker$color[6],lty=4)
lines(quarter_detect_3500, col=marker$color[7])
lines(quarter_detect_3500_filter, col=marker$color[7],lty=4)
legend('topright', legend = c('ICA, 100% detection (N=3500)','no ICA, 100% detection (N=3500)', 'ICA, 50% detection (N=3500)','no ICA, 50% detection (N=3500)',
                              'ICA, 25% detection (N=3500)', 'no ICA, 25% detection (N=3500)'),col=c(marker$color[5],marker$color[5],marker$color[6],marker$color[6],marker$color[7],marker$color[7]),lty=c(1,4,1,4,1,4),lwd=2, bty='n' )
dev.off()


png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_10000_response_paper.png'), width=1400, height=800, res=150)
plot(full_detect_10000, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Endemic potential statistic (probability of circulation)',xlim=c(0,10),col=marker$color[5])
lines(full_detect_10000_filter,type='l', col=marker$color[5],lty=4)
lines(half_detect_10000,type='l', col=marker$color[6])
lines(half_detect_10000_filter,type='l', col=marker$color[6],lty=4)
lines(quarter_detect_10000,type='l', col=marker$color[7])
lines(quarter_detect_10000_filter,type='l', col=marker$color[7],lty=4)
legend('topright', legend = c('ICA, 100% detection (N=10,000)','no ICA, 100% detection (N=10,000)', 'ICA, 50% detection (N=10,000)','no ICA, 50% detection (N=10,000)',
                              'ICA, 25% detection (N=10,000)', 'no ICA, 25% detection (N=10,000)'),col=c(marker$color[5],marker$color[5],marker$color[6],marker$color[6],marker$color[7],marker$color[7]),lty=c(1,4,1,4,1,4),lwd=2, bty='n' )

dev.off()

absRiskMat_20000_det_1 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_20000_det_1.csv'),header=FALSE)
absRiskMat_20000_det_50 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_20000_det_50.csv'),header=FALSE)
absRiskMat_20000_det_25 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_20000_det_25.csv'),header=FALSE)
absRiskMat_5000_det_1 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_5000.csv'),header=FALSE)
absRiskMat_5000_det_50 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_5000_det_50.csv'),header=FALSE)
absRiskMat_5000_det_25 = read.csv(file=paste0(dir,'absolute_risk_matrix_N_5000_det_25.csv'),header=FALSE)

full_detect_5000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_5000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
half_detect_5000 = read.table(paste0(dir,'e_and_d_value-0pcase_det_0.5_N_5000_beta_135_fast_migrate_0_response_paper.out'), col.names=c('time', 'E&D'))
quarter_detect_5000 = read.table(paste0(dir,'e_and_d_value-0pcase_det_0.25_N_5000_beta_135_fast_migrate_0_response_paper.out'), col.names=c('time', 'E&D'))


png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_20000_response_paper_absRisk.png'), width=1400, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(full_detect_20000, type='l',lty=1, xlab = '',ylab = 'Probability of circulation',xlim=c(0,20),xaxt='n',col=marker$color[5])
lines(half_detect_20000,type='l', col=marker$color[6])
lines(quarter_detect_20000,type='l', col=marker$color[7])
legend('topright', legend = c('100% detection (N=20000)','50% detection (N=20000)',
                              '25% detection (N=20000)'),col=marker$color[5:7],lwd=2, bty='n' )
plot(absRiskMat_20000_det_1[1:(nrow(absRiskMat_20000_det_1)-4),], type='l',lty=1, xlab="Time since last paralytic case (years)",ylab="Absolute risk of circulation",yaxt='n',xlim=c(0,20),col=marker$color[5],ylim=c(-.05,.02))
lines(absRiskMat_20000_det_50,type='l', col=marker$color[6])
lines(absRiskMat_20000_det_25[1:(nrow(absRiskMat_20000_det_25)-1),],type='l', col=marker$color[7])
ticks<-c(-.04,-.03,-.02,-.01,0,.01)
axis(2,at=ticks,labels=ticks)
abline(0,0)
legend('bottomright', legend = c('100% detection (N=20000)','50% detection (N=20000)',
                              '25% detection (N=20000)'),col=marker$color[5:7],lwd=2, bty='n' )
mtext("Time since last paralytic case (years)",side=1,line=2)
dev.off()

png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_5000_response_paper_absRisk.png'), width=1400, height=800, res=150)
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(3,3,1,1),mgp=c(2,1,0))
plot(full_detect_5000, type='l',lty=1, xlab = '',ylab = 'Probability of circulation',xlim=c(0,10),xaxt='n',col=marker$color[5])
lines(half_detect_5000,type='l', col=marker$color[6])
lines(quarter_detect_5000,type='l', col=marker$color[7])
legend('topright', legend = c('100% detection (N=5000)','50% detection (N=5000)',
                              '25% detection (N=5000)'),col=marker$color[5:7],lwd=2, bty='n' )
plot(absRiskMat_5000_det_1[1:(nrow(absRiskMat_5000_det_1)-3),], type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Absolute risk of circulation',xlim=c(0,10),col=marker$color[5],ylim=c(-.05,.05))
lines(absRiskMat_5000_det_50,type='l', col=marker$color[6])
lines(absRiskMat_5000_det_25,type='l', col=marker$color[7])
abline(0,0)
legend('topright', legend = c('100% detection (N=5000)','50% detection (N=5000)',
                              '25% detection (N=5000)'),col=marker$color[5:7],lwd=2, bty='n' )
mtext("Time since last paralytic case (years)",side=1,line=2)
dev.off()

