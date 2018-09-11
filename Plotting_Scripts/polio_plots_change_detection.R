dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

library(RColorBrewer)

full_detect_3500 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_3500.out'), col.names=c('time', 'E&D'))
half_detect_3500 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_50_N_3500_beta_135.out'), col.names=c('time', 'E&D'))
quarter_detect_3500 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_25_N_3500_beta_135.out'), col.names=c('time', 'E&D'))

full_detect_3500_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_3500.out'), col.names=c('time', 'E&D'))
half_detect_3500_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_50_N_3500_beta_135.out'), col.names=c('time', 'E&D'))
quarter_detect_3500_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_25_N_3500_beta_135.out'), col.names=c('time', 'E&D'))

full_detect_25000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
half_detect_25000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_50_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
quarter_detect_25000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_25_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))

full_detect_25000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
half_detect_25000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_50_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))
quarter_detect_25000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_25_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))

full_detect_20000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))
half_detect_20000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_50_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))
quarter_detect_20000 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_25_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))

full_detect_20000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))
half_detect_20000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_50_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))
quarter_detect_20000_filter = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_25_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))

marker = list(color = brewer.pal(8,"GnBu"))

png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_20000.png'), width=1400, height=800, res=150)
plot(full_detect_20000, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='E&D stat (probability of circulation)',xlim=c(0,10),col=marker$color[3])
lines(full_detect_20000_filter,type='l', col=marker$color[3],lty=4)
lines(half_detect_20000,type='l', col=marker$color[4])
lines(half_detect_20000_filter,type='l', col=marker$color[4],lty=4)
lines(quarter_detect_20000,type='l', col=marker$color[5])
lines(quarter_detect_20000_filter,type='l', col=marker$color[5],lty=4)
legend('topright', legend = c('ICA, 100% detection','no ICA, 100% detection', 'ICA, 50% detection','no ICA, 50% detection',
                              'ICA, 25% detection', 'no ICA, 25% detection'),col=c(marker$color[3],marker$color[3],marker$color[4],marker$color[4],marker$color[5],marker$color[5]),lty=c(1,4,1,4,1,4),lwd=2, bty='n' )

dev.off()

png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_3500.png'), width=1400, height=800, res=150)
plot(full_detect_3500, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='E&D stat (probability of circulation)',xlim=c(0,10),col=marker$color[3])
lines(full_detect_3500_filter, col=marker$color[3],lty=4)
lines(half_detect_3500, col=marker$color[4])
lines(half_detect_3500_filter, col=marker$color[4],lty=4)
lines(quarter_detect_3500, col=marker$color[5])
lines(quarter_detect_3500_filter, col=marker$color[5],lty=4)
legend('topright', legend = c('ICA, 100% detection','no ICA, 100% detection', 'ICA, 50% detection','no ICA, 50% detection',
                              'ICA, 25% detection', 'no ICA, 25% detection'),col=c(marker$color[3],marker$color[3],marker$color[4],marker$color[4],marker$color[5],marker$color[5]),lty=c(1,4,1,4,1,4),lwd=2, bty='n' )

dev.off()

png(paste0(dir,'Effect_of_detection_rate_on_ED_statistic_N_25000.png'), width=1400, height=800, res=150)
plot(full_detect_25000, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='E&D stat (probability of circulation)',xlim=c(0,10),col=marker$color[3])
lines(full_detect_25000_filter,type='l', col=marker$color[3],lty=4)
lines(half_detect_25000,type='l', col=marker$color[4])
lines(half_detect_25000_filter,type='l', col=marker$color[4],lty=4)
lines(quarter_detect_25000,type='l', col=marker$color[5])
lines(quarter_detect_25000_filter,type='l', col=marker$color[5],lty=4)
legend('topright', legend = c('ICA, 100% detection','no ICA, 100% detection', 'ICA, 50% detection','no ICA, 50% detection',
                              'ICA, 25% detection', 'no ICA, 25% detection'),col=c(marker$color[3],marker$color[3],marker$color[4],marker$color[4],marker$color[5],marker$color[5]),lty=c(1,4,1,4,1,4),lwd=2, bty='n' )

dev.off()