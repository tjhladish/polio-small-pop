dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

cens25 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast_trunc_response_paper.out'),col.names=c('time', 'E&D'))
cens20 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_trunc_response_paper.out'),col.names=c('time', 'E&D'))
cens10 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_trunc_response_paper.out'),col.names=c('time', 'E&D'))

filter_20000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
filter_25000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
filter_10000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))

png(paste0(dir,'ED_right_censoring_20000_response_paper.png'), width=1400, height=800, res=150)
plot(cens20, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Endemic potential statistic (probability of circulation)',xlim=c(0,6),ylim=c(0,1))
lines(filter_20000_no_pcase,col="red")
legend('topright', legend = c('Right-hand censoring (N=20,000)','No right-hand censoring (N=20,000)'),col=c('black','red'),lwd=2,bty='n' )
dev.off()

png(paste0(dir,'ED_right_censoring_25000_response_paper.png'), width=1400, height=800, res=150)
plot(cens25, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Endemic potential statistic (probability of circulation)',xlim=c(0,6),ylim=c(0,1))
lines(filter_25000_no_pcase,col="red")
legend('topright', legend = c('Right-hand censoring (N=25,000)','No right-hand censoring (N=25,000)'),col=c('black','red'),lwd=2,bty='n' )
dev.off()

png(paste0(dir,'ED_right_censoring_10000_response_paper.png'), width=1400, height=800, res=150)
plot(cens10, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='Endemic potential statistic (probability of circulation)',xlim=c(0,6),ylim=c(0,1))
lines(filter_10000_no_pcase,col="red")
legend('topright', legend = c('Right-hand censoring (N=10,000)','No right-hand censoring (N=10,000)'),col=c('black','red'),lwd=2,bty='n' )
dev.off()