dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

cens25 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_test_trunc_25000_beta_135_fast_det_1.out'),col.names=c('time', 'E&D'))
cens20 = read.table(paste0(dir,'e_and_d_values-0pcase_filter_test_trunc_20000_beta_135_fast_det_1.out'),col.names=c('time', 'E&D'))
filter_20000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))
filter_25000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast.out'), col.names=c('time', 'E&D'))

png(paste0(dir,'ED_censoring_20000.png'), width=1400, height=800, res=150)
plot(filter_20000_no_pcase, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='E&D stat (probability of circulation)',xlim=c(0,6),ylim=c(0,1))
lines(cens20,col="red")

png(paste0(dir,'ED_censoring_25000.png'), width=1400, height=800, res=150)
plot(filter_25000_no_pcase, type='l',lty=1, xlab='Time since last paralytic case (years)', ylab='E&D stat (probability of circulation)',xlim=c(0,6),ylim=c(0,1))
lines(cens25,col="red")

legend('topright', legend = c('No censoring','Censoring'),col=c('black','red'),lwd=2,bty='n' )
dev.off()