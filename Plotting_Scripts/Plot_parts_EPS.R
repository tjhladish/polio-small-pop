dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

MS_3500 = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond.out'), col.names=c('time', 'E&D'))
MS_3500_denom = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond_denom.out'), col.names=c('time', 'E&D'))
MS_3500_num = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond_numer.out'), col.names=c('time', 'E&D'))

#plot on top of each other -- change size of numerator/denominator graphs
par(mfrow=c(3,1))
plot(MS_3500_num,type='l',main="Numerator of EPS",xaxt='n')
plot(MS_3500_denom,type='l',main="Denominator of EPS",xaxt='n')
plot(MS_3500,type='l',main="EPS")

MS_20000 = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out'), col.names=c('time', 'E&D'))
MS_20000_denom = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond_denom.out'), col.names=c('time', 'E&D'))
MS_20000_num = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond_numer.out'), col.names=c('time', 'E&D'))

par(mfrow=c(3,1))
plot(MS_20000_num,type='l',main="Numerator of EPS",xaxt='n')
plot(MS_20000_denom,type='l',main="Denominator of EPS",xaxt='n')
plot(MS_20000,type='l',main="EPS")

EE_20000 = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_EE_starting_cond.out'), col.names=c('time', 'E&D'))
EE_20000_denom = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_EE_starting_cond_denom.out'), col.names=c('time', 'E&D'))
EE_20000_num = read.table(paste0(dir,'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_EE_starting_cond_numer.out'), col.names=c('time', 'E&D'))

png(paste0(dir,'EPS_breakdown_N_20000.png'), width=1200, height=800, res=150)
par(mfrow=c(2,1),mar=c(4,4,2,2))
plot(EE_20000_num,type='l',main="Parts of EPS - N=20000",xaxt='n',xlab='',ylab='Interval counts')
lines(EE_20000_denom,type='l',col='red')
legend('topright',legend=c('circulation intervals\n (numerator)','interparalytic case intervals\n (denominator)'),col=c('black','red'),lwd=2,bty='n')
#plot(EE_20000_denom,type='l',main="Denominator of EPS",xaxt='n')
plot(EE_20000,type='l',main="EPS - N=20000",xlab='Years since paralytic case',ylab='Probability of circulation',ylim=c(0,1))
dev.off()
