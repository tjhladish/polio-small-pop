dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

library(RColorBrewer)

filter_20000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))
filter_20000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast.out'), col.names=c('time', 'E&D'))

marker = list(color = brewer.pal(8,"GnBu"))
png(paste0(dir,'Population_size_on_ED_statistic_20000.png'), width=1200, height=800, res=150)
plot(filter_20000_no_pcase, type='l',lty=1, col = marker$color[5],xlab = 'Time since last paralytic case (years)',ylab = 'Endemic potential statistic (probability of circulation)',xlim=c(0,5))
lines(filter_20000_1_pcase,col=marker$color[5],lty=4)
legend('topright', legend = c('N=20000, ICA','N=20000, no ICA'),col=c(marker$color[5],marker$color[5]),lty=c(1,4),lwd=2, bty='n' )
dev.off()
