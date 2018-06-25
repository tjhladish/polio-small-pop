i1 = read.table("I1_10000.000000,beta_135.000000,detect_rate_1.000000rho_0.200000_ES_stat_multinomial_extinct_test.csv", col.names=1:7300, fill=TRUE)
ir = read.table("Ir_10000.000000,beta_135.000000,detect_rate_1.000000rho_0.200000_ES_stat_multinomial_extinct_test.csv", col.names=1:7300, fill=TRUE)
infecteds = ir+i1

#max(infecteds[1:60,], na.rm=T)
#pdf('Random_sample_of_replicates.pdf', width=11, height=8.5)
png('Random_sample_of_replicates.png', width=2100, height=1200, res=200)
par(mfrow=c(12,5), mar=c(.2,.2,.2,.2))
for (i in 1:60) {
    plot(as.numeric(infecteds[i,]), type='l', xlim=c(0,15*365), axes=F, ylim=c(0,250));
    box();
    abline(v=0:15*365, col='#88888888', lwd=0.5)
}
dev.off()

gt2000 = which(!is.na(infecteds[,2000]))
png('Random_sample_of_long_replicates.png', width=2100, height=1200, res=200)
#pdf('Random_sample_of_long_replicates.pdf', width=11, height=8.5)
par(mfrow=c(12,5), mar=c(.2,.2,.2,.2))
for (i in 1:60) {
    plot(as.numeric(infecteds[gt2000[i],]), type='l', xlim=c(0,15*365), axes=F, ylim=c(0,250));
    box();
    abline(v=0:15*365, col='#88888888', lwd=0.5)
}
dev.off()

ed = read.table('e_and_d_values.out', col.names=c('time', 'E&D'))
ed_filter = read.table('e_and_d_values-pcase_filter.out', col.names=c('time', 'E&D'))
ed_filter2 = read.table('e_and_d_values-2pcase_filter.out', col.names=c('time', 'E&D'))

#pdf('Effect_of_case_filter_on_ED_statistic.pdf', width=11, height=8.5)
png('Effect_of_case_filter_on_ED_statistic.png', width=1400, height=800, res=150)
par(mar=c(4.1,4.1,2,2))
plot(ed, type='l', xlab='Time since last paralytic case (years)', ylab='E&D stat (probability of circulation)')
lines(ed_filter, col='red')
lines(ed_filter2, col='blue')
legend('topright', legend=c('No filtering (N = 10000)', 'Having at least 1 paralytic case (N = 10000)', 'Having at least 2 paralytic cases (N = 10000)'), col=c('black','red','blue'), lwd=2, bty='n')
dev.off()
