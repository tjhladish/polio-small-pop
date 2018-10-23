dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'
  
library(RColorBrewer)

# replace this stuff:
# filter_3500_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_3500_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_10000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_10000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_5000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_5000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_5000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_5000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_20000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_20000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_25000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# filter_25000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_25000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
# with:

read_ed_data <- function(fn, dir) {
  res <- read.table(paste0(dir, fn), col.names=c('time', 'E&D'))
  N <- as.integer(gsub("^.+_N_(\\d+)_.+$","\\1",fn))
  ICA <- gsub("^.+-([01])pcase_.+$","\\1",fn) != "0" # TRUE when not "0"
  within(res,{
    N = N
    ICA = ICA
  })
}

# list.files(dir) # TODO make this approach possible

fnlist <- c(
  'e_and_d_values-0pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out',
  'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out',
  'e_and_d_values-0pcase_filter_det_1_N_5000_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_5000_beta_135_fast_response_paper.out',
  'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out',
  'e_and_d_values-0pcase_filter_det_1_N_25000_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_25000_beta_135_fast_response_paper.out'
)

thing <- lapply(fnlist, read_ed_data, dir=dir)
thing2 <- Reduce(rbind, thing)

marker = list(color = brewer.pal(8,"GnBu"))

png(paste0(dir,'Population_size_on_ED_statistic_response_paper.png'), width=1200, height=800, res=150)

require(ggplot2)

logit <- function(p) log(p/(1-p))

ggsave("ica_pop_size.png",
  ggplot(thing2) + theme_minimal() +
    aes(x=time, y=logit(E.D), linetype=ICA, color=log(N), group=interaction(ICA, N)) +
    geom_line() +
    coord_cartesian(xlim=c(0,5)),
  device = png(), height = 6, width = 6
)

for(i in 1:length(thing2)){
  subset(thing2, N==3500 & ICA==FALSE)
}
plot(filter_3500_no_pcase, type='l',lty=1, col = marker$color[3],xlab = 'Time since last paralytic case (years)',ylab = 'Endemic potential statistic (probability of circulation)',xlim=c(0,5))
  lines(filter_3500_1_pcase,col=marker$color[3],lty=4)
  lines(filter_5000_no_pcase,col=marker$color[4],lty=1)
  lines(filter_5000_1_pcase,col=marker$color[4],lty=4)
  lines(filter_10000_no_pcase,col=marker$color[5],lty=1)
  lines(filter_10000_1_pcase,col=marker$color[5],lty=4)
  lines(filter_20000_no_pcase,col=marker$color[6],lty=1)
  lines(filter_20000_1_pcase,col=marker$color[6],lty=4)
  lines(filter_25000_no_pcase,col=marker$color[7],lty=1)
  lines(filter_25000_1_pcase,col=marker$color[7],lty=4)
legend('topright', legend = c('N=3500, ICA','N=3500, no ICA', 'N=5000, ICA','N=5000, no ICA',
'N=10000, ICA','N=10000, no ICA',"N=20000, ICA","N=20000, no ICA",'N=25000, ICA','N=25000, no ICA'),col=c(marker$color[3],marker$color[3],marker$color[4],marker$color[4],marker$color[5],marker$color[5],marker$color[6],marker$color[6],marker$color[7],marker$color[7]),lty=c(1,4,1,4,1,4,1,4,1,4),lwd=2, bty='n' )
dev.off()

filter_22000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_22000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
filter_22000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_22000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
plot(filter_22000_no_pcase,type='l')
lines(filter_22000_1_pcase,col='red')

filter_19000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_19000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
filter_19000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_19000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
plot(filter_19000_no_pcase,type='l')
lines(filter_19000_1_pcase,col='red')

filter_21000_no_pcase = read.table(paste0(dir,'e_and_d_values-0pcase_filter_det_1_N_21000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
filter_21000_1_pcase = read.table(paste0(dir,'e_and_d_values-1pcase_filter_det_1_N_21000_beta_135_fast_response_paper.out'), col.names=c('time', 'E&D'))
plot(filter_21000_no_pcase,type='l')
lines(filter_21000_1_pcase,col='red')


