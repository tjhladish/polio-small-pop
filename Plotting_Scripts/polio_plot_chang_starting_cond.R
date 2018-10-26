dir = '/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/'

read_ed_data <- function(fn, dir) {
  res <- read.table(paste0(dir, fn), col.names=c('time', 'E&D'))
  N <- as.integer(gsub("^.+_N_(\\d+)_.+$","\\1",fn))
  ICA <- gsub("^.+-([01])pcase_.+$","\\1",fn) != "0" # TRUE when not "0"
  EE<-TRUE
  if(grepl("mostly_sus_starting_cond.out",fn)){
    EE<-FALSE #not EE if mostly sus starting cond
  }
  within(res,{
    N = N
    ICA = ICA
    EE = EE
  })
}

fnlist <- c(
  'e_and_d_values-0pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_3500_beta_135_fast_response_paper.out',
  'e_and_d_value-0pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
  'e_and_d_value-1pcase_det_1_N_3500_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
  'e_and_d_value-0pcase_det_1_N_10000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
  'e_and_d_value-1pcase_det_1_N_10000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
  'e_and_d_values-0pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_10000_beta_135_fast_response_paper.out',
  'e_and_d_values-0pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out',
  'e_and_d_values-1pcase_filter_det_1_N_20000_beta_135_fast_response_paper.out',
  'e_and_d_value-0pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out',
  'e_and_d_value-1pcase_det_1_N_20000_beta_135_fast_migrate_0_mostly_sus_starting_cond.out'
)

thing <- lapply(fnlist, read_ed_data, dir=dir)
thing2 <- Reduce(rbind, thing)
require(ggplot2)

for(k in 0:2){
  thing2 <- Reduce(rbind, thing[(4*k + 1):(4*k+4)])
  ggsave(paste0("change_starting_cond_N_",thing2$N[1],".png"),ggplot(thing2) + theme_minimal() +
  aes(x=time, y=E.D, linetype=ICA, color=EE, group=interaction(EE, ICA, N)) +
  geom_line() +
  coord_cartesian(xlim=c(0,5))+
    ggtitle(paste0("Comparing starting conditions, N = ",thing2$N[1]))+
    xlab("Years since paralytic case")+
    ylab("Probability of circulation"))
}

