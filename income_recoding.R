#convert income to midpoints of category
dat_sub_no_hisp %>%
  mutate(incimp1_mid=ifelse(incimp1==1, 2499.5,
                            ifelse(incimp1==2, 7499.5,
                                   ifelse(incimp1==3, 12499.5,
                                          ifelse(incimp1==4, 17499.5,
                                                 ifelse(incimp1==5, 22499.5,
                                                        ifelse(incimp1==10, 29999.5,
                                                               ifelse(incimp1==11, 27499.5,
                                                                      ifelse(incimp1==12, 32499.5,
                                                                             ifelse(incimp1==20, 39999.5,
                                                                                    ifelse(incimp1==21, 37499.5,
                                                                                           ifelse(incimp1==22, 42499.5,
                                                                                                  ifelse(incimp1==30, 49999.5,
                                                                                                         ifelse(incimp1==31, 47499.5,
                                                                                                                ifelse(incimp1==32, 52499.5,
                                                                                                                       ifelse(incimp1==40, 59999.5,
                                                                                                                              ifelse(incimp1==41, 57499.5,
                                                                                                                                     ifelse(incimp1==42, 62499.5,
                                                                                                                                            ifelse(incimp1==50, 69999.5,
                                                                                                                                                   ifelse(incimp1==51, 67499.5,
                                                                                                                                                          ifelse(incimp1==52, 72499.5,
                                                                                                                                                                 ifelse(incimp1==60, 75000,
                                                                                                                                                                        ifelse(incimp1==61, 77499.5,
                                                                                                                                                                               ifelse(incimp1==62, 82499.5,
                                                                                                                                                                                      ifelse(incimp1==63, 87499.5,
                                                                                                                                                                                             ifelse(incimp1==64, 92499.5,
                                                                                                                                                                                                    ifelse(incimp1==65, 97499.5,
                                                                                                                                                                                                           ifelse(incimp1==66, 100000,
                                                                                                                                                                                                                  ifelse(incimp1==67, 102499.5,
                                                                                                                                                                                                                         ifelse(incimp1==68, 107499.5,
                                                                                                                                                                                                                                ifelse(incimp1==69, 112499.5, 115000)))))))))))))))))))))))))))))),
         incimp1_mid_infl=incimp1_mid / (cpi/245.1)) -> del#, #inflation adjust midpoint to 2017 dollars (since 2018 NHIS respondents reported their prior-year incomes)
         incimp1_mid_infl_cat=ifelse(incimp1_mid_infl<25000, 1, #75000 is the top code in early years so we have to use that as the highest category; q1 is 29193, q2 is 56933
                                     ifelse(incimp1_mid_infl>=25000 & incimp1_mid_infl<50000, 2,
                                            ifelse(incimp1_mid_infl>=50000 & incimp1_mid_infl<75000, 3, 4)))) -> dat_sub_no_hisp 

dat_sub_no_hisp %>%
  filter(year>=2007) %>%
  mutate(sw_over=ipwpoint(exposure=class,
                          family="multinomial",
                          link="logit",
                          numerator=~1,
                          denominator=~rcs(age, 3) + sex + rcs(int_year_fin, 3), 
                          data = subset(dat_sub_no_hisp, year>=2007),
                          weights=mortwt_f,
                          trace=FALSE,
                          maxit=1000)$ipw.weights,
         sw_f_over=sw_over*mortwt_f)  -> dat_sub_no_hisp_2009

summary(dat_sub_no_hisp_2009$sw_over)


mod_sub <- coxph(Surv(time,dead)~class, robust=T, w=sw_f_over, data=dat_sub_no_hisp_2009)

summary(mod_sub)
