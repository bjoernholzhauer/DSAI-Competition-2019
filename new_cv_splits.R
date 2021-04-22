library(tidyverse)

adat <- read_csv(file="/files/feat_bjoern_trial.csv")
# fulldat <- read_csv("/shared_data/data/training_data/training_data_2015_split_on_outcome.csv") %>%
#   bind_rows(read_csv("/shared_data/data/test_data_full/testing_phase2_release.csv"))

# fulldat %>%
#   filter(intphaseendyear>1990) %>%
#   mutate(Oncology=1L*(str_detect(tolower(strtherapeuticarea),"oncology"))) %>%
#   dplyr::select(DrugKey, indicationkey, intphaseendyear, Oncology) %>%
#   distinct() %>%
#   group_by(Oncology, intphaseendyear) %>% 
#   summarize(n=n()) %>%
#   ungroup() %>%
#   ggplot(aes(x=intphaseendyear, y=n, col=as_factor(Oncology))) +
#   geom_line()

adat %>%
  filter(!is.na(outcome)) %>%
  group_by(years_since_1999) %>%
  summarize(n=n()) %>%
  arrange(years_since_1999) %>%
  ungroup() %>%
  mutate(cs=cumsum(n),
         csp = cs/sum(n))

set.seed(1234)

# adat %>%
#   mutate( x = dim(filter(.,years_since_1999==13))[1])
          
tmp1 <- adat %>%
  dplyr::select(DrugKey, indicationkey, years_since_1999, outcome) %>%
  filter(!is.na(outcome)) %>%
  group_by(outcome, years_since_1999) %>%
  mutate(cv_split1 = ifelse(years_since_1999<=11, 0, sample( rep(1:3, ceiling(n()/3)), size=n(), replace = F)),
         cv_split2 = ifelse(years_since_1999<=11, 0, sample( rep(1:3, ceiling(n()/3)), size=n(), replace = F)),
         
         cv_split3 = ifelse(years_since_1999<=13, 0, sample( rep(1:2, ceiling(n()/2)), size=n(), replace = F)),
         cv_split4 = ifelse(years_since_1999<=13, 0, sample( rep(1:2, ceiling(n()/2)), size=n(), replace = F)),
         cv_split5 = ifelse(years_since_1999<=13, 0, sample( rep(1:2, ceiling(n()/2)), size=n(), replace = F)),
         
         cv_split6 = ifelse(years_since_1999<=12, 0, sample( rep(1:2, ceiling(n()/2)), size=n(), replace = F)),
         cv_split7 = ifelse(years_since_1999<=12, 0, sample( rep(1:2, ceiling(n()/2)), size=n(), replace = F)),
         cv_split8 = ifelse(years_since_1999<=12, 0, sample( rep(1:2, ceiling(n()/2)), size=n(), replace = F)),
         
         cv_split9 = ifelse(years_since_1999<13, 0, 
                            ifelse(years_since_1999==13, sample( rep(0:1, ceiling(n()/2)), size=n(), replace = F),1)),
         cv_split10 = ifelse(years_since_1999<13, 0, 
                            ifelse(years_since_1999==13, sample( rep(0:1, ceiling(n()/2)), size=n(), replace = F), 1)),
         cv_split11 = ifelse(years_since_1999<13, 0, 
                            ifelse(years_since_1999==13, sample( rep(0:1, ceiling(n()/2)), size=n(), replace = F), 1)),
         
         cv_split12 = ifelse(years_since_1999<8, sample(c(rep(0, 9), 1:5), size=n(), replace = T), sample(1:5, size=n(), replace = T) )) %>% 
  ungroup() %>%
  group_by(DrugKey, indicationkey, years_since_1999, outcome) %>%
  summarize_all(max) %>%
  pivot_longer(cols = starts_with("cv_split"), names_to = "cv_split", values_to = "initfoldno") %>%
  ungroup() %>%
  mutate(splitno = as.numeric(str_extract(cv_split, "[0-9]+")))
  
#table(tmp1$splitno, tmp1$initfoldno)

new_cv_splits <- tibble(splitno = sort(unique(tmp1$splitno))) %>%
  mutate(valid = map(splitno, function(x) sort(unique(tmp1$initfoldno[tmp1$splitno==x & tmp1$initfoldno != 0] ))) ) %>%
  unnest(valid) %>%
  mutate(foldid = row_number()) %>%
  full_join(tmp1, by=c("splitno")) %>%
  ungroup() %>%
  #filter(valid == initfoldno | initfoldno==0) %>%
  mutate(set = ifelse(is.na(outcome), "test", ifelse(initfoldno!=valid, "train", "val"))) %>%
  # group_by(foldid, set) %>%
  # summarize(n=n(), sumdk=sum(DrugKey), sumik=sum(indicationkey), mo=mean(outcome)) %>%
  # print(n=100)
  dplyr::select(DrugKey, indicationkey, foldid, set)

write_csv(new_cv_splits, "/files/new_cv_splits.csv")