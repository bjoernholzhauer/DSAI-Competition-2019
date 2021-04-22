library(tidyverse)
#library(Amelia)

AICROWD_TRAIN_DATA_PATH <- Sys.getenv(c("AICROWD_TRAIN_DATA_PATH"), unset="/shared_data/data/training_data/training_data_2015_split_on_outcome.csv")
AICROWD_TEST_DATA_PATH <- Sys.getenv(c("AICROWD_TEST_DATA_PATH"), unset="/shared_data/data/test_data_full/testing_phase2_release.csv")
AICROWD_PREDICTIONS_OUTPUT_PATH <- Sys.getenv(c("AICROWD_PREDICTIONS_OUTPUT_PATH"), unset="random_prediction.csv")

create_trial_level_analysis_dataset <- function(){
fulldat <- read_csv(AICROWD_TRAIN_DATA_PATH) %>%
  bind_rows(read_csv(AICROWD_TEST_DATA_PATH)) %>%
  mutate( 
    phaseendyear = case_when(
      as.numeric(intphaseendyear)!=1900 & !is.na(intphaseendyear) ~ as.numeric(intphaseendyear),
      trimws(GenericName) == "9-aminocamptothecin" ~ 1999,
      trimws(GenericName) == "bevacizumab" ~ 2004,
      trimws(GenericName) == "CT-2584" ~ 2001,
      trimws(GenericName) == "levofloxacin" ~ 1997,
      trimws(GenericName) == "paclitaxel" ~ 1993,
      trimws(GenericName) == "raloxifene hydrochloride" ~ 1997,
      trimws(GenericName) == "ramosetron" ~ 1996,
      trimws(GenericName) == "heptaplatin" ~ 1999,
      trimws(GenericName) == "imatinib mesilate" ~ 2001,
      trimws(GenericName) == "SRP-299" ~ 2008,
      trimws(GenericName) == "gemcitabine hydrochloride" ~ 1995,
      trimws(GenericName) == "contusugene ladenovec" ~ 2009,
      TRUE ~ 1995)
  )

# Time since first outcome for the drug (i.e. early or late in life-cycle)
time_since_first = fulldat %>%
  dplyr::select(DrugKey, indicationkey, row_id, phaseendyear) %>%
  group_by(DrugKey) %>%
  mutate(time_since_first_outcome = phaseendyear - min(phaseendyear)) %>%
  ungroup() %>%
  dplyr::select(row_id, time_since_first_outcome)

# relative size
relsize = fulldat %>%
  dplyr::select(DrugKey, indicationkey, strDiseaseType, strtherapeuticarea, row_id, intactualaccrual, inttargetaccrual) %>%
  mutate(size = ifelse( !is.na(intactualaccrual), intactualaccrual,
                            ifelse( !is.na(inttargetaccrual), inttargetaccrual, NA_real_)),
         log_size = ifelse( !is.na(size), log(size+0.5), NA_real_),
         first_ta = str_remove(str_extract(strtherapeuticarea, "[A-Za-z]*\\|?"), "\\|")) %>%
  group_by(DrugKey, indicationkey) %>%
  mutate(phase2trials = n(),
         phase2size = ifelse(is.na(sum( size, na.rm=T  )), 0,sum( size, na.rm=T  )) ) %>%
  ungroup() %>%
  group_by(strDiseaseType) %>%
  mutate(n= sum(!is.na(log_size)),
         rel_log_size_dis = ifelse(n>1 & !is.na(log_size), log_size - (sum(log_size, na.rm=T)-log_size)/(n-1), NA_real_ ),
         ph2_pgm_denom = sum(1/phase2trials),
         rel_ph2_size_dis = log(phase2size+0.5) - log( sum(phase2size/phase2trials) + 0.5*ph2_pgm_denom) + log(ph2_pgm_denom) ) %>%
  group_by(first_ta) %>%
  mutate(n2= sum(!is.na(log_size)),
         rel_log_size_ta = ifelse(n>1 & !is.na(log_size), log_size - (sum(log_size, na.rm=T)-log_size)/(n-1), NA_real_ ),
         rel_log_size_dis=ifelse(n<10, rel_log_size_ta, rel_log_size_dis),
         rel_ph2_size_ta = log(phase2size+0.5) - log( sum(phase2size/phase2trials) + 0.5*sum(1/phase2trials)) + log(sum(1/phase2trials)),
         rel_ph2_size_dis = ifelse(ph2_pgm_denom<4, rel_ph2_size_ta, rel_ph2_size_dis) ) %>%
  ungroup() %>%
  dplyr::select(row_id, rel_log_size_dis, rel_log_size_ta, rel_ph2_size_dis, rel_ph2_size_ta)


# Get information on orphan indications
orphanind = fulldat %>%
  dplyr::select(DrugKey, indicationkey, strDiseaseType) %>%
  distinct() %>%
  left_join( read_rds("/files/orphmatched.rds"), by=c("DrugKey", "strDiseaseType") ) %>%
  ungroup() %>%
  mutate(
    any_orphan_designation = ifelse(is.na(any_orphan_designation), 0, any_orphan_designation),
    orphan_indication = ifelse(is.na(orphan_indication), 0, orphan_indication),
    any_oncology_orphan_designation = ifelse(is.na(any_oncology_orphan_designation), 0, any_oncology_orphan_designation),
    any_nononc_orphan_designation = ifelse(is.na(any_nononc_orphan_designation), 0, any_nononc_orphan_designation)) %>%
  group_by(DrugKey, indicationkey) %>%
  summarize(any_orphan_designation = max(any_orphan_designation),
            orphan_indication = max(orphan_indication),
            any_oncology_orphan_designation = max(any_oncology_orphan_designation),
            any_nononc_orphan_designation = max(any_nononc_orphan_designation)) %>%
  ungroup()
  

# Potential logitoffsets one could use for prediction to adjust for end of follow-up bias
logitoffsets <- fulldat %>%
  mutate(logitoffset = ( (phaseendyear==2017)*0.2546 + (phaseendyear>=2018)*0.8914 ) * (1-str_detect(tolower(fulldat$strtherapeuticarea), "oncology")) +
           ((phaseendyear==2017)*0.1335 + (phaseendyear>=2018)*1.3822) * str_detect(tolower(fulldat$strtherapeuticarea), "oncology") )  %>%
  dplyr::select(row_id, logitoffset)

# Insulin or flu vaccine?
insflu <- fulldat %>%
  mutate(insulin = 1L*str_detect(tolower(GenericName), "insulin"),
         fluvacc = 1L*(str_detect(tolower(GenericName), "influenz") & str_detect(tolower(GenericName), "vaccine"))) %>%
  dplyr::select(row_id, insulin, fluvacc)

# Trial termination score
terminScore = function(data){
  Scores = data %>%
    mutate(termreason = str_split(strTerminationReason, pattern="\\|")) %>%
    unnest(cols="termreason") %>%
    mutate(
      safety = ifelse(is.na(termreason),0,str_detect(tolower(termreason), "adverse eff")*1), # Was there a study termination due to safety?
      bad = case_when( # Rate trial outcomes as far as recorded from 0 = good to 4 = bad
        str_detect(tolower(termreason), "positive outcome") ~ 0,
        is.na(termreason) ~ 1, 
        str_detect(tolower(termreason), "terminated, other") ~ 2,
        str_detect(tolower(termreason), "indeterminat") ~ 2,
        str_detect(tolower(termreason), "enrollment") ~ 2,
        str_detect(tolower(termreason), "priorit") ~ 4,
        str_detect(tolower(termreason), "business") ~ 3,
        str_detect(tolower(termreason), "funding") ~ 3,
        str_detect(tolower(termreason), "lack of eff") ~ 4,
        str_detect(tolower(termreason), "negative outcome") ~ 4,
        str_detect(tolower(termreason), "adverse eff") ~ 4,
        TRUE ~ 5,
      ),
      pct_accrual = ifelse( is.na(intactualaccrual) | is.na(inttargetaccrual),  # Proportion of target accrual that is achieved
                            ifelse(bad>2, 0.1, 
                                   ifelse(bad==1, 0.75, 1)), 
                            pmin(1,intactualaccrual/inttargetaccrual))
    ) %>%
    group_by(row_id) %>% # Summarize on a trial level
    summarize(pct_accrual2=min(pct_accrual),
              safety=max(safety),
              mean_trialendscore=mean(bad),
              max_trialendscore=max(bad),
              trialendscore0=mean( (bad==0)),
              trialendscore1=mean( (bad==1)),
              trialendscore2=mean( (bad==2)),
              trialendscore3=mean( (bad==3)),
              trialendscore4=mean( (bad==4))) %>%
    mutate(accrued90 = (pct_accrual2>=0.90)*1,
           accrued50 = (pct_accrual2>=0.5)*1) %>%
    ungroup() 
  return(Scores)
}

trialcompl <- terminScore(fulldat)

drugtypes <- fulldat %>% 
  mutate(thdes = str_split(TherapyDescription, pattern="\\|")) %>%
  unnest(cols="thdes") %>%
  mutate(drugtype1 = 1L*(thdes %in% c("Anticancer, other", "Anticancer, vaccine")),
         drugtype2 = 1L*(thdes %in% c("Anticancer, immunological")),
         drugtype3 = 1L*(thdes %in% c("Anticancer, antimetabolite")),
         drugtype4 = 1L*(thdes %in% c("Anticancer, alkylating", "Anticancer, antibiotic")),
         drugtype5 = 1L*(thdes %in% c("Anticancer, hormonal", "Anticancer, interferon")),
         drugtype6 = 1L*(thdes %in% c("Immunosuppressant")),
         drugtype7 = 1L*(str_detect(tolower(thdes), "monoclonal") | str_detect(tolower(thdes), "antibody")),
         monoclonal = 1L* str_detect(tolower(thdes), "monoclonal"),
         drugtype8 = 1L*(thdes %in% c("Gene therapy")),
         drugtype9 = 1L*(thdes %in% c("Reformulation, fixed-dose combinations"))) %>%
  #dplyr::select(DrugKey, indicationkey, thdes, outcome) %>%
  group_by(DrugKey, indicationkey) %>%
  summarize(drugtype1=max(drugtype1),
            drugtype2=max(drugtype2),
            drugtype3=max(drugtype3),
            drugtype4=max(drugtype4),
            drugtype5=max(drugtype5),
            drugtype6=max(drugtype6),
            drugtype7=max(drugtype7),
            drugtype8=max(drugtype8),
            drugtype9=max(drugtype9),
            monoclonal=max(monoclonal)) %>%
  ungroup()

basics <- fulldat %>%  mutate(
         sponstype = ifelse(is.na(strSponsorType), NA_real_,
                            1*(str_detect(tolower(strSponsorType),"industry") | 
                          str_detect(tolower(strSponsorType),"contract research")) + 
           2*(str_detect(tolower(strSponsorType),"academic") | str_detect(tolower(strSponsorType),"cooperative") | 
                str_detect(tolower(strSponsorType),"government"))), # General sponsor type
         # One-hot-encoded (dummy variables) for specific sponsor types
         spons1=1*(str_detect(tolower(strSponsorType),"academic")),
         spons2=1*(str_detect(tolower(strSponsorType),"cooperativ")),
         spons3=1*(str_detect(tolower(strSponsorType),"government")),
         spons4=1*(str_detect(tolower(strSponsorType),"top 20")),
         spons5=1*(str_detect(tolower(strSponsorType),"generic")),
         spons6=1*(str_detect(tolower(strSponsorType),"diagnostics") | str_detect(tolower(strSponsorType),"device")),
         spons7=1*(str_detect(tolower(strSponsorType),"other pharma")),
         spons8=1*(str_detect(tolower(strSponsorType),"not for profi")),
         spons9=1*(str_detect(tolower(strSponsorType),"miscel")),
         
         termreason = case_when(
           str_detect(tolower(strTerminationReason), "lack of eff") ~ 4L,
           str_detect(tolower(strTerminationReason), "negative outcome") ~ 4L,
           str_detect(tolower(strTerminationReason), "adverse eff") ~ 4L,
           str_detect(tolower(strTerminationReason), "priorit") ~ 4L,
           str_detect(tolower(strTerminationReason), "business") ~ 3L,
           str_detect(tolower(strTerminationReason), "funding") ~ 3L,
           str_detect(tolower(strTerminationReason), "terminated, other") ~ 2L,
           str_detect(tolower(strTerminationReason), "indeterminat") ~ 2L,
           str_detect(tolower(strTerminationReason), "enrollment") ~ 2L,
           str_detect(tolower(strTerminationReason), "positive outcome") ~ 1L,
           TRUE ~ NA_integer_),
         
         pct_accrual = ifelse( is.na(intactualaccrual) | is.na(inttargetaccrual) | inttargetaccrual==0,  # Proportion of target accrual that is achieved
                               NA_real_,
                               intactualaccrual/inttargetaccrual),
         
         years_since_1999=ifelse(intphaseendyear>2000, # Years since 2000 (below 2000 also = 1)
                                 ifelse(intphaseendyear>2018, 2018-1999, # Truncate to 2018 (should not have data after)
                                        intphaseendyear-1999), 1),
         
         ta1=1*str_detect(strtherapeuticarea, "Autoimmune/Inflammation"),
         ta2=1*str_detect(strtherapeuticarea, "Cardiovascular"),
         ta3=1*str_detect(strtherapeuticarea, "CNS"),
         ta4=1*str_detect(strtherapeuticarea, "Genitourinary"),
         ta5=1*str_detect(strtherapeuticarea, "Infectious Disease"),
         ta6=1*str_detect(strtherapeuticarea, "Metabolic/Endocrinology"),
         ta7=1*str_detect(strtherapeuticarea, "Oncology"),
         ta8=1*str_detect(strtherapeuticarea, "Ophthalmology"),
         ta9=1*str_detect(strtherapeuticarea, "Vaccines"),
         is_a_generic = 1L*( !str_detect(GenericName, "stem cell") & !str_detect(GenericName, "vaccine") & !str_detect(GenericName, "insulin") & !str_detect(GenericName, "interferon") &
           ( str_detect(tolower(TherapyDescription), "biosimilar") | str_detect(tolower(TherapyDescription), "reformulation") ) &
             str_detect(GenericName, ",") ),
         unwilling_to_pay_12k = ifelse(is.na(GenericName), NA_integer_,
                                             1L*( ( str_starts(GenericName,"[A-Z]") & !str_starts(GenericName,"(HPV|MMR|DTP|HIV)")) | 
                                     str_detect(GenericName, "[A-Z][A-Z][A-Z][A-Z][A-Z]") | str_detect(GenericName, "[0-9][0-9]") )),
         is_a_combination = ifelse(is.na(GenericName), NA_integer_,1L*(str_detect(GenericName, "\\+") | (str_detect(GenericName, "[a-z]/[a-z]") & # Identify combination drugs
                                                                   !str_detect(GenericName, "alpha/beta") &
                                                                   !str_detect(GenericName, "beta/gamma") & 
                                                                   !str_detect(GenericName, "gamma/delta")))),
         is_mab = ifelse(is.na(GenericName), NA_integer_,
                               1L*(str_detect(GenericName, "mab$") | str_detect(GenericName, "mab ")  | # Is the drug a monoclonal antibody (ending in -mab or muromonab, which is the one exception to this rule)
                       str_detect(GenericName, "mab,") | str_detect(GenericName, "mab-") | str_detect(GenericName, "muromonab"))),
         
         prior_approval = ifelse(is.na(intpriorapproval), NA_integer_, 1L*(intpriorapproval=='[]') + 0*(intpriorapproval=='[0]'))
         
  ) %>%
  dplyr::select(-strSponsorType, -strTerminationReason, -strtherapeuticarea, -GenericName, -intpriorapproval,
                -X1, -intClinicalTrialID, -intphaseendyear,
                -originkey, -mediumdescription, -strregulatorystatus, -strdesignkeyword, -strDiseaseType,
                -DrugCountryName, -strMechanismOfAction, -DrugDeliveryRouteDescription, -DrugTarget, -strPatientPopulation, 
                -strPrimaryEndPoint, -strGender, -strExclusionCriteria, -decMinAge, -decMaxAge, -strMinAgeUnit, -strMaxAgeUnit,
                -strStudyDesign, -TherapyDescription, -strPatientSegment, -strSponsor, -strLocation)

dis1 <- fulldat %>%
  dplyr::select(DrugKey, indicationkey, row_id, strDiseaseType) %>% #strtherapeuticarea,  outcome
  mutate(strDiseaseType=str_remove_all(str_remove_all(str_remove_all(str_remove_all(strDiseaseType, "\\(N/A\\)\\|"),"\\|\\(N/A\\)"),"\\(\\)\\|"),"\\|\\(\\)")) %>%
  mutate(distyp = str_split(strDiseaseType, pattern="\\|")) %>%
  unnest(cols="distyp") %>%
  distinct() %>%
  mutate(distyp=str_remove_all(distyp, "[\\[\\]]"),
         lifestyle_aspect = case_when(
           str_squish(str_remove_all(tolower(distyp),' ')) == 'femalesexualdysfunction' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'endometriosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'contraception' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'menopausalsymptoms' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'infertility' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otherviralvaccines' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hepatitisvaccines' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'respiratoryvaccines' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'vector-bornediseasevaccines' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otherbacterialvaccines' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'influenzavaccines' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'overactivebladder' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'benignprostatichyperplasia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'urinaryincontinence' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'transplantation/gvhd' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'asthma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'chronicobstructivepulmonarydisease' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pulmonaryfibrosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pain(neuropathic)' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pain(nociceptive)' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == '(n/a)' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'supportivecare' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'na' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'renaldisease' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hepaticfibrosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'nafld' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'age-relatedmaculardegeneration' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'diabeticretinopathy' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'glaucoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'retinitispigmentosa' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'retinalveinocclusion' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'dryeyesyndrome' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'colorectal' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'softtissuesarcoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cns,glioblastoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'multiplemyeloma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'ovarian' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'prostate' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'renal' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'melanoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'esophageal' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'breast' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lung,non-smallcell' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'liver' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pancreas' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'head/neck' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'endometrial' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,acutemyelogenous' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "lymphoma,non-hodgkin's" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,acutelymphocytic' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gastric' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'unspecifiedsolidtumor' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'thyroid' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'mesothelioma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,chroniclymphocytic' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "lymphoma,hodgkin's" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bladder' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,chronicmyelogenous' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lung,smallcell' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cns,other' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cns,medulloblastoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'metastaticcancer' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'myeloproliferativeneoplasms' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gist' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'osteosarcoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'skin,basalcellcarcinoma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'unspecifiedhematologicalcancer' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'unspecifiedcancer' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'testicular' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cervical' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'neuroendocrine' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'insomnia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "alzheimer'sdisease" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'multiplesclerosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'anxiety' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'depression' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bipolardisorder' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'movementdisorders' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'attentiondeficithyperactivedisorder' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'migraine' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'amyotrophiclateralsclerosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "parkinson'sdisease" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'schizophrenia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'autism' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'epilepsy' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'restlesslegssyndrome' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "huntington'sdisease" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "dementia(non-alzheimer's)" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'diabeticcomplications' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'type2diabetes' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'obesity' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'dyslipidemia' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'type1diabetes' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hyponatremia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hyperuricemia/gout' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hbv' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hiv' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hcv' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hpv' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'rabies' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'westnilevirus(wnv)' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cytomegalovirusinfection(cmv)' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'sepsis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'respiratoryinfections' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bacterialskininfection' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otitismedia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'clostridiumdifficile' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'urinarytractinfections' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'onychomycosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'intra-abdominalinfections' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'myelodysplasticsyndrome' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hemostasis/hemophilia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'thalassemia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'anemia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'sicklecelldisease' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'ulcerativecolitis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gerd' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'irritablebowelsyndrome' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'uterinefibroids' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'constipation' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'functionaldyspepsia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gastroparesis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'psoriasis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'atopicdermatitis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'anti-aging(dermatology)' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cerebralpalsy' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'spinalmuscularatrophies' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'growthdisorders' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'neonatalbraininjury' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lysosomalstoragedisorders' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'infantrespiratorydistresssyndrome' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'acutecoronarysyndromes' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'coronaryarterydisease' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'peripheralarterialdisease' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'congestiveheartfailure' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hypertension' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'thromboticdisorders' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'arrhythmia' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cardiomyopathy' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'stroke(neuroprotection)' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'osteoporosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'osteoarthritis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otherinflammatoryarthritis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bonefracturehealing' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'allergicrhinitis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lupus' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "crohn'sdisease" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'rheumatoidarthritis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cysticfibrosis' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == "sjogren'ssyndrome" ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'scleroderma' ~ 0,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'smokingcessation' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'alcoholdependence' ~ 1,
           TRUE ~ NA_real_ ),
        new_disno = case_when(
           str_squish(str_remove_all(tolower(distyp),' ')) == 'femalesexualdysfunction' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'endometriosis' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'contraception' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'menopausalsymptoms' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'infertility' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otherviralvaccines' ~ 2,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hepatitisvaccines' ~ 2,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'respiratoryvaccines' ~ 2,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'vector-bornediseasevaccines' ~ 2,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otherbacterialvaccines' ~ 2,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'influenzavaccines' ~ 2,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'overactivebladder' ~ 3,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'benignprostatichyperplasia' ~ 3,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'urinaryincontinence' ~ 3,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'transplantation/gvhd' ~ 4,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'asthma' ~ 5,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'chronicobstructivepulmonarydisease' ~ 5,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pulmonaryfibrosis' ~ 5,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pain(neuropathic)' ~ 6,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pain(nociceptive)' ~ 6,
           str_squish(str_remove_all(tolower(distyp),' ')) == '(n/a)' ~ 7,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'supportivecare' ~ 7,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'na' ~ 7,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'renaldisease' ~ 8,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hepaticfibrosis' ~ 9,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'nafld' ~ 9,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'age-relatedmaculardegeneration' ~ 10,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'diabeticretinopathy' ~ 10,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'glaucoma' ~ 10,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'retinitispigmentosa' ~ 10,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'retinalveinocclusion' ~ 10,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'dryeyesyndrome' ~ 10,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'colorectal' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'softtissuesarcoma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cns,glioblastoma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'multiplemyeloma' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'ovarian' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'prostate' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'renal' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'melanoma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'esophageal' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'breast' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lung,non-smallcell' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'liver' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'pancreas' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'head/neck' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'endometrial' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,acutemyelogenous' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == "lymphoma,non-hodgkin's" ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,acutelymphocytic' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gastric' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'unspecifiedsolidtumor' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'thyroid' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'mesothelioma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,chroniclymphocytic' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == "lymphoma,hodgkin's" ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bladder' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'leukemia,chronicmyelogenous' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lung,smallcell' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cns,other' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cns,medulloblastoma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'metastaticcancer' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'myeloproliferativeneoplasms' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gist' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'osteosarcoma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'skin,basalcellcarcinoma' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'unspecifiedhematologicalcancer' ~ 12,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'unspecifiedcancer' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'testicular' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cervical' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'neuroendocrine' ~ 11,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'insomnia' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == "alzheimer'sdisease" ~ 14,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'multiplesclerosis' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'anxiety' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'depression' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bipolardisorder' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'movementdisorders' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'attentiondeficithyperactivedisorder' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'migraine' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'amyotrophiclateralsclerosis' ~ 14,
           str_squish(str_remove_all(tolower(distyp),' ')) == "parkinson'sdisease" ~ 14,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'schizophrenia' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'autism' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'epilepsy' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'restlesslegssyndrome' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == "huntington'sdisease" ~ 14,
           str_squish(str_remove_all(tolower(distyp),' ')) == "dementia(non-alzheimer's)" ~ 14,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'diabeticcomplications' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'type2diabetes' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'obesity' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'dyslipidemia' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'type1diabetes' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hyponatremia' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hyperuricemia/gout' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hbv' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hiv' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hcv' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hpv' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'rabies' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'westnilevirus(wnv)' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cytomegalovirusinfection(cmv)' ~ 16,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'sepsis' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'respiratoryinfections' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bacterialskininfection' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otitismedia' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'clostridiumdifficile' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'urinarytractinfections' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'onychomycosis' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'intra-abdominalinfections' ~ 17,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'myelodysplasticsyndrome' ~ 18,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hemostasis/hemophilia' ~ 18,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'thalassemia' ~ 18,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'anemia' ~ 18,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'sicklecelldisease' ~ 18,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'ulcerativecolitis' ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gerd' ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'irritablebowelsyndrome' ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'uterinefibroids' ~ 1,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'constipation' ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'functionaldyspepsia' ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'gastroparesis' ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'psoriasis' ~ 20,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'atopicdermatitis' ~ 20,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'anti-aging(dermatology)' ~ 20,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cerebralpalsy' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'spinalmuscularatrophies' ~ 14,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'growthdisorders' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'neonatalbraininjury' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lysosomalstoragedisorders' ~ 15,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'infantrespiratorydistresssyndrome' ~ 5,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'acutecoronarysyndromes' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'coronaryarterydisease' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'peripheralarterialdisease' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'congestiveheartfailure' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'hypertension' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'thromboticdisorders' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'arrhythmia' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cardiomyopathy' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'stroke(neuroprotection)' ~ 21,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'osteoporosis' ~ 22,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'osteoarthritis' ~ 23,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'otherinflammatoryarthritis' ~ 23,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'bonefracturehealing' ~ 22,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'allergicrhinitis' ~ 5,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'lupus' ~ 23,
           str_squish(str_remove_all(tolower(distyp),' ')) == "crohn'sdisease" ~ 19,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'rheumatoidarthritis' ~ 23,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'cysticfibrosis' ~ 5,
           str_squish(str_remove_all(tolower(distyp),' ')) == "sjogren'ssyndrome" ~ 23,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'scleroderma' ~ 20,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'smokingcessation' ~ 13,
           str_squish(str_remove_all(tolower(distyp),' ')) == 'alcoholdependence' ~ 13,
           TRUE ~ NA_real_ ),
         disno = case_when(
           trimws(tolower(distyp)) == trimws(tolower("Allergic Rhinitis")) ~1,
           trimws(tolower(distyp)) == trimws(tolower("Asthma")) ~2,
           trimws(tolower(distyp)) == trimws(tolower("Chronic Obstructive Pulmonary Disease")) ~3,
           trimws(tolower(distyp)) == trimws(tolower("Insomnia")) ~4,
           trimws(tolower(distyp)) == trimws(tolower("HBV")) ~5,
           trimws(tolower(distyp)) == trimws(tolower("Pain (neuropathic)")) ~6,
           trimws(tolower(distyp)) == trimws(tolower("Colorectal")) ~7,
           trimws(tolower(distyp)) == trimws(tolower("Soft Tissue Sarcoma")) ~8,
           trimws(tolower(distyp)) == trimws(tolower("CNS, Glioblastoma")) ~9,
           trimws(tolower(distyp)) == trimws(tolower("Multiple Myeloma")) ~10,
           trimws(tolower(distyp)) == trimws(tolower("(N/A)")) ~11,
           trimws(tolower(distyp)) == trimws(tolower("Ovarian")) ~12,
           trimws(tolower(distyp)) == trimws(tolower("Prostate")) ~13,
           trimws(tolower(distyp)) == trimws(tolower("Renal")) ~14,
           trimws(tolower(distyp)) == trimws(tolower("Melanoma")) ~15,
           trimws(tolower(distyp)) == trimws(tolower("Esophageal")) ~16,
           trimws(tolower(distyp)) == trimws(tolower("Breast")) ~17,
           trimws(tolower(distyp)) == trimws(tolower("Lung, Non-Small Cell")) ~18,
           trimws(tolower(distyp)) == trimws(tolower("Liver")) ~19,
           trimws(tolower(distyp)) == trimws(tolower("Pancreas")) ~20,
           trimws(tolower(distyp)) == trimws(tolower("Head/Neck")) ~21,
           trimws(tolower(distyp)) == trimws(tolower("Endometrial")) ~22,
           trimws(tolower(distyp)) == trimws(tolower("Leukemia, Acute Myelogenous")) ~23,
           trimws(tolower(distyp)) == trimws(tolower("Lymphoma, Non-Hodgkin's")) ~24,
           trimws(tolower(distyp)) == trimws(tolower("Myelodysplastic Syndrome")) ~25,
           trimws(tolower(distyp)) == trimws(tolower("Leukemia, Acute Lymphocytic")) ~26,
           trimws(tolower(distyp)) == trimws(tolower("Lupus")) ~27,
           trimws(tolower(distyp)) == trimws(tolower("Age-Related Macular Degeneration")) ~28,
           trimws(tolower(distyp)) == trimws(tolower("Gastric")) ~29,
           trimws(tolower(distyp)) == trimws(tolower("Alzheimer's Disease")) ~30,
           trimws(tolower(distyp)) == trimws(tolower("Diabetic Complications")) ~31,
           trimws(tolower(distyp)) == trimws(tolower("Acute Coronary Syndromes")) ~32,
           trimws(tolower(distyp)) == trimws(tolower("Coronary Artery Disease")) ~33,
           trimws(tolower(distyp)) == trimws(tolower("Peripheral Arterial Disease")) ~34,
           trimws(tolower(distyp)) == trimws(tolower("Psoriasis")) ~35,
           trimws(tolower(distyp)) == trimws(tolower("Unspecified Solid Tumor")) ~36,
           trimws(tolower(distyp)) == trimws(tolower("Multiple Sclerosis")) ~37,
           trimws(tolower(distyp)) == trimws(tolower("HIV")) ~38,
           trimws(tolower(distyp)) == trimws(tolower("Other Viral Vaccines")) ~39,
           trimws(tolower(distyp)) == trimws(tolower("Congestive Heart Failure")) ~40,
           trimws(tolower(distyp)) == trimws(tolower("Hypertension")) ~41,
           trimws(tolower(distyp)) == trimws(tolower("Crohn's Disease")) ~42,
           trimws(tolower(distyp)) == trimws(tolower("Osteoporosis")) ~43,
           trimws(tolower(distyp)) == trimws(tolower("Thyroid")) ~44,
           trimws(tolower(distyp)) == trimws(tolower("Mesothelioma")) ~45,
           trimws(tolower(distyp)) == trimws(tolower("Rheumatoid Arthritis")) ~46,
           trimws(tolower(distyp)) == trimws(tolower("Ulcerative Colitis")) ~47,
           trimws(tolower(distyp)) == trimws(tolower("Leukemia, Chronic Lymphocytic")) ~48,
           trimws(tolower(distyp)) == trimws(tolower("Thrombotic Disorders")) ~49,
           trimws(tolower(distyp)) == trimws(tolower("GERD")) ~50,
           trimws(tolower(distyp)) == trimws(tolower("Irritable Bowel Syndrome")) ~51,
           trimws(tolower(distyp)) == trimws(tolower("Supportive Care")) ~52,
           trimws(tolower(distyp)) == trimws(tolower("Type 2 Diabetes")) ~53,
           trimws(tolower(distyp)) == trimws(tolower("Obesity")) ~54,
           trimws(tolower(distyp)) == trimws(tolower("HCV")) ~55,
           trimws(tolower(distyp)) == trimws(tolower("Anxiety")) ~56,
           trimws(tolower(distyp)) == trimws(tolower("Depression")) ~57,
           trimws(tolower(distyp)) == trimws(tolower("Bipolar Disorder")) ~58,
           trimws(tolower(distyp)) == trimws(tolower("Atopic Dermatitis")) ~59,
           trimws(tolower(distyp)) == trimws(tolower("Cystic Fibrosis")) ~60,
           trimws(tolower(distyp)) == trimws(tolower("Lymphoma, Hodgkin's")) ~61,
           trimws(tolower(distyp)) == trimws(tolower("Bladder")) ~62,
           trimws(tolower(distyp)) == trimws(tolower("Leukemia, Chronic Myelogenous")) ~63,
           trimws(tolower(distyp)) == trimws(tolower("Lung, Small Cell")) ~64,
           trimws(tolower(distyp)) == trimws(tolower("Sepsis")) ~65,
           trimws(tolower(distyp)) == trimws(tolower("CNS, Other")) ~66,
           trimws(tolower(distyp)) == trimws(tolower("CNS, Medulloblastoma")) ~67,
           trimws(tolower(distyp)) == trimws(tolower("Movement Disorders")) ~68,
           trimws(tolower(distyp)) == trimws(tolower("Metastatic Cancer")) ~69,
           trimws(tolower(distyp)) == trimws(tolower("Uterine fibroids")) ~70,
           trimws(tolower(distyp)) == trimws(tolower("Myeloproliferative Neoplasms")) ~71,
           trimws(tolower(distyp)) == trimws(tolower("Attention Deficit Hyperactive Disorder")) ~72,
           trimws(tolower(distyp)) == trimws(tolower("Constipation")) ~73,
           trimws(tolower(distyp)) == trimws(tolower("Dyslipidemia")) ~74,
           trimws(tolower(distyp)) == trimws(tolower("Migraine")) ~75,
           trimws(tolower(distyp)) == trimws(tolower("Renal Disease")) ~76,
           trimws(tolower(distyp)) == trimws(tolower("Sjogren's Syndrome")) ~77,
           trimws(tolower(distyp)) == trimws(tolower("Overactive Bladder")) ~78,
           trimws(tolower(distyp)) == trimws(tolower("Arrhythmia")) ~79,
           trimws(tolower(distyp)) == trimws(tolower("Anti-aging (dermatology)")) ~80,
           trimws(tolower(distyp)) == trimws(tolower("Pain (nociceptive)")) ~81,
           trimws(tolower(distyp)) == trimws(tolower("Smoking Cessation")) ~82,
           trimws(tolower(distyp)) == trimws(tolower("Diabetic Retinopathy")) ~83,
           trimws(tolower(distyp)) == trimws(tolower("Female Sexual Dysfunction")) ~84,
           trimws(tolower(distyp)) == trimws(tolower("Amyotrophic Lateral Sclerosis")) ~85,
           trimws(tolower(distyp)) == trimws(tolower("Parkinson's Disease")) ~86,
           trimws(tolower(distyp)) == trimws(tolower("GIST")) ~87,
           trimws(tolower(distyp)) == trimws(tolower("Osteoarthritis")) ~88,
           trimws(tolower(distyp)) == trimws(tolower("Endometriosis")) ~89,
           trimws(tolower(distyp)) == trimws(tolower("Glaucoma")) ~90,
           trimws(tolower(distyp)) == trimws(tolower("Other Inflammatory Arthritis")) ~91,
           trimws(tolower(distyp)) == trimws(tolower("Transplantation/GVHD")) ~92,
           trimws(tolower(distyp)) == trimws(tolower("Hemostasis/Hemophilia")) ~93,
           trimws(tolower(distyp)) == trimws(tolower("Schizophrenia")) ~94,
           trimws(tolower(distyp)) == trimws(tolower("Osteosarcoma")) ~95,
           trimws(tolower(distyp)) == trimws(tolower("Hepatic Fibrosis")) ~96,
           trimws(tolower(distyp)) == trimws(tolower("NAFLD")) ~97,
           trimws(tolower(distyp)) == trimws(tolower("Cerebral Palsy")) ~98,
           trimws(tolower(distyp)) == trimws(tolower("Respiratory Infections")) ~99,
           trimws(tolower(distyp)) == trimws(tolower("Thalassemia")) ~100,
           trimws(tolower(distyp)) == trimws(tolower("Bacterial Skin Infection")) ~101,
           trimws(tolower(distyp)) == trimws(tolower("Benign Prostatic Hyperplasia")) ~102,
           trimws(tolower(distyp)) == trimws(tolower("Hepatitis Vaccines")) ~103,
           trimws(tolower(distyp)) == trimws(tolower("HPV")) ~104,
           trimws(tolower(distyp)) == trimws(tolower("Type 1 Diabetes")) ~105,
           trimws(tolower(distyp)) == trimws(tolower("Skin, Basal Cell Carcinoma")) ~106,
           trimws(tolower(distyp)) == trimws(tolower("Autism")) ~107,
           trimws(tolower(distyp)) == trimws(tolower("Retinitis Pigmentosa")) ~108,
           trimws(tolower(distyp)) == trimws(tolower("Unspecified Hematological Cancer")) ~109,
           trimws(tolower(distyp)) == trimws(tolower("Otitis Media")) ~110,
           trimws(tolower(distyp)) == trimws(tolower("Epilepsy")) ~111,
           trimws(tolower(distyp)) == trimws(tolower("Spinal Muscular Atrophies")) ~112,
           trimws(tolower(distyp)) == trimws(tolower("Anemia")) ~113,
           trimws(tolower(distyp)) == trimws(tolower("Retinal Vein Occlusion")) ~114,
           trimws(tolower(distyp)) == trimws(tolower("Hyponatremia")) ~115,
           trimws(tolower(distyp)) == trimws(tolower("Urinary Incontinence")) ~116,
           trimws(tolower(distyp)) == trimws(tolower("Alcohol Dependence")) ~117,
           trimws(tolower(distyp)) == trimws(tolower("Unspecified Cancer")) ~118,
           trimws(tolower(distyp)) == trimws(tolower("Restless Legs Syndrome")) ~119,
           trimws(tolower(distyp)) == trimws(tolower("Rabies")) ~120,
           trimws(tolower(distyp)) == trimws(tolower("Contraception")) ~121,
           trimws(tolower(distyp)) == trimws(tolower("Respiratory Vaccines")) ~122,
           trimws(tolower(distyp)) == trimws(tolower("Functional Dyspepsia")) ~123,
           trimws(tolower(distyp)) == trimws(tolower("Testicular")) ~124,
           trimws(tolower(distyp)) == trimws(tolower("Clostridium difficile")) ~125,
           trimws(tolower(distyp)) == trimws(tolower("Urinary Tract Infections")) ~126,
           trimws(tolower(distyp)) == trimws(tolower("Onychomycosis")) ~127,
           trimws(tolower(distyp)) == trimws(tolower("Cervical")) ~128,
           trimws(tolower(distyp)) == trimws(tolower("Vector-Borne Disease Vaccines")) ~129,
           trimws(tolower(distyp)) == trimws(tolower("West Nile Virus (WNV)")) ~130,
           trimws(tolower(distyp)) == trimws(tolower("Pulmonary Fibrosis")) ~131,
           trimws(tolower(distyp)) == trimws(tolower("Other Bacterial Vaccines")) ~132,
           trimws(tolower(distyp)) == trimws(tolower("Influenza Vaccines")) ~133,
           trimws(tolower(distyp)) == trimws(tolower("Growth Disorders")) ~134,
           trimws(tolower(distyp)) == trimws(tolower("Hyperuricemia/Gout")) ~135,
           trimws(tolower(distyp)) == trimws(tolower("Gastroparesis")) ~136,
           trimws(tolower(distyp)) == trimws(tolower("Scleroderma")) ~137,
           trimws(tolower(distyp)) == trimws(tolower("Dry Eye Syndrome")) ~138,
           trimws(tolower(distyp)) == trimws(tolower("Intra-abdominal Infections")) ~139,
           trimws(tolower(distyp)) == trimws(tolower("Menopausal Symptoms")) ~140,
           trimws(tolower(distyp)) == trimws(tolower("Huntington's Disease")) ~141,
           trimws(tolower(distyp)) == trimws(tolower("Cardiomyopathy")) ~142,
           trimws(tolower(distyp)) == trimws(tolower("Sickle Cell Disease")) ~143,
           trimws(tolower(distyp)) == trimws(tolower("Neuroendocrine")) ~144,
           trimws(tolower(distyp)) == trimws(tolower("Infertility")) ~145,
           trimws(tolower(distyp)) == trimws(tolower("Cytomegalovirus Infection (CMV)")) ~146,
           trimws(tolower(distyp)) == trimws(tolower("Dementia (non-Alzheimer's)")) ~147,
           trimws(tolower(distyp)) == trimws(tolower("Neonatal Brain Injury")) ~148,
           trimws(tolower(distyp)) == trimws(tolower("Stroke (neuroprotection)")) ~149,
           trimws(tolower(distyp)) == trimws(tolower("Bone Fracture Healing")) ~150,
           trimws(tolower(distyp)) == trimws(tolower("Lysosomal Storage Disorders")) ~151,
           trimws(tolower(distyp)) == trimws(tolower("Infant Respiratory Distress Syndrome")) ~152,
           TRUE ~ 11 )) %>%
  dplyr::select(row_id, disno, new_disno)

dis2 <- dis1 %>%
  dplyr::select(row_id, disno) %>%
  mutate(anytyp=1) %>%
  distinct() %>%
  group_by(row_id) %>%
  pivot_wider(names_from = disno, names_prefix = "disease", values_from = anytyp, id_cols = "row_id") %>%  
  replace(., is.na(.), 0)

dis2a <- dis1 %>%
  dplyr::select(row_id, new_disno) %>%
  arrange(new_disno) %>%
  mutate(anytyp=1) %>%
  distinct() %>%
  group_by(row_id) %>%
  pivot_wider(names_from = new_disno, names_prefix = "newta", values_from = anytyp, id_cols = "row_id") %>%  
  replace(., is.na(.), 0) %>%
  ungroup() %>%
  mutate(newtanos = rowSums(dplyr::select(.,starts_with("newta"))),
         newta=case_when(
           newtanos==1 ~ as.integer(max.col( dplyr::select(.,starts_with("newta")), 'first')),
           newtanos==2 & newta12==1 & newta18==1 ~ 12L, # hematology and liquid tumor = liquid tumor
           newta6==1 ~ 6L, # pain is pain
           newta2==1 ~ 2L, # vaccines primarily vaccines
           newta13==1 & newta14==1 ~ 14L, # Neuroscience + neurodegenration = neurodegeneration
           newta10==1 ~ 10L, # If it's ophtha, focus on that
           newta21==1 ~ 21L, # CV = primarily CV; even if combined with metabolism and/or renal and/or hematology or derm or tumors or resp
           newta11==1 & newta12==1 ~ 24L, # Extra category for liquid + solid tumor
           newta11==1 ~ 11L,
           newta12==1 ~ 12L,
           newta16==1 ~ 16L, # Infections are primarily infections
           newta17==1 ~ 17L, # Infections are primarily infections
           newta4==1 ~ 4L,
           newta8==1 ~ 8L,
           newta13==1 ~ 13L,
           newta5==1 ~ 5L,
           newta20==1 ~ 20L,
           newta19==1 ~ 19L,
           newta23==1 ~ 23L,
           newta18==1 ~ 18L,
           newta15==1 ~ 15L,
           newta9==1 ~ 9L,
           TRUE  ~ NA_integer_),
         newta = ifelse(newta==7, 23, ifelse(newta>7, newta-1, newta))) #%>%
  # filter(is.na(newta)) %>%
  # group_by(newta5, newta13, newta16, newta6, newta11, newta12, newta18, newta23, newta10, newta14, newta15, newta21, newta20, newta2, newta19, newta22, newta7, newta17, newta1, newta8, newta3, newta4, newta9) %>%
  # summarize(n=n()) %>%
  # arrange(desc(n)) %>%
  # print(n=60)


drugclass1 <- fulldat %>%
  dplyr::select(DrugKey, strMechanismOfAction, outcome, intpriorapproval, indicationkey, intphaseendyear, strDiseaseType, strtherapeuticarea) %>%
  mutate(moa = str_split(strMechanismOfAction, pattern="\\|"),
         prior_approval = 1*(intpriorapproval=='[]') + 0*(intpriorapproval=='[0]')) %>%
  unnest(cols="moa") %>%
  mutate(moa=str_remove_all(moa, "[\\[\\]]")) %>%
  ungroup() %>%
  dplyr::select(-intpriorapproval)

logit <- function(x){
  log(x)-log(1-x)
}

drugclass2 <- drugclass1 %>%
  left_join(drugclass1, by="moa") %>%
  filter(intphaseendyear.x>=intphaseendyear.y) %>%
  mutate(outcome.y=ifelse(DrugKey.x==DrugKey.y & indicationkey.x==indicationkey.y, NA, outcome.y),
         prior_approval.y=ifelse(DrugKey.x==DrugKey.y & indicationkey.x==indicationkey.y, NA, prior_approval.y)) %>%
  dplyr::select(DrugKey.x, indicationkey.x, moa, prior_approval.x, DrugKey.y, outcome.y, prior_approval.y) %>%
  distinct() %>%
  group_by(DrugKey.x, indicationkey.x, moa) %>%
  summarize(class_counts = ifelse(max(prior_approval.x)==1 & sum(0,outcome.y, na.rm=T)==0,n(),n()-1),
            class_approvals = max(prior_approval.x, sum(0,outcome.y, na.rm=T)),
            class_logit_lcl50 = logit(qbeta(0.25, 1/3+class_approvals, 
                                 2.7+(class_counts-class_approvals) )),
            class_logit_ucl50 = logit(qbeta(0.75, 1/3+class_approvals, 
                                     2.7+(class_counts-class_approvals) ))) %>%
  group_by(DrugKey.x, indicationkey.x) %>%
  summarize(moas=n(),
            meancll50=mean(class_logit_lcl50, na.rm=T),
            meanclu50=mean(class_logit_ucl50, na.rm=T),
            class_approvals = sum(class_approvals),
            class_counts=sum(class_counts)) %>%
  ungroup() %>%
  rename(DrugKey=DrugKey.x, indicationkey=indicationkey.x)

drugclass3 <- drugclass1 %>%
  left_join(drugclass1, by="moa") %>%
  filter(intphaseendyear.x>=intphaseendyear.y & strDiseaseType.x==strDiseaseType.y) %>%
  mutate(outcome.y=ifelse(DrugKey.x==DrugKey.y & indicationkey.x==indicationkey.y, NA, outcome.y),
         prior_approval.y=ifelse(DrugKey.x==DrugKey.y & indicationkey.x==indicationkey.y, NA, prior_approval.y)) %>%
  dplyr::select(DrugKey.x, indicationkey.x, moa, prior_approval.x, DrugKey.y, outcome.y, prior_approval.y) %>%
  distinct() %>%
  group_by(DrugKey.x, indicationkey.x, moa) %>%
  summarize(class_counts = ifelse(max(prior_approval.x)==1 & sum(0,outcome.y, na.rm=T)==0,n(),n()-1),
            class_approvals = max(prior_approval.x, sum(0,outcome.y, na.rm=T)),
            class_logit_lcl50 = logit(qbeta(0.25, 1/3+class_approvals, 
                                            2.7+(class_counts-class_approvals) )),
            class_logit_ucl50 = logit(qbeta(0.75, 1/3+class_approvals, 
                                            2.7+(class_counts-class_approvals) ))) %>%
  group_by(DrugKey.x, indicationkey.x) %>%
  summarize(dtmoas=n(),
            dtmeancll50=mean(class_logit_lcl50, na.rm=T),
            dtmeanclu50=mean(class_logit_ucl50, na.rm=T),
            dtclass_approvals = sum(class_approvals),
            dtclass_counts=sum(class_counts)) %>%
  ungroup() %>%
  rename(DrugKey=DrugKey.x, indicationkey=indicationkey.x)

drugclass4 <- drugclass1 %>%
  left_join(drugclass1, by="moa") %>%
  filter(intphaseendyear.x>=intphaseendyear.y & strtherapeuticarea.x==strtherapeuticarea.y) %>%
  mutate(outcome.y=ifelse(DrugKey.x==DrugKey.y & indicationkey.x==indicationkey.y, NA, outcome.y),
         prior_approval.y=ifelse(DrugKey.x==DrugKey.y & indicationkey.x==indicationkey.y, NA, prior_approval.y)) %>%
  dplyr::select(DrugKey.x, indicationkey.x, moa, prior_approval.x, DrugKey.y, outcome.y, prior_approval.y) %>%
  distinct() %>%
  group_by(DrugKey.x, indicationkey.x, moa) %>%
  summarize(class_counts = ifelse(max(prior_approval.x)==1 & sum(0,outcome.y, na.rm=T)==0,n(),n()-1),
            class_approvals = max(prior_approval.x, sum(0,outcome.y, na.rm=T)),
            class_logit_lcl50 = logit(qbeta(0.25, 1/3+class_approvals, 
                                            2.7+(class_counts-class_approvals) )),
            class_logit_ucl50 = logit(qbeta(0.75, 1/3+class_approvals, 
                                            2.7+(class_counts-class_approvals) ))) %>%
  group_by(DrugKey.x, indicationkey.x) %>%
  summarize(tamoas=n(),
            tameancll50=mean(class_logit_lcl50, na.rm=T),
            tameanclu50=mean(class_logit_ucl50, na.rm=T),
            taclass_approvals = sum(class_approvals),
            taclass_counts=sum(class_counts)) %>%
  ungroup() %>%
  rename(DrugKey=DrugKey.x, indicationkey=indicationkey.x)

tdesign <- fulldat %>%
  dplyr::select(row_id, strdesignkeyword) %>% #row_id,
  mutate(design = str_split(strdesignkeyword, pattern = "\\|")) %>%
  unnest(cols = "design") %>%
  dplyr::select(-strdesignkeyword) %>%
  mutate(design = str_remove_all(design, "[\\[\\]]")) %>%
  distinct() %>%
  mutate(
    designscore = case_when(
      design == "efficacy" ~ 1,
      design == "safety" ~ NA_real_,
      design == "open label" ~ -1,
      design == "randomizd" | design == "randomized" ~ 1,
      design == "multiple arm" ~ 1,
      design == "double blind/blinded" ~ 1,
      design == "single arm" ~ -1,
      design == "placebo control" ~ 1,
      design == "pharmacokinetics" ~ -1,
      design == "pharmacodynamics" ~ -1,
      design == "active comparator" ~ 1,
      design == "dose response" ~ 1,
      is.na(design) ~ NA_real_,
      design == "cross over" ~ 1,
      design == "immunogenicity" ~ -1,
      design == "fixed dose" ~ NA_real_,
      design == "adaptive" ~ NA_real_,
      design == "multiple ascending dose" ~ -1,
      design == "single ascending dose" ~ -1,
      design == "observational" ~ -1,
      design == "non inferiority" ~ 1,
      design == "superiority" ~ 1,
      design == "bioavailability" ~ -1,
      design == "basket" ~ 1,
      design == "drug-drug interaction" ~ -1,
      design == "umbrella" ~ 1,
      design == "bioequivalence" ~ -1,
      design == "non interventional" ~ -1,
      TRUE ~ NA_real_
    ),
    efficacy_assessed = 1*(design=="efficacy" | design=="dose response" |
                             design=="superiority" | design == "non inferiority"),
    randomized = 1*(design == "randomizd"  | design == "randomized"),
    blinded = 1*(design == "double blind/blinded"),
    controlled = 1*(design == "multiple arm" | design == "double blind/blinded" | 
                      design == "placebo control" | design == "cross over" | 
                      design == "active comparator" | design == "dose response" |
                      design == "non inferiority" | design == "superiority"),
    pk_aspects = 1*(design == "pharmacokinetics" | design == "immunogenicity" | 
                      design == "bioavailability" | design == "drug-drug interaction" |
                      design == "bioequivalence"),
    innovative_design = 1*(design == "basket" | design == "umbrella" | design == "adaptive")
  ) %>%
  group_by(row_id) %>%
  summarize(designscore=sum(designscore, na.rm=T),
            efficacy_assessed=ifelse( max(efficacy_assessed, na.rm=T)==-Inf, NA_real_, max(efficacy_assessed, na.rm=T)),
            randomized=ifelse( max(randomized, na.rm=T)==-Inf, NA_real_, max(randomized, na.rm=T)),
            blinded=ifelse( max(blinded, na.rm=T)==-Inf, NA_real_, max(blinded, na.rm=T)),
            controlled=ifelse( max(controlled, na.rm=T)==-Inf, NA_real_, max(controlled, na.rm=T)),
            pk_aspects=ifelse( max(pk_aspects, na.rm=T)==-Inf, NA_real_, max(pk_aspects, na.rm=T)),
            innovative_design=ifelse( max(innovative_design, na.rm=T)==-Inf, NA_real_, max(innovative_design, na.rm=T))) 

final <- basics %>%
  left_join(dis2, by="row_id") %>%
  left_join(dis2a, by="row_id") %>%
  left_join(tdesign, by="row_id") %>%
  left_join(trialcompl, by="row_id") %>%
  left_join(insflu, by="row_id") %>%  
  left_join(logitoffsets, by="row_id") %>%
  left_join(relsize, by="row_id") %>%
  left_join(time_since_first, by="row_id") %>%
  left_join(dplyr::select(fulldat, row_id, GenericName, strDiseaseType), by="row_id") %>%
  rename(supportive_care=disease52) %>%
  left_join(drugclass2, by=c("DrugKey", "indicationkey")) %>%
  left_join(drugclass3, by=c("DrugKey", "indicationkey")) %>%
  left_join(drugclass4, by=c("DrugKey", "indicationkey")) %>%
  left_join(drugtypes, by=c("DrugKey", "indicationkey")) %>%
  left_join(orphanind, by=c("DrugKey", "indicationkey")) %>%
  ungroup() %>%
  mutate( detrended_trialendscore = mean( mean_trialendscore-1.14671-years_since_1999*0.02393),
          after2007 = (years_since_1999>8)*1L,
          monoclonal = ifelse(is.na(monoclonal), 0, monoclonal),
          monoclonal_no_inn = 1L*(monoclonal==1 & is_mab==0) )   # based on lm(data=adat, mean_trialendscore ~ years_since_1999)
  
casewgts <- 1/pmax(1, glm(data=final, outcome ~ years_since_1999, family=binomial("logit")) %>% predict(object=., newdata=tibble(years_since_1999=1:19), type="response") / 
                     mean(final$outcome[final$years_since_1999>=12], na.rm=T))

adat <- final %>%
  mutate(casewgt1 = ifelse(is.na(outcome), 1, ifelse(outcome==1, casewgts[years_since_1999], 1)),
         casewgt2 = ifelse(is.na(outcome), 1, casewgts[years_since_1999]*casewgt1),
         drugtype1=ifelse(is.na(drugtype1), 0, drugtype1),
         drugtype2=ifelse(is.na(drugtype2), 0, drugtype2),
         drugtype3=ifelse(is.na(drugtype3), 0, drugtype3),
         drugtype4=ifelse(is.na(drugtype4), 0, drugtype4),
         drugtype5=ifelse(is.na(drugtype5), 0, drugtype5),
         drugtype6=ifelse(is.na(drugtype6), 0, drugtype6),
         drugtype7=ifelse(is.na(drugtype7), 0, drugtype7),
         drugtype8=ifelse(is.na(drugtype8), 0, drugtype8),
         drugtype9=ifelse(is.na(drugtype9), 0, drugtype9),
         lcm_onc = ta7*prior_approval) %>%
  group_by(DrugKey, indicationkey) %>%
  mutate(predgroup = group_indices()) %>% 
  ungroup()
return(adat)
}

adat <- create_trial_level_analysis_dataset()

# final %>% group_by(years_since_1999) %>% summarize(mo=weighted.mean(outcome, casewgt1, na.rm=T)) %>% ggplot(aes(x=years_since_1999, y=mo)) +geom_line()
# 
# mean(final$casewgt1)
# mean(final$casewgt2)

if (AICROWD_TRAIN_DATA_PATH == "/shared_data/data/training_data/training_data_2015_split_on_outcome.csv") {
  write_csv(adat, "/files/feat_bjoern_trial.csv")  
}

    # mutate(newta2310=case_when(
  #   newta2==1 ~ 1,
  #   newta3==1 ~ 2,
  #   newta10==1 ~ 3,
  #   TRUE ~ 0)) %>%
  # dplyr::select(-newta2, -newta3, -newta10)

# final %>% 
#   dplyr::select(starts_with("newta")) %>%
#   mutate(xxx=rowSums(data.frame(.))) %>%
#   filter(xxx>1 ) %>%
#   group_by(newta5, newta13, newta16, newta6, newta11, newta12, newta18, newta23, newta10, newta14, newta15, newta21, newta20, newta2, newta19, newta22, newta7, newta17, newta1, newta8, newta3, newta4, newta9) %>%
#   summarize(n=n()) %>%
#   arrange(desc(n)) %>%
#   print(n=60)


  # dplyr::select(xxx) %>% 
  # summarize(min=min(xxx), mean=mean(xxx), median=median(xxx), p80=quantile(xxx, p=0.8), max=max(xxx))
  # #mutate(ntsm=sum(vars(starts_with("newta"))))

# final %>%
#   group_by(newta) %>%
#   summarize(n=n(),
#             t1=sum(ta1),
#             t2=sum(ta2),
#             t3=sum(ta3),
#             t4=sum(ta4),
#             t5=sum(ta5),
#             t6=sum(ta6),
#             t7=sum(ta7),
#             t8=sum(ta8),
#             t9=sum(ta9)) %>% print(n=25)

