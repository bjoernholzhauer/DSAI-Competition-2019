library(tidyverse)
library(stringdist)
library(stringi)
library(furrr)

AICROWD_TRAIN_DATA_PATH <- Sys.getenv(c("AICROWD_TRAIN_DATA_PATH"), unset="/shared_data/data/training_data/training_data_2015_split_on_outcome.csv")
AICROWD_TEST_DATA_PATH <- Sys.getenv(c("AICROWD_TEST_DATA_PATH"), unset="/shared_data/data/test_data_full/testing_phase2_release.csv")

fulldat <- read_csv(AICROWD_TRAIN_DATA_PATH) %>%
  bind_rows(read_csv(AICROWD_TEST_DATA_PATH))

#countries = "australia|the usa|the us|europe|and the usa|and europe|and the us|eu and switzerland|the eu and the us|the eu and us|unspecified countries|australia and the us|and the us|and the eu|the eu|the us and eu||the usa|japan|australia|australia|australia|australian|china|europe|japan|th us|tha us|ths us|the eu|the us|japan|mexico|s korea|s korea|south korea|swizerland|s korea|australia|unspecified|usa|japan|switzerland|mexico|usa|taiwan|korea"

countries = c("europe and the us,", "the eu and the us;", "the eu and the us,", "the us and eu;", "australia and the us,", "europe and the us;", "the eu and switzerland;", "the eu and us;",
              ", the eu", "; the eu", "; the us", ", the us", "the us;", "the eu;", "the us,", "s korea,",   "australia;", "japan;", "the eu,",  "the usa;", "s korea;",  "japan,", "south korea;", 
              "switzerland;",  "australia,", "europe;", "eu,", "australian;",  "switzerland,", "swizerland;", "taiwan,",  " the us$",
              "th us;", "tha us;", "the eu ",  "ths us;", "australian", "australia", "mexico;",  "eu;","usa,",  "usa;","us;", "china;", "europe,", "korea;", "korea", "japan")


replace_list = list(
  list(a = "Cancer", b = "Cancer (Oncology)"),
  list(a = "CML", b = "Chronic myeloid leukaemia"),
  list(a = "MAS", b = "Macrophage activation syndrome"),
  list(a = "APL", b = "Acute promyelocytic leukemia"),
  list(a = "CTL", b = "Cutaneous T cell lymphoma"),
  list(a = "GHIS", b = "Growth hormone insensitivity syndrome"),
  list(a = "AML", b = "Acute myeloid leukemia"),
  list(a = "CLL", b = "Chronic lymphocytic leukemia"),
  list(a = "AA", b = "Aortic Aneurysm"),
  list(a = "SMA", b = "Spinal muscular atrophy"),
  list(a = "AMD", b = "Age-related macular degeneration (AMD)"),
  list(a = "ALS", b = "Amyotrophic lateral sclerosis"),
  list(a = "LBL", b = "Lymphoblastic lymphoma"),
  list(a = "ALL", b = "Acute lymphocytic leukemia"),
  list(a = "CTCL", b ="Cutaneous T cell lymphoma"),
  list(a = "nsclc", b = "Non-small-cell lung carcinoma"),
  list(a = "MELAS", b ="Mitochondrial encephalomyopathy, lactic acidosis, and stroke-like episodes"),
  list(a = "NOH", b = "Neurogenic orthostatic hypotension"),
  list(a ="nOH", b = "Neurogenic orthostatic hypotension"),
  list(a = "PAF", b = "Pure Autonomic Failure"),
  list(a ="ALCL", b = "Anaplastic large cell lymphoma"),
  list(a = "NHL", b = "Non-Hodgkin Lymphoma"),
  list(a ="under 34wk old", b = "pediatric neonatal"),
  list(a = "CF", b = "cystic fibrosis"),
  list(a ="GvHD", b = "graft versus host disease"),
  list(a ="GVHD", b = "graft versus host disease"),
  list(a = "IPF", b = "Idiopathic pulmonary fibrosis"),
  list(a ="GIST", b = "Gastrointestinal stromal tumors (GIST)"),
  list(a = "MALT", b = "mucosa-associated lymphoid tissue (MALT) lymphoma"),
  list(a ="ZES", b = "Zollinger Ellison syndrome"),
  list(a = "Weight loss in AIDS", b = "Weight loss in AIDS/HIV"),
  list(a ="AIDS-related Kaposi's sarcoma", b = "AIDS/HIV-related Kaposi's sarcoma"),
  list(a = "AIDS", b = "AIDS/HIV"),
  list(a = "Neuralgia", b = "neuralgia (Pain (neuropathic))"),
  list(a = "neuralgia", b = "neuralgia (Pain (neuropathic))"),
  list(a = "Pain,", b = "pain (Pain (neuropathic))"),
  list(a = "pain,", b = "pain (Pain (neuropathic))"),
  list(a = "Palsy", b = "palsy (Movement Disorders)"),
  list(a = "palsy", b = "palsy (Movement Disorders)"),
  list(a = "Spinal cord injury", b = "spinal cord injury (Movement Disorders)"),
  list(a = "spinal cord injury", b = "spinal cord injury (Movement Disorders)"),
  list(a = "Spasticity", b = "spasticity (Movement Disorders)"),
  list(a = "spasticity", b = "spasticity (Movement Disorders)"),
  list(a = "Dyskinesia", b = "dyskinesia (Movement Disorders)"),
  list(a = "dyskinesia", b = "dyskinesia (Movement Disorders)"),
  list(a = "tardive", b = "Tardive (Movement Disorders)"),
  list(a = "tardive", b = "tardive (Movement Disorders)"),
  list(a = "Dystonia,", b = "dystonia (Movement Disorders)"),
  list(a = "dystonia", b = "dystonia (Movement Disorders)"),
  list(a = "Dyskinesia", b = "Dyskinesia (Movement Disorders)"),
  list(a = "dyskinesia", b = "dyskinesia (Movement Disorders)"),
  list(a = "Improvement of walking", b = "improvement of walking (Movement Disorders)"),
  list(a = "improvement of walking", b = "improvement of walking (Movement Disorders)"),
  list(a = "Sleep,", b = "sleep (Insomnia)"),
  list(a = " sleep,", b = " sleep (Insomnia)"),
  list(a = "Radiation-induced", b = "radiation-induced (Supportive Care)"),
  list(a = "radiation-induced", b = "radiation-induced (Supportive Care)"),
  list(a = "reperfusion injury", b = "reperfusion injury (Acute Coronary Syndromes)"),
  list(a = "Reperfusion injury", b = "reperfusion injury (Acute Coronary Syndromes)"),
  list(a = "infarction, myocardial", b = "infarction, myocardial (Acute Coronary Syndromes)"),
  list(a = "Infarction, myocardial", b = "infarction, myocardial (Acute Coronary Syndromes)"),
  list(a = "Pulmonary tuberculosis", b = "Pulmonary tuberculosis infection"),
  list(a = "Buerger's", b = "buerger's (thromboangiitis obliterans/Peripheral Arterial Disease)"),
  list(a = "thromboangiitis obliterans", b = "buerger's (thromboangiitis obliterans/Peripheral Arterial Disease)"),
  list(a = "HIV", b = "HIV (human immunodeficiency virus)"),
  list(a = "HCV", b = "HCV (Hepatitis C virus)"),
  list(a = "HBV", b = "HBV (Hepatitis B virus)"),
  list(a ="HSCT", b = "Hematopoietic stem cell transplantation"),
  list(a = "RDS", b = "Respiratory Distress Syndrome"),
  list(a ="ANCA", b = "Antineutrophil cytoplasmic antibodies"),
  list(a ="[CNS", b="[CNS (non-Hodgkin lymphoma)"),
  list(a ="hepatocellular", b="hepatocellular (liver)"),
  list(a ="Hepatocellular", b="hepatocellular (liver)"),
  list(a ="hypercholesterolaemia", b = "hypercholesterolaemia (Dyslipidemia)"),
  list(a ="hypercholesterolemia", b = "hypercholesterolaemia (Dyslipidemia)"),
  list(a ="hypertriglyceridaemia", b = "hypertriglyceridaemia (Dyslipidemia)"),
  list(a ="Hypercholesterolaemia", b = "hypercholesterolaemia (Dyslipidemia)"),
  list(a ="Hypercholesterolemia", b = "hypercholesterolaemia (Dyslipidemia)"),
  list(a ="Hypertriglyceridaemia", b = "hypertriglyceridaemia (Dyslipidemia)"),
  list(a ="Hyperlipoproteinemia", b = "hyperlipoproteinemia (Dyslipidemia)"),
  list(a ="hyperlipoproteinemia", b = "hyperlipoproteinemia (Dyslipidemia)"),
  list(a ="stomach", b = "stomach (Gastric)"),
  list(a ="Stomach", b = "stomach (Gastric)"),
  list(a ="neuroblastoma", b = "neuroblastoma (CNS)"),
  list(a ="Neuroblastoma", b = "neuroblastoma (CNS)"),
  list(a ="short stature", b = "short stature (growth disorders)"),
  list(a ="Short stature", b = "short stature (growth disorders)"),
  list(a ="delay of growth", b = "delay of growth (growth disorders)"),
  list(a ="Delay of growth", b = "delay of growth (growth disorders)"),
  list(a ="acromegaly", b = "acromegaly (growth disorders)"),
  list(a ="Acromegaly", b = "acromegaly (growth disorders)"),
  list(a ="Carnitine deficiency", b = "carnitine deficiency (growth disorders/Spinal Muscular Atrophies)"),
  list(a ="carnitine deficiency", b = "carnitine deficiency (growth disorders/Spinal Muscular Atrophies)"),
  list(a ="Progeria", b = "progeria (growth disorders/congenital conditions/aging)"),
  list(a ="progeria", b = "progeria (growth disorders/congenital conditions/aging)"),
  list(a ="Anemia", b = "Anemia (insufficient healthy red blood cells)"),
  list(a ="anemia", b = "anemia (insufficient healthy red blood cells)"), 
  list(a ="Osteogenesis imperfecta", b = "osteogenesis imperfecta (osteoporosis)"), 
  list(a ="osteogenesis imperfecta", b = "osteogenesis imperfecta (osteoporosis)"), 
  list(a ="duchenne's", b = "duchenne's (growth disorders)"),
  list(a ="Duchenne's", b = "Duchenne's (growth disorders)"),
  list(a ="insulin resistance", b = "insulin resistance (diabetes)"),
  list(a ="Insulin resistance", b = "insulin resistance (diabetes)"),
  list(a ="biliary cirrhosis", b = "biliary cirrhosis (hepatic disorder)"),
  list(a ="Biliary cirrhosis", b = "biliary cirrhosis (hepatic disorder)"),
  list(a ="Encephalitis", b = "encephalitis (Vaccines (Infectious Disease))"),
  list(a ="encephalitis ", b = "encephalitis (Vaccines (Infectious Disease))"),
  list(a ="duchenne ", b = "duchenne (growth disorders)"),
  list(a ="Duchenne ", b = "Duchenne (growth disorders)"),
  list(a ="Chronic pain", b = "chronic pain (nociceptive/neuropathic)"),
  list(a ="chronic pain", b = "chronic pain (nociceptive/neuropathic)"),
  list(a ="Purpura", b = "purpura (skin hemorrhages)"),
  list(a ="purpura", b = "purpura (skin hemorrhages)"),
  list(a ="Turner syndrome", b = "turner syndrome (growth disorders)"),
  list(a ="Turner Syndrome", b = "turner syndrome (growth disorders)"),
  list(a ="nephritis", b = "nephritis (renal disease/inflammation)"),
  list(a ="Nephritis", b = "nephritis (renal disease/inflammation)"),
  list(a ="Lymphoma", b = "Lymphoma (lymphocytic leukaemia)"),
  list(a ="GERD", b = "Gastroesophageal reflux disease (GERD)"),
  list(a ="graft ", b = "graft (organ transplant)"),
  list(a ="Graft ", b = "graft (organ transplant)"), 
  list(a ="graft-", b = "graft (organ transplant)"),
  list(a ="Graft-", b = "graft (organ transplant)"),
  list(a ="veno-occlusive", b = "veno-occlusive (vascular)"),
  list(a ="veno occlusive", b = "veno occlusive (vascular)"),
  list(a ="veno-occlusion", b = "veno-occlusion (vascular)"),
  list(a ="veno occlusion", b = "veno occlusion (vascular)"),
  list(a ="hepatic", b = "hepatic (liver)"), 
  list(a ="Uveitis", b = "uveitis (Ophthalmology)"), 
  list(a ="uveitis", b = "uveitis (Ophthalmology)"),
  list(a ="kidney disease", b = "kidney disease (renal Disease)"), 
  list(a ="Kidney disease", b = "kidney disease (renal Disease)"), 
  list(a ="congenital muscular dystrophy", b = "congenital muscular dystrophy (growth disorders)"), 
  list(a ="Congenital muscular dystrophy", b = "congenital muscular dystrophy (growth disorders)"), 
  list(a ="acute rejection", b = "acute rejection (graft loss)"),
  list(a ="Acute rejection", b = "acute rejection (graft loss)"),
  list(a ="transplant rejection", b = "transplant rejection (graft loss)"),
  list(a ="Transplant rejection", b = "transplant rejection (graft loss)"),
  list(a ="circadian rhythm disorders", b = "circadian rhythm disorders (Insomnia)"),
  list(a ="Circadian rhythm disorders", b = "circadian rhythm disorders (Insomnia)"),
  list(a ="lennox-gastaut syndrome", b = "lennox-gastaut syndrome (epilepsy)"),
  list(a ="lennox gastaut syndrome", b = "lennox-gastaut syndrome (epilepsy)"),
  list(a ="Lennox-Gastaut syndrome", b = "lennox-gastaut syndrome (epilepsy)"),
  list(a ="Lennox-Gastaut Syndrome", b = "lennox-gastaut syndrome (epilepsy)"),
  list(a ="seizures", b = "seizures (epilepsy)"),
  list(a ="Seizures", b = "seizures (epilepsy)"),
  list(a ="organ rejection", b = "organ rejection (Transplantation/graft versus host)"),
  list(a ="Organ rejection", b = "organ rejection (Transplantation/graft versus host)"),
  list(a ="Acute rejection", b = "Acute rejection (Transplantation/graft versus host)"),
  list(a ="acute rejection", b = "acute rejection (Transplantation/graft versus host)"),
  list(a = "ESRD", b = "end stage renal disease (ESRD)"),
  list(a = "esophagus", b = "esophagus (esophageal)"),
  list(a = "nephrotic", b = "nephrotic (kidney/renal"),
  list(a = "Nephrotic", b = "nephrotic (kidney/renal"),
  list(a = "Renal", b = "Renal (kidney)"),
  list(a = "renal", b = "renal (kidney)"),
  list(a = "myelofibrosis", b = "myelofibrosis (leukemia)"),
  list(a = "Myelofibrosis", b = "myelofibrosis (leukemia)"),
  list(a = "myeloproliferative", b = "myeloproliferative (leukemia)"),
  list(a = "Myeloproliferative", b = "myeloproliferative (leukemia)"),
  list(a = "Melanoma", b = "melanoma (skin)"),
  list(a = "melanoma", b = "melanoma (skin)"),
  list(a = "short bowel syndrome", b = "short bowel syndrome (e.g. Crohn's disease/necrotising enterocolitis/congenital)"),
  list(a = "Short Bowel Syndrome", b = "short bowel syndrome (e.g. Crohn's disease/necrotising enterocolitis/congenital)"),
  list(a = "Glioma", b = "glioma (brain/spinal cord/CNS cancer)"),
  list(a = "glioma", b = "glioma (brain/spinal cord/CNS cancer)"),
  list(a = "Inhalation Anthrax", b = "inhalation anthrax (Respiratory Infections)"),
  list(a = "Inhalation anthrax", b = "inhalation anthrax (Respiratory Infections)"),
  list(a = "inhalation anthrax", b = "inhalation anthrax (Respiratory Infections)"),
  list(a = "inhalation Anthrax", b = "inhalation anthrax (Respiratory Infections)"),
  list(a = "esrd", b = "end stage renal disease (ESRD)"),
  list(a = "wilson's disease", b = "wilson's disease (metabolic/liver disorder)"),
  list(a = "Wilson's disease", b = "wilson's disease (metabolic/liver disorder)"),
  list(a = "MDS", b = "Myelodysplastic syndromes (MDS)"),
  list(a = "Leukaemia", b = "leukemia"),
  list(a = "leukaemia", b = "leukemia"))

remove_countries <- function(x){
    
    future_map(x, function(y){
      for (j in 1:length(countries)){
        y = str_remove_all(tolower(y), countries[j])
      }
      return(y)
      }, .options=future_options(globals=c("countries"), packages=c("stringr"))) %>% 
      as_vector() %>%
      trimws()
    #for (i in 1:length(tmp)) {
      
    #}
    
    #return(trimws(tmp))
}

term_replacement <- function(x){
    #tmp = x
    future_map(x, function(y) {
      for (j in 1:length(replace_list)){
        y = str_replace_all(y, coll(replace_list[[j]]$a), coll(replace_list[[j]]$b) )
      }
      return(y)
    }, .options=future_options(globals=c("replace_list"), packages=c("stringr"))) %>% 
      as_vector() %>%
      trimws()
    
    # for (i in 1:length(tmp)) {
    #   for (j in 1:length(replace_list)){
    #     tmp[i] = str_replace_all(tmp[i], coll(replace_list[[j]]$a), coll(replace_list[[j]]$b) )
    #   }
    # }
    # 
    # return(trimws(tmp))
}

plan(multiprocess) #plan(tweak(multiprocess, workers=8))

events <- read_delim(file="/shared_data/data/raw/MIT_data_full_190723/pharmaprojects_Events.csv",
           delim="|") %>%
  filter(HCode==117) %>%
  # mutate(Information = tolower(trimws(str_extract(Information, "[A-Z][A-Za-z ]++[;,]")))) %>%
  # group_by(Information) %>% summarize(n=n()) %>% arrange(desc(n)) %>%
  # .$Information
  # mutate(Information = trimws(tolower( str_replace_all(Information, "[^[:alnum:]]", " ")))) %>%
  mutate(Information = ifelse(is.na(Information), "unspecified", ifelse(trimws(Information)=="", "unspecified", Information))) %>%
  mutate(Information=term_replacement(Information)) %>%
  mutate(Information=remove_countries(Information)) %>%
  #group_by(Information) %>% summarize(n=n()) %>% arrange(desc(n)) %>% print(n=300)
  #mutate(Information = trimws(str_remove(Information, "[A-Z][A-Za-z ]++[;,]"))) %>%
  left_join(dplyr::select(read_csv("/shared_data/data/training_data/training_data_2015_split_on_outcome.csv"), DrugKey, indicationkey), by="DrugKey") %>%
  filter(!is.na(indicationkey)) 

# events %>%
#   dplyr::select(Information) %>%
#   distinct() %>%
#   print(n=865)

  
# events %>%
#   group_by(Information) %>%
#   summarize(n=n()) %>%
#   arrange(desc(n)) %>%
#   print(n=100)

# We need just not the longest substring, but also second longest.
#longest_string <- function(s){return(s[which.max(nchar(s))])}

longest_common_substring <- function(a,b){
  if (is.na(a) | is.na(b) | a=="" | b==""){
    return("")
  } else {
    ## get all forward substrings of 'b'
    sb <- as_vector( map(1:nchar(b), function(x) stri_sub(b, x, x:nchar(b))) )
    ## extract them from 'a' if they exist
    sstr <- unique( trimws(na.omit(str_extract(a, coll(sb)))) )
    ## match the longest one
    return( sstr[which.max(nchar(sstr))] )
  }
}

# longest_common_substring("abc", "def")
# longest_common_substring("oncology  soft tissue sarcoma cns  glioblastoma ", "cancer (oncology)  myeloma")
# a = "oncology  soft tissue sarcoma cns  glioblastoma "
# b = "cancer (oncology)  myeloma"
# as_vector( map(1:nchar(b), function(x) stri_sub(b, x, x:nchar(b))) )

# Example shows that longest substring just does not work!
# longest_common_substring(a="oncology  soft tissue sarcoma cns glioblastoma b",
#                          b="b glioblastoma multiforme")



lcsbstr_no_lib <- function(a,b) { 
  if (is.na(a) | is.na(b)) {
    return(0)
  } else {
    # print(str_replace_all(a, "[\\[\\]\\|\\,\\.]", " "))
    # print(str_replace_all(b, "[\\[\\]\\|\\,\\.]", " "))
    # print(longest_common_substring(str_replace_all(a, "[\\[\\]\\|\\,\\.]", " "),
    #                                str_replace_all(b, "[\\[\\]\\|\\,\\.]", " ")))
               
    longest_cmn_sbstr = longest_common_substring(str_replace_all(a, "[\\[\\]\\|\\,\\.]|oncology|cancer|cardiovascular|infectious disease| disease", " "),
                                        str_replace_all(b, "[\\[\\]\\|\\,\\.]|oncology|cancer|cardiovascular|infectious disease| disease", " "))
    if (identical(longest_cmn_sbstr,character(0))) {
      return(0)
    } else {
      return(str_length(longest_cmn_sbstr))
    }
  }
}


# lcsbstr_no_lib(tolower("Bla Glioblastoma blee "),
#                tolower("Blub glioblastoma multiforme"))
# 
# longest_common_substring(str_replace_all(tolower("Oncology [Prostate]"), "[\\[\\]\\|\\,\\.oncology]", " "),
#                          str_replace_all(tolower("cancer (oncology), myeloma"), "[\\[\\]\\|\\,\\.oncology]", " "))
# 
# lcsbstr_no_lib(tolower("Oncology [Prostate]"),
#               tolower("cancer (oncology), myeloma"))
#lcsbstr_no_lib(tolower("Oncology [Renal]"), tolower("cancer (oncology), myeloma"))



lcsbstr2_no_lib <- function(a,b) { 
  if (is.na(a) | is.na(b)) {
    return(0)
  } else {
    a2 = str_replace_all(a, "[\\[\\]\\|\\,\\.]|oncology|cancer|cardiovascular|infectious disease| disease", " ")
    b2 = str_replace_all(b, "[\\[\\]\\|\\,\\.]|oncology|cancer|cardiovascular|infectious disease| disease", " ")
    longest_cmn_sbstr  <- longest_common_substring(a2,b2)
    if (identical(longest_cmn_sbstr,character(0))) {
      return(0)
    } else {
      #print(longest_cmn_sbstr)
      a2 = str_remove(a2, coll(longest_cmn_sbstr))
      b2 = str_remove(b2, coll(longest_cmn_sbstr))
      longest_cmn_sbstr  <- longest_common_substring(a2,b2)
      if (identical(longest_cmn_sbstr,character(0))) {
        return(0)
      } else {
        return(str_length(longest_cmn_sbstr))
      }
    }
  }
}

#lcsbstr2_no_lib("I only wish to understand what is going on.","My understandinging is incomplete.")
# lcsbstr2_no_lib("Truth","Truth")
# lcsbstr_no_lib("Truth","Truth")

# lcsbstr2_no_lib(tolower("Oncology [Renal]"), tolower("cancer (oncology), myeloma"))
# lcsbstr2_no_lib(tolower("Oncology [Multiple Myeloma]"), tolower("cancer (oncology), myeloma"))

# read_csv("/shared_data/data/training_data/training_data_2015_split_on_outcome.csv") %>% 
#   filter(GenericName=="selexipag") %>% #str_detect(strDiseaseType, "Pulmonar")
#   dplyr::select(GenericName, strDiseaseType) %>%
#   distinct()


# tibble(x=c(1,1,1,1,2,3,3,3,4),y=1:length(x)) %>%
#   rowwise() %>%
#   mutate(xx= map(tibble(.), function(x) {
#     print(x) 
#     return(rename(x, xxx=x, yyy=y))
#     }
#     )) %>%
#   unnest(xx)

# read_csv("/shared_data/data/training_data/training_data_2015_split_on_outcome.csv") %>%
#   left_join(dplyr::select(events, DrugKey, HCode, Information), by="DrugKey") %>%
#   ungroup() %>% mutate() %>% dplyr::select(strtherapeuticarea, oncology) %>%
#   group_by(strtherapeuticarea) %>%
#   summarize(sum(oncology), mean(oncology)) %>% print(n=300)

orphmatches <- fulldat %>%
  dplyr::select(DrugKey, strtherapeuticarea, strDiseaseType) %>%
  left_join(dplyr::select(events, DrugKey, HCode, Information), by="DrugKey") %>%
  ungroup() %>%
  mutate(oncology = str_detect(strtherapeuticarea, "Oncology"),
         strDiseaseType=term_replacement(strDiseaseType),
         tadis = paste(strtherapeuticarea, strDiseaseType)) %>%
  dplyr::select(DrugKey, tadis, strDiseaseType, oncology, HCode, Information) %>% #strther, 
  distinct() %>%
  group_by(DrugKey) %>%
  mutate(anyorph=max(!is.na(HCode))) %>%
  filter(anyorph==1) %>%
  ungroup()%>%
  mutate(
    oncorph = str_detect(tolower(Information), 
              "cancer|leukaemia|leukemia|lymphoma|myelodysplastic|carcinoma|radiation-induced|chemotherapy|myeloma|neoplasm|melanoma|malignant|metasta|lymphoblast|acute myeloid|hematopoietic stem|myeloprolif|myelofibro|malignancy|blastoma|sarcoma| tumour|mesothelioma| tumor")) %>% 
  mutate(    bothlupus = ( str_detect(tolower(Information), "lupus") & str_detect(tolower(tadis), "lupus") ) ) %>% 
  mutate(
    oncmismatch = (!oncology & oncorph),
    mismatch2 = ( (tadis=="Autoimmune/Inflammation [Chronic Obstructive Pulmonary Disease]" & str_detect(Information, "ophthalmology|bronchopulmonary dysplasia|hypertension\\, pulmonary")) |
                    ( oncorph & str_detect(Information, "purpura")) |
                    ( oncorph & str_detect(Information, "progeria")) |
                    ( str_detect(tadis, "Arthritis|arthritis") & str_detect(Information, "myopathy|myopathies|purpura|familial cold|neonatal onset multisystem|anaemia")) |
                    ( str_detect(tadis, "Transplantation|Ulcerative Colitis|Acute Coronary Syndromes|Restless Legs Syndrome") & 
                      str_detect(Information, "demyelinating polyneuropathy|inflammation, muscle|respiratory distress|purpura|tourette")) |
                    ( str_detect(tadis, "Skin Infection") & str_detect(Information, "lung infection")) |
                    ( str_detect(tadis, "Respiratory Infections") & str_detect(Information, "cryptosporidiosis")) |
                    ( str_detect(tadis, "Hepatic Fibrosis") & str_detect(Information, "fibrosis\\, pulmonary")) |
                    ( str_detect(tadis, "Glioblastoma") & str_detect(Information, coll("stomach (gastric)")))) ) %>% 
  mutate(
    bonus = ((str_detect(tadis, "Breast") & str_detect(Information, "breast")) |
               (str_detect(tadis, "Ovarian") & str_detect(Information, "ovarian")) |
               (str_detect(tadis, "Bladder") & str_detect(Information, "bladder")) ) ) %>% 
  mutate(
    painmatch = (str_detect(tadis, "Pain (neuropathic) pain") & str_detect(Information, coll("post-herpetic "))),
    sepsismatch = (str_detect(tadis, "[Ss]epsis") & str_detect(Information, "sepsis")) ) %>% 
  mutate(
    delta= stringdist(a=tolower(tadis), 
                      b=tolower(Information), method="jw") + 
      stringdist(a=tolower(strDiseaseType), 
                 b=tolower(Information), method="jaccard") + 
      stringdist(a=tolower(strDiseaseType), 
                 b=tolower(Information), method="cosine")) %>%
  mutate(
    #lcs = (str_length(strDiseaseType) + str_length(Information) - stringdist(a=tolower(strDiseaseType), b=tolower(Information), method="lcs") )/2,
    #qg= stringdist(a=tolower(strDiseaseType), b=tolower(Information), method="qgram"),
    # lcs2 = future_map2(tolower(tadis), tolower(Information), lcsbstr_no_lib,
    #                    .options=future_options(globals=c("lcsbstr_no_lib", "longest_common_substring"), 
    #                                            packages=c("tidyverse", "stringi", "stringr", "stringdist"))))


                #                            .options=future_options(globals=c("lcsbstr_no_lib", "lcsbstr2_no_lib", "longest_common_substring"), 
                #                                                    packages=c("tidyverse", "stringi", "stringr", "stringdist")) )) 
                
    lcs2 = ifelse(oncmismatch | mismatch2, 0,
                  lcsbstr_no_lib(tolower(tadis),
                                 tolower(Information))),
    lcs3 = ifelse(oncmismatch | mismatch2, 0, 
                  lcsbstr2_no_lib(tolower(tadis),
                                  tolower(Information))),
    crit2 =  (lcs2 >0.275 | delta<0.9),
    sumcrit = 2.0*bothlupus + lcs2/10 + lcs3/10 - delta + bonus*1) %>% 
  # mutate( scomp = future_map(tibble(.), string_comparisons, 
  #                            .options=future_options(globals=c("lcsbstr_no_lib", "lcsbstr2_no_lib", "longest_common_substring"), 
  #                                                    packages=c("tidyverse", "stringi", "stringr", "stringdist")) )) 
  ungroup() %>%
  group_by(DrugKey, Information) %>%
  arrange(DrugKey, Information, desc(lcs2), desc(lcs3), delta) %>%
  mutate(ordno = row_number(),
         bestmatch=ifelse(ordno==1, ifelse(( sepsismatch | painmatch | bothlupus | (lcs2==5 & sumcrit>-0.265 & (oncology & !oncmismatch)) | (lcs2>5 & ( sumcrit>-0.3095 | (oncology & !oncmismatch & sumcrit>-0.5))) |
                                               (lcs2>6 & sumcrit>-0.5) | (lcs2>7 & sumcrit>-0.85) ), 1, 0), 0) #, #bestmatch=ifelse(sum(bestmatch)==0, 0.5, bestmatch)
         ) %>%  #arrange(DrugKey, Information, desc(sumcrit)) #%>%  #filter(bestmatch>0) %>%  print(n=100)
  ungroup()

write_rds(orphmatches, "/files/orphmatches1.rds", compress = "bz2")
#check sumcrit<0, perhaps allow >-0.5, if first word match>5?
# great = sumcrit>0.5

# All matches with sumcrit>1 seem great
# between 0.5 and 1:Irritable Bowel Syndrome != short bowel syndrome
#  -? [Diabetic Complications|Diabetic Retinopathy] = Uveitis ???

# meematches %>%
#   group_by(DrugKey, Information) %>%
#   mutate(bestscore = max(sumcrit)) %>%
#   filter(bestscore>2.1 & ordno<3 & max(bestmatch)==1) %>%
#   #arrange(bestscore, DrugKey, Information, ordno) %>%
#   dplyr::select(DrugKey, Information, bestscore, ordno, tadis) %>%
#   pivot_wider(values_from = tadis, names_from=ordno, names_prefix = "tadis") %>% print(n=30)
#   write_csv("/files/orphan.csv")

orphmatched <- orphmatches %>%
  mutate(orphan_indication = (bestmatch==1 & sumcrit>=0.7)*1.0 + (bestmatch==1 & sumcrit<0.7)*0.75,
         any_oncology_orphan_designation = anyorph*oncorph,
         any_nononc_orphan_designation = anyorph*(1-oncorph),
         ) %>%
  group_by(DrugKey,strDiseaseType) %>%
  summarize(any_orphan_designation = max(anyorph,na.rm=T),
            orphan_indication = max(orphan_indication,na.rm=T),
            any_oncology_orphan_designation=max(any_oncology_orphan_designation,na.rm=T),
            any_nononc_orphan_designation = max(any_nononc_orphan_designation, na.rm=T))

write_rds(orphmatched, "/files/orphmatched.rds", compress = "bz2")
# %>%
#   filter(oncology == 1 & any_oncology_orphan_designation == 0) %>% print(n=600)


# meematches %>%
#   filter(bestmatch==1 & sumcrit<0) %>%
#   dplyr::select(DrugKey, tadis, Information, delta, lcs2, lcs3, sumcrit) %>%
#   arrange(tadis,desc(sumcrit)) %>%
#   print(n=500)
# 
# meematches %>%
#   filter(ordno==1 & bestmatch!=1) %>%
#   dplyr::select(DrugKey, tadis, Information, delta, lcs2, lcs3, sumcrit) %>%
#   arrange(tadis,desc(sumcrit)) %>%
#   print(n=500)
# 
# 
# meematches %>% 
#   filter(DrugKey %in% c(5510,   14862,  17878,  19039,  23267,  28103,  29888,  30206,  30740,  30816,  33085,  36354 ,  37697) ) %>%
#   dplyr::select(DrugKey, tadis, Information, bestmatch, ordno, delta, lcs2, lcs3, sumcrit) %>%
#   print(n=500)
#   
# lcsbstr_no_lib(tolower("Oncology Soft Tissue Sarcoma CNS Glioblastoma"),
#                tolower("glioblastoma multiforme"))
# 
# lcsbstr_no_lib(tolower("Infection-associated haemolytic uraemic syndrome"), tolower("[Age-Related Macular Degeneration]"))
# 
# str_length("Infection-associated haemolytic uraemic syndrome")
# str_length("[Age-Related Macular Degeneration]")
# stringdist(a=tolower("Infection-associated haemolytic uraemic syndrome"), b=tolower("[Age-Related Macular Degeneration]"), method="qgram")
# 
# read_csv("/shared_data/data/test_data_full/testing_phase2_release.csv") %>%
#   dplyr::select(decMinAge) %>%
#   filter(str_length(decMinAge)>3)
# 
# 
# /shared_data/data/raw/MIT_data_full_190723/pharmaprojects_MarketSizeRating.csv
# 
# /shared_data/data/raw/MIT_data_full_190723/sitetrove_tblTrialResults.csv