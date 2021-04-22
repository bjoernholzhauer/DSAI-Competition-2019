library(tidyverse)
library(uwot)

# Auxilliary functions

# Function to add column
fncols <- function(data, cname, what_entry=0) {
  add <-cname[!cname%in%names(data)]
  
  if(length(add)!=0) data[add] <- what_entry
  data
}

# Function to calculate geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


################################################
# UMAP embedding for indications
################################################

indumap <- function(fulldat){
  
  message("-----  Indication UMAP  -----")
  # Check whether we are on evaluation sever or training server, pre-pend path accordingly
  if (Sys.getenv("AICROWD_TRAIN_DATA_PATH")==""){
    filepath = paste0("/home/desktop1/files/", "umapind.rds")
  } else {
    filepath = "umapind.rds"
  }
  
  if (file.exists(filepath)){
    ctime = file.info(filepath)$ctime
    message( paste("Found existing umap-embedding on disk, using that. Created at ", ctime, ", if that is too old, delete it.") )
    rval <- read_rds(filepath)
  } else {
    message("No existing umap embedding on disk, creating a new one - will not reproduce old results (this is random, make sure to read in the umapind.rds to obtain consistent results.")
    
    dis1 <- fulldat %>%
      dplyr::select(indicationkey, strDiseaseType) %>% #strtherapeuticarea,  outcome
      mutate(strDiseaseType=str_remove_all(str_remove_all(str_remove_all(str_remove_all(strDiseaseType, "\\(N/A\\)\\|"),"\\|\\(N/A\\)"),"\\(\\)\\|"),"\\|\\(\\)")) %>%
      mutate(distyp = str_split(strDiseaseType, pattern="\\|")) %>%
      unnest(cols="distyp") %>%
      distinct() %>%
      mutate(distyp=str_remove_all(distyp, "[\\[\\]]"),
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
               TRUE ~ 11)) %>%
      dplyr::select(indicationkey, disno) %>%
      mutate(anytyp=1) %>%
      distinct() %>%
      group_by(indicationkey) %>%
      pivot_wider(names_from = disno, names_prefix = "disease", values_from = anytyp, id_cols = c("indicationkey")) %>%  
      replace(., is.na(.), 0)
    
    tastuff <- fulldat %>%
      dplyr::select(indicationkey, strtherapeuticarea, DrugKey) %>%
      group_by(indicationkey) %>% 
      mutate(sta = str_split(strtherapeuticarea, pattern="\\|"), # Therapeutic areas
             number_of_projects=length(unique(DrugKey))) %>% 
      unnest(cols="sta") %>%
      ungroup() %>%
      mutate(
        intsta = case_when(
          sta=="Autoimmune/Inflammation" ~ 1,
          sta=="Cardiovascular" ~ 2,
          sta=="CNS" ~ 3,
          sta=="Genitourinary" ~ 4,
          sta=="Infectious Disease" ~ 5,
          sta=="Metabolic/Endocrinology" ~ 6,
          sta=="Oncology" ~ 7,
          sta=="Ophthalmology" ~ 8,
          sta=="Vaccines (Infectious Disease)" ~ 9,
          TRUE ~ 10),
        flg=2) %>%
      dplyr::select(indicationkey, intsta, number_of_projects, flg) %>%
      distinct() %>%
      group_by(indicationkey, number_of_projects) %>%
      pivot_wider(names_from=intsta, names_prefix="ta", values_from = flg,
                  values_fill=list(flg=0)) %>%
      ungroup() %>%
      fncols(paste0("ta", 1:9)) %>%
      mutate(number_of_projects=sqrt(number_of_projects))#,sqrt_sumt_patno=(sumt_patno)^(1/3))
    
    ages <- fulldat %>%
      dplyr::select(indicationkey, strMinAgeUnit, strMaxAgeUnit, decMinAge, decMaxAge) %>%
      mutate(veryyoung = 1 *(trimws(tolower(strMinAgeUnit)) %in% c("days","weeks","months") | 
                               trimws(tolower(strMaxAgeUnit))  %in% c("days", "weeks", "months")),
             peds = 1*(veryyoung==1 | ifelse(is.na(decMinAge), FALSE, (decMinAge<18)) | ifelse(is.na(decMaxAge), FALSE, (decMaxAge<18))),
             above65 = 1*( ( (trimws(tolower(strMinAgeUnit))=="years" | is.na(strMinAgeUnit)) & ifelse(is.na(decMinAge), FALSE, (decMinAge>65))) |
                             ((trimws(tolower(strMaxAgeUnit))=="years" | is.na(strMaxAgeUnit)) & ifelse(is.na(decMaxAge), FALSE, (decMaxAge>65))) ),
             above75 = 1*( ((trimws(tolower(strMinAgeUnit))=="years" | is.na(strMinAgeUnit)) & ifelse(is.na(decMinAge), FALSE, (decMinAge>75))) |
                             ((trimws(tolower(strMaxAgeUnit))=="years" | is.na(strMaxAgeUnit)) & ifelse(is.na(decMaxAge), FALSE, (decMaxAge>75))) ),
             noextreme = 1*(veryyoung==0 & peds==0 & above65==0 & above75==0)
      ) %>%
      group_by(indicationkey) %>%
      summarize(veryyoung=mean(veryyoung),
                peds=mean(peds),
                above65=mean(above65),
                above75=mean(above75),
                noextreme=mean(noextreme)) %>%
      ungroup()
    
    ddrdat <- fulldat %>% 
      dplyr::select(indicationkey, DrugDeliveryRouteDescription) %>%
      mutate(ddrd = str_split(DrugDeliveryRouteDescription, pattern="\\|")) %>%
      unnest(cols="ddrd") %>%
      mutate(ddrd=str_remove_all(ddrd, "[\\[\\]]"),
             ddrdg = case_when(
               str_detect(tolower(ddrd),"subcu") | str_detect(tolower(ddrd),"intraven") | str_detect(tolower(ddrd),"intra-articul") |
                 str_detect(tolower(ddrd),"intradermal") | str_detect(tolower(ddrd),"intra-arterial") | str_detect(tolower(ddrd),"intraarterial") |
                 str_detect(tolower(ddrd),"intramusc") | str_detect(tolower(ddrd),"subcu") | trimws(tolower(ddrd))=="injectable" | trimws(tolower(ddrd))=="injected" ~ 1,
               str_detect(tolower(ddrd),"intravesical") | str_detect(tolower(ddrd),"intralymphatic") | str_detect(tolower(ddrd),"intracardiac") |
                 str_detect(tolower(ddrd),"intraperitoneal") | str_detect(tolower(ddrd),"intracavity") | str_detect(tolower(ddrd),"intratumoral")  ~ 2,
               str_detect(tolower(ddrd),"intracerebral") | str_detect(tolower(ddrd),"intraspinal") | 
                 str_detect(tolower(ddrd),"intracranial") | str_detect(tolower(ddrd),"intrathecal") ~ 3,
               str_detect(tolower(ddrd),"ocular") | str_detect(tolower(ddrd),"intravitreal") | str_detect(tolower(ddrd),"ophthalmo") ~ 4,
               str_detect(tolower(ddrd),"oral") ~ 5,
               str_detect(tolower(ddrd),"inhal") ~ 6,
               str_detect(tolower(ddrd),"topical") | str_detect(tolower(ddrd),"transderm") |  
                 str_detect(tolower(ddrd),"topical, mucosal") | str_detect(tolower(ddrd),"otic") ~ 7,
               str_detect(tolower(ddrd),"depot") | str_detect(tolower(ddrd),"implant") ~ 8,
               str_detect(tolower(ddrd),"multiple") ~ 9,
               str_detect(tolower(ddrd),"rectal") | str_detect(tolower(ddrd),"vaginal") ~ 10,
               TRUE ~ 0
             ),
             flg=1) %>%
      group_by(indicationkey, ddrdg) %>%
      summarize(flags=sum(flg)) %>%
      ungroup() %>%
      group_by(indicationkey) %>%
      pivot_wider(names_from=ddrdg, values_from = flags, names_prefix = "ddr") %>%
      replace(., is.na(.), 0) %>%
      fncols(paste0("ddr", 0:10)) %>%
      mutate(thesum=sum(ddr0, ddr1, ddr5, ddr7, ddr8, ddr6, ddr3, ddr2, ddr9, ddr10, ddr4),
             ddr0=ddr0/thesum, 
             ddr1=ddr1/thesum, 
             ddr5=ddr5/thesum, 
             ddr7=ddr7/thesum, 
             ddr8=ddr8/thesum, 
             ddr6=ddr6/thesum, 
             ddr3=ddr3/thesum, 
             ddr2=ddr2/thesum, 
             ddr9=ddr9/thesum, 
             ddr10=ddr10/thesum, 
             ddr4=ddr4/thesum) %>%
      ungroup() %>%
      dplyr::select(-thesum)
    
    eptype <- fulldat %>%
      dplyr::select(indicationkey, strPrimaryEndPoint) %>%
      mutate(mortality = 1*(str_detect(tolower(strPrimaryEndPoint), "death") | str_detect(tolower(strPrimaryEndPoint), "mortality") | str_detect(tolower(strPrimaryEndPoint), "survival")),
             transpl = 1*(str_detect(tolower(strPrimaryEndPoint), "reject") | str_detect(tolower(strPrimaryEndPoint), "graft") | 
                            str_detect(tolower(strPrimaryEndPoint), "bpar") | str_detect(tolower(strPrimaryEndPoint), "transplant")),
             safety = 1*(str_detect(tolower(strPrimaryEndPoint), "safety") | str_detect(tolower(strPrimaryEndPoint), "adverse")  | str_detect(tolower(strPrimaryEndPoint), "dose-limit") | 
                           str_detect(tolower(strPrimaryEndPoint), "toxicity") | str_detect(tolower(strPrimaryEndPoint), "dlt") 
                         | str_detect(tolower(strPrimaryEndPoint), "grade") | str_detect(tolower(strPrimaryEndPoint), "mtd") | str_detect(tolower(strPrimaryEndPoint), "tolerab") ),
             progression = 1*(str_detect(tolower(strPrimaryEndPoint), "progression") | str_detect(tolower(strPrimaryEndPoint), "onset") 
                              | str_detect(tolower(strPrimaryEndPoint), "new diagno") | str_detect(tolower(strPrimaryEndPoint), "(cr)") | str_detect(tolower(strPrimaryEndPoint), "recist") | 
                                str_detect(tolower(strPrimaryEndPoint), "objective resp") | str_detect(tolower(strPrimaryEndPoint), "complete resp")),
             visual = 1*(str_detect(tolower(strPrimaryEndPoint), "visual") | str_detect(tolower(strPrimaryEndPoint), "acuity") ),
             bp = 1*(str_detect(tolower(strPrimaryEndPoint), "pressure") | str_detect(tolower(strPrimaryEndPoint), "sbp") | str_detect(tolower(strPrimaryEndPoint), "dbp") ),
             resp = 1*(str_detect(tolower(strPrimaryEndPoint), "respon") | str_detect(tolower(strPrimaryEndPoint), "threshold") | str_detect(tolower(strPrimaryEndPoint), "achiev") ),
             bleed = 1*(str_detect(tolower(strPrimaryEndPoint), "bleed") | str_detect(tolower(strPrimaryEndPoint), "hemorr") | str_detect(tolower(strPrimaryEndPoint), "haemorr") ),
             pro = 1*(str_detect(tolower(strPrimaryEndPoint), "patient report") | str_detect(tolower(strPrimaryEndPoint), "patient-report") | 
                        str_detect(tolower(strPrimaryEndPoint),"e-diary") | str_detect(tolower(strPrimaryEndPoint),"ediary") | str_detect(tolower(strPrimaryEndPoint),"qol") | 
                        str_detect(tolower(strPrimaryEndPoint),"quality of") | str_detect(tolower(strPrimaryEndPoint),"quality-of") | str_detect(tolower(strPrimaryEndPoint),"symptom score") | 
                        str_detect(tolower(strPrimaryEndPoint),"questionnair") | str_detect(tolower(strPrimaryEndPoint),"pain") | str_detect(tolower(strPrimaryEndPoint),"rating") ),
             occur = 1*(str_detect(tolower(strPrimaryEndPoint), "occur") | str_detect(tolower(strPrimaryEndPoint), "event") | str_detect(tolower(strPrimaryEndPoint), "observ") ),
             pain = 1*str_detect(tolower(strPrimaryEndPoint),"pain"),
             efficacy = 1*(str_detect(tolower(strPrimaryEndPoint),"efficac") | str_detect(tolower(strPrimaryEndPoint),"effectiv")),
             spiro = 1*(str_detect(tolower(strPrimaryEndPoint),"fev1")|str_detect(tolower(strPrimaryEndPoint),"pef") |str_detect(tolower(strPrimaryEndPoint),"fvc") |
                          str_detect(tolower(strPrimaryEndPoint),"spirometr") |str_detect(tolower(strPrimaryEndPoint),"expira")),
             brain = 1*(str_detect(tolower(strPrimaryEndPoint),"lesions")|str_detect(tolower(strPrimaryEndPoint),"brain mri")| str_detect(tolower(strPrimaryEndPoint)," brain ")),
             mri = 1*(str_detect(tolower(strPrimaryEndPoint),"magnetic resonan") | str_detect(tolower(strPrimaryEndPoint)," mri") | str_detect(tolower(strPrimaryEndPoint),"\\(mri")),
             gluc = 1*(str_detect(tolower(strPrimaryEndPoint), "hba1c") | str_detect(tolower(strPrimaryEndPoint), "glucose") | 
                         str_detect(tolower(strPrimaryEndPoint), "fpg") | str_detect(tolower(strPrimaryEndPoint), "ogtt")
                       | str_detect(tolower(strPrimaryEndPoint), "hypogly") | str_detect(tolower(strPrimaryEndPoint), "insulin"))) %>%
      replace(., is.na(.), 0) %>%
      group_by(indicationkey) %>%
      summarize(mort=qbeta(p=0.5,1+sum(mortality),1+sum(1-mortality)),
                safety=qbeta(p=0.5,1+sum(safety),1+sum(1-safety)),
                prog=qbeta(p=0.5,1+sum(progression),1+sum(1-progression)),
                visual=qbeta(p=0.5,1+sum(visual),1+sum(1-visual)),
                bp=qbeta(p=0.5,1+sum(bp),1+sum(1-bp)),
                transpl=qbeta(p=0.5,1+sum(transpl),1+sum(1-transpl)),
                resp=qbeta(p=0.5,1+sum(resp),1+sum(1-resp)),
                bleed=qbeta(p=0.5,1+sum(bleed),1+sum(1-bleed)),
                pro=qbeta(p=0.5,1+sum(pro),1+sum(1-pro)),
                occur=qbeta(p=0.5,1+sum(occur),1+sum(1-occur)),
                pain=qbeta(p=0.5,1+sum(pain),1+sum(1-pain)),
                efficacy=qbeta(p=0.5,1+sum(efficacy),1+sum(1-efficacy)),
                spiro=qbeta(p=0.5,1+sum(spiro),1+sum(1-spiro)),
                brain=qbeta(p=0.5,1+sum(brain),1+sum(1-brain)),
                mri=qbeta(p=0.5,1+sum(mri),1+sum(1-mri)),
                gluc=qbeta(p=0.5,1+sum(gluc),1+sum(1-gluc)))
    
    prepro <- fulldat %>%
      dplyr::select(indicationkey, strSponsor, strSponsorType, strdesignkeyword) %>%
      rowwise() %>%
      mutate(nspons = length(str_split(strSponsor, pattern="\\|")[[1]]), 
             nsponst = length(str_split(strSponsorType, pattern="\\|")[[1]]),
             spons = list( str_split(strSponsor, pattern="\\|")[[1]][ c(1:nspons, rep(nspons, max(nspons,nsponst)-nspons)) ] ), 
             sponst = list( str_split(strSponsorType, pattern="\\|")[[1]][ c(1:nsponst, rep(nsponst, max(nspons,nsponst)-nsponst)) ] ),
             rando = 1*str_detect(strdesignkeyword, "random"),
             efficacy = 1*str_detect(strdesignkeyword, "efficacy"),
             safety = 1*str_detect(strdesignkeyword, "safety"),
             pk = 1*str_detect(strdesignkeyword, "pharmacokin"),
             pd = 1*str_detect(strdesignkeyword, "pharmacodyn"),
             drf = 1*str_detect(strdesignkeyword, "dose resp"),
             rubbish = 1*(str_detect(strdesignkeyword, "open")) + 1*(str_detect(strdesignkeyword, "single")) + 1*(str_detect(strdesignkeyword, "ascend")),
             sophis = 1*str_detect(strdesignkeyword, "random") + 1*str_detect(strdesignkeyword, "multiple arm") + 1*str_detect(strdesignkeyword, "umbrella") + 
               1*str_detect(strdesignkeyword, "basket") + 1*str_detect(strdesignkeyword, "efficacy") + 1*str_detect(strdesignkeyword, "super") + 1*str_detect(strdesignkeyword, "infer") +
               1*str_detect(strdesignkeyword, "immunogen") + 1*str_detect(strdesignkeyword, "dose resp") + 1*str_detect(strdesignkeyword, "active") + 
               1*str_detect(strdesignkeyword, "double")  + 1*str_detect(strdesignkeyword, "placebo")) %>%
      ungroup() %>%
      replace(., is.na(.), 0) %>%
      group_by(indicationkey) %>%
      summarize(nspons=log(mean(nspons)+0.5),
                nsponst=log(mean(nsponst)+0.5),
                rando=qbeta(0.5, 1+sum(rando), 1+sum(1-rando)),
                efficacy=qbeta(0.5, 1+sum(efficacy), 1+sum(1-efficacy)),
                safety=qbeta(0.5, 1+sum(safety), 1+sum(1-safety)),
                pk=qbeta(0.5, 1+sum(pk), 1+sum(1-pk)),
                pd=qbeta(0.5, 1+sum(pd), 1+sum(1-pd)),
                drf=qbeta(0.5, 1+sum(drf), 1+sum(1-drf)),
                rubbish=log(mean(rubbish)+0.5),
                sophis=log(mean(sophis)+0.5)) %>%
      ungroup()
    
    
    inddis <- dis1 %>%
      left_join(tastuff, by="indicationkey") %>%
      left_join(ages, by="indicationkey") %>%
      left_join(eptype, by="indicationkey") %>%
      left_join(ddrdat, by="indicationkey") %>%
      left_join(prepro, by="indicationkey") %>%
      ungroup() %>%
      replace(., is.na(.), 0)
    
    # inddis %>% print(n=100)
    # colnames(inddis)[colSums(is.na(inddis)) > 0]
    
    indumap <- umap(dplyr::select(inddis, -indicationkey), 
                    n_neighbors = dim(inddis)[1], 
                    n_epochs = 2500,learning_rate = 0.5, init = "random")
    
    rval <- inddis %>%
      ungroup() %>%
      mutate(ind_umap1 = indumap[,1],
             ind_umap2 = indumap[,2]) %>%
      dplyr::select(indicationkey, ind_umap1, ind_umap2)
    
    write_rds(x=rval, path=filepath)
  }
  
  return( rval )
}

################################################
# UMAP embedding for sponsors
################################################

umap_sponsor <- function(fulldat){
  message("----- Sponsor UMAP ----- ")
  # Check whether we are on evaluation sever or training server, pre-pend path accordingly
  if (Sys.getenv("AICROWD_TRAIN_DATA_PATH")==""){
    filepath = paste0("/home/desktop1/files/", "umapspo.rds")
  } else {
    filepath = "umapspo.rds"
  }
  
  if (file.exists(filepath)){
    ctime = file.info(filepath)$ctime
    message( paste("Found existing umap-embedding on disk, using that. Created at ", ctime, ", if that is too old, delete it.") )
    rval <- read_rds(filepath)
  } else {
    message("No existing umap embedding on disk, creating a new one.")
    
    message("Start generating sponsor embedding")
    # Some TA based features for sponsors 
    tas <- fulldat %>%
      dplyr::select(row_id, strtherapeuticarea) %>%
      mutate(sta = str_split(strtherapeuticarea, pattern="\\|")) %>% 
      unnest(cols="sta") %>%
      group_by(row_id, sta) %>%
      mutate(
        intsta = case_when(
          sta=="Autoimmune/Inflammation" ~ 1,
          sta=="Cardiovascular" ~ 2,
          sta=="CNS" ~ 3,
          sta=="Genitourinary" ~ 4,
          sta=="Infectious Disease" ~ 5,
          sta=="Metabolic/Endocrinology" ~ 6,
          sta=="Oncology" ~ 7,
          sta=="Ophthalmology" ~ 8,
          sta=="Vaccines (Infectious Disease)" ~ 9,
          TRUE ~ 10),
        #mutate(intsta=group_indices(), # Create group incdices - may want to redo this to avoid issues with new groups
        flg=1) %>% #group_by(sta, sta) %>%summarize(n=n())
      ungroup() %>%
      dplyr::select(-strtherapeuticarea, -sta) %>%
      distinct() %>%
      group_by(row_id) %>%
      pivot_wider(names_from=intsta, names_prefix="ta", values_from = flg,
                  values_fill=list(flg=0)) %>%
      ungroup() %>%
      fncols(paste0("ta", 1:9))
    
    # Get UMAP embedding for each sponsor
    message("Sponsor embedding: start pre-prossing on basic sponsor features")
    prepro <- fulldat %>%
      left_join(tas, by="row_id") %>%
      ungroup() %>%
      rowwise() %>%
      mutate(nspons = length(str_split(strSponsor, pattern="\\|")[[1]]), 
             nsponst = length(str_split(strSponsorType, pattern="\\|")[[1]]),
             spons = list( str_split(strSponsor, pattern="\\|")[[1]][ c(1:nspons, rep(nspons, max(nspons,nsponst)-nspons)) ] ), 
             sponst = list( str_split(strSponsorType, pattern="\\|")[[1]][ c(1:nsponst, rep(nsponst, max(nspons,nsponst)-nsponst)) ] ),
             rando = 1*str_detect(strdesignkeyword, "random"),
             pk = 1*str_detect(strdesignkeyword, "pharmacokin"),
             pd = 1*str_detect(strdesignkeyword, "pharmacodyn"),
             drf = 1*str_detect(strdesignkeyword, "dose resp"),
             rubbish = (str_detect(strdesignkeyword, "open"))*1 + 1*(str_detect(strdesignkeyword, "single")) + 1*(str_detect(strdesignkeyword, "ascend")),
             sophis = 1*str_detect(strdesignkeyword, "random") + 1*str_detect(strdesignkeyword, "multiple arm") + 1*str_detect(strdesignkeyword, "umbrella") + 
               1*str_detect(strdesignkeyword, "basket") + 1*str_detect(strdesignkeyword, "efficacy") + 1*str_detect(strdesignkeyword, "super") + 1*str_detect(strdesignkeyword, "infer") +
               1*str_detect(strdesignkeyword, "immunogen") + 1*str_detect(strdesignkeyword, "dose resp") + 1*str_detect(strdesignkeyword, "active") - rubbish + 
               1*str_detect(strdesignkeyword, "double")  + 1*str_detect(strdesignkeyword, "placebo") 
      ) %>%
      ungroup() %>%
      unnest(cols=c("spons","sponst")) %>%
      group_by(spons, sponst) %>%
      summarize(log_nrec=log(n()+0.5),
                ta1=log(sum(ta1)+0.5) - log_nrec,
                ta2=log(sum(ta2)+0.5) - log_nrec,
                ta3=log(sum(ta3)+0.5) - log_nrec,
                ta4=log(sum(ta4)+0.5) - log_nrec,
                ta5=log(sum(ta5)+0.5) - log_nrec,
                ta6=log(sum(ta6)+0.5) - log_nrec,
                ta7=log(sum(ta7)+0.5) - log_nrec,
                ta8=log(sum(ta8)+0.5) - log_nrec,
                ta9=log(sum(ta9)+0.5) - log_nrec,
                nspons = log( length(nspons) + 0.5),
                log_totalpats = log( sum( ifelse(is.na(intactualaccrual), ifelse(is.na(inttargetaccrual), 0, inttargetaccrual), intactualaccrual),na.rm=T) + 0.5),
                log_totalpats_pt = log( sum( ifelse(is.na(intactualaccrual), ifelse(is.na(inttargetaccrual), 0, inttargetaccrual), intactualaccrual),na.rm=T) + 0.5) - log_nrec,
                log_totalna = log( sum(is.na(inttargetaccrual)) + sum(is.na(intactualaccrual)) + sum(is.na(intduration)) + sum(is.na(strdesignkeyword)) + sum(is.na(strTerminationReason)) + 0.5),
                log_totalna_pt = log( sum(is.na(inttargetaccrual)) + sum(is.na(intactualaccrual)) + sum(is.na(intduration)) + sum(is.na(strdesignkeyword)) + sum(is.na(strTerminationReason)) + 0.5) - log_nrec,
                lognarate = log_totalna-log_nrec,
                meanyear=mean(ifelse(intphaseendyear<1990, 2009, ifelse(intphaseendyear>2014, 2014, intphaseendyear) ))-2009,
                log_drugs = log( length(unique(DrugKey)) + 0.5),
                log_indications = log( length(unique(indicationkey)) + 0.5),
                log_indications_pt = log( length(unique(indicationkey)) + 0.5) - log_nrec,
                log_indications_pd = log_indications - log_drugs,
                log_stas = log( length(unique(strtherapeuticarea)) + 0.5),
                log_stas_pt =log_stas - log_nrec,
                log_stas_pt =log_stas - log_drugs) %>%
      ungroup() %>%
      # mutate(odd = ifelse(tolower(sponst) %in% c("academic","contract research organization","cooperative group","government",
      #                                            "industry, all other pharma","industry, generic","industry, top 20 pharma"),0,1),
      #        toppharma = ifelse(tolower(sponst) == "industry, top 20 pharma", 1, 0),
      #        nonindustry = ifelse( tolower(sponst) %in% c("academic","cooperative group","government","not for profit funding entity "), 1, 0),
      #        bucket = 1*str_detect(prepro$spons, "\\(Oth"),
      #        bucket2 = rnorm(n=n(), bucket*10, 0.5),
      #        sponstype = 1*(str_detect(tolower(spons),"industry") | 
      #                         str_detect(tolower(spons),"contract research")) + 
      #          2*(str_detect(tolower(spons),"academic") | str_detect(tolower(spons),"cooperative") | 
      #               str_detect(tolower(spons),"government")),
      #        spons_nobrack = trimws(str_remove_all(spons, "\\{.*\\}")),
    #        spons_inbrack = trimws(str_extract(paste0(spons,"{none}"), "\\{.*\\}"))) %>%
    # ungroup()
    mutate(odd = ifelse(tolower(sponst) %in% c("academic","contract research organization","cooperative group","government",
                                               "industry, all other pharma","industry, generic","industry, top 20 pharma"),0,1),
           toppharma = ifelse(tolower(sponst) == "industry, top 20 pharma", 1, 0),
           nonindustry = ifelse( tolower(sponst) %in% c("academic","cooperative group","government","not for profit funding entity "), 1, 0),
           bucket = 1*str_detect(spons, "\\(Oth"),
           bucket2 = rnorm(n=n(), bucket*10, 0.5),
           sponstype = 1*(str_detect(tolower(spons),"industry") | 
                            str_detect(tolower(spons),"contract research")) + 
             2*(str_detect(tolower(spons),"academic") | str_detect(tolower(spons),"cooperative") | 
                  str_detect(tolower(spons),"government")),
           clean_spons = trimws(str_remove_all(spons, "\\{.*\\}")),
           spons_nobrack = str_detect(spons, "\\{.*\\}"),
           spons_inbrack = str_remove_all( ifelse(spons_nobrack, trimws(str_extract(spons, "\\{.*\\}")), "none"), "[\\{\\}]"))
    
    message("Sponsor embedding: get country information")
    tmp <- fulldat %>% 
      dplyr::select(row_id, strLocation) %>%
      mutate(location = str_split(strLocation, pattern="\\|")) %>%
      unnest(location) %>%
      mutate(location = str_remove_all(location, "[\\[\\]]"),
             country_us = 1*str_detect(tolower(location),"united states"),
             country_ge = 1*str_detect(tolower(location),"germany"),
             country_uk = 1*str_detect(tolower(location),"united king"),
             country_ca = 1*str_detect(tolower(location),"canada"),
             country_fr = 1*str_detect(tolower(location),"france"),
             country_it = 1*str_detect(tolower(location),"italy"),
             country_be = 1*str_detect(tolower(location),"belgium"),
             country_sp = 1*str_detect(tolower(location),"spain"),
             country_nl = 1*str_detect(tolower(location),"netherland"),
             country_pl = 1*str_detect(tolower(location),"poland"),
             country_au = 1*str_detect(tolower(location),"australia"),
             country_ru = 1*str_detect(tolower(location),"russia"),
             country_cz = 1*str_detect(tolower(location),"czech"),
             country_se = 1*str_detect(tolower(location),"sweden"),
             country_hu = 1*str_detect(tolower(location),"hungary"),
             country_dk = 1*str_detect(tolower(location),"denmark"),
             country_aa = 1*str_detect(tolower(location),"austria"),
             country_ch = 1*str_detect(tolower(location),"switzer"),
             
             country_jp = 1*str_detect(tolower(location),"japan"),
             country_ii = 1*str_detect(tolower(location),"israel"),
             country_ro = 1*str_detect(tolower(location),"romania"),
             country_sa = 1*str_detect(tolower(location),"south africa"),
             country_in = 1*str_detect(tolower(location),"india"),
             country_mx = 1*str_detect(tolower(location),"mexico"),
             country_rk = 1*str_detect(tolower(location),"south korea"),
             country_ar = 1*str_detect(tolower(location),"argentin"),
             country_fi = 1*str_detect(tolower(location),"finland"),
             country_bu = 1*str_detect(tolower(location),"bulga"),
             country_ur = 1*str_detect(tolower(location),"ukrai"),
             country_sl = 1*str_detect(tolower(location),"slovak"),
             country_no = 1*str_detect(tolower(location),"norway"),
             country_nz = 1*str_detect(tolower(location),"zealand"),
             country_br = 1*str_detect(tolower(location),"brazil"),
             country_pr = 1*str_detect(tolower(location),"puerto rico"),
             country_tw = 1*str_detect(tolower(location),"taiwan"),
             country_pe = 1*str_detect(tolower(location),"peru"),
             country_cl = 1*str_detect(tolower(location),"chile"),
             country_sp = 1*str_detect(tolower(location),"singap"),
             country_cp = 1*str_detect(tolower(location),"china"),
             country_hk = 1*str_detect(tolower(location),"kong"),
             region_yu = 1*(str_detect(tolower(location),"serbia") | str_detect(tolower(location),"yugosl") | str_detect(tolower(location),"macedon") | 
                              str_detect(tolower(location),"bosnia") | str_detect(tolower(location),"croatia") + str_detect(tolower(location),"montenegro") > 0),
             
             region_eu = 1*(country_ge + country_uk + country_fr + country_be + country_sp + country_nl + country_pl + country_cz + country_se + country_hu + country_dk +
                              country_aa + country_ch + country_ro + country_fi + country_bu + country_ur + country_sl + country_no + region_yu +
                              str_detect(tolower(location),"malta") + str_detect(tolower(location),"luxem") + str_detect(tolower(location),"moldov") + str_detect(tolower(location),"iceland") + 
                              + str_detect(tolower(location),"europe") + str_detect(tolower(location),"portugal") + str_detect(tolower(location),"latvia") + 
                              str_detect(tolower(location),"ireland") + str_detect(tolower(location),"lithuan") + str_detect(tolower(location),"estonia") + str_detect(tolower(location),"greece") + 
                              str_detect(tolower(location),"norway") + str_detect(tolower(location),"malta") >0), 
             region_ee = 1*(country_ru + country_pl + country_cz  + country_hu + country_ro + country_bu + country_ur + country_sl + region_yu + str_detect(tolower(location),"moldov") +
                              str_detect(tolower(location),"latvia") + str_detect(tolower(location),"lithuan") + str_detect(tolower(location),"estonia") + 
                              str_detect(tolower(location),"belarus") + str_detect(tolower(location),"eastern europe") + str_detect(tolower(location),"armenia") +str_detect(tolower(location),"georgia") >0),
             region_we = 1*(country_ge + country_uk + country_fr + country_be + country_sp + country_nl + 
                              country_aa + country_ch + str_detect(tolower(location),"malta") + str_detect(tolower(location),"luxem") + 
                              str_detect(tolower(location),"portugal") + str_detect(tolower(location),"ireland") >0),
             region_ne = 1*(country_se + country_dk + country_fi +  country_no + str_detect(tolower(location),"iceland") >0),
             region_usca = 1*(country_us +  country_ca>0),
             region_af_nsa = 1*(str_detect(tolower(location),"egypt") + str_detect(tolower(location),"africa") & !str_detect(tolower(location),"south africa") +  str_detect(tolower(location),"morocc") + 
                                  str_detect(tolower(location),"tunisia") + str_detect(tolower(location),"ghana") + str_detect(tolower(location),"gambia") + str_detect(tolower(location),"mali") + 
                                  str_detect(tolower(location),"tanza") + str_detect(tolower(location),"zambia") + str_detect(tolower(location),"kenya") >0),
             region_un = 1*(region_af_nsa + str_detect(tolower(location),"arabia") + str_detect(tolower(location),"iran") + str_detect(tolower(location),"middle east") + 
                              str_detect(tolower(location),"arab") + str_detect(tolower(location),"egypt") + str_detect(tolower(location),"cuba") + str_detect(tolower(location),"north korea") >0),
             region_as = 1*(str_detect(tolower(location),"asia") + str_detect(tolower(location),"vietna") + str_detect(tolower(location),"indones") +
                              str_detect(tolower(location),"indones") + str_detect(tolower(location),"thail") + str_detect(tolower(location),"malaysia") + str_detect(tolower(location),"philip") + 
                              country_hk + country_cp + country_rk + country_sp + country_tw >0),
             region_is = 1*(str_detect(tolower(location),"sri lanka") + str_detect(tolower(location),"india") + str_detect(tolower(location),"pakista") + str_detect(tolower(location),"bangla") +
                              str_detect(tolower(location),"afghan") + str_detect(tolower(location),"nepal")>0), 
             region_ans = 1*(country_nz + country_au + country_sa>0)
      ) %>% 
      dplyr::select(-strLocation, -location) %>%
      group_by(row_id) %>%
      summarize_all(sum)
    
    message("Sponsor embedding: Merge in country")
    countries <- fulldat %>%
      dplyr::select(row_id, strSponsor) %>%
      mutate(spons = str_split(strSponsor, pattern="\\|")) %>%
      unnest(cols="spons") %>%
      left_join(tmp, by="row_id") %>%
      ungroup() %>%
      dplyr::select(-row_id, -strSponsor) %>%
      group_by(spons) %>%
      summarize_all(mean) %>% 
      ungroup()
    
    message("Sponsor embedding: merge country back in together")
    prepro <- prepro %>% 
      left_join(countries, by="spons") %>%
      ungroup()
    
    # message("Summarize missing before UMAP")
    # prepro %>%
    #   summarise_all(funs(sum(is.na(.)))) %>%
    #   gather(key="variable", value="value") %>%
    #   filter(value>0) %>%
    #   print(n=200)
    
    message("Sponsor embedding pre-processing done - start UMAP")
    #spons_umap <- umap(, n_neighbors = min(dim(prepro)[1], 500), learning_rate = 0.5, init = "random")
    spons_umap <- umap(dplyr::select(prepro, -log_totalna), 
                       n_neighbors = dim(prepro)[1], 
                       n_epochs = 2500,learning_rate = 0.5, init = "random")
    
    message("Sponsor embedding pre-processing done - UMAP done")
    prepro2 <- prepro %>%
      mutate(spons_umap1=spons_umap[,1],
             spons_umap2=spons_umap[,2]) 
    
    message("Sponsor embedding: Finished, generate return tibble")
    rval <- fulldat %>%
      dplyr::select(DrugKey, indicationkey, row_id, strSponsor) %>%
      mutate(spons = str_split(strSponsor, pattern="\\|")) %>%
      unnest(cols="spons") %>%
      dplyr::select(DrugKey, indicationkey, row_id,spons) %>%
      left_join(dplyr::select(prepro2, spons, sponst, spons_umap1, spons_umap2), 
                by="spons") 
    
    write_rds(x=rval, path=filepath)
  }
  
  #%>%
  # group_by(DrugKey, indicationkey) %>%
  # summarize(n=n()) %>%
  # ungroup() %>%
  # summarize(min=min(n), max=max(n), mean=mean(n), q3=quantile(n,0.75) , p95=quantile(n,0.95))
  
  # group_by(DrugKey, indicationkey) %>%
  # mutate(num_of_trials = length(unique(row_id))) %>%
  # ungroup() %>%
  # group_by(DrugKey, indicationkey, spons, sponstype) %>%
  # summarize(num_of_trials = mean(num_of_trials),
  #   trial_prop = length(unique(row_id))/num_of_trials) %>%
  # ungroup() %>%
  
  #group_by(DrugKey,indicationkey) %>%
  #   group_by(DrugKey, indicationkey) %>%
  #   arrange(desc(num_of_trials), DrugKey,indicationkey, desc(trial_prop)) %>%
  #   mutate(rowno=row_number()) %>%
  #   filter(rowno==1) %>% print(n=100)
  # # %>%
  #   mutate(rowno=row_number(),
  #          xn=n()) %>%
  #   ungroup()
  # # %>%
  # group_by(row_id) %>%
  # summarize(n=n()) %>%
  # ungroup() %>%
  # group_by(n) %>%
  # summarize(hm=n())
  
  return( rval )
}

################################################
# Function to get the features derived by Bjoern
################################################

get_features_bjoern <- function(file1="/shared_data/data/training_data/training_data_2015_split_on_outcome.csv", file2=NA){
  
  # Read in the files, both need to be either character strings or tibbles/dataframes (mixing the two not foreseen)
  if (!is.na(file2)){
    if (is.character(file1)){
      fulldat <- read_csv(file1) %>%
        bind_rows(read_csv(file2))
    } else {
      fulldat <- as_tibble(file1) %>% 
        bind_rows(as_tibble(file2))
    }
  } else {
    if (is.character(file1)){
      fulldat <- read_csv(file1) 
    } else {
      fulldat <- as_tibble(file1)
    }
  }
  
  if (!("outcome" %in% colnames(fulldat))) {
    message("'outcome' is not present adding missing")
    fulldat <- fulldat %>%
      mutate(outcome=NA_integer_)
  }
  
  # Get sponsor information
  sponsor <- fulldat %>%
    mutate(spons = str_split(strSponsorType, pattern="\\|")) %>%
    unnest(cols="spons") %>%
    mutate(spons_partners = str_count(strSponsor,"\\|"),
           spons_curlybr = str_count(strSponsor,"[\\{\\}]"),
           sponstype = 1*(str_detect(tolower(spons),"industry") | 
                            str_detect(tolower(spons),"contract research")) + 
             2*(str_detect(tolower(spons),"academic") | str_detect(tolower(spons),"cooperative") | 
                  str_detect(tolower(spons),"government")), # General sponsor type
           # One-hot-encoded (dummy variables) for specific sponsor types
           spons1=1*(str_detect(tolower(spons),"academic")),
           spons2=1*(str_detect(tolower(spons),"cooperativ")),
           spons3=1*(str_detect(tolower(spons),"government")),
           spons4=1*(str_detect(tolower(spons),"top 20")),
           spons5=1*(str_detect(tolower(spons),"generic")),
           spons6=1*(str_detect(tolower(spons),"diagnostics") | str_detect(tolower(spons),"device")),
           spons7=1*(str_detect(tolower(spons),"other pharma")),
           spons8=1*(str_detect(tolower(spons),"not for profi")),
           spons9=1*(str_detect(tolower(spons),"miscel")),
    ) %>%
    group_by(DrugKey, indicationkey, row_id) %>% # Take maximum for each trial 
    summarize(sponstype=max(sponstype), # highest "quality" sponsor for each trial
              spons_partners = mean(spons_partners),
              spons_curlybr = mean(spons_curlybr),
              spons1=max(spons1),
              spons2=max(spons2),
              spons3=max(spons3),
              spons4=max(spons4),
              spons5=max(spons5),
              spons6=max(spons6),
              spons7=max(spons7),
              spons8=max(spons8)) %>%
    ungroup() %>%
    group_by(DrugKey, indicationkey) %>% # Take maximum for the drug
    summarize(sponstype=max(sponstype),
              spons_partners = mean(spons_partners),
              spons_curlybr = mean(spons_curlybr),
              spons_partners_sd = ifelse(is.na(sd(spons_partners)), 0, sd(spons_partners)),
              spons_curlybr_sd = ifelse(is.na(sd(spons_curlybr)), 0, sd(spons_curlybr)),
              spons1=max(spons1),
              spons2=max(spons2),
              spons3=max(spons3),
              spons4=max(spons4),
              spons5=max(spons5),
              spons6=max(spons6),
              spons7=max(spons7),
              spons8=max(spons8)) 
  
  indiumap <- indumap(fulldat)
  
  ris <- umap_sponsor(fulldat)
  
  ris2 <- ris %>% 
    group_by(DrugKey, indicationkey, spons, sponst, spons_umap1, spons_umap2 ) %>%
    summarize(count=length(unique(row_id))) %>%
    group_by(DrugKey, indicationkey) %>%
    mutate(srecords=n(),
           industry=str_detect(tolower(sponst), "industry")) %>%
    arrange(DrugKey, indicationkey,desc(count)) %>%
    mutate(srank=row_number()) %>%
    group_by(DrugKey, indicationkey, industry) %>%
    arrange(DrugKey, indicationkey, industry, desc(count)) %>%
    mutate(srank2=row_number(),
           flag1=(srank==1),
           flag2=(srank2==1)) %>%
    ungroup()
  # arrange(DrugKey, indicationkey,desc(count)) %>%
  # filter(srecords>5) %>% print(n=30)
  
  ris2oth <- ris2 %>%
    filter(!flag1 & !flag2) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize(avgoth_umap1=mean(spons_umap1),
              avgoth_umap2=mean(spons_umap2))
  
  ris3 <- dplyr::select(filter(ris2, flag1), DrugKey, indicationkey, spons_umap1, spons_umap2) %>%
    left_join( rename( dplyr::select(filter(ris2, flag2 & industry==T), DrugKey, indicationkey, spons_umap1, spons_umap2), 
                       spons_ind_umap1=spons_umap1, spons_ind_umap2=spons_umap2),
               by=c("DrugKey","indicationkey")) %>%
    left_join( rename( dplyr::select(filter(ris2, flag2 & industry==F), DrugKey, indicationkey, spons_umap1, spons_umap2),
                       spons_nind_umap1=spons_umap1, spons_nind_umap2=spons_umap2),
               by=c("DrugKey","indicationkey")) %>%
    left_join(ris2oth, by=c("DrugKey","indicationkey") ) %>%
    ungroup()
  
  
  # Information for trial completion
  trialcompl <- fulldat %>%
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
    group_by(DrugKey, indicationkey, row_id) %>% # Summarize on a trial level
    summarize(pct_accrual=min(pct_accrual),
              bad=max(bad),
              safety=max(safety)) %>%
    mutate(accrued90 = (pct_accrual>=0.90)*1,
           accrued50 = (pct_accrual>=0.5)*1) %>%
    ungroup() %>%
    group_by(DrugKey, indicationkey) %>% # summarize on a drug/indication level
    summarize(pct_accrual=gm_mean(pct_accrual),
              accrued50=mean(accrued50),
              accrued90=mean(accrued90),
              mean_trialendscore=mean(bad),
              max_trialendscore=max(bad),
              trialendscore0=mean( (bad==0)),
              trialendscore1=mean( (bad==1)),
              trialendscore2=mean( (bad==2)),
              trialendscore3=mean( (bad==3)),
              trialendscore4=mean( (bad==4)),
              safety=max(safety)) %>%
    ungroup()
  
  # Look for safety stopping of any trial of the drug
  drug_compl <- trialcompl %>% 
    group_by(DrugKey) %>%
    summarize(drug_safety=mean(safety),
              drug_max_trialendscore=max(max_trialendscore)) %>%
    ungroup()
  trialcompl <- trialcompl %>%
    left_join(drug_compl, by="DrugKey") %>%
    ungroup()
  
  # Some further simple features (on disease area): 
  simpdat1 <- fulldat %>%
    mutate(patno = ifelse(is.na(intactualaccrual), ifelse(is.na(inttargetaccrual), # Number of patients
                                                          0, inttargetaccrual), intactualaccrual)) %>%
    group_by(DrugKey, indicationkey) %>% #, intphaseendyear, strtherapeuticarea) %>%
    mutate(number_of_trials=n(), # Number of trials
           total_pats=sum(patno)) %>% # Total number of patients
    mutate(total_pats=ifelse(total_pats==0, 1, total_pats)) %>% 
    ungroup() %>%
    dplyr::select(DrugKey, indicationkey, intphaseendyear, strtherapeuticarea, number_of_trials, total_pats) %>%
    mutate(sta = str_split(strtherapeuticarea, pattern="\\|"), # Therapeutic areas
           years_since_1999=ifelse(intphaseendyear>2000, # Years since 2000 (below 2000 also = 1)
                                   ifelse(intphaseendyear>2018, 2018-1999, # Truncate to 2018 (should not have data after)
                                          intphaseendyear-1999), 1)) %>% 
    unnest(cols="sta") %>%
    group_by(sta) %>%
    mutate(
      intsta = case_when(
        sta=="Autoimmune/Inflammation" ~ 1,
        sta=="Cardiovascular" ~ 2,
        sta=="CNS" ~ 3,
        sta=="Genitourinary" ~ 4,
        sta=="Infectious Disease" ~ 5,
        sta=="Metabolic/Endocrinology" ~ 6,
        sta=="Oncology" ~ 7,
        sta=="Ophthalmology" ~ 8,
        sta=="Vaccines (Infectious Disease)" ~ 9,
        TRUE ~ 10),
      #mutate(intsta=group_indices(), # Create group incdices - may want to redo this to avoid issues with new groups
      flg=1) %>% #group_by(sta, sta) %>%summarize(n=n())
    ungroup() %>%
    dplyr::select(-strtherapeuticarea, -sta, -intphaseendyear) %>%
    distinct() %>%
    group_by(DrugKey, indicationkey, years_since_1999, number_of_trials, total_pats) %>%
    pivot_wider(names_from=intsta, names_prefix="ta", values_from = flg,
                values_fill=list(flg=0)) %>%
    ungroup() %>%
    fncols(paste0("ta", 1:9))
  
  if ( length(unique(fulldat$intpriorapproval)) == 2 ) {
    if ( min(sort(unique(fulldat$intpriorapproval))==c("['0']","[]"))!=1 ) {
      warning("Help! They've changed the data on intpriorapproval - populated values don't match")
    }
  } else {
    warning("Help! They've changed the data on intpriorapproval - more than two possible values now")
  }
  
  # Determine prior approval status overall y/n, by number and TA
  print("lcm1")
  lcm1 <- fulldat %>% 
    mutate(prior_approval = 1*(intpriorapproval=='[]') + 0*(intpriorapproval=='[0]')) %>% # Prior approvals?
    group_by(DrugKey, indicationkey, outcome, intphaseendyear, strtherapeuticarea) %>% 
    summarize(prior_approval=max(prior_approval)) %>%
    ungroup() %>% 
    dplyr::select(DrugKey, indicationkey, outcome, intphaseendyear, prior_approval, strtherapeuticarea) %>%
    distinct()
  
  
  # Number of prior approvals
  print("lcm2")
  lcm2 <- lcm1 %>% 
    dplyr::select(DrugKey, intphaseendyear, prior_approval, outcome) %>%
    distinct() %>% 
    group_by(DrugKey, intphaseendyear) %>%
    summarize(approvals_in_year=sum(outcome, na.rm=T),
              prior_approval=min(prior_approval)) %>% 
    ungroup() %>%
    group_by(DrugKey) %>%
    arrange(DrugKey, intphaseendyear) %>%
    mutate(#totalapprovals=mean(approvals_in_year, na.rm=T),
      fpa = first(prior_approval),
      dcum= replace_na( lag(cumsum(approvals_in_year), 1), 0),
      past_approvals = pmax( 1*(cumsum(prior_approval)>0), fpa+dcum)) %>%
    dplyr::select(DrugKey, intphaseendyear, past_approvals) %>%
    distinct()
  
  # Number of prior approvals within TA
  print("lcm3")
  lcm3 <- lcm1 %>% 
    group_by(DrugKey) %>%
    mutate(number_of_tas=length(unique(strtherapeuticarea))) %>%
    ungroup() %>% #dplyr::select(DrugKey, strtherapeuticarea, intphaseendyear) %>% distinct() %>% dim()
    group_by(DrugKey, strtherapeuticarea, intphaseendyear, number_of_tas) %>%
    summarize(approvals_in_year=sum(outcome, na.rm=T),
              prior_approval = min(prior_approval)) %>%
    ungroup() %>%
    group_by(DrugKey, strtherapeuticarea) %>%
    arrange(DrugKey, strtherapeuticarea, intphaseendyear) %>%
    mutate(#totalapprovals=mean(approvals_in_year, na.rm=T), 
      fpa = ifelse(number_of_tas==1, first(prior_approval),0),
      dcum= replace_na( lag(cumsum(approvals_in_year), 1), 0),
      past_approvals_ta = pmax( 1*((number_of_tas==1)*cumsum(prior_approval)>0), fpa+dcum)) %>%
    ungroup() %>% #dplyr::select(DrugKey, strtherapeuticarea, intphaseendyear) %>% 
    #arrange(desc(totalapprovals), DrugKey, intphaseendyear) %>%
    dplyr::select(DrugKey, strtherapeuticarea, intphaseendyear, past_approvals_ta)
  
  print("lcm-merge")
  lcm <- distinct(dplyr::select(lcm1, DrugKey, indicationkey, strtherapeuticarea, intphaseendyear, prior_approval)) %>% 
    left_join(lcm2, by=c("DrugKey", "intphaseendyear")) %>% # distinct() %>% dim() dplyr::select(DrugKey, indicationkey, intphaseendyear) %>% distinct()
    left_join(lcm3, by=c("DrugKey", "strtherapeuticarea", "intphaseendyear")) %>%
    dplyr::select(DrugKey, indicationkey, prior_approval, past_approvals, past_approvals_ta) %>%
    group_by(DrugKey, indicationkey, prior_approval, past_approvals) %>%
    summarize(past_approvals_ta=max(past_approvals_ta)) %>%
    distinct()
  
  print("Based on GenericName")
  # Things based on the generic name and life-cycle/prior approval:
  wiian <- fulldat %>% # INN costs 9000 dollars pre 2014, 12000 afterwards, see: https://extranet.who.int/tools/inn_online_application/
    mutate(unwilling_to_pay_12k = 1*( ( str_starts(GenericName,"[A-Z]") & !str_starts(GenericName,"(HPV|MMR|DTP|HIV)")) | 
                                        str_detect(GenericName, "[A-Z][A-Z][A-Z][A-Z][A-Z]") | str_detect(GenericName, "[0-9][0-9]") ),
           is_a_combination = 1*(str_detect(GenericName, "\\+") | (str_detect(GenericName, "[a-z]/[a-z]") & # Identify combination drugs
                                                                     !str_detect(GenericName, "alpha/beta") &
                                                                     !str_detect(GenericName, "beta/gamma") & 
                                                                     !str_detect(GenericName, "gamma/delta"))),
           is_mab = 1*(str_detect(GenericName, "mab$") | str_detect(GenericName, "mab ")  | # Is the drug a monoclonal antibody (ending in -mab or muromonab, which is the one exception to this rule)
                         str_detect(GenericName, "mab,") | str_detect(GenericName, "mab-") | str_detect(GenericName, "muromonab"))) %>%
    dplyr::select(DrugKey, indicationkey, unwilling_to_pay_12k, is_a_combination, is_mab) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize(unwilling_to_pay_12k=max(unwilling_to_pay_12k),
              is_a_combination=max(is_a_combination),
              is_mab=max(is_mab)) # ,prior_approval=max(prior_approval)
  
  print("Disease type")
  # Using Disease Type variable
  dis1 <- fulldat %>%
    dplyr::select(DrugKey, indicationkey, strDiseaseType) %>% #strtherapeuticarea,  outcome
    mutate(strDiseaseType=str_remove_all(str_remove_all(str_remove_all(str_remove_all(strDiseaseType, "\\(N/A\\)\\|"),"\\|\\(N/A\\)"),"\\(\\)\\|"),"\\|\\(\\)")) %>%
    mutate(distyp = str_split(strDiseaseType, pattern="\\|")) %>%
    unnest(cols="distyp") %>%
    distinct() %>%
    mutate(distyp=str_remove_all(distyp, "[\\[\\]]"),
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
             TRUE ~ 11 ),
           keytyp = 1*(disno %in% c(111,48,63,78,24,26,133,99,53,83,5,29,25,10,23,71,61,49,92,3,16,20,19,74,2,1,15,35,41,54,76,
                                    81,86,46,38,37,31,17,40,66,14,14,18,68,13,8,9,21,6,12, 62, 64, 30, 34, 22, 12, 62, 6, 56, 97, 15,
                                    113, 59, 72,101,60,55,91,105,52,7)),
           vkeytyp = 1*(disno %in% c(30,2,62,17,9,7,31,16,29,5,21,38,133,26,23,38,63,19,18,24,15,10,25,12,6,81,20,34,13,14,99,46,8,53))
    ) %>%
    dplyr::select(DrugKey, indicationkey, disno, new_disno, keytyp, vkeytyp)
  
  print("Flags for diseases")
  dis2 <- dis1 %>%
    dplyr::select(DrugKey, indicationkey, disno) %>%
    mutate(anytyp=1) %>%
    distinct() %>%
    group_by(DrugKey, indicationkey) %>%
    pivot_wider(names_from = disno, names_prefix = "disease", values_from = anytyp, id_cols = c("DrugKey", "indicationkey")) %>%  
    replace(., is.na(.), 0)
  
  dis2a <- dis1 %>%
    dplyr::select(DrugKey, indicationkey, new_disno) %>%
    mutate(anytyp=1) %>%
    distinct() %>%
    group_by(DrugKey, indicationkey) %>%
    pivot_wider(names_from = new_disno, names_prefix = "newta", values_from = anytyp, id_cols = c("DrugKey", "indicationkey")) %>%  
    replace(., is.na(.), 0)
  
  dis3 <- dis1 %>%
    filter(keytyp==1) %>%
    dplyr::select(DrugKey, indicationkey, disno) %>%
    mutate(anytyp=1) %>%
    distinct() %>%
    group_by(DrugKey, indicationkey) %>%
    pivot_wider(names_from = disno, names_prefix = "kdisease", values_from = anytyp, id_cols = c("DrugKey", "indicationkey")) %>%  
    replace(., is.na(.), 0)
  
  dis4 <- dis1 %>%
    filter(vkeytyp==1) %>%
    dplyr::select(DrugKey, indicationkey, disno) %>%
    mutate(anytyp=1) %>%
    distinct() %>%
    group_by(DrugKey, indicationkey) %>%
    pivot_wider(names_from = disno, names_prefix = "vkdisease", values_from = anytyp, id_cols = c("DrugKey", "indicationkey")) %>%  
    replace(., is.na(.), 0)
  
  dises <- dis2 %>%
    full_join(dis2a, by=c("DrugKey", "indicationkey")) %>%
    full_join(dis3, by=c("DrugKey", "indicationkey")) %>%
    full_join(dis4, by=c("DrugKey", "indicationkey")) %>%
    replace(., is.na(.), 0)
  
  # Features for class effects
  drugclass1 <- fulldat %>%
    dplyr::select(DrugKey, strMechanismOfAction, outcome, intpriorapproval, indicationkey) %>%
    mutate(moa = str_split(strMechanismOfAction, pattern="\\|"),
           prior_approval = 1*(intpriorapproval=='[]') + 0*(intpriorapproval=='[0]')) %>%
    unnest(cols="moa") %>%
    mutate(moa=str_remove_all(moa, "[\\[\\]]")) %>%
    group_by(DrugKey) %>%
    mutate(ever_approved=max(outcome, prior_approval, na.rm=T),
           indications=length(unique(indicationkey))) %>%
    ungroup() %>%
    dplyr::select(DrugKey, moa, ever_approved, indications) %>%
    distinct() 
  
  drugclass1a <- drugclass1 %>%
    dplyr::select(-ever_approved,-indications) %>%
    left_join(drugclass1, by="moa") %>%
    filter(DrugKey.x != DrugKey.y) %>%
    group_by(DrugKey.x, moa) %>%
    summarize(class_any_approval_yn = 1*(sum(ever_approved)>0),
              class_approval_rate=qbeta(0.5,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_approval_rate_lcl80=qbeta(0.1,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_approval_rate_ucl80=qbeta(0.9,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_indications=mean(indications),
              class_approval_per_ind=class_approval_rate/class_indications,
              class_approval_lcl_per_ind=class_approval_rate_lcl80/class_indications,
              class_approval_ucl_per_ind=class_approval_rate_ucl80/class_indications) %>%
    ungroup() %>%
    rename(DrugKey=DrugKey.x)
  
  
  drugclass2 <- fulldat %>%
    dplyr::select(DrugKey, strMechanismOfAction, outcome, indicationkey) %>%
    mutate(moa = str_split(strMechanismOfAction, pattern="\\|")) %>%
    unnest(cols="moa") %>%
    mutate(moa=str_remove_all(moa, "[\\[\\]]")) %>%
    group_by(moa, DrugKey, indicationkey) %>%
    mutate(ever_approved=max(0, outcome, na.rm=T)) %>%
    ungroup() %>%
    dplyr::select(DrugKey, indicationkey, moa, ever_approved) %>%
    distinct() 
  
  drugclass2a <- drugclass2 %>%
    dplyr::select(-ever_approved) %>%
    left_join(drugclass2, by=c("moa","indicationkey")) %>%
    filter(DrugKey.x != DrugKey.y) %>%
    group_by(DrugKey.x, indicationkey, moa) %>%
    summarize(class_ind_approval_yn = 1*(sum(ever_approved)>0),
              class_ind_approval_rate=qbeta(0.5,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_ind_approval_rate_lcl80=qbeta(0.1,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_ind_approval_rate_ucl80=qbeta(0.9,0.12+sum(ever_approved), 0.88+sum(1-ever_approved))) %>%
    ungroup() %>%
    rename(DrugKey=DrugKey.x)
  
  drugclass3 <- fulldat %>%
    dplyr::select(DrugKey, strMechanismOfAction, outcome, strDiseaseType, indicationkey) %>%
    mutate(moa = str_split(strMechanismOfAction, pattern="\\|")) %>%
    unnest(cols="moa") %>%
    mutate(moa=str_remove_all(moa, "[\\[\\]]"),
           strDiseaseType=str_remove_all(str_remove_all(str_remove_all(str_remove_all(strDiseaseType, "\\(N/A\\)\\|"),"\\|\\(N/A\\)"),"\\(\\)\\|"),"\\|\\(\\)")) %>%
    mutate(distyp = str_split(strDiseaseType, pattern="\\|")) %>%
    unnest(cols="distyp") %>%
    mutate(distyp=str_remove_all(distyp, "[\\[\\]]")) %>%
    dplyr::select(-strMechanismOfAction, -strDiseaseType) %>%
    distinct() %>%
    group_by(moa, DrugKey, distyp) %>%
    mutate(ever_approved=max(0, outcome, na.rm=T)) %>%
    ungroup() %>%
    dplyr::select(DrugKey, indicationkey, distyp, moa, ever_approved) %>%
    distinct() 
  
  drugclass3a <- drugclass3 %>%
    dplyr::select(-ever_approved) %>%
    left_join(dplyr::select(drugclass3,-indicationkey), by=c("moa","distyp")) %>%
    filter(DrugKey.x != DrugKey.y) %>%
    group_by(DrugKey.x, indicationkey, moa) %>%
    summarize(class_dt_approval_yn = 1*(sum(ever_approved)>0),
              class_dt_approval_rate=qbeta(0.5,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_dt_approval_rate_lcl80=qbeta(0.1,0.12+sum(ever_approved), 0.88+sum(1-ever_approved)),
              class_dt_approval_rate_ucl80=qbeta(0.9,0.12+sum(ever_approved), 0.88+sum(1-ever_approved))) %>%
    ungroup() %>%
    rename(DrugKey=DrugKey.x)
  
  drugclass <- fulldat %>%
    dplyr::select(DrugKey, indicationkey) %>%
    distinct() %>%
    left_join(drugclass3a, by=c("DrugKey", "indicationkey")) %>%
    left_join(drugclass2a, by=c("DrugKey", "indicationkey", "moa")) %>%
    left_join(drugclass1a, by=c("DrugKey", "moa")) %>%
    ungroup() %>%
    dplyr::select(-moa) %>%
    mutate(class_dt_approval_yn=ifelse(is.na(class_dt_approval_yn),0 , class_dt_approval_yn),
           class_dt_approval_rate = ifelse(is.na(class_dt_approval_rate), qbeta(0.5,0.12, 0.88) , class_dt_approval_rate),
           class_dt_approval_rate_lcl80 = ifelse(is.na(class_dt_approval_rate_lcl80), qbeta(0.5,0.12, 0.88) , class_dt_approval_rate_lcl80),
           class_dt_approval_rate_ucl80 = ifelse(is.na(class_dt_approval_rate_ucl80), qbeta(0.5,0.12, 0.88) , class_dt_approval_rate_ucl80),
           class_ind_approval_yn = ifelse(is.na(class_ind_approval_yn),0 , class_ind_approval_yn),
           class_ind_approval_rate= ifelse(is.na(class_ind_approval_rate), qbeta(0.5,0.12, 0.88) , class_ind_approval_rate),
           class_ind_approval_rate_lcl80= ifelse(is.na(class_ind_approval_rate_lcl80), qbeta(0.1,0.12, 0.88) , class_ind_approval_rate_lcl80),
           class_ind_approval_rate_ucl80= ifelse(is.na(class_ind_approval_rate_ucl80), qbeta(0.9,0.12, 0.88) , class_ind_approval_rate_ucl80),
           class_approval_rate = ifelse(is.na(class_approval_rate), qbeta(0.5,0.12, 0.88), class_approval_rate),
           class_approval_rate_lcl80= ifelse(is.na(class_approval_rate_lcl80), qbeta(0.1,0.12, 0.88), class_approval_rate_lcl80),
           class_approval_rate_ucl80= ifelse(is.na(class_approval_rate_ucl80), qbeta(0.9,0.12, 0.88), class_approval_rate_ucl80),
           class_indications=ifelse(is.na(class_indications),0 , class_indications),         
           class_any_approval_yn=ifelse(is.na(class_any_approval_yn),0 , class_any_approval_yn),
           class_approval_per_ind=ifelse(is.na(class_approval_per_ind),0 , class_approval_per_ind),                   
           class_approval_lcl_per_ind=ifelse(is.na(class_approval_lcl_per_ind),0 , class_approval_lcl_per_ind),                   
           class_approval_ucl_per_ind=ifelse(is.na(class_approval_lcl_per_ind),0 , class_approval_lcl_per_ind)
    ) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize_all(function(x) mean(x, na.rm=T)) %>%
    ungroup()
  
  missunkmoa <- fulldat %>%
    dplyr::select(DrugKey, indicationkey, strMechanismOfAction) %>%
    mutate(moa = str_split(strMechanismOfAction, pattern="\\|")) %>%
    unnest(cols="moa") %>%
    mutate(moa=str_remove_all(moa, "[\\[\\]]")) %>%
    dplyr::select(DrugKey, indicationkey, moa) %>%
    distinct() %>%
    mutate( miss_unk_moa = 1*(is.na(moa) | str_detect(tolower(moa),"undisclosed")  | str_detect(tolower(moa),"undefined")  | str_detect(tolower(moa),"unknown"))) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize(miss_unk_moa=mean(miss_unk_moa)) %>%
    ungroup()
  
  outcomedat <- fulldat %>%
    dplyr::select(DrugKey, indicationkey, outcome) %>%
    distinct()
  
  message("Merging all data")
  adat <- simpdat1 %>%
    full_join(outcomedat, by=c("DrugKey", "indicationkey")) %>%
    full_join(trialcompl, by=c("DrugKey", "indicationkey")) %>%
    full_join(sponsor, by=c("DrugKey", "indicationkey")) %>%
    full_join(wiian, by=c("DrugKey", "indicationkey")) %>%
    full_join(lcm, by=c("DrugKey", "indicationkey")) %>%
    full_join(ris3, by=c("DrugKey", "indicationkey")) %>%
    full_join(dises, by=c("DrugKey", "indicationkey")) %>%
    left_join(missunkmoa, by=c("DrugKey", "indicationkey")) %>%
    full_join(drugclass, by=c("DrugKey", "indicationkey")) %>%
    ungroup() %>%
    left_join(indiumap, by=c("indicationkey")) %>%
    ungroup() %>%
    mutate(dataset=ifelse(is.na(outcome), "test", "training"))
  # %>% fncols(c("vkdisease2", "vkdisease5", "vkdisease6", "vkdisease7", "vkdisease8", "vkdisease9", "vkdisease10", "vkdisease12", 
  #            "vkdisease13", "vkdisease14", "vkdisease15", "vkdisease16", "vkdisease17", "vkdisease18", "vkdisease19", 
  #            "vkdisease20", "vkdisease21", "vkdisease23", "vkdisease24", "vkdisease25", "vkdisease26", "vkdisease29", 
  #            "vkdisease30", "vkdisease31", "vkdisease34", "vkdisease38", "vkdisease46", "vkdisease53", "vkdisease62", "vkdisease63",
  #            "vkdisease81", "vkdisease99", "vkdisease133", "kdisease1",  "kdisease2",  "kdisease3", "kdisease5",  "kdisease6",  
  #            "kdisease7",  "kdisease8",  "kdisease9",  "kdisease10","kdisease12", "kdisease13", "kdisease14", "kdisease15",
  #            "kdisease16", "kdisease17", "kdisease18", "kdisease19", "kdisease20", "kdisease21", "kdisease22", "kdisease23",                   
  #            "kdisease24", "kdisease25", "kdisease26", "kdisease29", "kdisease30", "kdisease31", 
  #            "kdisease34", "kdisease35", "kdisease37", "kdisease38", "kdisease40", "kdisease41", 
  #            "kdisease46", "kdisease48", "kdisease49", "kdisease52", "kdisease53", "kdisease54",
  #            "kdisease55", "kdisease56", "kdisease59", "kdisease60", "kdisease61", "kdisease62",
  #            "kdisease63", "kdisease64", "kdisease66", "kdisease68", "kdisease71", "kdisease72",
  #            "kdisease74", "kdisease76", "kdisease78", "kdisease81", "kdisease83", "kdisease86",
  #            "kdisease91", "kdisease92", "kdisease97", "kdisease99", "kdisease101", "kdisease105",
  #            "kdisease111", "kdisease113", "kdisease133", "disease1",   "disease2",   "disease3",  
  #            "disease4",   "disease5",   "disease6",   "disease7",   "disease8",   "disease9",  
  #            "disease10",  "disease12",  "disease13",  "disease14",  "disease15",  "disease16", 
  #            "disease17",  "disease18",  "disease19",  "disease20",  "disease21",  "disease22", 
  #            "disease23",  "disease24",  "disease25",  "disease26",  "disease27",  "disease28", 
  #            "disease29",  "disease30",  "disease31",  "disease32",  "disease33",  "disease34", 
  #            "disease35",  "disease36",  "disease37",  "disease38",  "disease39",  "disease40", 
  #            "disease41",  "disease42",  "disease43",  "disease44",  "disease45",  "disease46", 
  #            "disease47",  "disease48",  "disease49",  "disease50",  "disease51",  "disease52", 
  #            "disease53",  "disease54",  "disease55",  "disease56",  "disease57",  "disease58", 
  #            "disease59",  "disease60",  "disease61",  "disease62",  "disease63",  "disease64", 
  #            "disease65",  "disease66",  "disease67",  "disease68",  "disease69",  "disease70", 
  #            "disease71",  "disease72",  "disease73",  "disease74",  "disease75",  "disease76", 
  #            "disease77",  "disease78",  "disease79",  "disease80",  "disease81",  "disease82", 
  #            "disease83",  "disease84",  "disease85",  "disease86",  "disease87",  "disease88", 
  #            "disease89",  "disease90",  "disease91",  "disease92",  "disease93",  "disease94", 
  #            "disease95",  "disease96",  "disease97",  "disease98",  "disease99",  "disease100", 
  #            "disease101", "disease102", "disease103", "disease104", "disease105", "disease106",
  #            "disease107", "disease108", "disease109", "disease110", "disease111", "disease112",
  #            "disease113", "disease114", "disease115", "disease116", "disease117", "disease118",
  #            "disease119", "disease120", "disease121", "disease122", "disease123", "disease124",
  #            "disease125", "disease126", "disease127", "disease128", "disease129", "disease130",
  #            "disease131", "disease132", "disease133", "disease134", "disease135", "disease136",
  #            "disease137", "disease138", "disease139", "disease140", "disease141", "disease142",
  #            "disease143", "disease144", "disease145", "disease146", "disease147", "disease148",
  #            "disease149", "disease150"))
  # 
  # message("Now checking whether outcome is in file")
  # 
  # if ("outcome" %in% colnames(fulldat)) {
  #   message("'outcome' is present merging it in")
  #   
  #   outcomes <- fulldat %>%
  #     group_by(DrugKey, indicationkey) %>%
  #     summarize(outcome=max(outcome)) %>%
  #     ungroup()
  #   adat <- adat %>%
  #     left_join(outcomes, by=c("DrugKey", "indicationkey")) %>%
  #     ungroup()
  # } else {
  #   message("'outcome' is not present will not be in returned tibble")
  # }
  
  drugno1 <- adat %>% 
    filter(dataset=="training") %>%
    dplyr::select(DrugKey) %>%
    distinct() %>%
    arrange(DrugKey) %>%
    mutate(drugno=row_number())
  drugno2 <- adat %>% 
    filter(dataset=="test") %>%
    dplyr::select(DrugKey) %>%
    distinct() %>% 
    left_join(drugno1, by="DrugKey") %>%
    filter(is.na(drugno)) %>%
    arrange(DrugKey) %>%
    mutate(drugno=max(drugno1$drugno) + row_number())
  drugnos <- drugno1 %>%
    bind_rows(drugno2)
  
  indno1 <- adat %>% 
    filter(dataset=="training") %>%
    dplyr::select(indicationkey) %>%
    distinct() %>%
    arrange(indicationkey) %>%
    mutate(indno=row_number())
  indno2 <- adat %>% 
    filter(dataset=="test") %>%
    dplyr::select(indicationkey) %>%
    distinct() %>% 
    left_join(indno1, by="indicationkey") %>%
    filter(is.na(indno)) %>%
    arrange(indicationkey) %>%
    mutate(indno=max(indno1$indno) + row_number())
  indnos <- indno1 %>%
    bind_rows(indno2)
  
  adat <- adat %>%
    left_join(drugnos, by="DrugKey") %>%
    left_join(indnos, by="indicationkey") %>%
    ungroup()
  
  
  if (file1 == "/shared_data/data/training_data/training_data_2015_split_on_outcome.csv" & is.na(file2)) {
    if (dim(adat)[1] != 4903 ) {
      message( paste("Warning: there should be 4903 records, but there are ", dim(adat)[1], ". Did something go wrong?"))
    } else {
      message( "Great: Expected number of records in generated dataset.")
    }
  }
  
  return(adat)
}

adat <- get_features_bjoern(file1="/shared_data/data/training_data/training_data_2015_split_on_outcome.csv",
                            file2="/shared_data/data/test_data_full/testing_phase2_release.csv")

adat %>%
  group_by(dataset) %>%
  summarise_all(funs(sum(is.na(.)))) %>%
  gather(key="variable", value="value", -dataset) %>%
  filter(value>0) %>%
  arrange(dataset, desc(value)) %>%
  print(n=200)

write_csv(x=adat,path="/files/feat_bjoern.csv")
