# Code shared by Misa Graff

lipid_data <- read.table("cchc_meds_pheno_0323.txt", sep = "|", header = TRUE)

# Convert all the variable names to lowercase
colnames(lipid_data) = tolower(colnames(lipid_data))


count <- sum(lipid_data$med1 != "" | lipid_data$med2 != "" | lipid_data$med3 != "" | lipid_data$med4 != "" | lipid_data$med5 != "" | 
             lipid_data$med6 != "" | lipid_data$med7 != "" | lipid_data$med8 != "" | lipid_data$med9 != "" | lipid_data$med10 != "")
# 3,869

for (var in c("med1", "med2", "med3", "med4", "med5", "med6", "med7", "med8", "med9", "med10")) {
  lipid_data[lipid_data[[var]] == "", var] <- "NA"
}

##Define Statin use

lipid_data$statin <- 0

for (num in 1:10) {
  lipid_data$statin[lipid_data$med[num] %in% c("ATORVASTATIN", "ATROVASTATIN", "FLUVASTATIN", "LOVASTATIN", "PITAVASTATIN", "PRAVASTATIN", "ROSUVASTATIN", "SIMVASTATIN", "STATIN", "LIPITOR", "LESCOL XL", "LESCOL", "ALTOPREV", "LIVALO", "PRAVACHOL", "CRESTOR", "EZALLOR", "ZOCOR", "FLOLIPID")] <- 1
}

##Define Niacinamide use

lipid_data$niacin <- 0

for (num in 1:10) {
  lipid_data$niacin[lipid_data$med[num] %in% c("NIACIN", "VITAMIN B3", "NIASPAN", "NIACOR", "SIMCOR", "NICOTINIC ACID", "NIACIN SR")] <- 1
}


##Define Bile acid sequestrant use

lipid_data$bileacid_sequ <- 0

for (num in 1:10) {
  lipid_data$bileacid_sequ[lipid_data$med[num] %in% c("WELCHOL", "QUESTRAN LIGHT", "QUESTRAN", "PREVALITE", "CHOLESTYRAMINE LIGHT", "COLESTID FLAVORED", "COLESTID", "LOCHOLEST")] <- 1
  lipid_data$bileacid_sequ[lipid_data$med[num] %in% c("COLESEVELAM", "CHOLESTYRAMINE", "COLESTIPOL", "CHOLES")] <- 1
}

##Define Fibrate use

lipid_data$fibrates <- 0

for (num in 1:10) {
  lipid_data$fibrates[lipid_data$med[num] %in% c("GEMFIBROZIL", "LOPID", "FENOFIBRATE", "TRICOR", "FIBRICOR", "LOFRIBRA", "CLOFIBRATE", "ATROMID-S", "FIBRATE")] <- 1
}

##Define cholesterol inhibator use

lipid_data$cholabs_inh <- 0

for (num in 1:10) {
  lipid_data$cholabs_inh[lipid_data$med[num] %in% c("PHYTOSTEROLS", "SOLUBLE FIBERS", "PHOSPHOLIPIDS", "STEARIC ACIDS", "EZETIMIBE", "ZETIA", "SCH-48461", "2-AZETIDINONE")] <- 1
}

##Define metformin use (not as lipid med, but for added PRS descriptive analyses)

lipid_data$metformin <- 0

for (num in 1:10) {
  lipid_data$metformin[lipid_data$med[num] %in% c("METFORMIN", "FORTAMET", "GLUCOPHAGE", "GLUCOPHAGE XR", "GLUMETZA", "RIOMET", "RIOMET ER", "BIGUANIDE")] <- 1
}

##Define NSAID use (not as lipid med, but for added PRS descriptive analyses)

lipid_data$nsaid <- 0

for (num in 1:10) {
  lipid_data$nsaid[lipid_data$med[num] %in% c("ASPIRIN", "ADVIL", "IBUPROFEN", "NAPROXEN", "EXCEDRIN", "MOTRIN", "ALEVE", "ANACIN", "CELECOXIB", "CELEBREX", "DICLFENAC", "VOLTAREN", "INDOMETHACIN", "KETOROLAC", "TORADOL", "MECLOFENAMATE", "DIFLUNISAL", "TOLMETIN", "KETOPROFEN", "FLURBIPROFEN")] <- 1
}

lipid_data <- lipid_data[, !(grepl("*freq", names(lipid_data)))]
lipid <- lipid_data[order(lipid_data$labid),]


# Load "cchcsummarymed28Jan22.txt" dataset in R
data <- read.table("cchcsummarymed28Jan22.txt", sep = "|", header = TRUE)

# Select relevant variables
data <- data[, c("med_lipid", "cholmedcode", "rrid", "labid")]
data <- data[order(data$labid), ]

# Merge datasets
data <- merge(data, lipid, by = "labid", all.x = TRUE)

# Recode med_lipid variable based on other variables
data$med_lipid[data$statin == 1] <- 1
data$med_lipid[data$fibrates == 1] <- 1
data$med_lipid[data$cholabs_inh == 1] <- 1
data$med_lipid[data$bileacid_sequ == 1] <- 1
data$med_lipid[data$niacin == 1] <- 1
data$med_lipid[is.na(data$med_lipid)] <- 0

# Recode other variables
data$statin[data$cholmedcode == "statin"] <- 1
data$statin[is.na(data$statin)] <- 0
data$fibrates[is.na(data$fibrates)] <- 0
data$cholabs_inh[is.na(data$cholabs_inh)] <- 0
data$bileacid_sequ[is.na(data$bileacid_sequ)] <- 0
data$niacin[is.na(data$niacin)] <- 0
data$metformin[is.na(data$metformin)] <- 0

# Replace empty cholmedcode with "NA"
data$cholmedcode[data$cholmedcode == ""] <- "NA"

# Remove unnecessary variables
data <- data[order(data$labid), ]
data <- data[, !grepl("^med[0-9]+$", names(data))]

# Save "temp" dataset
write.table(data, "temp.csv", sep = ",", row.names = FALSE, quote = FALSE)

# Save final output
write.table(data, "CCHC_MED1-10_Lipids_classes_NSAIDs_2023.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Remove med1-med10 variables
data <- data[, !grepl("^med[0-9]+$", names(data))]

# Save final output without med1-med10 variables
write.table(data, "CCHC_Lipid_classes_2023.txt", sep = "\t", row.names = FALSE, quote = FALSE)

