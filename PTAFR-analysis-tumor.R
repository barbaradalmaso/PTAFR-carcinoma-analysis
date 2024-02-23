# -------------------------------------------- Part 1: PTAFR expression on global tumor samples and select high and low PAFR groups -------------------------------------------------------
# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(ggrepel)
library(jsonlite)
##### First, I will extract metadata from the JSON file to create a list of file_names, case_id, tissue, and type (tumor or normal tissue) for further analysis. This part was performed using Linux system at DSMZ.
# Case_id and file_name information
raw.metadata <- fromJSON("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata.cart.2023-08-09.json")
raw.metadata_bladder <- fromJSON("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/metadata.cart.2023-09-06.json")
raw.metadata_breast <- fromJSON("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/metadata.cart.2023-09-06.json")
raw.metadata <- rbind(raw.metadata[,3:4], raw.metadata_bladder[,3:4], raw.metadata_breast[,3:4])

case_ids <- character(nrow(raw.metadata))
for (i in 1:nrow(raw.metadata)) { # Loop for extract case_id information from JSON files

	case_ids[i] <- raw.metadata[[1]][[i]]$case_id
}

case_ids <- unlist(case_ids)
raw.metadata$case_id <- case_ids

# Tissue information 
clinical <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/clinical_data/clinical.tsv", header = TRUE, sep = "\t") 
clinical_bladder <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/clinical.tsv", header = TRUE, sep = "\t") 
clinical_breast <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/clinical.tsv", header = TRUE, sep = "\t") 

clinical <- rbind(clinical[,c(1,123)], clinical_bladder[,c(1,123)], clinical_breast[,c(1,123)])

# Type information
sample <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/sample_data/sample.tsv", header = TRUE, sep = "\t")
sample_bladder <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/sample.tsv", header = TRUE, sep = "\t")
sample_breast <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/sample.tsv", header = TRUE, sep = "\t")

sample <- rbind(sample[,c(2,35)], sample_bladder[,c(2,35)], sample_breast[,c(2,35)])

# Merge all, clean tissue types and save file
metadata <- merge(clinical, raw.metadata[,2:3], by = "case_id")
metadata <- merge(metadata, sample, by = "case_id")
metadata_all <- unique(metadata)

metadata_tumor <- metadata[metadata$tissue_type != "Normal", ]

# Clean tissues types
Bladder = c("Anterior wall of bladder",
		  "Bladder neck",
		  "Bladder, NOS",
		  "Dome of bladder",
		  "Gallbladder",
		  "Lateral wall of bladder",
		  "Posterior wall of bladder",
		  "Trigone of bladder ",
		  "Ureteric orifice")

Colorectal = c("Ascending colon",
			 "Cecum",
			 "Colon, NOS",
			 "Connective, subcutaneous and other soft tissues of abdomen",
			 "Descending colon",
			 "Hepatic flexure of colon",
			 "Rectosigmoid junction",
			 "Rectum, NOS",
			 "Sigmoid colon",
			 "Small intestine, NOS",
			 "Splenic flexure of colon",
			 "Transverse colon")

Lung = c("Lower lobe, lung",
	    "Lung, NOS",
	    "Middle lobe, lung",
	    "Overlapping lesion of lung",
	    "Upper lobe, lung",
	    "Main bronchus")

Stomach = c("Body of stomach",
		  "Cardia, NOS",
		  "Esophagus, NOS",
		  "Fundus of stomach",
		  "Gastric antrum",
		  "Lesser curvature of stomach, NOS",
		  "Lower third of esophagus",
		  "Middle third of esophagus",
		  "Pylorus",
		  "Stomach, NOS",
		  "Thoracic esophagus")

Uterus = c("Cervix uteri",
		 "Corpus uteri",
		 "Endometrium",
		 "Fundus uteri",
		 "Isthmus uteri",
		 "Uterus, NOS")

Pancreas = c("Ampulla of Vater",
		   "Body of pancreas",
		   "Head of pancreas",
		   "Overlapping lesion of pancreas",
		   "Pancreas, NOS",
		   "Tail of pancreas",
		   "Extrahepatic bile duct") 

Breast = c("Breast, NOS",
		 "Lower-inner quadrant of breast",
		 "Lower-outer quadrant of breast",
		 "Overlapping lesion of breast",
		 "Upper-inner quadrant of breast",
		 "Upper-outer quadrant of breast")

metadata_tumor <- metadata_tumor %>%
	mutate(tissue_or_organ_of_origin = case_when(
		tissue_or_organ_of_origin %in% Bladder ~ "Bladder",
		tissue_or_organ_of_origin %in% Colorectal ~ "Colorectal",
		tissue_or_organ_of_origin %in% Lung ~ "Lung",
		tissue_or_organ_of_origin %in% Stomach ~ "Stomach",
		tissue_or_organ_of_origin %in% Uterus ~ "Uterus",
		tissue_or_organ_of_origin %in% Pancreas ~ "Pancreas",
		tissue_or_organ_of_origin %in% Breast ~ "Breast",
		tissue_or_organ_of_origin == "Kidney, NOS" ~ "Kidney",
		tissue_or_organ_of_origin == "Prostate gland" ~ "Prostate",
		tissue_or_organ_of_origin == "Thyroid gland" ~ "Thyroid",
		TRUE ~ NA
	))

# Save cleaned metadata table
write.csv(metadata_tumor, file = "~/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata.csv", row.names = FALSE)

### Now, With the organized metadata, I will use the column "file_names" to export RNA-seq files per tissue for PAFR expression analysis.
metadata <- read.delim("/Volumes/Extreme SSD/DSMZ/files_processed/metadata.csv", header = TRUE, sep = ",")
metadata <- na.omit(metadata)
colnames(metadata)[2] <- "tissue"
metadata <- unique(metadata)

# Download a separate table as an example
table <- read.table("/Volumes/Extreme SSD/DSMZ/05ea028e-6b8a-4352-88aa-057494cf5e80.rna_seq.augmented_star_gene_counts.tsv", header = T, sep ="\t", skip = 1)
table <- table[,c(1:3, 7)]
colnames(table) <- c("gene_id", "gene_name", "gene_type", "tpm" )
table <- subset(table, gene_type == "protein_coding")

# Set the directory and get unique tissue names
directory <- "/Volumes/Extreme SSD/DSMZ/files_raw" 
tissues <- unique(metadata$tissue)

# Function to create a table for each tissue
table_for_tissue <- function(x, metadata, directory) {
    metadata_files <- metadata %>%
        filter(tissue == x)
    
    # List to store tables
    tables_for_tissue <- list()
    
    # Loop on metadata for each tissue
    for (file_name in metadata_files$file_name) {
        file_path <- file.path(directory, file_name)
        
        # Check if the file exists before attempting to read it
        if (file.exists(file_path)) {
            file <- read.table(file_path, header = TRUE, sep = "\t")
            
            # Extract relevant columns (gene_name, gene_type, tpm_unstranded) and add to the list
            filtered_data <- file %>%
                select(gene_name, tpm_unstranded) %>%
                filter(gene_name == "PTAFR")
            
            tables_for_tissue[[file_name]] <- filtered_data
            cat("File successfully downloaded:", file_path, "\n")
        } else {
            cat("File not found:", file_path, "\n")
        }
    }
    
    return(tables_for_tissue)
}

# List to store tables for each tissue
table_list <- list()

# Loop over the unique tissues to create separate tables
for (tissue_name in tissues) {
    table_list[[tissue_name]] <- table_for_tissue(tissue_name, metadata, directory)
}

# Clean table
# Create a loop for extracting data from the list, naming PTAFR expression according to file_name
tissues <- names(table_list)
tissues.df <- data.frame()

for (i in tissues) {
	tissue.list <- table_list[[i]]
	for (j in names(tissue.list)) {
		sample.df <- tissue.list[[j]]
		colnames(sample.df)[2] <- j
		if (ncol(tissues.df) == 0) {
			tissues.df <- sample.df
		} else {
			tissues.df <- cbind(tissues.df, sample.df[, 2, drop = FALSE])
			colnames(tissues.df)[ncol(tissues.df)] <- j
		}
	}
}

PTAFR_counts <- t(tissues.df[,2:ncol(tissues.df)])
colnames(PTAFR_counts) <- "PTAFR"

PTAFR_counts <- data.frame(PTAFR = PTAFR_counts,
					  file_name = rownames(PTAFR_counts))

PTAFR_counts_data <- merge(PTAFR_counts, metadata, by = "file_name")
PTAFR_counts_data$PTAFR <- as.numeric(PTAFR_counts_data$PTAFR)

# Quartile calculation
q1 <- quantile(PTAFR_counts_data$PTAFR, 0.25, na.rm = TRUE)
q2 <- quantile(PTAFR_counts_data$PTAFR, 0.50, na.rm = TRUE)
q3 <- quantile(PTAFR_counts_data$PTAFR, 0.75, na.rm = TRUE)

# Histogram
# Configure PTAFR expression table
PTAFR_neg <- subset(PTAFR_counts_data, PTAFR <= q1)
PTAFR_neg$type <- c(rep("Low", nrow(PTAFR_neg)))

PTAFR_pos <- subset(PTAFR_counts_data, PTAFR >= q3)
PTAFR_pos$type <- c(rep("High", nrow(PTAFR_pos)))

PTAFR_neutral <- subset(PTAFR_counts_data, PTAFR <= q3 & PTAFR >= q1)
PTAFR_neutral$type <- c(rep("Control", nrow(PTAFR_neutral)))

PTAFR_expression <- rbind(PTAFR_neg, PTAFR_pos, PTAFR_neutral)
PTAFR_expression$tissue <- str_to_title(PTAFR_expression$tissue)
order <- c("Uterus", "Thyroid", "Stomach", "Prostate", "Pancreas", "Lung", "Kidney", "Colorectal", "Breast", "Bladder")
PTAFR_expression$tissue <- factor(PTAFR_expression$tissue, levels = order)

# Plote o gráfico
ggplot(PTAFR_expression, aes(x = PTAFR, y = tissue, fill = tissue)) +
	geom_density_ridges(
		jittered_points = TRUE,
		position = position_points_jitter(width = 0.05, height = 0),
		point_shape = '|', point_size = 1.5, point_alpha = 0.4, point_color = "black", alpha = 0.8, color = NA
	) +
	geom_vline(xintercept = q2, color = "black", size = 0.25, linetype = "dashed") +
	geom_vline(xintercept = q1, color = "black", size = 0.25) +
	geom_vline(xintercept = q3, color = "black", size = 0.25) +
	scale_fill_manual(values = cores_tissue) +
	labs(x = "PTAFR Expression Log???(TPM)") +
	theme(
		axis.title.x = element_text(hjust = 0.5),
		axis.title.y = element_blank(),
		axis.text.y = element_text(size = 9),
		legend.position = "none",
		panel.border = element_rect(color = "black", size = 1.5, fill = NA),
		panel.background = element_rect(fill = "white")
	) +
	coord_cartesian(xlim = c(-2, 100))

# Histogram for all tissues
cores_tissue <- c("#001F3F", "#003366", "#004080", "#005A8D", "#0072BB", "#0088CC", "#3299CC", "#66B3CC", "#99CCFF", "#B3E0FF")

######## Basic histogram
ggplot(PTAFR_expression, aes(x = PTAFR, fill = tissue)) +
	geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.8) +
	geom_vline(xintercept = q2, color = "black", size = 0.3, linetype = "dashed") + # Line median
	geom_vline(xintercept = q1, color = "black", size = 0.3) + # line quartile 1
	geom_vline(xintercept = q3, color = "black", size = 0.3) + # line quartile 2
	labs(x = NULL, fill = "Type") +
	scale_fill_manual(values = cores_tissue) +  # Definir a paleta de cores manualmente
	theme(
		axis.title.y = element_blank(),
		axis.text.y = element_text(size = 10),
		axis.text.x = element_blank(),
		legend.position = "none",
		panel.border = element_rect(color = "black", size = 1.5, fill = NA),
		panel.background = element_rect(fill = "white")) +
	coord_cartesian(xlim = c(-2,100))

PTAFR_expression$log2_transformed <- log2(PTAFR_expression$PTAFR)
PTAFR_expression$log2_transformed <- ifelse(is.finite(PTAFR_expression$log2_transformed), PTAFR_expression$log2_transformed, NA)

# Individual histogram for each type of tumor
	ggplot(PTAFR_expression, aes(x = PTAFR, y = tissue, fill = tissue)) +
	geom_density_ridges(
		jittered_points = TRUE,
		position = position_points_jitter(width = 0.05, height = 0),
		point_shape = '|', point_size = 1.5, point_alpha = 0.4, point_color = "black", alpha = 0.8, color = NA) +
	geom_vline(xintercept = q2, color = "black", size = 0.25, linetype = "dashed") +
	geom_vline(xintercept = q1, color = "black", size = 0.25) +
	geom_vline(xintercept = q3, color = "black", size = 0.25) +
	scale_fill_manual(values = cores_tissue) +
	labs(x = "PTAFR Expression Log2(TPM)") +
		theme(
		axis.title.x = element_text(hjust = 0.5),
		axis.title.y = element_blank(),
		axis.text.y = element_text(size = 9),
		legend.position = "none",
		panel.border = element_rect(color = "black", size = 1.5, fill = NA),
		panel.background = element_rect(fill = "white")
	) +
	coord_cartesian(xlim = c(-2,100))

# PTAFR postive and negative for each type of tumor
count <- PTAFR_expression %>%
	filter(type != "Control")

count <- count %>%
	group_by(tissue, type) %>% 
	summarise(count = n())

count <- count %>%
	group_by(tissue) %>%
	mutate(sum = sum(count))

count_percentage <- count
count_percentage$perc <- count_percentage$count / count_percentage$sum

count_percentage <- count_percentage %>%
	select(-count, -sum)

count_percentage$tissue <- as.factor(count_percentage$tissue)
count_percentage$type <- as.factor(count_percentage$type)
count_percentage$perc <- count_percentage$perc * 100
count_percentage$perc <- round(count_percentage$perc)
count_percentage$type <- gsub("Negative", "Low", count_percentage$type)
count_percentage$type <- gsub("Positive", "High", count_percentage$type)

# Plot
ggplot(count_percentage, aes(x = tissue, y = perc, fill = tissue, alpha = type)) +
	geom_col(aes(position = "stack")) +
	labs(x = NULL, y = "Patients (%)", alpha = "PTAFR", fill = NULL) +
	theme(axis.text.x = element_text(size = 8),
		 axis.text.y = element_blank(),
		 axis.title.x = element_text(size = 9),
		 axis.line.x = element_line(lineend = "butt"),
		 panel.grid.major.y = element_blank(),
		 panel.grid.minor.y = element_blank(),
		 panel.grid.major.x = element_blank(),
		 panel.background = element_blank()) +
	coord_flip() +
	scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
	scale_fill_manual(values = cores_tissue) +
	scale_alpha_manual(values = c(High = 0.8, Low = 0.4))


# Plot for figure on PTAFR protein expression downloaded from The Human Protein Atlas https://www.proteinatlas.org/ENSG00000169403-PTAFR/pathology/ovarian+cancer#Quantity
protein_PTAFR <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/protein_expression.csv", header = TRUE, sep = ",")

protein_PTAFR$Percentage <- factor(protein_PTAFR$Percentage, 
							levels = c("<75%", "75%-25%", "<25%", "None"))

protein_PTAFR$Tissue <- factor(protein_PTAFR$Tissue, 
							levels = order)

ggplot(protein_PTAFR, aes(x = Tissue, y = Counting, fill = Percentage)) +
	geom_col(aes(position = "stack")) +
	labs(x = NULL, y = "Patients (%)", alpha = "PTAFR", fill = "Positive Tumor Cells") +
	theme(axis.text.x = element_text(size = 10),
		 axis.text.y = element_text(size = 10),
		 axis.title.x = element_text(size = 12),
		 axis.line.x = element_line(lineend = "butt"),
		 panel.grid.major.y = element_blank(),
		 panel.grid.minor.y = element_blank(),
		 panel.grid.major.x = element_blank(),
		 panel.background = element_blank()) +
	coord_flip() +
	scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
	scale_fill_manual(values = c("#1b9e77", "#d95f02", "#001F3F", "gray"))

# Save metadata table
write.table(PTAFR_expression, file = "~/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata_ptafr.csv", row.names = FALSE)

# -------------------------------------------- Part 2: Metadata analysis: SEER Survival, % Tumor size, and Metadata T and N  -------------------------------------------------------
# Survival rates analysis using SEER-NIH 5-year relative survival 
# Load packages
	library(ggplot2)
	library(dplyr)
	library(ggrepel)
	library(stringr)
	
PAFR_expression <- read.table("/Volumes/Extreme SSD/DSMZ/files_processed/metadata_ptafr.csv", header = TRUE, sep = " ") # Updated to SSD (15.12.23)
	
# Create a data frame based on the values obtained on the SEER website
survival <- data.frame(tissue = c("Breast", "Uterus", "Colorectal", "Stomach", "Kidney", "Lung", "Pancreas", "Prostate", "Thyroid", "Bladder"),
				   survival = c(87.9, 81.4, 62.9, 24.9, 67.5, 17.8, 5.8,
				   		   98.0, 96.9, 74.3))

sample_size <- PAFR_expression %>% # take off control group
	group_by(tissue, type) %>%
	summarize(count = n())

sample_size <- sample_size %>%
	group_by(tissue) %>%
	mutate(sample = sum(count))

sample_size <- sample_size[,c(-2,-3)]
sample_size <- unique(sample_size)
	
ptafr_rates <- PAFR_expression %>%
	group_by(tissue) %>%
	summarize(mean = mean(PTAFR),
			sd_value = sd(PTAFR))

ptafr_rates <- merge(ptafr_rates, sample_size, by = "tissue")
survival <- merge(survival, ptafr_rates, by = "tissue")
survival <- unique(survival)
correlation_test <- cor.test(survival$survival, survival$mean, method = "pearson")

# SEER survival plot
cores_tissue <- c("#001F3F", "#003366", "#004080", "#005A8D", "#0072BB", "#0088CC", "#3299CC", "#66B3CC", "#99CCFF", "#B3E0FF")

ggplot(data = survival, aes(x = survival, y = log2(mean), color = tissue)) +
	geom_point(aes(size = sample), shape = 16) +  
	geom_smooth(method = "lm", color = "black", alpha = 0.15) + 
	labs(x = "Tumor 5-Year Relative Survival (%)", y = "PTAFR Expression Log2(TPM)", size = "Sample size") +
	scale_color_manual(values = cores_tissue) +
	theme(axis.title.y = element_text(size = 12),
	axis.text.y = element_text(size = 12),
	axis.text.x = element_text(size = ),
	panel.border = element_rect(color = "black", size = 1.5, fill = NA),
	panel.background = element_rect(fill = "white"))

ggplot(data = survival, aes(x = survival, y = log2(mean), color = tissue)) +
	geom_point(aes(size = sample), shape = 16) +
	geom_smooth(method = "lm", color = "black", alpha = 0.15) +
	labs(x = "Time (months)", y = "Patients (%)") +
	scale_color_manual(values = cores_tissue) +
	theme(axis.title.y = element_text(size = 12),
		 axis.text.y = element_text(size = 12),
		 axis.text.x = element_text(size = 12),  # Defina o tamanho do texto no eixo x
		 panel.border = element_rect(color = "black", size = 1.5, fill = NA),
		 panel.background = element_rect(fill = "white")) +
	xlim(0, 300) +  # Define os limites do eixo x
	ylim(0, 150)    # Define os limites do eixo y


# Metadata with tumor size and tumor stage
###### Survival analysis based on clinical data (days to death)
# Tumor size
metadata_ptafr <- read.delim2("/Volumes/Extreme SSD/DSMZ/files_processed/metadata_ptafr.csv", sep = " ", header = TRUE)
clinical <- read.delim("/Volumes/Extreme SSD/DSMZ/metadata/pathology_detail.tsv", sep = " ", header = TRUE)

clinical <- clinical[,c(1,17,48,49)]
metadata_clinical <- merge(clinical, metadata_ptafr, by = "case_id")
metadata_clinical <- metadata_clinical %>% 
	distinct(case_id, .keep_all = TRUE)
metadata_clinical$gross_tumor_weight <- gsub("'--", NA, metadata_clinical$gross_tumor_weight)
metadata_clinical$tumor_largest_dimension_diameter <- gsub("'--", NA, metadata_clinical$tumor_largest_dimension_diameter)
metadata_clinical$tumor_thickness <- gsub("'--", NA, metadata_clinical$tumor_thickness)
metadata_clinical <- metadata_clinical[,c(1,3,9)]
metadata_clinical$tumor_largest_dimension_diameter <- as.numeric( metadata_clinical$tumor_largest_dimension_diameter)
metadata_clinical <- metadata_clinical %>%
	filter(type != "Control")
high <- metadata_clinical %>%
	filter(type != "Low")
low <- metadata_clinical %>%
	filter(type != "High")
low <- low[c(-10, -22, -27, -78),] # Remove outliers
test <- wilcox.test(high$tumor_largest_dimension_diameter, low$tumor_largest_dimension_diameter)

metadata_clinical <- rbind(high, low)
write.table(metadata_clinical, file = "metadata_clinical.tsv", sep = "\t", row.names = FALSE)
######## Finished part tumor tickness

# Tumor stage
clinical_data <- read.delim("/Volumes/Extreme SSD/DSMZ/metadata/clinical.tsv", sep = " ", header = TRUE)
clinical_data <- clinical_data[,c(1,26:29)]

# Clean data
clinical_data <- clinical_data %>% mutate_all(~ifelse(. %in% c("'--", "Unknown", "Not Reported"), NA, .)) # Add NA on unknown data
clinical_data <- clinical_data %>% 
	distinct(case_id, .keep_all = TRUE)

## Pathological cleaning 
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_m = ifelse(ajcc_pathologic_m %in% c("cM0 (i+)"), "M1", ajcc_pathologic_m)) # pathologic m
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IA", "Stage IB"), "Stage I", ajcc_pathologic_stage)) # pathologic stage I
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IIA", "Stage IIB", "Stage IIC"), "Stage II", ajcc_pathologic_stage)) # pathologic stage II
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IIIA", "Stage IIIA1", "Stage IIIB", "Stage IIIC", "Stage IIIC1", "Stage IIIC2"), "Stage III", ajcc_pathologic_stage)) # pathologic stage II
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IVA", "Stage IVB", "Stage IVC"), "Stage IV", ajcc_pathologic_stage)) # pathologic stage IV
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_t = ifelse(ajcc_pathologic_t %in% c("T1a", "T1b", "T1b1", "T1b2", "T1c"), "T1", ajcc_pathologic_t)) # clinical t 
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_t = ifelse(ajcc_pathologic_t %in% c("T2a", "T2a1", "T2a2", "T2b", "T2c"), "T2", ajcc_pathologic_t)) # clinical t 
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_t = ifelse(ajcc_pathologic_t %in% c("T3a", "T3b", "T3c"), "T3", ajcc_pathologic_t)) # clinical t 
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_t = ifelse(ajcc_pathologic_t %in% c("T4a", "T4b", "T4c", "T4d", "Tis"), "T4", ajcc_pathologic_t)) # clinical t 
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_n = ifelse(ajcc_pathologic_n %in% c("N1a", "N1b", "N1c", "N1mi"), "N1", ajcc_pathologic_n)) # clinical n
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_n = ifelse(ajcc_pathologic_n %in% c("N2a", "N2b"), "N2", ajcc_pathologic_n)) # clinical n
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_n = ifelse(ajcc_pathologic_n %in% c("N3a", "N3b", "N3c"), "N3", ajcc_pathologic_n)) # clinical n
clinical_data <- clinical_data %>% mutate(ajcc_pathologic_n = ifelse(ajcc_pathologic_n %in% c("N0 (i-)", "N0 (i+)"), "N0", ajcc_pathologic_n)) # clinical n

clinical_data <- merge(clinical_data, PAFR_expression, by = "case_id")
clinical_data <- clinical_data %>%
	filter(type != "Control") # take off control

clinical_data <- clinical_data %>% 
	distinct(case_id, .keep_all = TRUE)

write.table(clinical_data, file = "clinical_data", sep = "\t", row.names = FALSE)

#------------------------------------ Part 3: Preparing samples for DESEQ Analysis -------------------------------------------------------
# Create separated data frames based on PAFR low and high expression and tissue type
# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)

# Load selected data from metadata_ptafr from different tumor types
metadata_ptafr <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata_ptafr.csv", header = TRUE, sep = " ")

metadata <- metadata_ptafr[metadata_ptafr$type != "Control", ]

# Set the directory and get unique tissue names
directory <- "~/Desktop/Patients_data/files_counts"
tissues <- unique(metadata$tissue)

# Function to create a table for each tissue
table_for_tissue <- function(x, metadata, directory) {
	metadata_files <- metadata %>%
		filter(tissue == x)
	
	# List to store tables
	tables_for_tissue <- list()
	
	# Loop on metadata for each tissue
	for (file_name in metadata_files$file_name) {
		file_path <- file.path(directory, file_name)
		
		# Check if the file exists before attempting to read it
		if (file.exists(file_path)) {
			file <- read.table(file_path, header = TRUE, sep = "\t", skip = 1)
			
			# Extract relevant columns (gene_type, tpm_unstranded) and add to the list
			filtered_data <- file %>%
				select(gene_type, tpm_unstranded) %>%
				filter(!is.na(tpm_unstranded)) %>%
				filter(gene_type == "protein_coding") %>%
				select(tpm_unstranded)
				
			tables_for_tissue[[file_name]] <- filtered_data
			cat("File sucessfully found:", file_path, "\n")
		} else {
			cat("File not found:", file_path, "\n")
		}
	}
	
	return(tables_for_tissue)
}

# List to store tables for each tissue and loop
table_list <- list() 
for (tissue_name in tissues) {
	table_list[[tissue_name]] <- table_for_tissue(tissue_name, metadata, directory)
}

# First, I'll try to rename tpm_unstranded according to file name
# Loop over the tissues in the table_list and rename 'tpm_unstranded'
table_list <- lapply(table_list, function(tissue_df) {
	lapply(names(tissue_df), function(df_name) {
		df <- tissue_df[[df_name]]
		colnames(df) <- df_name
		return(df)
	})
})

# Inicialize combined_dataframe with same line number as individual dataframes
combined_dataframe <- data.frame(matrix(NA, nrow = nrow(table_list[[1]][[1]]), ncol = 0))

# Loop over dataframes on list of lists
for (tissue_df in table_list) {
	# Loop over dataframes on tissue df
	for (df in tissue_df) {
		# Add current dataframe with new column on combined dataframe
		combined_dataframe <- cbind(combined_dataframe, df)
		cat("File successfully downloaded:")
	}
}

## Create a separated metadata table for each tissue type
for (i in metadata$tissue) {
	filtered_df <- metadata %>%
		filter(tissue == i) %>%
		select(file_name, case_id, tissue, type)
	
	assign(paste0("metadata.", i), filtered_df, envir = .GlobalEnv)
}

# Create a separated counts table for each tissue type
tissues <- unique(metadata$tissue)
metadata_names <- paste("metadata.", tissues, sep = "")

for (i in metadata_names) {
	current_df <- get(i)
	col_name <- as.character(current_df$file_name)
	selected_col <- combined_dataframe[col_name]
	assign(paste0("cts_", i), selected_col, envir = .GlobalEnv)
}

### Save new dataframes
## counts tables
tissues <- unique(metadata$tissue)
cts_names <- paste0("cts_metadata.", tissues)
dir_path <- "/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/"

for (i in cts_names) {
	table <- get(i)
	write.table(table, file = paste0 (dir_path, i, ".tsv"), sep = "\t", row.names = FALSE)
}

## Metadata tables
## counts tables
tissues <- unique(metadata$tissue)
metadata_names <- paste0("metadata.", tissues)
dir_path <- "/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/"

for (i in metadata_names) {
	table <- get(i)
	write.table(table, file = paste0 (dir_path, i, ".tsv"), sep = "\t", row.names = FALSE)
}

# -------------------------------------------------- Part 4: DESEQ Analysis --------------------------------------------------------
# Perform DESEQ analysis over different tumor types according to PTAFR low and high
# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(DESeq2)

# Load data
# Linux
metadata_ptafr <- read.delim("/Volumes/Extreme SSD/DSMZ/files_processed/metadata_ptafr.csv", header = TRUE, sep = " ")

# Cts files
tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("cts_metadata.", tissues, ".tsv")

for (i in file_names) {
	table <- read.delim(paste0("/Volumes/Extreme SSD/DSMZ/files_processed/", i), header = TRUE, sep = "\t")
	assign(gsub(".tsv", "", i), table, envir = .GlobalEnv)
}

for (i in file_names) {
	i <- gsub(".tsv", "", i)
	table <- get(i)
	colnames(table) <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", colnames(table))
	colnames(table) <- gsub("X", "", colnames(table))
	colnames(table) <- gsub("\\.", "-", colnames(table))
	assign(i, table, envir = .GlobalEnv)
}

# Metadata files
tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("metadata.", tissues, ".tsv")
for (i in file_names) {
	table <- read.delim(paste0("/Volumes/Extreme SSD/DSMZ/files_processed/", i), header = TRUE, sep = "\t")
	table$file_name <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", table$file_name)
	assign(gsub(".tsv", "", i), table, envir = .GlobalEnv)
}

# I will downlod one table as an example
table <- read.table("/Volumes/Extreme SSD/DSMZ/05ea028e-6b8a-4352-88aa-057494cf5e80.rna_seq.augmented_star_gene_counts.tsv", header = T, sep ="\t", skip = 1)
table <- table[,c(1:3, 7)]
colnames(table) <- c("gene_id", "gene_name", "gene_type", "tpm" )
table <- subset(table, gene_type == "protein_coding")

# Transform cts tables on matrix
tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("cts_metadata.", tissues)

for (i in file_names) {
	file <- get(i)
	file <- as.matrix(file) 
	rownames(file) <- table$gene_name
	assign(i, file, envir = .GlobalEnv)
}

# Get coldata
rownames(metadata.Colorectal) <- colnames(cts_metadata.Colorectal)
rownames(metadata.Kidney) <- colnames(cts_metadata.Kidney)
rownames(metadata.Lung) <- colnames(cts_metadata.Lung)
rownames(metadata.Pancreas) <- colnames(cts_metadata.Pancreas)
rownames(metadata.Prostate) <- colnames(cts_metadata.Prostate)
rownames(metadata.Stomach) <- colnames(cts_metadata.Stomach)
rownames(metadata.Thyroid) <- colnames(cts_metadata.Thyroid)
rownames(metadata.Uterus) <- colnames(cts_metadata.Uterus)

# Deseq analysis
metadata_names <- paste0("metadata.", tissues)
cts_names <- paste0("cts_metadata.", tissues)
metadata_names <- sort(metadata_names)
cts_names <- sort(cts_names)
for (i in seq_along(cts_names)) {
	cts <- get(cts_names[i])
	metadata <- get(metadata_names[i])

		if (!is.null(cts) && !is.null(metadata)) {
		res <- DESeqDataSetFromMatrix(countData = round(cts),
								colData = metadata,
								design = ~type)
		deseq_name <- paste0("deseq.", tissues[i])
		assign(deseq_name, res, envir = .GlobalEnv)
		cat("DESeqDataSet criado para", tissues[i], "\n")
	} else {
		cat("Erro: no 'cts.' or 'metadata.' for", tissues[i], "\n")
	}
}

file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	file <- get(i)
	res <- DESeq(file)
	res <- results(res, contrast = c("type", "High", "Low"))
	res <- data.frame(res)
	res$diffexpressed <- "no"
	res$diffexpressed[res$log2FoldChange > 1 & res$padj < 0.05] <- "up"
	res$diffexpressed[res$log2FoldChange < -1 & res$padj < 0.05] <- "down"
	assign(i, res, envir = .GlobalEnv)
	}

# Volcano-plot (general overview)
file_names <- paste0("deseq.", tissues)
file_names <- sort(file_names)
tissues <- sort(tissues)

for (variable in file_names) {

	file <- get(variable) 
	file <- file[row.names(file) != "PTAFR", ]
	
	plot <- ggplot(data = file,
				aes(x = log2FoldChange,
				    y = -log10(padj),
				    col = diffexpressed)) +
		geom_point() +
		labs(title = variable) +
		scale_color_manual(values = c("#1E90FF", "#a9a9a9", "#FF0000")) +
		theme(panel.background = element_blank(),
			 plot.background = element_blank(),
			 panel.border = element_rect(color = "black", fill = NA, size = 2),
			 legend.position = "none",
			 axis.text.x = element_text(size = 18), 
			 axis.text.y = element_text(size = 18),
			 axis.title.y = element_blank(),
			 axis.title.x = element_blank())
	
	print(plot)
	}

# Count % DEG
file_names 
df <- data.frame()
for (i in file_names) {
	file <- get(i)
	gene_counts <- data.frame(table(file$diffexpressed))
		if (ncol(df) == 0) {
		df <- gene_counts
		colnames(df) <- c("Genes", i)
	} else {
		df[[i]] <- gene_counts$Freq
	}
}
rownames(df) <- df$Genes
df <- df[,-1]
DEG <- data.frame(t(df))
total <- rowSums(DEG)
DEG$total <- total

DEG <- DEG %>%
	mutate(perc.up = 100*up/total,
		  perc.down = 100*down/total,
		  perc.no = 100*no/total)

DEG <- DEG[,c(-1,-2,-3,-4)]
DEG$tissue <- c("Breast", "Pancreas",
			 "Thyroid", "Bladder", "Kidney", "Prostate", 
			 "Colorectal", "Lung", "Uterus", "Stomach")

DEG <- DEG %>%
	arrange(tissue)

DEG_all <- data.frame(perc = c(DEG$perc.up, DEG$perc.down, DEG$perc.no),
				  tissue = c(DEG$tissue, DEG$tissue, DEG$tissue),
				  type = c(rep.int("Up", 10), rep.int("Down", 10),
				  	    rep.int("No", 10)))

# Plot
DEG_all$tissue <- factor(DEG_all$tissue, levels = rev(unique(DEG_all$tissue)))

DEG_all$type <- factor(DEG_all$type, levels = c("No", "Up", "Down"))

ggplot(DEG_all, aes(x = tissue, y = perc, fill = type)) +
	geom_col(position = "stack") +
	labs(x = NULL, y = "Genes (%)", fill = "Differentially Expressed") +
	theme(axis.text.x = element_text(size = 10),
		 axis.title.x = element_text(size = 12),
		 axis.line.x = element_line(lineend = "butt"),
		 axis.text.y = element_text(size = 12),
		 panel.grid.major.y = element_blank(),
		 panel.grid.minor.y = element_blank(),
		 panel.grid.major.x = element_blank(),
		 panel.background = element_blank()) +
	scale_fill_manual(values = c(Down = "#001F3F", 
						    No = "#e0e0e0", 
						    Up = "#d95f02")) +
		coord_flip() +
	scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
	

# PCA-plot (optional?)
metadata_names <- paste0("metadata.", tissues)
cts_names <- paste0("cts.", tissues)
metadata_names <- sort(metadata_names)
cts_names <- sort(cts_names)

for (i in cts_names) {
	for (o in metadata_names) {
		cts <- get(i)
		coldata <- get(o)
		
		if (ncol(cts) == nrow(coldata)) { 
			
			# Create DESeqDataSetFromMatrix
			res_data <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~type)
			
			# Apply VST transformation
			res_data <- vst(res_data, blind = FALSE)
			
			# Create PCA plot
			plot <- plotPCA(res_data, intgroup = "type") +
				geom_point(aes(color = type)) + 
				labs(title = i) +
				scale_color_manual(values = c("#00DAE0", "#FF9289")) +
				theme(panel.background = element_blank(),
					 plot.background = element_blank(),
					 panel.border = element_rect(color = "black", fill = NA, size = 2),
					 legend.position = "none",
					 axis.text.x = element_text(size = 18), 
					 axis.text.y = element_text(size = 12),
					 axis.title.y = element_text(size = 18),
					 axis.title.x = element_text(size = 18))
			
			print(plot)
			
		} else {
			print(paste0(i, " nao tem o mesmo numero de colunas que ", o))
		}
	}
}

# After getting the DEG across the different tumor types, the next step is to evaluate the biological pathways related to these genes in each tumor types
# Ora analysis for each tumor type
# Load packages
library(clusterProfiler)
library(circlize)
library(org.Hs.eg.db)

file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	data <- get(i)
	data <- as.data.frame(data)
	data <- data %>%
		filter(padj <= 0.05, log2FoldChange >= 1 | log2FoldChange < -1)
	data_genes <- rownames(data)
	assign(paste0("sig.", i), data_genes, envir = .GlobalEnv)
}

file_names <- paste0("deseq.", tissues)
file_names <- paste0("sig.", file_names)
for (i in file_names) {
	file <- get(i)
	genes <- enrichGO(gene = file,
				   OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
				   ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
	
	genes <- data.frame(genes)
	assign(paste0("ora", i), genes, envir = .GlobalEnv)
}

file_names <- paste0("ora", file_names)
df <- data.frame(Description = character(0), padjust = numeric(0), Tissue = character(0))
for (i in file_names) {
	file <- get(i)
	file <- file[1:3, c(3, 7)]
	file$Tissue <- rep(gsub("orasig.deseq.", "", i), nrow(file))
	if (nrow(df) == 0) {
		df <- file
	} else {
		df <- rbind(df, file)
	}
}

df <- df %>% 
	mutate(Description = gsub("extracellular structure organization", 
						 "extracellular matrix organization", Description),
		  Description = gsub("leukocyte chemotaxis", 
		  			    "leukocyte migration", Description),
		  Description = gsub("positive regulation of cell adhesion", 
		  			    "leukocyte cell-cell adhesion", Description),
		  Description = gsub("regulation of cell-cell adhesion", 
		  			    "leukocyte cell-cell adhesion", Description),
		  Description = gsub("organic acid catabolic process", 
		  			    "small molecule catabolic process", Description),
		  p.adjust = -log2(p.adjust))

# Chord diagram
chord_data <- df %>%
	group_by(Tissue, Description) %>%
	summarise(p.adjust_sum = sum(p.adjust))
chord_matrix <- chord_data %>%
	spread(key = Description, value = p.adjust_sum, fill = 0) %>%
	as.data.frame()
rownames(chord_matrix) <- chord_matrix$Tissue
chord_matrix <- chord_matrix[,-1]
chord_matrix <- as.matrix(chord_matrix)
chord_matrix <- t(chord_matrix)

# save chord
chordDiagram(chord_matrix, annotationTrack = c("name", "grid"), 
		   directional = 1, transparency = 0.4, 
		   big.gap = 50, small.gap = 1)

# Now I analyzed  the significant pathways regulated in each type of tumor. Now I will analyze the similar pathways and differently expressed genes across the different tumor types
# Load packages
library(radsets)
library(circlize)
library(gplots)
library(purrr)
library(RVenn)
library(ggplot2)

# Get sig up genes 
file_names <- paste0("deseq.", tissues)
for (i in file_names) {
	data <- get(i)
	data <- as.data.frame(data)
	data <- data %>%
		filter(padj <= 0.05, log2FoldChange >= 1)
	data_genes <- rownames(data)
	assign(paste0("up.", i), data_genes, envir = .GlobalEnv)
}

### Vennn
file_names <- paste0("up.", file_names)
new_list <- list()
for (i in file_names) {
	list <- list(get(i))
	tissue <- gsub("up.deseq.", "", i)
	names(list) <- tissue
	new_list <- append(new_list,list)
}

toy <- Venn(new_list)
setmap(toy, element_clustering = T, set_clustering = T)
sig_toy <- overlap(toy) # Overlap genes

# Heatmap with similar DEG across the tumor types
# Bind all cts files
cts <- cbind(cts_metadata.Bladder, cts_metadata.Breast, cts_metadata.Colorectal,
		   cts_metadata.Kidney, cts_metadata.Lung, 
		   cts_metadata.Pancreas, cts_metadata.Prostate, cts_metadata.Stomach, 
		   cts_metadata.Thyroid, cts_metadata.Uterus, cts_metadata.Kidney)

# Select just genes DEG in all samples
# Cts files
cts_deg <- subset(cts, rownames(cts) %in% sig_toy)
cts_deg <- log10(cts_deg)
cts_deg <- replace(cts_deg, is.infinite(cts_deg), NA)
sample_names <- unique(colnames(cts_deg))
cts_deg <- cts_deg[, colnames(cts_deg) %in% sample_names]

# Annotation with metadata 
annotation <- subset(metadata_ptafr, !(tissue == "Endometrium"))
annotation <- subset(annotation, !(type == "Control"))
annotation$file_name <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", annotation$file_name)
annotation <- annotation[!duplicated(annotation$file_name), ]
annotation <- subset(annotation, file_name %in% sample_names)
rownames(annotation) <- annotation$file_name
annotation <- annotation[,c(-1, -2, -3, -5)]

annotation <- annotation %>%
	arrange(tissue, type)
annotation <- annotation %>%
	arrange(type)

colnames(annotation) <- c("Tissue", "PTAFR")

common_colnames <- intersect(colnames(cts_deg), annotation$file_name)
cts_deg <- cts_deg[, common_colnames]
colnames(cts_deg) <- common_colnames

order_colnames <- rownames(annotation)
cts_deg <- cts_deg[, order_colnames]

# Define os nomes das colunas em cts_deg
colnames(cts_deg) <- order_colnames

# create heatmap
# annotation colors
cores_tissue <- c("#001F3F", "#003366", "#004080", "#005A8D", "#0072BB", "#0088CC", "#3299CC", "#66B3CC", "#99CCFF", "#B3D9FF")
cores_type <- c("#D3D3D3", "#A9A9A9")
cores_annotation <- list(Tissue = setNames(cores_tissue,  sort(tissues)),
					PTAFR = setNames(cores_type, c("Low", "High")))

# Delete duplicated genes
genes <- unique(rownames(cts_deg))
cts_deg <- cts_deg[genes,]


library(pheatmap)
pheatmap(cts_deg, 
	    cluster_rows = FALSE,  # Agrupar as linhas
	    cluster_cols = FALSE, scale = "row", # Agrupar as colunas
	    color = colorRampPalette(c("white", "#FBF2F2", "#AF0000"))(20),  # Esquema de cores
	    show_rownames = TRUE,  # N??o mostrar os nomes das linhas
	    show_colnames = FALSE,
	    clustering_distance_cols = "canberra",
	    fontsize_row = 4,
	    annotation_col = annotation,
	    annotation_colors = cores_annotation)

# I dindt find any downregulated DEG similar across all tumor types, so I`ll not put the code here.
# Then, I`ll compare the ORA analysis for the specific DEG in all tumor samples
# Load packages
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(org.Hs.eg.db)
library(stringr)

## ORA Analysis for all
genes <- enrichGO(gene = sig_toy,
			   OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = TRUE,
			   ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

genes@result[["Description"]] <- str_to_title(genes@result[["Description"]])

# Cnetplot DEG
# Log2FoldChange between samples
file_names <- paste0("deseq.", tissues)
fc_all <- data.frame()

for (i in file_names) {
	file <- get(i)
	if (length(fc_all) == 0) {
		fc_all <- file[, 1:2]
		colnames(fc_all) <- c("mean", gsub("deseq.", "", i))
	} else {
		colname <- gsub("deseq.", "", i) # Corrigido: extrair o nome da coluna corretamente
		fc_all[, colname] <- file$log2FoldChange
		cat(i, "adicionado\n") # Adicionado "\n" para uma nova linha na mensagem de sa??da
	}
}

fc_all <- fc_all[,-1]

fd_siggenes <- fc_all[rownames(fc_all) %in% sig_toy, ]
fd_siggenes$mean <- rowMeans(fd_siggenes)

# Tirar ORA redundantes
genes@result[["Description"]] <- gsub("Regulation Of Leukocyte Cell-Cell Adhesion", "Leukocyte Cell-Cell Adhesion", genes@result[["Description"]])
genes@result[["Description"]] <- gsub("Regulation Of Cell-Cell Adhesion", "Leukocyte Cell-Cell Adhesion", genes@result[["Description"]])


gene_list <- fd_siggenes$mean
names(gene_list) <- rownames(fd_siggenes)
cnetplot(genes, circular = FALSE,showCategory = 7, colorEdge = T, foldChange=gene_list) +
	scale_colour_gradient2(name = "Log2 Fold-Change", low = "blue", mid = "green", high = "#f33119",
					   node_label = "none")


cnetplot(genes, node_label="gene", 
	    showCategory = 7, foldChange=gene_list, colorEdge = T) +
	scale_colour_gradient2(name = "Log2 Fold-Change", low = "blue", mid = "lightgreen", high = "#f33119")

# -------------------------------------------------- Part 5: Tumor cell composition between PTAFR Low and High --------------------------------------------------------
# Metadata analysis for loot at tumor micro-environment between PTAFR low and high groups
# Load packages
library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(tidyr)
library(ggpattern)
library(ggridges)
library(DESeq2) 
library(stats)

# Load complete metadata
clinical <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/clinical_data/clinical.tsv", header = TRUE, sep = "\t")
pathology <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/clinical_data/pathology_detail.tsv", header = TRUE, sep = "\t")
slide <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/sample_data/slide.tsv", header = TRUE, sep = "\t")
sample <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/gdc_sample_sheet.2023-08-09.tsv", header = TRUE, sep = "\t")

# Load breast e bladder metadata
# Breast 
clinical_breast <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/clinical.tsv", header = TRUE, sep = "\t")
pathology_breast <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/pathology_detail.tsv", header = TRUE, sep = "\t")
slide_breast <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/slide.tsv", header = TRUE, sep = "\t")
sample_breast <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_breast/gdc_sample_sheet.2023-09-06.tsv", header = TRUE, sep = "\t")

# Bladder
clinical_bladder <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/clinical.tsv", header = TRUE, sep = "\t")
pathology_bladder <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/pathology_detail.tsv", header = TRUE, sep = "\t")
slide_bladder <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/slide.tsv", header = TRUE, sep = "\t")
sample_bladder <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/metadata_bladder/gdc_sample_sheet.2023-09-06.tsv", header = TRUE, sep = "\t")

# Merge all metadata 
clinical <- rbind(clinical, clinical_bladder, clinical_breast)
pathology <- rbind(pathology, pathology_bladder, pathology_breast)
slide <- rbind(slide, slide_bladder, slide_breast)
sample <- rbind(sample, sample_bladder, sample_breast)

# Metadata A - with immune cell data information
metadata_c <- read.delim("/home/bad23/Documents/Data_Analysis/Bioinformatic-Analysis-DSMZ/Cancer/patient_counts_processed/metadata_ptafr.csv", header = TRUE, sep = " ")
metadata_a <- merge(metadata_c, slide[,c(2,16,17,19)], by = "case_id")
pathology_a <- pathology[,c(1,21)]
clinical_a <- clinical[,c(1,27)]
slide_a <- slide[,c(2,16,17,19)]
metadata_a <- merge(metadata_c, pathology_a, by = "case_id", all = TRUE)
metadata_a <- merge(metadata_a, clinical_a, by = "case_id", all = TRUE)
metadata_a <- merge(metadata_a, slide_a, by = "case_id", all = TRUE)
metadata_a <- metadata_a[,-1]
metadata_a <- metadata_a[!duplicated(metadata_a),]
metadata_a <- metadata_a[!is.na(metadata_a$file_name), ]

metadata_a$ajcc_pathologic_n <- gsub("'--", NA, metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N0 \\(i-\\)|N0 \\(i\\+\\)|N0 \\(mol\\+\\)", "N0", metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N1a|N1b|N1c|N1mi", "N1", metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N2a|N2b", "N2", metadata_a$ajcc_pathologic_n)
metadata_a$ajcc_pathologic_n <- gsub("N3a|N3b|N3c", "N3", metadata_a$ajcc_pathologic_n)

metadata_a$lymph_nodes_positive <- gsub("'--", NA, metadata_a$lymph_nodes_positive)
metadata_a$lymph_nodes_positive <- as.numeric(metadata_a$lymph_nodes_positive)

metadata_a$percent_lymphocyte_infiltration <- gsub("'--", NA, metadata_a$percent_lymphocyte_infiltration)
metadata_a$percent_monocyte_infiltration <- gsub("'--", NA, metadata_a$percent_monocyte_infiltration)
metadata_a$percent_neutrophil_infiltration <- gsub("'--", NA, metadata_a$percent_neutrophil_infiltration)
metadata_a$percent_lymphocyte_infiltration <- as.numeric(metadata_a$percent_lymphocyte_infiltration)
metadata_a$percent_monocyte_infiltration <- as.numeric(metadata_a$percent_monocyte_infiltration)
metadata_a$percent_neutrophil_infiltration <- as.numeric(metadata_a$percent_neutrophil_infiltration)

# Create plot with mean and sd
metadata_a <- metadata_a[metadata_a$type != "Control", ]

cores <- c("Low" = "#bcd5f5", "High" = "#3b648b")
# Lymphocyte infiltration
# Lymphocyte infiltration
ggplot(data = metadata_a, aes(x = type, y = percent_lymphocyte_infiltration, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Lymphocyte Infiltration (%)",
		fill = "Type") +
	scale_fill_manual(values = cores) +
	theme_classic() +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.border = element_rect(color = "black", fill = NA, size = 1))

# T test
grupo_negative <- metadata_a$percent_lymphocyte_infiltration[metadata_a$type == "Low"]
grupo_positive <- metadata_a$percent_lymphocyte_infiltration[metadata_a$type == "High"]
t.test(grupo_negative, grupo_positive)

# Graph for monocyte infiltration
ggplot(data = metadata_a, aes(x = type, y = percent_monocyte_infiltration, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Monocyte Infiltration (%)",
		fill = "Type") +
	scale_fill_manual(values = cores) +
	theme_classic() +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.border = element_rect(color = "black", fill = NA, size = 1))

# T test 
grupo_negative <- metadata_a$percent_monocyte_infiltration[metadata_a$type == "Low"]
grupo_positive <- metadata_a$percent_monocyte_infiltration[metadata_a$type == "High"]
t.test(grupo_negative, grupo_positive)

# Graph for neutrophil infiltration
ggplot(data = metadata_a, aes(x = type, y = percent_neutrophil_infiltration, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Neutrophil Infiltration (%)",
		fill = "Type") +
	scale_fill_manual(values = cores) +
	theme_classic() +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.border = element_rect(color = "black", fill = NA, size = 1))

# T test
grupo_negative <- metadata_a$percent_neutrophil_infiltration[metadata_a$type == "Low"]
grupo_positive <- metadata_a$percent_neutrophil_infiltration[metadata_a$type == "High"]
t.test(grupo_negative, grupo_positive)

# Metadata B - with tumor cell data information
metadata_b <- merge(metadata_c, slide[,c(2,18,20,23,24,25)], by = "case_id")
pathology_b <- pathology[,c(1,48)]
clinical_b <- clinical[,c(1,26,28,29)]
metadata_b <- merge(metadata_b, pathology_b, by = "case_id", all = TRUE)
metadata_b <- merge(metadata_b, clinical_b, by = "case_id", all = TRUE)
metadata_b <- metadata_b[,-1]
metadata_b <- metadata_b[!duplicated(metadata_b$file_name),]
metadata_b <- metadata_b[!is.na(metadata_b$file_name), ]

metadata_b$percent_necrosis <- gsub("'--", NA, metadata_b$percent_necrosis)
metadata_b$percent_normal_cells <- gsub("'--", NA, metadata_b$percent_normal_cells)
metadata_b$percent_stromal_cells <- gsub("'--", NA, metadata_b$percent_stromal_cells)
metadata_b$percent_tumor_cells <- gsub("'--", NA, metadata_b$percent_tumor_cells)
metadata_b$percent_tumor_nuclei <- gsub("'--", NA, metadata_b$percent_tumor_nuclei)
metadata_b$tumor_largest_dimension_diameter <- gsub("'--", NA, metadata_b$tumor_largest_dimension_diameter)
metadata_b$ajcc_pathologic_m <- gsub("'--", NA, metadata_b$ajcc_pathologic_m)
metadata_b$ajcc_pathologic_stage <- gsub("'--", NA, metadata_b$ajcc_pathologic_stage)
metadata_b$ajcc_pathologic_t <- gsub("'--", NA, metadata_b$ajcc_pathologic_t)
metadata_b <- metadata_b[metadata_b$type != "Control", ]
metadata_b[, 6:11] <- lapply(metadata_b[, 6:11], as.numeric)

# Graph with percentage of tumor cells
cores <- c("Low" = "#DEEBF7", "High" = "#3282BD")
ggplot(data = metadata_b, aes(x = type, y = percent_tumor_cells, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Tumor cells (%)",
		fill = "Type") +
	scale_fill_manual(values = cores) +
	theme_classic() +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.border = element_rect(color = "black", fill = NA, size = 1))

# T test
grupo_negative <- metadata_b$percent_tumor_cells[metadata_b$type == "Low"]
grupo_positive <-  metadata_b$percent_tumor_cells[metadata_b$type == "High"]
t.test(grupo_negative, grupo_positive)

# Graph for percent of normal cells
ggplot(data = metadata_b, aes(x = type, y = percent_normal_cells, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Normal cells (%)",
		fill = "Type") +
	scale_fill_manual(values = cores) +
	theme_classic() +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.border = element_rect(color = "black", fill = NA, size = 1))

# Graph for percentage of stromal cells
ggplot(data = metadata_b, aes(x = type, y = percent_stromal_cells, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Stromal cells (%)",
		fill = "PTAFR") +
	scale_fill_manual(values = cores) +
	theme_classic() +
	theme(
		axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.border = element_rect(color = "black", fill = NA, size = 1))

# Graph for percent of necrosis
ggplot(data = metadata_b, aes(x = type, y = percent_necrosis, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Necrosis (%)") +
	scale_fill_manual(values = cores) +
	theme_classic()
# T test
grupo_negative <- metadata_b$percent_necrosis[metadata_b$type == "Low"]
grupo_positive <-  metadata_b$percent_necrosis[metadata_b$type == "High"]
t.test(grupo_negative, grupo_positive)

# Graph for percent of tumor nuclei
ggplot(data = metadata_b, aes(x = type, y = percent_tumor_nuclei, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.2, size = 1) +
	labs(x = "Type", y = "Necrosis (%)") +
	scale_fill_manual(values = cores) +  # Define as cores manualmente
	theme(panel.background = element_blank(),
		 plot.background = element_blank(),
		 panel.border = element_rect(color = "black", fill = NA, size = 2),
		 legend.position = "none",
		 axis.text.x = element_text(size = 10), 
		 axis.text.y = element_text(size = 10),
		 axis.title.y = element_text(size = 13),
		 axis.title.x = element_blank())
# T test
grupo_negative <- metadata_b$percent_tumor_nuclei[metadata_b$type == "Low"]
grupo_positive <-  metadata_b$percent_tumor_nuclei[metadata_b$type == "High"]
t.test(grupo_negative, grupo_positive)

# Graph for tumor diameter
ggplot(data = metadata_b, aes(x = type, y = tumor_largest_dimension_diameter, fill = type)) +
	geom_bar(stat = "summary", fun = "mean", position = "dodge") +
	geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.1, size = 0.5) +
	labs(x = NULL, y = "Tumor diameter (cm2)") +
	scale_fill_manual(values = cores) +
	theme_classic()
# T test
grupo_negative <- metadata_b$tumor_largest_dimension_diameter[metadata_b$type == "Low"]
grupo_positive <-  metadata_b$tumor_largest_dimension_diameter[metadata_b$type == "High"]
t.test(grupo_negative, grupo_positive)

# -------------------------------------------------- Part 6: Deconvolution Analysis --------------------------------------------------------
setwd("~/Desktop/immunodeconv.linux")
# Set packages
library(dplyr)

metadata_ptafr <- read.delim("~/Desktop/immunodeconv.linux/files_processed/metadata_ptafr.csv", header = TRUE, sep = " ")

positive_tumor <- metadata_ptafr %>%
	filter(type == "High")

negative_tumor <- metadata_ptafr %>%
	filter(type == "Low")

tissues <- unique(metadata_ptafr$tissue)
file_names <- paste0("cts_metadata.", tissues, ".tsv")
positive_files <- positive_tumor$file_name
negative_files <- negative_tumor$file_name
negative_files <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", negative_files)
positive_files <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", positive_files)
positive_cts <- data.frame()
negative_cts <- data.frame()

# Export data files
for (i in file_names) {
	table <- read.delim(paste0("~/Desktop/immunodeconv.linux/files_processed/", i), header = TRUE, sep = "\t")
	colnames(table) <- gsub("X", "", colnames(table))
	colnames(table) <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", colnames(table))
	colnames(table) <- gsub("\\.", "-", colnames(table))
	selected_cols <- table %>%
		select(any_of(negative_files))
	if (nrow(negative_cts) == 0) {
		negative_cts <- selected_cols
	} else {
		negative_cts <- bind_cols(negative_cts, selected_cols)
	}
}

for (i in file_names) {
	table <- read.delim(paste0("~/Desktop/immunodeconv.linux/files_processed/", i), header = TRUE, sep = "\t")
	colnames(table) <- gsub("X", "", colnames(table))
	colnames(table) <- gsub(".rna_seq.augmented_star_gene_counts.tsv", "", colnames(table))
	colnames(table) <- gsub("\\.", "-", colnames(table))
	selected_cols <- table %>%
		select(any_of(positive_files))
	if (nrow(positive_cts) == 0) {
		positive_cts <- selected_cols
	} else {
		positive_cts <- bind_cols(positive_cts, selected_cols)
	}
}

# Organize tables (as.matrix) with gene ID as rownames and file.name as colnames
table <- read.table("~/Desktop/immunodeconv.linux/05ea028e-6b8a-4352-88aa-057494cf5e80.rna_seq.augmented_star_gene_counts.tsv", header = T, sep ="\t", skip = 1)
table <- table[,c(1:3)]
table <- subset(table, gene_type == "protein_coding")

negative_cts <- as.matrix(negative_cts)
positive_cts <- as.matrix(positive_cts)

rownames(negative_cts) <- table$gene_name
rownames(positive_cts) <- table$gene_name

negative_cell_score <- immunedeconv::deconvolute(negative_cts, "quantiseq")
negative_tumor_score <- immunedeconv::deconvolute_estimate(negative_cts)

positive_cell_score <- immunedeconv::deconvolute(positive_cts, "quantiseq")
positive_tumor_score <- immunedeconv::deconvolute_estimate(positive_cts)

# Now that I have the data, I'll process the results
matrix_negative_cell_score <- as.matrix(negative_cell_score[,2:ncol(negative_cell_score)])
rownames(matrix_negative_cell_score) <- negative_cell_score$cell_type

matrix_positive_cell_score <- as.matrix(positive_cell_score[,2:ncol(positive_cell_score)])
rownames(matrix_positive_cell_score) <- positive_cell_score$cell_type

# Statistical test
wilcox.test(matrix_negative_cell_score[1,], matrix_positive_cell_score[1,])

# Positive processing
mean_pos <- data.frame(Mean = rowMeans(matrix_positive_cell_score))
mean_pos$cell <- rownames(mean_pos)
mean_pos$PTAFR <- "High"

# Negative processing
mean_neg <- data.frame(Mean = rowMeans(matrix_negative_cell_score))
mean_neg$cell <- rownames(mean_neg)
mean_neg$PTAFR <- "Low"

mean_cell <- rbind(mean_neg, mean_pos)
rownames(mean_cell) <- NULL

# Take off uncharacterized cell
mean_cell <- mean_cell %>%
	filter(cell != "uncharacterized cell")

mean_cell$cell_type <- c("lymphoid", "myeloid", "myeloid", "myeloid", 
					"myeloid", "lymphoid", "lymphoid",
					"lymphoid", "lymphoid", "myeloid")

mean_cell <- mean_cell %>%
	arrange(cell_type, cell)

mean_cell_lymphoid <- mean_cell[1:10,]
mean_cell_myeloid <- mean_cell[11:20,]

write.table(mean_cell_lymphoid, file = "~/Desktop/immunodeconv.linux/lymphoid_cells.tsv", sep = "\t", row.names = FALSE)
write.table(mean_cell_myeloid, file = "~/Desktop/immunodeconv.linux/myeloid_cells.tsv", sep = "\t", row.names = FALSE)

# Graph
library(ggplot2)
ggplot(mean_cell_myeloid, aes(x = PTAFR, y = cell, size = Mean, color = Mean)) +
	geom_point(alpha = 0.7) +  
	scale_size_continuous(range = c(1, 10)) +
	scale_color_gradient(
		low = "blue",
		high = "red"
	) +
	labs(
		x = NULL,
		size = "Cell proportion (%)",
		color = "p.value"
	) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "white"),
		panel.border = element_rect(color = "black", fill = NA, size = 1))


# Opcao de deconvolucao 2
negative_cell_score <- immunedeconv::deconvolute(negative_cts, "mcp_counter")
cells <- negative_cell_score$cell_type
negative_cell_score <- as.matrix(negative_cell_score[,2:ncol(negative_cell_score)])
rownames(negative_cell_score) <- cells

positive_cell_score <- immunedeconv::deconvolute(positive_cts, "mcp_counter")
positive_cell_score <- as.matrix(positive_cell_score[,2:ncol(positive_cell_score)])
rownames(positive_cell_score) <- cells

tumor_composition <- data.frame(Positive = rowMeans(positive_cell_score),
						  Negative = rowMeans(negative_cell_score))

write.table(tumor_composition, file = "~/Desktop/immunodeconv.linux/tumor_cells.tsv", sep = "\t", row.names = FALSE)

