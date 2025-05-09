
#Generating all plots for Optimising detection of viral outbreaks via wastewater surveillance, Bellekom et al

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(gridExtra)



#creating plot for sample vs log copy number per L
ggplot(metagen, aes(x=sample_comb1, y=log10_quantity_L, fill=sample_full_name)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',stackratio=1.5, dotsize=1.2, position = position_dodge(0.2)) +
  xlab("Wastewater Batch") + ylab("Log10 Covid abundance per L") +
  labs(fill = "Sample")


# plot for band presence vs copies per L

metaplot<-ggplot(metagen, aes(x=Band, y=log10_quantity_L, fill=sample_full_name)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',stackratio=1.5, dotsize=1.2, position = position_dodge(0.2)) +
  xlab("PCR band presence") + ylab("Log10 Covid abundance per L") +
  labs(fill = "Sample")

# Plot the binomial band presence vs copies per L
bionomialpcrband<-ggplot(metagen, mapping=aes(x=log10_quantity_L, y=Band_num)) + 
  stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial))+
  geom_point(aes(color = sample_full_name))+ 
  xlab("Log10 Covid abundance per L")+
  ylab("PCR band presence")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size=13),
        axis.ticks.length = unit(0.25, "cm") ) 




# Creating stacked bar charts for read distribution, averaged across replicates within an aliquot


# Create a compiled species column by coalescing taxonomy levels
barcodedata <- barcodedata %>%
  mutate(compiled_species = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain))

# Calculate total reads for each species across all samples
total_reads <- barcodedata %>%
  rowwise() %>%
  mutate(Total_Reads = sum(c_across(starts_with("Sample")))) %>%
  ungroup()

# Identify the top 20 species based on total read counts
top_species <- total_reads %>%
  arrange(desc(Total_Reads)) %>%
  slice(1:20)

# Mark other species as "Other"
barcodedata <- barcodedata %>%
  mutate(compiled_species = if_else(compiled_species %in% top_species$compiled_species, compiled_species, "Other"))

# Reshape data and average replicates
plot_data <- barcodedata %>%
  pivot_longer(cols = starts_with("Sample"), names_to = "Sample", values_to = "Reads") %>%
  mutate(Sample_Group = case_when(
    grepl("Sample 1", Sample) ~ "Sample 1",
    grepl("Sample 2", Sample) ~ "Sample 2",
    grepl("Sample 3", Sample) ~ "Sample 3"
  )) %>%
  group_by(Sample_Group, compiled_species) %>%
  summarise(Average_Reads = mean(Reads), .groups = "drop") %>%
  group_by(Sample_Group) %>%
  mutate(Total_Reads_Per_Sample = sum(Average_Reads),
         Percentage = Average_Reads / Total_Reads_Per_Sample * 100) %>%
  ungroup()

# Plot the averaged stacked bar chart
ggplot(plot_data, aes(x = Sample_Group, y = Percentage, fill = compiled_species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis_d(option = "magma", name = "Top 20 contributors") +
  labs(x = "Sample", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.spacing.x = unit(0.1, "cm"),
    legend.position = "right",
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(0),
    legend.key.height = unit(0.8, "lines"),
    legend.key.width = unit(0.8, "lines")
  )


grid.arrange(, nrow = 2)




# Creating stacked bar charts for read distribution for all wastewater batches across all aliquots and replicates 

# coalesce the taxnomic data to make a single column with the best classification for that row, replace NA in species with genus etc
barcodedata <- barcodedata %>%
  mutate(compiled_species = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain))


# Calculate total reads for each species across all sample
total_reads <- barcodedata %>%
  rowwise() %>%
  mutate(Total_Reads = sum(c_across(starts_with("Sample")))) %>%
  ungroup()

# Identify the top 10 species based on global total read counts
top_species <- total_reads %>%
  arrange(desc(Total_Reads)) %>%
  slice(1:20)

# Mark the remaining species as 'Other'
barcodedata <- barcodedata %>%
  mutate(compiled_species = if_else(compiled_species %in% top_species$compiled_species, compiled_species, "Other"))



# Reshape data for plotting
plot_data <- barcodedata %>%
  pivot_longer(cols = starts_with("Sample"), names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(Total_Reads_Per_Barcode = sum(Reads),
         Percentage = Reads / Total_Reads_Per_Barcode * 100) %>%
  ungroup()

# Create a stacked bar chart with species contributions as percentages
ggplot(plot_data, aes(x = Sample, y = Percentage, fill = compiled_species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis_d(option = "magma", name = "Top 20 contributors") +  # Use the viridis color palette and set legend title
  labs(x = "Sample", y = "Percentage of Total Reads") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.5, "lines"),        # Adjust the size of the legend key
    legend.text = element_text(size = 8),        # Adjust the text size in the legend
    legend.title = element_text(size = 10),      # Adjust the title size in the legend
    legend.spacing.x = unit(0.1, "cm"),          # Adjust the horizontal spacing between legend key and text
    legend.position = "right",                   # Position the legend to the right
    legend.box.spacing = unit(0.1, "cm"),        # Adjust spacing within the legend box
    legend.margin = margin(0),                   # Remove margin around the legend
    legend.key.height = unit(0.8, "lines"),      # Adjust the height of the legend keys
    legend.key.width = unit(0.8, "lines"),       # Adjust the width of the legend keys
    #legend.box.background = element_rect(colour = "grey50", size = 0.5)  # Add a border to the legend box
  )



# Create heat maps of read distributions by Taxon for each wastewater batch, by alqituot and replicate 

# Aggregate data by Family, summing read counts for each sample
family_data <- barcodedata %>%
  group_by(Domain_all) %>%
  summarize(across(starts_with("Sample"), sum, .names = "Total_{col}")) %>%
  ungroup()

# Reshape data for plotting
plot_data <- family_data %>%
  pivot_longer(cols = starts_with("Total_Sample"), names_to = "Sample", values_to = "Reads") %>%
  mutate(Sample = sub("Total_", "", Sample)) %>%  # Remove "Total_" prefix for cleaner sample names
  group_by(Sample) %>%
  mutate(Total_Reads_Per_Sample = sum(Reads),
         Percentage = Reads / Total_Reads_Per_Sample * 100) %>%
  ungroup()

# Plot heatmap of read percentages by domain
ggplot(plot_data, aes(x = Sample, y = Domain_all, fill = Percentage)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", name = "Read Percentage") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 10)
  )






# Generate stacked barchats with spiked SARS-CoV-2 by sample and replicate, highlighting SARs-CoV-2 reads in red


# Coalesce taxonomy columns into one
barcodedata <- barcodedata %>%
  mutate(compiled_species = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain))

# Calculate total reads per species
total_reads <- barcodedata %>%
  rowwise() %>%
  mutate(Total_Reads = sum(c_across(starts_with("Sample")))) %>%
  ungroup()

# Identify top 20 species
top_species <- total_reads %>%
  arrange(desc(Total_Reads)) %>%
  slice(1:20)

# Ensure SARS-CoV-2 is included
top_species_names <- unique(c(top_species$compiled_species, "SARS-CoV-2"))

# Replace species not in top 20 + SARS-CoV-2 with 'Other'
barcodedata <- barcodedata %>%
  mutate(compiled_species = if_else(compiled_species %in% top_species_names, compiled_species, "Other"))

# Reshape for plotting
plot_data <- barcodedata %>%
  pivot_longer(cols = starts_with("Sample"), names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(
    Total_Reads_Per_Barcode = sum(Reads),
    Percentage = Reads / Total_Reads_Per_Barcode * 100
  ) %>%
  ungroup()

# Set order for x-axis
plot_data <- plot_data %>%
  mutate(Sample = factor(Sample, levels = c(
    "Sample Neat replica 1", "Sample Neat replica 2", "Sample Neat replica 3", 
    "Sample -1 replica 1", "Sample -1 replica 2", "Sample -1 replica 3", 
    "Sample -2 replica 1", "Sample -2 replica 2", "Sample -2 replica 3", 
    "Sample -3 replica 1", "Sample -3 replica 2", "Sample -3 replica 3", 
    "Sample -4 replica 1", "Sample -4 replica 2", "Sample -4 replica 3"
  )),
  Sample_num = as.numeric(Sample))


library(scales)
default_palette <- viridis_pal(option = "magma")(length(unique(plot_data$compiled_species)))

# Get species levels in the order ggplot would assign
species_levels <- levels(factor(plot_data$compiled_species))

# Map colors, replacing SARS-CoV-2 with red
custom_colors <- setNames(default_palette, species_levels)
custom_colors["SARS-CoV-2"] <- "red"

# Plot
plot1 <- ggplot(plot_data, aes(x = Sample, y = Percentage, fill = compiled_species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors, name = "Top 20 contributors") +
  labs(x = "Sample", y = "Percentage of Total Reads") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 10),
    legend.spacing.x = unit(0.1, "cm"),
    legend.position = "right",
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(0),
    legend.key.height = unit(0.8, "lines"),
    legend.key.width = unit(0.8, "lines")
  )

plot2 <- ggplot(plot_data, aes(x = Sample, y = Percentage, fill = compiled_species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors, name = "Top 20 contributors") +
  labs(x = "Sample", y = "Percentage of Total Reads") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 10),
    legend.spacing.x = unit(0.1, "cm"),
    legend.position = "right",
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(0),
    legend.key.height = unit(0.8, "lines"),
    legend.key.width = unit(0.8, "lines")
  )



grid.arrange(plot1, plot2, nrow = 2)





# Create NDMS figure and calculate stress SMART-9N community composition data by wastewater batch and sample


#NMDS to show dissimilarity between batches
# Calculate Bray-Curtis dissimilarity matrix
dissimilarity_matrix <- vegdist(adonismerged [,c(4:16300)], method = "bray")

nmds_result <- metaMDS(dissimilarity_matrix, k = 2, trymax = 100)

# Check stress value (indicates goodness of fit; below 0.2 is ideal)
print(nmds_result$stress)

# Extract NMDS coordinates and add group information
nmds_data <- as.data.frame(nmds_result$points)
nmds_data$Batch <- factor(rep(c("Batch 4", "Batch 3", "Batch 2", "Batch 1"), each = 9))

# Calculate Convex Hulls for each group
convex_hulls <- do.call(rbind, lapply(levels(nmds_data$Batch), function(group) {
  group_data <- subset(nmds_data, Batch == group)
  chull_indices <- chull(group_data$MDS1, group_data$MDS2)  # Calculate convex hull indices
  hull_points <- group_data[chull_indices, ]  # Get points forming the hull
  hull_points$Group <- group  # Add group information to the hull points
  return(hull_points)  # Return the convex hull points
}))

# Plot NMDS result with convex hulls and colored areas for each group

stress_value <- nmds_result$stress

ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Batch)) +
  geom_point(size = 3) +  # Add points for each sample
  geom_polygon(data = convex_hulls, aes(x = MDS1, y = MDS2, fill = Batch, group = Batch), 
               color = "black", alpha = 0.3) +  # Add colored areas (convex hulls)
  labs(x = "NMDS Dimension 1", y = "NMDS Dimension 2") +
  annotate("text", x = Inf, y = Inf, label = paste("Stress = ", round(stress_value, 3)), 
           hjust = 1, vjust = 1, size = 5, color = "black") +  # Add stress value to the plot
  theme_minimal()





##################

# Combined figure for the log10 maximum coverage and log 10 percentage of reads that are SARS-CoV-2 for enriched/non-enriched samples

enrichedboxcov<-ggplot(meta_enriched_comp, aes(x = as.factor(Enriched), 
                                               y = log10(Max_Coverage + ifelse(Max_Coverage == 0, runif(nrow(meta_enriched_comp), 1e-4, 1e-3), 0)))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(aes(color = Sample_type), width = 0.2, size = 2, alpha = 0.8) +  ype
  scale_x_discrete(labels = c("0" = "Non-enriched", "1" = "Enriched")) +
  labs(x = "", y = "Log10 maximum Coverage", color = "Sample Type") +
  theme_minimal()

enrichboxperc<-ggplot(meta_enriched_comp, aes(x = as.factor(Enriched), 
                                              y = log10(Percent_covid + ifelse(Percent_covid == 0, 
                                                                               runif(nrow(meta_enriched_comp), 1e-4, 1e-3), 0)))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(aes(color = Sample_type), width = 0.2, size = 2, alpha = 0.8) +  
  scale_x_discrete(labels = c("0" = "Non-enriched", "1" = "Enriched")) +
  labs(x = "", y = "Log10 percent SARS-CoV-2", color = "Sample Type") +
  theme_minimal()


combined_plot <- (enrichedboxcov | enrichboxperc) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")  
    )
  )
