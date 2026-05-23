library(ggplot2)
library(tidyverse)
library(readr)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(gridtext)

# Figure 4d
# Similar for plots S13.1, S13.2, S13.3, S13.4
imp_var <- read_excel("final_pgx_data.xlsx")

imp_var <- imp_var |>
  mutate(populations = as.character(populations)) |>
  mutate(pop_num = suppressWarnings(as.numeric(populations))) |>
  arrange(is.na(pop_num), pop_num) |>
  mutate(populations = factor(populations, levels = c(as.character(1:82), "SAS", "AFR", "AMR", "EAS", "EUR")),
         type = factor(type, levels = c("Dosage", "Toxicity", "Efficacy", "Other")),
         variant_gene = paste0(variant, " (<i>", gene, "</i>)"),
         chemicals = str_to_title(chemicals),
         drug_category = paste0(chemicals, " (", category, ")"),
         Linguistic_group_Tribe = factor(Linguistic_group_Tribe,
                                         levels = c("AA_T", "DR_T", "DR_NT",
                                                    "IE_T", "IE_NT", "TB_T",
                                                    "NT", "Global"))) |>
  arrange(type, Testing, variant_gene) |> 
  group_by(variant_gene) |>
  mutate(norm_frequency = as.numeric(scale(frequency))) |>
  ungroup()

mat <- imp_var |>
  select(populations, variant_gene, norm_frequency) |>
  pivot_wider(names_from = populations, values_from = norm_frequency) |>
  column_to_rownames("variant_gene") |>
  as.matrix()

testing_colors <- c(
  "Testing Required" = "#B30568",
  "Testing Recommended" = "#0568B3",
  "Actionable PGx" = "#69B305"
)
linguistic_colors <- c(
  "AA_T" = "#4682b4", "DR_T" = "#000080", "DR_NT" = "#8A8AF4",
  "IE_T" = "#cd3333", "IE_NT" = "#FAC4C4", "TB_T" = "#7ccd7c",
  "NT" = "#C7F2C7", "Global" = "#000000"
)
type_colors <- c("Toxicity"="#c6caca", "Dosage"="#c6caca", "Efficacy"="#c6caca", "Other"="#c6caca")

min_val <- min(imp_var$norm_frequency, na.rm = TRUE)
max_val <- max(imp_var$norm_frequency, na.rm = TRUE)

col_fun <- colorRamp2(
  c(min_val, -2.5, 0, 2.5, 5, max_val),
  c("#2E627AFF", "#1283c0ff", "#F5F5F5FF", "#D6604DFF", "#B2182BFF", "#67001FFF")
)

row_meta <- imp_var %>%
  distinct(variant_gene, drug_category, type, Testing)
row_meta <- row_meta[match(rownames(mat), row_meta$variant_gene), ]

col_meta <- imp_var %>%
  distinct(populations, Linguistic_group_Tribe)
col_meta <- col_meta[match(colnames(mat), col_meta$populations), ]

col_ha <- HeatmapAnnotation(
  Linguistic_Group = col_meta$Linguistic_group_Tribe,
  col = list(Linguistic_Group = linguistic_colors),
  spacer = anno_empty(height = unit(1, "mm"), border = FALSE),
  annotation_name_gp = gpar(col = NA),
  show_legend = FALSE
)

right_ha <- rowAnnotation(
  Drug = anno_text(row_meta$drug_category,
                   gp = gpar(fontsize = 28, fontface = "italic"),
                   just = "left",
                   location = 0.01)
)

left_ha <- rowAnnotation(
  Type = row_meta$type,
  Variant = anno_text(gt_render(rownames(mat)), 
                      gp = gpar(fontsize = 32, fontface = "bold",
                                col = testing_colors[row_meta$Testing]),
                      just = "right",
                      location = 0.99),
  col = list(Type = type_colors),
  gap = unit(4, "mm"), 
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(
  mat,
  name = "Z-score Frequency",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_split = row_meta$type,
  row_title_gp = gpar(fontsize = 30, fontface = "bold"),
  column_split = col_meta$Linguistic_group_Tribe,
  column_title_gp = gpar(fontsize = 28, fontface = "bold"),
  column_names_gp = gpar(fontsize = 30, fontface = "bold"),
  row_gap = unit(3, "mm"),
  column_gap = unit(3, "mm"),
  top_annotation = col_ha,
  right_annotation = right_ha,
  left_annotation = left_ha,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  show_heatmap_legend = FALSE
)

lg_z <- Legend(
  col_fun = col_fun,
  title = "Z-score Frequency",
  at = seq(floor(min_val), ceiling(max_val), by = 2),
  labels = seq(floor(min_val), ceiling(max_val), by = 2),
  direction = "horizontal",
  title_gp = gpar(fontsize = 28, fontface = "bold"),
  labels_gp = gpar(fontsize = 22),
  legend_width = unit(16, "cm"),
  grid_height = unit(8, "mm"),
  title_gap = unit(5, "mm")
)

lg_testing <- Legend(
  labels = names(testing_colors),
  title = "Testing",
  legend_gp = gpar(fill = testing_colors),
  direction = "horizontal",
  title_gp = gpar(fontface = "bold", fontsize = 28),
  labels_gp = gpar(fontsize = 26),
  title_gap = unit(5, "mm"),
  row_gap = unit(2, "mm")
)

combined_legends <- packLegend(
  lg_z, lg_testing,
  direction = "horizontal",
  gap = unit(6, "cm")
)

png("final_heatmap.png", width = 60, height = 15, units = "in", res = 600)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = FALSE,
  heatmap_legend_list = list(combined_legends)
)

dev.off()

# Figure 4e, 4f
cyp_genes <- c("cyp2b6", "cyp2c9", "cyp3a5", "cyp2c19", "cyp4f2", "cyp2d6")
non_cyp_genes <- c("cftr", "slco1b1", "dpyd", "tpmt", "ifnl3", "ugt1a1", "nudt15", "vkorc1")

pop_info <- read_excel("final_pop_info.xlsx") %>%
  select(FID, Linguistic_group_Tribe) %>% 
  distinct()

cyp2d6 <- read_excel("cyp2d6_merged_genotypes.xlsx") %>%
  rename(name = FID, phenotype = Metabolizer_Status) %>% 
  select(name, Linguistic_group_Tribe, phenotype) %>% 
  mutate(gene = "cyp2d6")

cyrius_samples <- cyp2d6$name

merged_rows_df <- data.frame()

for (gene in c(cyp_genes, non_cyp_genes)) {
  star_allele <- file.path("stargazer/", paste0("output-", gene, ".sg-genotype.txt"))
  
  if (!file.exists(star_allele)) {
    message("Gene not genotyped: ", gene)
    next
  }
  
  df <- read_tsv(star_allele, col_types = cols(.default = "c")) %>%
    select(name, phenotype) %>% 
    mutate(gene = gene)
  
  if(gene == "cyp2d6") {
  df <- df %>%
    filter(!name %in% cyrius_samples)
  }
  
  merged_df <- merge(df, pop_info, by.x = "name", by.y = "FID") %>%
    filter(Linguistic_group_Tribe != "")
  
  merged_rows_df <- bind_rows(merged_rows_df, merged_df)
}

final_df <- bind_rows(merged_rows_df, cyp2d6)

pop_gene_phenotype <- final_df %>%
  group_by(Linguistic_group_Tribe, gene, phenotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Linguistic_group_Tribe, gene) %>%
  mutate(freq = (count / sum(count)) * 100)

pop_gene_phenotype <- pop_gene_phenotype %>%
  mutate(phenotype = recode(phenotype,
                            "normal_metabolizer" = "Normal",
                            "Normal Metabolizer" = "Normal",
                            "normal_function" = "Normal",
                            "favorable_response" = "Normal",
                            "intermediate_metabolizer" = "Intermediate",
                            "Intermediate Metabolizer" = "Intermediate",
                            "poor_metabolizer" = "Poor",
                            "Poor Metabolizer" = "Poor",
                            "poor_function" = "Poor",
                            "decreased_function" = "Poor",
                            "rapid_metabolizer" = "Rapid",
                            "increased_function" = "Rapid",
                            "ultrarapid_metabolizer" = "Rapid",
                            "Ultrarapid Metabolizer" = "Rapid",
                            "unknown_metabolizer" = "Unknown",
                            "Unknown Metabolizer" = "Unknown",
                            "unknown_function" = "Unknown",
                            "Indeterminate" = "Unknown"),
         Linguistic_group_Tribe = recode(Linguistic_group_Tribe, 
                            "AA_Tribal" = "AA_T",
                            "DR_Non_Tribal" = "DR_NT",
                            "DR_Tribal" = "DR_T",
                            "IE_Non_Tribal" = "IE_NT",
                            "IE_Tribal" = "IE_T",
                            "TB_Non_Tribal" = "TB_NT",
                            "TB_Tribal" = "TB_T",
                            "CAO" = "CAO"),
          Linguistic_group_Tribe = factor(Linguistic_group_Tribe, levels = rev(c("AA_T", "DR_NT", "DR_T", "IE_NT",
                                                                            "IE_T", "TB_NT", "TB_T", "CAO"))),
         gene = toupper(as.character(gene)))

gene_order <- c(
  "CYP2B6",
  "CYP2C9",
  "CYP3A5",
  "CYP2C19",
  "CYP4F2",
  "CYP2D6",
  "CFTR",
  "SLCO1B1",
  "DPYD",
  "TPMT",
  "IFNL3",
  "UGT1A1",
  "NUDT15",
  "VKORC1"
)

pop_gene_phenotype$gene <- factor(
  pop_gene_phenotype$gene,
  levels = gene_order
)

phenotype_colors <- c(
  "Normal" = "#1f77b4",
  "Intermediate" = "#ff7f0e",
  "Poor" = "#d62728",
  "Rapid" = "#2ca02c",
  "Unknown" = "grey70"
)

p <- ggplot(pop_gene_phenotype, aes(x = freq, y = forcats::fct_rev(gene), fill = phenotype)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = phenotype_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Proportion", y = NULL, fill = "Metabolizer") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold.italic", size = 6),
    axis.title.x = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 6),
    panel.spacing.x = unit(0.08, "cm"),
    panel.spacing.y = unit(-0.02, "cm")
  )

ggsave(
  filename = "gene_matrix_main.png",
  plot = q,
  device = "png",
  width = 3,
  height = 2.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

pop_gene_phenotype <- pop_gene_phenotype %>% filter(gene %in% c("CYP2B6", "CYP2C19", "CYP2D6", "CYP3A5", "CYP4F2", "SLCO1B1", "UGT1A1", "VKORC1"))

q <- ggplot(pop_gene_phenotype, aes(x = freq, y = Linguistic_group_Tribe, fill = phenotype)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phenotype_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Proportion", y = NULL, fill = "Metabolizer") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 6),
    axis.title.x = element_text(size = 8),
    strip.background = element_rect(fill = "grey85", color = "grey50"),
    strip.text.x = element_text(face = "bold.italic", size = 6),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 6),
    panel.spacing.x = unit(0.08, "cm"),
    panel.spacing.y = unit(-0.02, "cm")
  )

ggsave(
  filename = "pop_gene_matrix_main.png",
  plot = r,
  device = "png",
  width = 9,
  height = 2.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Figure S13.4
actionable <- read_excel("actionable_per_sample.xlsx", sheet = "Total_Actionable")
actionable <- actionable |> filter(Population != "Siddi")

# Part A

(each_sample <- ggplot(actionable, aes(x = as.factor(total_actionable))) +
  geom_bar(fill = "#8da0cb", color = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 5) +
  labs(x = "Number of actionable variants", y = "Count of individuals") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold")))

ggsave("actionable_variants.jpeg",
       plot = each_sample, device = "jpeg", height = 5, width = 6, units = "in", dpi = 600)

# Part B

(each_group <- ggplot(actionable, aes(x = as.factor(total_actionable), fill = Linguistic_tribal)) +
  geom_bar(color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            vjust = -0.5, size = 3) +
  labs(x = "Number of actionable variants", 
       y = "Count of individuals") +
  facet_wrap(~ Linguistic_tribal, nrow = 1) +
  scale_fill_manual(values = c("#4682b4", "#8A8AF4", "#000080", "#FAC4C4", "#cd3333", "#C7F2C7", "#7ccd7c")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 10, face = "bold", color = "white"),
    strip.background = element_rect(fill = "#4d4d4d", color = NA),
    legend.position = "none"
  ))

ggsave("actionable_variants_by_ling_tribal.jpeg",
       plot = each_group, device = "jpeg", height = 7, width = 15, units = "in", dpi = 600)
