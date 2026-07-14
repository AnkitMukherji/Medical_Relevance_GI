library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(gridtext)
library(readr)
library(readxl)
library(scales)
library(tidyverse)

################################ Figure 4d ################################
# Similar for plots S13.1, S13.2, S13.3, S13.4

imp_var <- read_excel("final_pgx_data.xlsx")

imp_var <- imp_var |>
  mutate(populations = as.character(populations)) |>
  mutate(pop_num = suppressWarnings(as.numeric(populations))) |>
  arrange(is.na(pop_num), pop_num) |>
  mutate(populations = factor(populations, levels = c(as.character(1:82), "SAS", "AFR", "AMR", "EAS", "EUR")),
         type = factor(type, levels = c("Dosage", "Toxicity", "Efficacy", "Other")),
         Testing = factor(Testing, levels = c("Testing Required", "Testing Recommended", "Actionable PGx")),
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

test_symbols <- c(
  "Testing Required" = "<span style='font-size:50pt;'>\u25CF</span>",
  "Testing Recommended" = "<span style='font-size:35pt;'>\u25B2</span>",
  "Actionable PGx" = "<span style='font-size:40pt;'>\u25A0</span>"
)

variant_labels <- paste0(test_symbols[row_meta$Testing], " ", rownames(mat))

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
  
  Variant = anno_text(
    gt_render(variant_labels),
    gp = gpar(
      fontsize = 32,
      fontface = "bold",
      col = "black"
    ),
    just = "right",
    location = 0.99
  ),
  
  col = list(Type = type_colors),
  
  gap = unit(8, "mm"),
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
  row_title_gp = gpar(fontsize = 28, fontface = "bold"),
  column_split = col_meta$Linguistic_group_Tribe,
  column_title_gp = gpar(fontsize = 32, fontface = "bold"),
  column_names_gp = gpar(fontsize = 32, fontface = "bold"),
  row_gap = unit(4, "mm"),
  column_gap = unit(4, "mm"),
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
  title_gp = gpar(fontsize = 32, fontface = "bold"),
  labels_gp = gpar(fontsize = 28),
  legend_width = unit(16, "cm"),
  grid_height = unit(8, "mm"),
  title_gap = unit(5, "mm")
)

lg_testing <- Legend(
  labels = c(
    " Testing Required",
    " Testing Recommended",
    " Actionable PGx"
  ),

  title = "Testing",
  type = "points",
  pch = c(16, 17, 15),
  size = unit(8, "mm"),
  legend_gp = gpar(
    col = "black"
  ),
  title_gp = gpar(
    fontface = "bold",
    fontsize = 32
  ),
  labels_gp = gpar(
    fontsize = 32
  ),
  direction = "horizontal",
  title_gap = unit(5, "mm"),
  row_gap = unit(2, "mm")
)

combined_legends <- packLegend(
  lg_testing, lg_z, 
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

################################ Figure 4e, 4f ################################

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

################################ Figure S13.5 ################################

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

################################ Extended Data Figure 8c ################################

actionable <- read_excel("actionable_per_sample.xlsx", sheet = "Total_Actionable")

actionable <- actionable %>%
  mutate(
    Group = recode(
      Linguistic_tribal,
      "AA_Tribal"      = "AA_T",
      "DR_Tribal"      = "DR_T",
      "DR_Non_Tribal"  = "DR_NT",
      "IE_Tribal"      = "IE_T",
      "IE_Non_Tribal"  = "IE_NT",
      "TB_Tribal"      = "TB_T",
      "TB_Non_Tribal"  = "TB_NT",
      "CAO"            = "CAO"
    )
  )

overall_gi <- actionable %>% mutate(Group = "Overall_GI")

GI_data_combined <- bind_rows(
  overall_gi[, c("total_actionable", "Group")],
  actionable[, c("total_actionable", "Group")]
)

label_df <- actionable %>% count(Group)

group_levels <- c("Overall_GI", "AA_T", "DR_NT", "DR_T", "IE_NT", "IE_T", "TB_NT", "TB_T", "CAO")

group_labels <- setNames(
  paste0(label_df$Group, "\n(n=", label_df$n, ")"),
  label_df$Group
)

group_labels["Overall_GI"] <- paste0(
  "Overall_GI\n(n=", nrow(actionable), ")"
)

GI_data_combined$Group <- factor(
  GI_data_combined$Group,
  levels = rev(group_levels)
)

group_colors <- c(
  "Overall_GI" = "orange",
  "AA_T" = "#4682b4",
  "DR_T" = "#000080",
  "DR_NT" = "#8A8AF4",
  "IE_T" = "#cd3333",
  "IE_NT" = "#FAC4C4",
  "TB_T" = "#7ccd7c",
  "TB_NT" = "#C7F2C7",
  "CAO" = "grey60"
)

group_medians <- GI_data_combined %>%
  group_by(Group) %>%
  summarise(
    median_actionable = median(total_actionable),
    .groups = "drop"
  )

overall_median <- median(actionable$total_actionable)

ridge_plot <- ggplot(GI_data_combined, aes(x = total_actionable, y = Group, fill = Group)) +
  geom_density_ridges(
    stat = "binline",
    binwidth = 1,
    scale = 0.95,
    alpha = 0.85,
    colour = "black",
    linewidth = 0.4,
    draw_baseline = FALSE
  ) +
  geom_vline(
    xintercept = overall_median,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_point(
    data = group_medians,
    aes(
      x = median_actionable,
      y = Group
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "black",
    size = 2.5,
    stroke = 0.8
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(labels = group_labels, expand = expansion(add = c(0.1, 0.1))) +
  scale_x_continuous(
    breaks = 0:9,
    limits = c(-0.5, 9.5),
    expand = c(0, 0)
  ) +
  labs(
    x = "Number of actionable variants",
    y = NULL
  ) +
  theme_ridges(font_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    # panel.grid = element_blank(),
    axis.line.x = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.ticks.x = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    ),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  "actionable_variants_ridge_plot_binline.png",
  plot = ridge_plot,
  device = "png",
  height = 6,
  width = 8,
  units = "in",
  dpi = 600
)

################################ Extended Data Figure 8d ################################

df <- read_excel("final_pop_info.xlsx")

pop_info <- df |>
  select(FID, Population_Order) |> 
  distinct()|>
  mutate(Population_Order = as.numeric(Population_Order))

ling_info <- df |> 
  select(Population_Order, Linguistic_group_Tribe) |> 
  distinct() |>
  mutate(
    Linguistic_group_Tribe = recode(
      Linguistic_group_Tribe,
      "IE_Non_Tribal" = "IE_NT",
      "IE_Tribal" = "IE_T",
      "DR_Non_Tribal" = "DR_NT",
      "DR_Tribal" = "DR_T",
      "AA_Tribal" = "AA_T",
      "TB_Tribal" = "TB_T",
      "TB_Non_Tribal" = "TB_NT"
    )
  )

star_allele <- file.path("stargazer/", paste0("output-", "cyp2c19", ".sg-genotype.txt"))

df <- read_tsv(star_allele, col_types = cols(.default = "c")) |>
  select(name, phenotype) |> 
  mutate(gene = "CYP2C19")

overall_poor_met_freq <- df |>
  filter(phenotype %in% c("poor_metabolizer", "poor_function", "decreased_function")) |>
  summarise(overall_freq = n() / nrow(df)) |>
  pull(overall_freq)

merged_df <- merge(df, pop_info, by.x = "name", by.y = "FID") |>
  filter(Population_Order != "")

pop_phenotype <- merged_df |>
  group_by(Population_Order, phenotype) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(Population_Order) |>
  mutate(freq = count / sum(count)) |>
  ungroup() |>
  arrange(Population_Order) |>
  left_join(ling_info, by = "Population_Order") |> 
  mutate(
    Population_Order = factor(
      Population_Order,
      levels = unique(Population_Order)
    ),
    phenotype = recode(
      phenotype,
      "normal_metabolizer" = "Normal Metabolizer",
      "normal_function" = "Normal Metabolizer",
      "favorable_response" = "Normal Metabolizer",
      "intermediate_metabolizer" = "Intermediate Metabolizer",
      "poor_metabolizer" = "Poor Metabolizer",
      "poor_function" = "Poor Metabolizer",
      "decreased_function" = "Poor Metabolizer",
      "rapid_metabolizer" = "Rapid Metabolizer",
      "increased_function" = "Rapid Metabolizer",
      "ultrarapid_metabolizer" = "Ultrarapid Metabolizer",
      "unknown_metabolizer" = "Unknown Metabolizer",
      "unknown_function" = "Unknown Metabolizer",
      "Indeterminate" = "Unknown Metabolizer"
    )
  )

phenotype_colors <- c(
  "Normal Metabolizer" = "#1f77b4",
  "Intermediate Metabolizer" = "#ff7f0e",
  "Poor Metabolizer" = "#d62728",
  "Rapid Metabolizer" = "#2ca02c",
  "Ultrarapid Metabolizer" = "#9467bd",
  "Unknown Metabolizer" = "grey70"
)

linguistic_colors <- c("IE_T"="#cd3333", "DR_T"="#000080", "AA_T" ="#4682b4", "TB_T"="#7ccd7c", 
                      "IE_NT"="#FAC4C4", "DR_NT"="#8A8AF4", "TB_NT"="#C7F2C7")

pop_phenotype <- pop_phenotype |> 
  filter(phenotype == "Poor Metabolizer")

clopidogrel_cyp2c19 <- ggplot(pop_phenotype,
       aes(x = Population_Order,
           y = freq,
           fill = Linguistic_group_Tribe)) +

  geom_bar(stat = "identity") +
  
  ## mean line
  geom_hline(yintercept = overall_poor_met_freq,
             linetype = "dotted",
             linewidth = 1) +

  ## percentage labels on bars
  geom_label_repel(
  aes(label = scales::percent(freq, accuracy = 0.1)),
  size = 5,
  fill = alpha("white", 0.8),
  label.size = 0,
  fontface = "bold",
  direction = "y",       # only move vertically
  nudge_y = 0.002,       # start just above bar
  force = 0.1,           # weak repulsion
  force_pull = 10,       # strongly pull back to bar
  box.padding = 0.05,
  point.padding = 0.05,
  min.segment.length = 0) +

  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.1))
  ) +

  scale_fill_manual(values = linguistic_colors) +

  theme_minimal() +
  labs(x = "Population", y = "Percentage of CYP2C19 Poor Metabolizers") +
  theme(
    axis.title = element_text(
      size = 17,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 16,
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_blank(),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.spacing.x = unit(0.2, "cm"),
    legend.key.width = unit(0.8, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("clopidogrel_cyp2c19_ling_trib.png", clopidogrel_cyp2c19, device = "png", height = 8, width = 27, dpi = 600)