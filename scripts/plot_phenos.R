# ============================================================================
# Script: Phenotypic trait analysis and visualization
#
# Description:
#   Analyzes and visualizes two phenotypic traits (development time and
#   fecundity) across different mitochondrial backgrounds (Bei, Zim, Yak)
#   and population cages. Performs mixed-model ANOVA and generates estimated
#   marginal means tables and boxplots with individual data points.
#
# Author: Leah Darwin
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(lme4)
library(lmerTest)
library(knitr)
library(kableExtra)
library(emmeans)

# Color palette for mitotypes
color_palette <- c(
  "Bei"   = "#1C448E",
  "Zim"  = "#52C2BA",
  "Yak"  = "#FCAB10"
)

mts = names(color_palette)

# Load and prepare development time data
dev= read.csv("data/PupalDev_2.0.csv") %>%
  na.omit() %>%
  mutate(MitoCage = paste0(Mito,Cage)) %>%
  uncount(Count)

dev$Mito = factor(dev$Mito, levels=mts)

# Load and prepare fecundity data
fec = read.csv("data/Feucndity2.0.csv") %>%
  mutate(MitoCage = paste0(Mito,Cage))

fec$Mito = factor(fec$Mito, levels=mts)

# Fit mixed-effects models with mitotype as fixed effect and cage/replicate as random effects
dev_lm = lmer(Time ~ Mito + (1|MitoCage/Replicate), data=dev)
dev_avg = dev %>%
  summarise(Time = mean(Time), .by=c(Mito,Cage,Replicate,MitoCage))

fec_lm = lmer(Count ~ Mito + (1|MitoCage), data = fec)

# Calculate estimated marginal means for each mitotype
dev_emmeans = emmeans(dev_lm, ~Mito)
fec_emmeans = emmeans(fec_lm, ~Mito)

anova_table = bind_rows(
  as.data.frame(anova(dev_lm)) %>% mutate(Trait = "Development Time"),
  as.data.frame(anova(fec_lm)) %>% mutate(Trait = "Fecundity")
) %>%
  select(Trait, everything()) %>%
  rename(`Sum Sq` = `Sum Sq`, `Mean Sq` = `Mean Sq`,
         `Num DF` = NumDF, `Den DF` = DenDF,
         `F` = `F value`, `p` = `Pr(>F)`)

kable(anova_table, format = "latex", booktabs = TRUE, digits = 3,
      caption = "Mixed model ANOVA results for development time and fecundity.") %>%
  collapse_rows(columns = 1, latex_hline = "major")

# Combined estimated marginal means table
emmeans_table = bind_rows(
  as.data.frame(dev_emmeans) %>%
    select(Mito, emmean, SE, lower.CL, upper.CL) %>%
    mutate(Trait = "Development Time"),
  as.data.frame(fec_emmeans) %>%
    select(Mito, emmean, SE, lower.CL, upper.CL) %>%
    mutate(Trait = "Fecundity")
) %>%
  select(Trait, everything()) %>%
  rename(`Mitotype` = Mito, `Estimate` = emmean, `SE` = SE,
         `Lower 95% CI` = lower.CL, `Upper 95% CI` = upper.CL)

emmeans_kable = kable(emmeans_table, format = "latex", booktabs = TRUE, digits = 3,
                      caption = "Estimated marginal means by trait and mitotype.") %>%
  collapse_rows(columns = 1, latex_hline = "major")

emmeans_kable

# Save to tex file
writeLines(as.character(emmeans_kable), "output/emmeans_table.tex") 

# Create boxplot with overlaid points for development time
p1 = ggplot(dev_avg, aes(x = MitoCage, y = Time, fill=Mito)) +
  geom_boxplot(alpha=0.6) +
  geom_point(aes(color = Mito),
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  facet_wrap(~Mito, scales="free_x") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_linedraw() +
  labs(x = "Population", y = "Development Time (hours)") +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.grid = element_line(color = "gray90", linewidth = 0.5))

# Create boxplot with overlaid points for fecundity
p2 = ggplot(fec, aes(x = MitoCage, y = Count, fill=Mito))+
  geom_boxplot(alpha=0.6) +
  geom_point(aes(color = Mito),
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
  facet_wrap(~Mito, scales="free_x") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_linedraw() +
  labs(x = "Population", y = "Female Fecundity (# offspring)") +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.grid = element_line(color = "gray90", linewidth = 0.5))

# Combine plots and save
final = p2 + p1 + plot_layout(widths = c(2,1), guides="collect") + plot_annotation(tag_levels = "a")
final
ggsave("output/pheno_plot.pdf", final, width = 7, height = 3)
