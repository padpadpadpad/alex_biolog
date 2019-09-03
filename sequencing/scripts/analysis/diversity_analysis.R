# look at alpha diversity and pielous evenness

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(gridExtra)
library(lme4)
library(patchwork)
library(flextable)
library(webshot)
library(rmarkdown)
library(officer)

# figure path
path_fig <- 'plots'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('sequencing/data/output/20171024_17:18/ps_no_NA_phyla.rds')

# replace metadata with new metadata
# when wanting to add columns to metadata, it is better to edit metadata_creation and overwrite the metadata file as then it can be overwritten in all future files
meta_new <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 28,000. Woof.

# alpha diversity estimates ####

# prune OTUs that are not present in any of the samples
ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)

# metadata ###
m <- sample_data(ps_sub) %>%
  select(., SampleID, sample_name, treatment, evolution, preadapt_pop, n_clones) %>%
  data.frame() %>%
  janitor::clean_names()

# calculate diversity measures of each sample ####
a_div <- estimate_richness(ps, measures = c('Shannon', 'Observed')) %>%
  mutate(., SampleID = row.names(.)) %>%
  mutate(., pielou = Shannon / log(Observed)) %>%
  janitor::clean_names() %>%
  merge(., m, by = 'sample_id') %>%
  mutate_at(., c('sample_id', 'treatment', 'sample_name', 'evolution'), as.character) 

a_div <- select(a_div, sample_id, treatment, evolution, preadapt_pop, observed, pielou, n_clones) %>%
  filter(! treatment %in% c('wt_ancestor', 'nmc_t0', 'negative_control')) %>%
  mutate(., n_clones = paste('C_', n_clones, sep = ''))

# first plot of observed OTUs and evenness ####
gather(a_div, 'metric', 'value', c(observed, pielou)) %>%
  ggplot(., aes(treatment, value)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  ggtitle('Diversity and evenness across treatments') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~metric, scale = 'free_y')

ggsave(file.path(path_fig, 'alpha_diversity.png'), last_plot(), height = 5, width = 10)

# look at individual clones only, does evolution with and without community impact diversity and evenness ####

# run a test on evenness and diversity across evolution lines
head(a_div)

a_div_preadapt <- filter(a_div, ! treatment %in% c('lacz_ancestor', '4_unrelated_clones')) %>%
  mutate(n_clones =forcats::fct_relevel(n_clones,
                                "C_1", "C_4", "C_24"))

facet <- c(C_1 = '(a) single clone', C_4 = '(b) 4 clones', C_24 = '(c) 24 clones')


# plot alpha diversity across pre-adaptation treatments
plot_div1 <- ggplot(a_div_preadapt, aes(evolution, observed, col = evolution, fill = evolution)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  scale_x_discrete(labels = c('pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  facet_wrap(~n_clones, labeller = labeller(n_clones = facet)) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.title.x = element_blank()) +
  ylim(c(0, 1400)) +
  ylab('Alpha diversity\n(observed ASVs)') +
  xlab('') +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

# plot evenness across pre-adaptation treatments
plot_even1 <- ggplot(a_div_preadapt, aes(evolution, pielou, col = evolution, fill = evolution)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  scale_x_discrete(labels = c('pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  facet_wrap(~n_clones, labeller = labeller(n_clones = facet)) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.title.x = element_blank()) +
  ylim(c(0, 1)) +
  ylab("Pielou's evenness") +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

# model these differences
d_mods <- group_by(a_div_preadapt, n_clones) %>%
  nest() %>%
  mutate(mod_alpha = purrr::map(data, ~lm(observed ~ evolution, .x)),
         mod_evenness = purrr::map(data, ~lm(pielou ~ evolution, .x)))

d_mod_alpha <- d_mods %>%
  unnest(mod_alpha %>% purrr::map(broom::tidy)) %>%
  filter(term == 'evolutionwithout_community') %>%
  mutate(., model = 'diversity',
         p_adj = p.adjust(p.value, method = 'fdr'))
d_mod_pielou <- d_mods %>%
  unnest(mod_evenness %>% purrr::map(broom::tidy)) %>%
  filter(term == 'evolutionwithout_community') %>%
  mutate(., model = 'evenness',
         p_adj = p.adjust(p.value, method = 'fdr'))
d_mods <- bind_rows(d_mod_alpha, d_mod_pielou)

#--------------------------------------#
# analysis over levels of diversity ####
#--------------------------------------#

a_div <- mutate(a_div, n_clones = ifelse(treatment == 'lacz_ancestor', 'lacz_ancestor', n_clones))

# plot
plot_div2 <- ggplot(a_div, aes(forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')), observed, col = forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')), fill = forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')))) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = 'none',
        plot.title = element_text(size = 14)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones')) +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, 1400)) +
  ylab('Alpha diversity\n(observed ASVs)') +
  xlab(c('')) +
  ggtitle('(d)')

plot_even2 <- ggplot(a_div, aes(forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')), pielou, col = forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')), fill = forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')))) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = 'none',
        plot.title = element_text(size = 14)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones')) +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, 1)) +
  ylab("Pielou's evenness") +
  ggtitle('(d)') +
  xlab('')

plot_div <- plot_div1 + plot_div2 + plot_layout(ncol = 1, heights = c(0.4, 0.6))
plot_even <- plot_even1 + plot_even2 + plot_layout(ncol = 1, heights = c(0.4, 0.6))

ggsave(file.path(path_fig, 'alpha_diversity.png'), plot_div, height = 8, width = 9)
ggsave(file.path(path_fig, 'pielous_evenness.png'), plot_even, height = 8, width = 9)
ggsave(file.path(path_fig, 'alpha_diversity.pdf'), plot_div, height = 8, width = 9)
ggsave(file.path(path_fig, 'pielous_evenness.pdf'), plot_even, height = 8, width = 9)

mod_even <- lm(pielou ~ n_clones, a_div)
mod_even2 <- lm(pielou ~ 1, a_div)
anova(mod_even, mod_even2)
emmeans::emmeans(mod_even, pairwise ~ n_clones)

mod_div <- lm(observed ~ n_clones, a_div)
mod_div2 <- lm(observed ~ 1, a_div)
anova(mod_div, mod_div2)
emmeans::emmeans(mod_div, pairwise ~ n_clones)

# create table of all contrasts

levels <- c(`single clone vs. LacZ ancestor` = 'C_1 - lacz_ancestor',
            `24 clones vs. LacZ ancestor` = 'C_24 - lacz_ancestor',
            `4 clones vs. LacZ ancestor` = 'C_4 - lacz_ancestor',
            `24 clones vs. 4 clones` = 'C_24 - C_4',
            `single clone vs. 24 clones` = 'C_1 - C_24',
            `single clone vs. 4 clones` = 'C_1 - C_4')

contrasts_div <- emmeans::emmeans(mod_div, pairwise ~ n_clones) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(contrast = forcats::fct_recode(contrast, !!!levels),
         metric = 'Alpha diversity') %>%
  mutate_if(is.numeric, function(x)round(x, 2)) %>%
  mutate_all(as.character) %>%
  select(., metric, everything())

contrasts_even <- emmeans::emmeans(mod_even, pairwise ~ n_clones) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(contrast = forcats::fct_recode(contrast, !!!levels),
         metric = "Pielou's evenness") %>%
  mutate_if(is.numeric, function(x)round(x, 2)) %>%
  mutate_all(as.character) %>%
  select(., metric, everything())

contrasts <- bind_rows(contrasts_div, contrasts_even)

# super_fp
super_fp <- fp_text(vertical.align = "superscript", font.size = 8, font.family = 'Times')

# italics_fp
italic_fp <- fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

table <- flextable(contrasts) %>%
  merge_v(., j = 'metric') %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', j = 'contrast', part = 'all') %>%
  rotate(j = 'metric', rotation = 'lrtb', align = 'top') %>%
  align(align = 'left', j = 'metric', part = 'all') %>%
  set_header_labels(t.ratio = "t ratio",
                    p.value = "p value") %>%
  add_footer_row(., colwidths = 7, values = '') %>%
  compose(., j = "df", part = "header", 
          value = as_paragraph(as_chunk("d.f.", props = italic_fp))) %>%
  compose(., j = "metric", part = "footer", 
          value = as_paragraph('P values were adjusted using the Tukey method for comparing a family of 4 estimates')) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  padding(padding.top = 5, part = 'footer') %>%
  hline(., j = 'metric', border = fp_border(color="black", width = 2), part = "body") %>%
  hline(., i = 6, border = fp_border(color="black", width = 2), part = "body")
  
# save as a png ####
# create an Rmd file
rmd_name <- tempfile(fileext = ".Rmd")
cat("```{r echo=FALSE}\ntable\n```", file = rmd_name)

# render as an html file ----
html_name <- tempfile(fileext = ".html")
render(rmd_name, output_format = "html_document", output_file = html_name )

# get a png from the html file with webshot ----
webshot(html_name, zoom = 2, file = "table_diversity.png", 
        selector = "body > div.container-fluid.main-container > div.tabwid > table")


#--------------------------------------------------------------------#
# look at diversity of pseudomonads across samples and treatments ####
#--------------------------------------------------------------------#

SBW25 = "ACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAA"

ps2 <- subset_taxa(ps, Genus == 'Pseudomonas')
ps2 <- subset_taxa(ps2, ! rownames(tax_table(ps2)) %in% c(SBW25))
sample_sums(ps2)

a_div_pseu <- #subset_taxa(ps2, rownames(tax_table(ps)) %in% c(SBW25)) %>%
  estimate_richness(ps2, measures = c('Shannon', 'Observed')) %>%
  mutate(., SampleID = row.names(.)) %>%
  mutate(., pielou = Shannon / log(Observed)) %>%
  janitor::clean_names() %>%
  merge(., m, by = 'sample_id') %>%
  mutate_at(., c('sample_id', 'treatment', 'sample_name', 'evolution'), as.character) 

a_div_pseu <- select(a_div_pseu, sample_id, treatment, evolution, preadapt_pop, observed, pielou) %>%
  filter(treatment != 'wt_ancestor')

gather(a_div_pseu, 'metric', 'value', c(observed, pielou)) %>%
  ggplot(., aes(treatment, value)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15, height = 0), size = 3) +
  ggtitle('Diversity and evenness across treatments') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~metric, scale = 'free_y')

# calculate max proportion of ASVs per sample
d_pseu <- psmelt(ps2) %>%
  janitor::clean_names() %>%
  filter(treatment != 'wt_ancestor') %>%
  mutate(SBW25 = ifelse(otu == SBW25, 'yes', 'no')) %>%
  group_by(sample, SBW25, treatment) %>%
  summarise(., abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(., SBW25, abundance) %>%
  mutate(., prop = yes / (no + yes))

select(d_pseu, otu) %>%
  do(tibble(otu = paste('pseudomonad', 1:length(unique(.$otu)), sep = '_'), seq = unique(.$otu))) %>%
  write.csv(., 'sequencing/data/output/pseudomonads.csv', row.names = FALSE)
  
# plot
group_by(d_pseu, sample) %>%
  ggplot(., aes(treatment, prop)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  ggtitle('Proportion of SBW25 relative to the other Pseudomonads across treatments') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('proportion of reads assigned to SBW25') 

ggsave(file.path(path_fig, 'SBW25_success.pdf'), last_plot(), height = 5, width = 7)
ggsave(file.path(path_fig, 'SBW25_prop.png'), last_plot(), height = 5, width = 7)


# per OTU presence in samples
d_pres <- filter(d_pseu, abundance > 0) %>%
  group_by(otu) %>%
  summarise(., n_samples = n(),
            ave_prop = mean(prop))
