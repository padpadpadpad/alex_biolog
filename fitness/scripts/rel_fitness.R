#-----------------------------------#
# relative fitness data analysis ####
#-----------------------------------#

# load packages ####
library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)
library(lme4)
library(MicrobioUoE)

# plot path
fig_path <- 'plots/phenotype'

# function for facet labels
label_facets <- function(string){
  len <- length(string)
  string = paste('(', letters[1:len], ') long-term ', string, sep = '')
  return(string)
}

# load in fitness assay data ####
d_fit <- readxl::read_excel('fitness/data/Div-fitness_malthusian_fitness.xlsx', range = 'A4:N52', .name_repair = 'minimal') %>%
  janitor::clean_names() %>%
  pivot_longer(-c(clone, initial_density), names_to = c('.value', 'replicate'), names_sep = '_')

# load data
d <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  filter(sample_name %in% as.character(1:48)) %>%
  select(., clone = sample_name, evolution, preadapt_pop)

d_fit <- merge(d_fit, d, by = 'clone')

# calculate malthusian params
d_fit <- mutate(d_fit, mal_evo = log(evo/initial_density),
                mal_comp = log(anc/initial_density),
                rel_fit = mal_evo/mal_comp,
                sel_rate = mal_evo - mal_comp,
                sel_rate2 = exp(sel_rate))

# remove if no evo numbers were found
d_fit <- filter(d_fit, evo > 0)

# plot these
ggplot(d_fit, aes(evolution, sel_rate)) +
  MicrobioUoE::geom_pretty_boxplot(aes(col = evolution, fill = evolution)) +
  geom_point(aes(col = evolution), fill = 'white', shape = 21, position  = position_jitter(width = 0.1), size = 4) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(size = 16, color = 'black'),
        legend.position = 'none') +
  xlab('') +
  ylab('Relative fitness') +
  scale_x_discrete(labels = c('pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black')) +
  #facet_wrap(~clone) +
  NULL

ggsave(file.path(fig_path, 'relative_fitness.pdf'), last_plot(), height = 4, width = 8)
ggsave(file.path(fig_path, 'relative_fitness.png'), last_plot(), height = 4, width = 8)

d_fit <- mutate(d_fit, pop = group_indices(d_fit, long_term_mixed, rep))

# lets do an analysis on these
d_fit_mod <- lmer(sel_rate2 ~ environment*long_term_mixed + (1|pop), d_fit)
d_fit_mod2 <- lmer(sel_rate2 ~ environment+long_term_mixed + (1|pop), d_fit)
anova(d_fit_mod, d_fit_mod2)

# look at pairwise differences
em_mod <- emmeans::emmeans(d_fit_mod, ~ long_term_mixed*environment)

contr_mat <- coef(pairs(em_mod))[, c('c.2', 'c.3', 'c.5')]

emmeans::emmeans(d_fit_mod, pairwise ~ long_term_mixed*environment)
emmeans::emmeans(d_fit_mod, pairwise ~ long_term_mixed)

emmeans::emmeans(d_fit_mod, ~ long_term_mixed*environment, contr = contr_mat, adjust = 'Holm')

# try and make output table

d_table <- emmeans::emmeans(d_fit_mod, ~ long_term_mixed*environment, contr = contr_mat, adjust = 'Hochberg')$contrasts %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(., p.value = round(p.value, 4)) %>%
  mutate_at(., c('estimate', 'SE', 'df', 't.ratio'), function(x){round(x, 2)}) %>%
  rename(`p value` = p.value, `t-ratio` = t.ratio, `d.f.` = df) %>%
  mutate(., contrast = case_when(contrast == 'c.2' ~ 'mixed,mixed - mixed,static',
                                 contrast == 'c.3' ~ 'mixed,mixed - static,static',
                                 contrast == 'c.5' ~ 'static,mixed - static,static'))

table <- gt(d_table) %>%
  cols_align('center') %>%
  tab_source_note(
    source_note = "P value adjustment: Holm-Bonferroni method for 3 tests."
  ) %>%
  tab_source_note(
    source_note = "Degrees-of-freedom method: kenward-roger."
  ) %>%
  tab_style(
    style = cells_styles(
      text_weight = "bold"),
    locations = cells_data(
      rows = `p value` < 0.05)
  ) %>%
  tab_spanner(label = 'contrast', columns = 'contrast') %>%
  cols_label(contrast = 'long-term,short-term') %>%
  gt:::as.tags.gt_tbl()

# change font
table[[1]]$children[[1]] <- gsub(
  "font-family: [[:print:]]*\n",
  "font-famuly: 'Times New Roman';\n",
  table[[1]]$children[[1]]
)

html_print(table)


