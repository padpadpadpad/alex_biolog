# look at metadata across treatments

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(flextable)
library(webshot)
library(rmarkdown)
library(officer)


# figure path
path_fig <- 'plots'

# load data
d <- read.csv('sequencing/data/metadata.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names()

d_sub <- filter(d, !treatment %in% c('wt_ancestor', 'nmc_t0', 'negative_control')) %>%
  mutate(., log_density = log10(density_cfu_g)) %>%
  mutate_all(., function(x)ifelse(x == 0|is.infinite(x)|is.nan(x), NA, x))

# plot
ggplot(d_sub, aes(treatment, log_density)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, position  = position_jitter(width = 0.1), size = 3) +
  theme_bw() +
  ggtitle('Abundance of P. fluorescens at the end of the experiment') +
  xlab('Treatment') +
  ylab('Density (CFU / g)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_fig, 'pseudomonas_abundance.png'), last_plot(), height = 6, width = 8)

mod1 <- lm(log10(density_cfu_g) ~ treatment, d_sub)

# pairwise of pre_adaptation treatments
d_preadapt <- filter(d_sub, !treatment %in% c('4_unrelated_clones', 'lacz_ancestor')) %>%
  select(log_density, evolution, n_clones) %>%
  filter(!is.na(log_density))

ggplot(d_preadapt, aes(evolution, log_density)) +
  geom_point() +
  facet_wrap(~ n_clones)

mod_preadapt <- d_preadapt %>%
  nest(data = c(log_density, evolution)) %>%
  mutate(., model = purrr::map(data, ~lm(log_density ~ evolution, .x)),
         tidy_model = purrr::map(model, broom::tidy)) %>%
  select(-c(data, model)) %>%
  unnest(tidy_model) %>%
  filter(term == 'evolutionwithout_community') %>%
  mutate(., p.value = round(p.value, 3))

mod_preadapt

# pool all these together to look at the other effect
d_nclones <- select(d_sub, log_density, treatment, n_clones) %>%
  mutate(., n_clones = ifelse(treatment == 'lacz_ancestor', 'lacz_ancestor', n_clones)) %>%
  filter(!is.na(log_density))

ggplot(d_nclones, aes(n_clones, log_density)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, size = 3)

# problem with doing unbalanced anova is usually homogeneity of variances assumption
group_by(d_nclones, n_clones) %>%
  summarise(var = var(log_density),
            sd = sd(log_density))
# even though lacZ has lower n, has higher variance than individual clone where sample size is much bigger
# so can continue with normal lm and normal anovas

d_nclones <- mutate(d_nclones, preadapt = ifelse(n_clones == 'lacz_ancestor', 'none', 'yes'))

mod_treat <- lm(log_density ~ n_clones, d_nclones)
mod_treat2 <- lm(log_density ~ 1, d_nclones)
anova(mod_treat, mod_treat2)

levels <- c(`single clone vs. LacZ ancestor` = '1 - lacz_ancestor',
            `24 clones vs. LacZ ancestor` = '24 - lacz_ancestor',
            `4 clones vs. LacZ ancestor` = '4 - lacz_ancestor',
            `24 clones vs. 4 clones` = '24 - 4',
            `single clone vs. 24 clones` = '1 - 24',
            `single clone vs. 4 clones` = '1 - 4')

contrasts = emmeans::emmeans(mod_treat, pairwise ~ n_clones) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(contrast = forcats::fct_recode(contrast, !!!levels)) %>%
  mutate_if(is.numeric, function(x)round(x, 2)) %>%
  mutate_all(as.character)

mod_overall <- lm(log_density ~ preadapt, d_nclones)

emmeans::emmeans(mod_overall, pairwise ~ preadapt)

# create plot
ggplot(d_nclones, aes(forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')), log_density, col = forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')), fill = forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')))) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab(expression(log[10]~Density~(cfu~g^-1~soil))) +
  xlab('') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones', 'negative\ncontrol')) +
  ylim(c(4, 7.5))

ggsave(file.path(path_fig, 'abundance.png'), last_plot(), height = 5, width = 6)
ggsave(file.path(path_fig, 'abundance.pdf'), last_plot(), height = 5, width = 6)

# super_fp
super_fp <- fp_text(vertical.align = "superscript", font.size = 8, font.family = 'Times')

# italics_fp
italic_fp <- fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

table <- flextable(contrasts) %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', j = 'contrast', part = 'all') %>%
  set_header_labels(t.ratio = "t ratio",
                    p.value = "p value") %>%
  add_footer_row(., colwidths = 6, values = '') %>%
  compose(., j = "df", part = "header", 
          value = as_paragraph(as_chunk("d.f.", props = italic_fp))) %>%
  compose(., j = "contrast", part = "footer", 
          value = as_paragraph('P values were adjusted using the Tukey method for comparing a family of 4 estimates')) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  padding(padding.top = 5, part = 'footer')

# save as a png ####
# create an Rmd file
rmd_name <- tempfile(fileext = ".Rmd")
cat("```{r echo=FALSE}\ntable\n```", file = rmd_name)

# render as an html file ----
html_name <- tempfile(fileext = ".html")
render(rmd_name, output_format = "html_document", output_file = html_name )

# get a png from the html file with webshot ----
webshot(html_name, zoom = 2, file = "table_s1.png", 
        selector = "body > div.container-fluid.main-container > div.tabwid > table")

