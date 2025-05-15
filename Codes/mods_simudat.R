rm(list=ls())

# Models ------------------------------------------------------------------
# Aconteceu um acidente e apaguei os modelos. Depois que lembrar como fiz o sorteio
# (e tiver paciência pra fazer de novo), adiciono-os aqui

dim(HGEI[[1]]$HHHG$simudf)
aux = HGEI[[1]]$HHHG$simudf
aa = expand.grid(unique(aux$id), unique(aux$env))

# Summary -----------------------------------------------------------------
HHHG = readRDS(file = 'saves/HHHG_res.rds')
LHHG = readRDS(file = 'saves/LHHG_res.rds')
HHLG = readRDS(file = 'saves/HHLG_res.rds')
LHLG = readRDS(file = 'saves/LHLG_res.rds')

summa = rbind(
  do.call(rbind, HHHG) |>
    dplyr::mutate(scenario = "HHHG", rep = rep(1:50, each = 13),
                  mislev = as.numeric(gsub("mislev_", "", mislev))) |> 
    pivot_longer(`EM-SVD`:Pedigree, names_to = "strategy", values_to = "corr"),
  do.call(rbind, LHHG) |>
    dplyr::mutate(scenario = "LHHG", rep = rep(1:50, each = 13),
                  mislev = as.numeric(gsub("mislev_", "", mislev))) |> 
    pivot_longer(`EM-SVD`:Pedigree, names_to = "strategy", values_to = "corr"),
  do.call(rbind, LHLG) |>
    dplyr::mutate(scenario = "LHLG", rep = rep(1:50, each = 13),
                  mislev = as.numeric(gsub("mislev_", "", mislev))) |> 
    pivot_longer(`EM-SVD`:Pedigree, names_to = "strategy", values_to = "corr"),
  do.call(rbind, HHLG) |>
    dplyr::mutate(scenario = "HHLG", rep = rep(1:50, each = 13),
                  mislev = as.numeric(gsub("mislev_", "", mislev))) |> 
    pivot_longer(`EM-SVD`:Pedigree, names_to = "strategy", values_to = "corr")
)

ggplot(data = subset(summa, subset = strategy != "Mean"), 
       aes(x = factor(mislev), y = corr, color = strategy))  +
  geom_boxplot(alpha = .8) +
  geom_jitter(alpha = .5, size = .4) +
  facet_grid(strategy~scenario) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        strip.text = element_text(size = 12)) +
  labs(y = "Correlation between predicted and true value", x = "Percentage of missing data") + 
  scale_colour_viridis_d(option = 'turbo') +
  scale_x_discrete(labels = seq(20, 80, 5))

ggsave(filename = "simu_result.pdf", device = 'pdf', height = 5, width = 7, path = "outputs")

# temp = summa[which(summa$mislev == 0.2), ]
# temp$mislev = 0; temp$corr = 0
summa = summa |>  dplyr::mutate(wrap = paste(rep, strategy, sep = "_"))

ggplot(data = subset(summa, subset = strategy != "Mean"), 
       aes(x = mislev, y = corr, color = strategy)) + 
  geom_line(aes(group = wrap), alpha = .1)  +
  geom_line(data = subset(summa, subset = strategy != "Mean") |> 
              reframe(corr = mean(corr, na.rm = TRUE), .by = c(mislev, scenario, strategy)),
            aes(x = mislev, y = corr, color = strategy, group = strategy), linewidth = 1.3, 
            linetype = "solid") + 
  facet_grid(strategy~scenario) + 
  theme_bw() + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) + 
  scale_colour_viridis_d(option = 'turbo') +
  scale_x_continuous(labels = seq(20, 80, 5), breaks = seq(0.2, 0.8, 0.05)) +
  labs(y = "Correlation between predicted and true value", x = "Percentage of missing data") 

ggsave(filename = "simu_result_line.pdf", device = 'pdf', height = 5, width = 7, path = "outputs")


# Simulation overview -----------------------------------------------------

HGEI = read.csv(file = "Output_PhenoHigh.csv")
LGEI = read.csv(file = "Output_PhenoLow.csv")

summa = rbind(HGEI, LGEI)


summa |> filter(stage == "YT4") |> 
  ggplot(aes(x = year, y = mean_gv_met_female)) + 
  # geom_boxplot()
  geom_line(aes(color = GEI)) 





summa |> filter(stage != "F1") |> 
  reframe(across(mean_gv_tpe_male:met_tpe, \(x) mean(x, na.rm = TRUE)), 
          .by = c(year, scenario, GEI, stage)) |> 
  ggplot(aes(x = as.factor(year), y = var_gv_met_male)) + 
  geom_point() +
  theme_bw() + 
  theme() + 
  facet_grid(stage~GEI, scales = "free_y")

rbind(HGEI, LGEI) |> filter(stage != "F1") |> 
  ggplot(aes(x = as.factor(year), y = mean_gv_tpe_male, fill = GEI)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme() + 
  facet_grid(stage~GEI, scales = "free_y")

rbind(HGEI, LGEI) |> filter(!stage %in% c("F1", "Parents")) |> 
  ggplot(aes(x = as.factor(year), y = mean_gv_met_male, fill = GEI)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme() + 
  facet_grid(stage~GEI, scales = "free_y")

rbind(HGEI, LGEI) |> filter(!stage %in% c("F1", "Parents")) |> 
  ggplot(aes(x = as.factor(year), y = acc_met, fill = GEI)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme() + 
  facet_grid(stage~GEI, scales = "free_y")




