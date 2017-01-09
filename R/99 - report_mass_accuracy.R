# Function plotting mass error against mass from PQI ions output

plot_mass_ions <- function(ions_file) {
  ions <- read.ions(ions_file)

  ppm_sd <- sd(ions$Mass_error_ppm)
  ppm_mean <- mean(ions$Mass_error_ppm)

  ggplot(ions, aes(x = mz, y = Mass_error_ppm)) +
    geom_point(shape = 1) +
    geom_hline(colour = c("green", "orange", "orange", "red", "red"),
               linetype = "longdash",
               yintercept = ppm_mean + c(0, 1, -1, 2, -2)*ppm_sd)
}

# Function plotting mass error against mass from Mascot file
# plot_mass_mascot <- function(mascot_file) {
#   mascot <- read.mascot(mascot_file)
#   ggplot(mascot, aes(x = pep_exp_mz, y = 10^6*pep_delta/pep_exp_mr)) +
#     geom_point(shape = 1)
#     #scale_colour_hue(l=50) # Use a slightly darker palette than normal
# }

