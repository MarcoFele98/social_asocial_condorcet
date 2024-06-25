#________________________________________________________________________________________________________________________#
#___ Wisdom of the crowds, social and asocial decision-making, Roger's paradox and evolutionary social dilemmas _________#
#___ https://iopscience.iop.org/article/10.1088/1742-5468/2006/03/P03013 ________________________________________________#
#___ 11/2023 - Marco Fele _______________________________________________________________________________________________#
#________________________________________________________________________________________________________________________#

library(mosaic) # for findZeros
library(tidyverse)
library(pracma) # for erfc
library(numDeriv) # for finding derivatives

rm(list = ls())
load("social and asocial condorcet/data.RData")

# Condorcet jury theorem ____________________________________________________________________________________________________________________________________----
## Functions ----
calculate_jury_correct <- function(group_size, prob_individual_correct, quorum) {
  dbinom((ceiling((group_size + 1) * quorum)):group_size,
         size = group_size,
         prob = prob_individual_correct) |>
    sum()
}

## Results ----
condorcet_data <- expand.grid(jury_size = seq(1, 100, by = 6),
                              prob_individual_correct = seq(0, 1, by = 0.01)) |>
  mutate(prob_jury_correct = map2(jury_size, 
                                  prob_individual_correct,
                                  .f = calculate_jury_correct,
                                  quorum = 1/2) |> 
           unlist())

ggplot(condorcet_data) +
  geom_line(aes(prob_individual_correct, prob_jury_correct, 
                color = as.factor(jury_size))) +
  scale_color_discrete(name = "Jury size") +
  ylab("Probabaility jury (group) correct") +
  xlab("Probability judge (individual) correct") +
  ggtitle("Condorcet jury theorem")

ggsave("social and asocial condorcet/figures/condorcet.png",
       height = 5,
       width = 7,
       bg = "white")

# Condorcet jury theorem with social information ___________________________________________________________________________________________________________________________________-----
## Functions ----
prob_social_correct <- function(previous_q, # average probability of being correct in previous time step
                                eta, # proportion of social learners
                                prob_asocial_is_correct, # probability of being correct for asocial learners
                                group_size, 
                                quorum) {
  # calculate probability the members of the group seen by the focus individual are correct
  avg_prob_group_member_correct <- (1 - eta) * prob_asocial_is_correct + eta * previous_q   # this is a weighted average
  # calculate probability focus individual is correct
  calculate_jury_correct(group_size = group_size, 
                         prob_individual_correct = avg_prob_group_member_correct, 
                         quorum = quorum)
}
# this is the same but I need to name an argument "x" for finding derivative (R is shit)
prob_social_correct2 <- function(x, # average probability of being correct in previous time step
                                 eta, # proportion of social learners
                                 prob_asocial_is_correct, # probability of being correct for asocial learners
                                 group_size, 
                                 quorum) {
  # calculate probability the members of the group seen by the focus individual are correct
  avg_prob_group_member_correct <- (1 - eta) * prob_asocial_is_correct + eta * x   # this is a weighted average
  # calculate probability focus individual is correct
  calculate_jury_correct(group_size = group_size, 
                         prob_individual_correct = avg_prob_group_member_correct, 
                         quorum = quorum)
}
# vectorized
prob_social_correct_v <- Vectorize(prob_social_correct)
prob_social_correct2_v <- Vectorize(prob_social_correct)
# for evolutionary stuff later
calc_fitness_social <- function(x, expected_average_fitness, p) {
  expected_average_fitness / (expected_average_fitness * x + p * (1 - x))
}
calc_fitness_asocial <- function(x, expected_average_fitness, p) {
  p / (expected_average_fitness * x + p * (1 - x))
}

## Parameters ----
p = 0.55 # probability of being correct for asocial learners
group_size = 200 # number of individuals
neighbour_size <- 11 # number of individuals
quorum <- 1/2

# Example to visualize dynamics
eta = 0.7 # proportion of social learners
previous_q <- seq(0, 1, l = 1000) # previous state
next_q <- prob_social_correct_v(previous_q = previous_q, # average probability of being correct in previous time step
                                eta = eta, # proportion of social learners
                                prob_asocial_is_correct = p, # probability of being correct for asocial learners
                                group_size = neighbour_size, 
                                quorum = quorum)

plot(previous_q, next_q, xlim = c(0,1), ylim = c(0,1)) # this is the recurrence equation to solve
lines(previous_q, previous_q)

## Find equilibria ----
parameter_space <- seq(0, 1, l = 1000)
results <- data.frame(eta = parameter_space, equ_1 = NA, equ_2 = NA, equ_3 = NA)
for(i in 1:length(parameter_space)) { 
  # should rewrite using set() !!!! should be much faster
  eta <- parameter_space[i]
  
  print(eta)
  zeros <- findZeros(prob_social_correct(previous_q = mean_q, 
                                         eta = eta, 
                                         prob_asocial_is_correct = p, 
                                         group_size = neighbour_size,
                                         quorum = quorum) - mean_q ~ mean_q, 
                     xlim = c(0, 1))[[1]]
  results[i, "eta"] <- eta
  results[i, 1 + 1:length(zeros)] <- zeros
}

## Find stability of equilibria ----
results_l <- results |>
  pivot_longer(cols = contains("equ"),
               names_to = "equilibrium",
               values_to = "value") |>
  filter(!is.na(value) & eta < 0.995) |> # for making grad function work and be sure all equilibria have been found
  mutate(derivative = numDeriv::grad(func = prob_social_correct2_v, # also in pracma there is a grad function
                                     x = value,
                                     eta = eta, 
                                     prob_asocial_is_correct = p, 
                                     group_size = neighbour_size, 
                                     quorum = quorum),
         stable = ifelse(derivative < 1 & derivative > -1, T, F))

ggplot(results_l) + # fold bifurcation (blue-sky bifurcation)
  geom_point(aes(eta, value, color = stable), size = 4)

## Evolutionary analysis ----
results_evolution <- results_l |> 
  group_by(eta) |>
  mutate(mulitstabilty = ifelse(n() > 1, T, F),
         probability_wrong_basin = ifelse(stable == F, 
                                    erfc(sqrt(eta * group_size / 2) * (1 - 2 * value)) / 2, # a bit on error functions: https://math.stackexchange.com/questions/4537946/error-function-for-a-different-standard-deviation#:~:text=for%20a%20normally%20distributed%20random,%5B%E2%88%92x%2Cx%5D.
                                    NA),
         probability_wrong_basin = ifelse(mulitstabilty == T,
                                    unique(na.omit(probability_wrong_basin)),
                                    NA),
         expected_average_fitness = ifelse(mulitstabilty == T,
                                           probability_wrong_basin * min(value) + (1 - probability_wrong_basin) * max(value),
                                           value),
         pi = eta * expected_average_fitness + (1 - eta) * p) |>
  ungroup() |>
  mutate(eta_nash = eta[which.min(abs(expected_average_fitness - p))],
         # Add well defined fitness function (with demographic control) and selection gradient
         fitness_social = expected_average_fitness / (expected_average_fitness * eta + p * (1 - eta)),
         fitness_asocial = p / (expected_average_fitness * eta + p * (1 - eta))
         #fitness_group = fitness_social * eta + fitness_asocial * (1 - eta),
         # selection_gradient_social = numDeriv::grad(func = calc_fitness_social,
         #                                            x = eta,
         #                                            expected_average_fitness = expected_average_fitness,
         #                                            p = p),
         # selection_gradient_asocial = numDeriv::grad(func = calc_fitness_asocial,
         #                                            x = 1 - eta, # I think this is correct
         #                                            expected_average_fitness = expected_average_fitness,
         #                                            p = p)
         )

# Evolutionary analysis
ggplot(results_evolution) +
  geom_line(aes(eta, fitness_social), color = "red", linewidth = 2) +
  geom_line(aes(eta, fitness_asocial), color = "blue", linewidth = 2) +
  geom_vline(aes(xintercept = eta_nash), color = "purple", size = 2) +
  geom_hline(aes(yintercept = 1), lty = "dotted", linewidth = 2)  +
  geom_point(aes(eta_nash, 1), color = "purple", size = 8) +
  geom_text(aes(0.5, 1.1, label = "social"), color = "red") +
  geom_text(aes(0.5, 0.9, label = "asocial"), color = "blue") +
  geom_segment(aes(x = 0.5, y = 1.05, xend = 0.75, yend = 1.05),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red") +
  geom_segment(aes(x = 0.5, y = 0.95, xend = 0.75, yend = 0.95),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "blue") +
  geom_segment(aes(x = 1.1, y = 1.05, xend = 1, yend = 1.05),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red") +
  geom_segment(aes(x = 1.1, y = 0.95, xend = 1, yend = 0.95),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "blue") +
  ylab("Fitness")

ggsave("social and asocial condorcet/figures/fitness.png",
       height = 3,
       width = 5,
       bg = "white")

# ggplot(results_evolution) +
#   geom_line(aes(eta, selection_gradient_social), color = "red", linewidth = 2) +
#   geom_line(aes(eta, selection_gradient_asocial), color = "blue", linewidth = 2) +
#   geom_vline(aes(xintercept = eta_nash), color = "purple", size = 2) +
#   geom_hline(aes(yintercept = 0), lty = "dotted", linewidth = 2) +
#   geom_point(aes(eta_nash, 0), color = "purple", size = 8) +
#   ylab("Selection gradient")
# equilibria -> selection gradient == 0
# stability -> derivative of selection gradient (second derivative of fitness) evaluated at equilibrium: positive unstable, negative stable

# ggplot(results_evolution) +
#   geom_line(aes(eta, fitness_social), color = "red", linewidth = 2) +
#   geom_line(aes(eta, fitness_asocial), color = "blue", linewidth = 2) +
#   geom_line(aes(eta, selection_gradient_social), color = "red", linewidth = 2) +
#   geom_line(aes(eta, selection_gradient_asocial), color = "blue", linewidth = 2) 

# ggsave("social and asocial condorcet/figures/selection_gradient.png",
#        height = 3,
#        width = 5,
#        bg = "white")

# test <- data.frame(x = seq(0, 5, l = 100),
#                    f = sin(seq(0, 5, l = 100))) |>
#   mutate(dfdx = grad(func = sin,
#                      x = x))
# 
# test <- data.frame(x = seq(0, 1, l = 1000),
#                    f = sin(seq(0, 5, l = 100))) |>
#   mutate(dfdx = numDeriv::grad(func = calc_fitness_social,
#                                 x = eta,
#                                 expected_average_fitness = expected_average_fitness,
#                                 p = p))
# ggplot(test) +
#   geom_line(aes(x, f)) +
#   geom_line(aes(x, dfdx), color = "red") 

# Replicate figure ----
ggplot(results_evolution) +
  geom_line(aes(eta, pi), color = "orange", linewidth = 2) +
  geom_point(aes(eta, value, color = stable), size = 4) +
  geom_line(aes(eta, expected_average_fitness), size = 2) +
  geom_line(aes(eta, p), size = 2, lty = "dashed") +
  #geom_vline(aes(xintercept = eta_nash), color = "purple", size = 2) +
  geom_point(aes(eta_nash, p), color = "purple", size = 8)

ggsave("social and asocial condorcet/figures/bifurcation_social_dilemma.png",
       height = 4,
       width = 6,
       bg = "white")
