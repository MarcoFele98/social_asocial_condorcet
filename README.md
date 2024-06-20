![fitness](https://github.com/MarcoFele98/social_asocial_condorcet/assets/122376407/fa9e230a-a171-490c-8472-a45c4f40cb47)# social_asocial_condorcet
Here is the replication of: "Phase coexistence in a forecasting game", Philippe Curty and Matteo Marsili J.Stat.Mech. (2006)

## Backgorund

This model is an expansion of Condorcet Jury theorem, a basic model to explore the wisdom of the crowds effect where a group is faced with a decision task of choosing between two discrete options. The model is expanded by considering both social and asocial decision makers, and a typical "blue-sky" (or fold) bifurcation is present when varying the proportion of social learners. This type of bifurcation is widespread in collective-decion making, and shows that when opinions are not independent two stable-states can coexist, one where the group is "very correct", and one where the group is "very wrong" (misinformation cascades).Furthermore, by adding evolution, the model shows how frequency-dependency theoretically predicts social learning to have no advantage compared to asocial learning, an example of Rogers' paradox. Lastly, it could be possible to consider this a (particular) example of a social dilemma, where evolution does not maximise the outcome at the level of the group. The main tools used to analyse the model belong to the qualitative analysis of discrete time one dimensional dynamical systems.
![bifurcation_social_dilemma](https://github.com/MarcoFele98/social_asocial_condorcet/assets/122376407/d3e62455-d9e7-4721-9593-75ef25fdb1d9)


## Code
First I show Condorect Jury theorem.
![condorcet](https://github.com/MarcoFele98/social_asocial_condorcet/assets/122376407/ed0ed3c6-b8ea-496d-a6b7-726d7a099fab)

I then drop the assumtion of independency of opinions and include both social and asocial decision makers and show how these individuals update their opinion. In the paper there was no explict dynamical equation for tracking how belives change in the population, but the update rules can be seen as a recurrence equation which operates at the time scale of the population size, i.e, in the next time step on average everyone updated its opinion once. Despite the update rules are not presented in a closed form, the equilibrium can still be solved with R optimization tools that find the zeros of a function, while the stability criteria can be found calculating the first derivative of the recurrence equation at the equilibrium in respect to the frequency of social learners. I expanded the paper by doing a simple evolutionary analysis where I defined a fitness function that includes population regulation through negative density-dependency (such that at the evolutionary equilibirium individual fitness is equal to one).
![fitness](https://github.com/MarcoFele98/social_asocial_condorcet/assets/122376407/56cff993-e5b5-498a-b211-c8e1fe34eef2)

I then calculate the selection gradient to find the equilibria and check the second derviative of the fitness function for evolutionary stability.
![selection_gradient](https://github.com/MarcoFele98/social_asocial_condorcet/assets/122376407/f4e94429-9e47-4098-96f4-c5c45f7f6364)
