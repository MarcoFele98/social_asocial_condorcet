# social_asocial_condorcet
Here is the replication of: "Phase coexistence in a forecasting game", Philippe Curty and Matteo Marsili J.Stat.Mech. (2006)

## Backgorund

This model is an expansion of Condorcet Jury theorem, a basic model to explore the wisdom of the crowds effect where a group is faced with the decision task of deciding between two discrete options. The model is expanded by simulating both social and asocial decision makers, and a typical "blue-sky" bifucation comes out when varying the proportion of social information. This type of bifurcaiton is widspread in collective-decion making, and shows that when opinions are not independent two stable-states can coexist, one where the group is very correct, and one where the group is very wrong (misinformation cascades). Furthermore, by adding evolution, the model shows how frequency-dependency theoretically predicts social learning to have no advantage compared to asocial learning (an example of Rogers' paradox). Lastly, it could be possible to consider this a (particular) example of a social dilemma, where evolution does not maximise the outcome at the level of the group (so long for capitalism and the invisible hand). The main tools used to anlayse the model belong to the qualitative analysis of discrete time one dimensional dynamical systems.

## Code
First I show Condorect Jury theorem, then I drops the assumtion of independency of opinions and include both social and asocial decision makers and show how these individuals update their opinion. In the paper there was no explict dynamical equation for tracking how belives change in the population, but the update rules can be seen as a recurrence equation which operates at the time scale of the population size, i.e, in the next time step on average everyone updated its opinion once. The equilibrium of the recurrence equation can be solved with R optimization tools that find the zeros of a function, while the stability criteria can be found calculating the first derivative of the recurrence equation at the equilibrium. I expanded the paper by doing a simple evolutionary analysis where I defined a fitenss function that includes population regulation through negative density-dependency (such that at the evolutionary equilibirium individual fitness is equal to one), calculate the selection gradient to find the equilibria, and check the second derviative of the fitness function for evolutionary stability. 
