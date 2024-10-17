Lung adenocarcinoma, a prevalent subtype of lung cancer, represents one of the most lethal human malignancies. Despite
substantial efforts to elucidate its biological underpinnings, the underlying mechanisms governing lung adenocarcinoma
remain enigmatic. Modeling and comprehending the dynamics of gene regulatory networks are crucial for unraveling the
fundamental mechanisms of lung adenocarcinoma. Conventionally, the cancer is modeled as an equilibrium process based
on a time-invariant gene regulatory network to investigate stable cell states. However, the cancer is a nonequilibrium
process and the gene regulatory network should be regarded as time-varying in actual. Therefore, a feasible framework
was developed to explore the formation and progression of lung adenocarcinoma. On one hand, to delve into the underlying
mechanisms of lung adenocarcinoma formation, the time-invariant gene regulatory network for lung adenocarcinoma was
initially undertaken, and the composition of stable cell states was elucidated based on landscape theory. Furthermore, the
plasticity of different states was quantified using energy landscape decomposition theory by incorporating cell proliferation.
And transition probabilities between different states were defined to elucidate the transition between stable cell states.
Additionally, the global sensitivity analysis was performed and a total of 3 genes and 3 regulations were identified to be
more critical for the formation lung adenocarcinoma, offering a novel strategy for designing network-based therapies for
its treatment. On the other hand, the time-invariant gene regulatory network is extended as time-varying to delve into
the underlying mechanisms of lung adenocarcinoma progression. The lung adenocarcinoma progression was characterized
as four different disease stages based on the mixed states of cell population and the evolutionary direction. And the
progressionary mechanism of transition between stages was expounded by evaluating their dynamical transport, with the
dynamical transport cost between different stages quantified using Wasserstein metrics.

To unravel the fundamental mechanisms of lung adenocarcinoma, 
a feasible framework was developed to explore the formation and progression of lung adenocarcinoma.
We model the equilibrium process based on time-invariant gene regulatory network, 
and "fmin_mean_withoutBRD.m" and "fmian_sigma_draw_phenotypic_Landscape.m" is the approximation of phenotypic 
landscape based on constrained nonlinear optimization problem.
And we model the nonequilibrium process based on time-invariant gene regulatory network, 
and "fmin_mean_withBRD.m" and "fmian_sigma_draw_plasticity_Landscape.m" is the approximation of plasticity 
landscape based on constrained nonlinear optimization problem.
