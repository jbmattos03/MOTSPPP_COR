# Mixed-integer linear programming model and algorithms for the multi-objective traveling salesman problem with profits and passengers 

This repository contains the implementation, benchmark instances, and experimental data for the paper:  

> **"Mixed-integer linear programming model and algorithms for the multi-objective traveling salesman problem with profits and passengers"**  
> *Submitted to Expert System With Applications*  

## 📄 Abstract  

Ridesharing systems have emerged as effective mechanisms to enhance urban mobility by reducing costs, congestion, and environmental impacts. Despite the growing interest in optimization models for ridesharing, most of them are nonlinear and single-objective, focusing on minimizing travel costs while overlooking the conflicting relationships among cost, time, and bonus collection. This paper investigates the Multi-objective Traveling Salesman Problem with Profits and Passengers (MoTSPPP), an optimization problem that integrates passenger constraints with bonus collection mechanisms, capturing the trade-offs between minimizing route cost and travel time while maximizing collected bonuses. A multi-objective mixed-integer linear programming model is proposed, along with a proof of NP-hardness. Nine algorithms are investigated: an exact solver, three constructive heuristics, and five state-of-the-art metaheuristics (NSGA-II, MOEA/D, IBEA, SPEA2, MPLS). A comprehensive experimental study is conducted on 252 instances, comprising symmetric and asymmetric graphs with varying edge-weight correlations. Performance is assessed regarding processing time, solution quality, and solution diversity. Results supported by statistical tests confirm the computational difficulty of the MoTSPPP and demonstrate that metaheuristic approaches yield significantly better results than constructive heuristics.


## 📂 Repository Structure  
```
.
├── /Code/               # Source code (C/C++)
├── /Instances/          # Benchmark datasets
├── /results/            # Raw experimental and processed results
```

---
