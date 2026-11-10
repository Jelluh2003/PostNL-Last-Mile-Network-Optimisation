# PostNL Last-Mile Network Optimisation: Service Point & Parcel Locker Placement in Maastricht

A mixed-integer linear program (R) that redesigns PostNL's pickup network in Maastricht. It decides which service points to keep, close or open, where to place automated parcel lockers (APLs), and how much storage capacity each facility needs, at minimum total annual cost and with all demand still covered.

Group project for **Elab Case II**, Maastricht University School of Business and Economics (submitted June 2026).

**Authors (Group 13):** Jelle Hermans, Theodoros Peroulakis, Jooa Havia, Jasper Keck, Ege Demirelli, Arseniy Spisivyy

## Problem

PostNL's Maastricht network of 35 service points grew organically over decades. Each service point carries a fixed annual cost regardless of how many parcels it handles, and demand is highly uneven: the busiest location handles around 457 packages a day, the smallest around 3. Pickups are much cheaper for PostNL than home delivery, so where service points sit also drives the cost mix.

**Research question:** How can PostNL optimally locate service points and automated parcel lockers in Maastricht to minimise costs, and under what conditions are lockers cost-effective?

## Results

| | Current network | Optimised network |
|---|---|---|
| Service points | 35 | 24 |
| Parcel lockers (APLs) | 0 | 2 |
| Setup cost | €1.75M | €1.02M |
| Capacity cost | €0.30M | €0.32M |
| Delivery cost | €0.49M | €0.82M |
| **Total annual cost** | **€2.53M** | **€2.15M** |

- Saves about **€380,000 per year (~15%)** with every demand node still assigned to a pickup and a delivery facility.
- 16 existing service points are kept, 19 are closed, and 10 new sites are opened.
- Closed locations averaged 48 parcels a day versus 215 for retained ones, at €6.89 vs €1.07 fixed cost per package. Pickup rates were similar (~65%), so volume, not customer behaviour, drove the closures.
- Savings come from lower fixed costs; capacity and delivery costs rise slightly because fewer, larger facilities mean longer truck distances.
- **Sensitivity analysis:** lockers are not worth opening at a €20,000 annual cost under any scenario. At €10,000 they only pay off when adoption is high, and at €3,000-5,000 they are usually worth it. This analysis (varying APL setup cost, adoption rate, home delivery cost and demand intensity) is not included in the script.

## Method

1. **Road networks.** The OpenStreetMap graph of Maastricht is cleaned (unclassified roads are classified by maximum speed) and split into a directed **truck network** that respects one-way streets and an undirected **pedestrian network**.
2. **Demand grid.** Each CBS demographic square is split into a 10×10 grid. Counts are divided across sub-squares and rates are copied, then each sub-square is snapped to the nearest pedestrian-safe node. Shortest walking paths to the existing service points give each sub-square a catchment.
3. **Pickup model.** A distance-decay curve `f(d) = a·exp(−b·d)` is fitted to observed pickup rates per service point (population-weighted, via `optim`). A second-stage OLS regression on the residuals adds a demographic correction using catchment variables with |correlation| > 0.15.
4. **Costs.** Service point setup costs are proxied by local CBS home values, scaled so existing service points average €50,000. APLs cost a fraction of that in the script. Capacity costs €0.10 per slot per day (€36.50 a year) and delivery costs €1.50 per km.
5. **Candidate locations.** About 8,000 safe nodes are reduced to 100 demand-weighted k-means candidates (taking the highest-demand node per cluster), and the 35 existing locations are added, giving 133 candidates.
6. **Mixed-integer linear program.** Binary variables open a service point or APL at each candidate and assign demand nodes to a pickup and a delivery facility; continuous variables set capacity and overflow. Constraints:
   - a candidate is a service point or an APL, not both
   - pickups only at open facilities, deliveries only from open service points
   - at most one pickup facility per node, exactly one delivery facility per node
   - capacity covers peak-day load, with at most 2% overflow per facility and 1% city-wide

   Pickup savings enter the objective as negative delivery cost, which keeps the model linear. It is solved with GLPK.

## Running it

Requires R and these packages:

```r
install.packages(c("igraph", "readxl", "gridExtra", "dplyr", "ggplot2", "corrplot",
                   "ggnewscale", "ompr", "ompr.roi", "ROI.plugin.glpk", "tidyr",
                   "broom", "scales", "Rglpk", "Matrix"))
```

The script reads one Excel workbook, `data_Maastricht_2025.xlsx`, with the sheets `Nodes`, `Edges`, `Service Point Locations`, `CBS Squares` and `Daily Activity`. It is course/case data and is **not included in this repository**. Place it in `data/` (or update the path at the top of the script) and run the script top to bottom. Steps are numbered in the script, and the ILP is built in Step 13.

## Repository contents

| File | Description |
|---|---|
| `PostNL-Last-Mile-Delivery-Optimisation.R` | Full pipeline: road networks, demand grid, pickup model, candidate selection, MILP, results and plots |

## Limitations

The model ignores traffic and route optimisation (deliveries are point to point), competitors, parcel size and weight, and customer behaviour beyond distance and demographics. Setup costs are proxies based on home values, and APL costs and capacities come from academic sources rather than PostNL data.
