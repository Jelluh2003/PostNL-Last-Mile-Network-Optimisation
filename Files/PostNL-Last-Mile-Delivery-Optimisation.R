# ==============================================================================
# STEP 0: LOAD LIBRARIES & DATA
# ==============================================================================
library(igraph)
library(readxl)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(corrplot)
library(ggnewscale)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(tidyr)
library(broom)
library(scales)
library(Rglpk)
library(Matrix)

# Load the Nodes and Edges datasets
nodes <- read_excel("data_Maastricht_2025.xlsx", 
                    sheet = "Nodes")
edges <- read_excel("data_Maastricht_2025.xlsx", 
                    sheet = "Edges")

# Ensure MAX_SPEED is numeric so our mathematical comparisons work correctly
edges$MAX_SPEED <- as.numeric(edges$MAX_SPEED)

# ==============================================================================
# STEP 1: TRANSFORM UNCLASSIFIED ROADS
# ==============================================================================]

# We overwrite "unclassified" Types based on their MAX_SPEED limit
edges$TYPE <- ifelse(edges$TYPE == "unclassified" & edges$MAX_SPEED <= 15, "living_street", edges$TYPE)
edges$TYPE <- ifelse(edges$TYPE == "unclassified" & edges$MAX_SPEED == 30, "residential",   edges$TYPE)
edges$TYPE <- ifelse(edges$TYPE == "unclassified" & edges$MAX_SPEED %in% c(50, 60), "tertiary", edges$TYPE)
edges$TYPE <- ifelse(edges$TYPE == "unclassified" & edges$MAX_SPEED %in% c(70, 80), "secondary", edges$TYPE)

# ==============================================================================
# STEP 2: CREATE TRUCK AND PEDESTRIAN ALLOWANCE VARIABLES
# ==============================================================================

# Initialize both as FALSE (Default assumption is no entry unless specified)
edges$truck_allowed <- FALSE
edges$ped_allowed   <- FALSE

# --- Rule 1: Motorway & Trunk ---
# Expressways / provincial bypasses: Trucks YES, Peds NO
idx_motor <- edges$TYPE %in% c("motorway", "motorway_link", "trunk", "trunk_link")
edges$truck_allowed[idx_motor] <- TRUE

# --- Rule 2: Primary Roads ---
idx_primary <- edges$TYPE %in% c("primary", "primary_link")
edges$truck_allowed[idx_primary] <- TRUE
edges$ped_allowed[idx_primary & edges$MAX_SPEED <= 50] <- TRUE

# --- Rule 3: Secondary & Tertiary Roads ---
idx_sec_tert <- edges$TYPE %in% c("secondary", "tertiary")
edges$truck_allowed[idx_sec_tert] <- TRUE
edges$ped_allowed[idx_sec_tert & edges$MAX_SPEED <= 60] <- TRUE

# --- Rule 4: Residential, Service, & Living Streets ---
idx_res <- edges$TYPE %in% c("residential", "service", "living_street")
edges$truck_allowed[idx_res] <- TRUE
edges$ped_allowed[idx_res]   <- TRUE

# --- Rule 5: Pedestrian Paths & Steps ---
idx_ped <- edges$TYPE %in% c("pedestrian", "steps")
edges$ped_allowed[idx_ped] <- TRUE

# --- Rule 6: Transit, Construction, Services, Bridleway, We initialized everything as FALSE


# ==============================================================================
# STEP 3: CREATE THE IGRAPH NETWORK
# ==============================================================================

# 1. Force V1 and V2 to be the first two columns (igraph requirement)
edges_graph <- edges[, c("V1", "V2", names(edges)[!names(edges) %in% c("V1", "V2")])]

# 2. Duplicate only the two-way streets so cars can drive both ways
edges_twoway <- edges_graph[edges_graph$ONE_WAY == FALSE, ]
edges_twoway_rev <- edges_twoway
edges_twoway_rev$V1 <- edges_twoway$V2  # Reverse the From node
edges_twoway_rev$V2 <- edges_twoway$V1  # Reverse the To node

# Combine One-way, Two-way, and Reversed Two-way into one final list
edges_final <- rbind(edges_graph[edges_graph$ONE_WAY == TRUE, ], edges_twoway, edges_twoway_rev)

# 3. Build the specific graphs
truck_graph <- graph_from_data_frame(d = edges_final, vertices = nodes, directed = TRUE)
E(truck_graph)$weight <- ifelse(E(truck_graph)$truck_allowed, E(truck_graph)$DIST, Inf)

# Collapse the truck graph into an undirected graph for pedestrians
ped_graph <- as.undirected(truck_graph, mode = "collapse", edge.attr.comb = "first")
E(ped_graph)$weight <- ifelse(E(ped_graph)$ped_allowed, E(ped_graph)$DIST, Inf)

# Extract exact X, Y spatial coordinates from your nodes data
city_layout <- as.matrix(nodes[, c("X", "Y")])

# Setup the plotting window to show both maps side-by-side
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

# Plot 1: Truck Network (Blue with arrows)
plot(truck_graph, 
     layout = city_layout, 
     vertex.size = 0.5,        # Tiny dots for intersections
     vertex.label = NA,        # Turn off messy text labels
     edge.arrow.size = 0.1,    # Tiny arrows to show one-way directions
     edge.color = "royalblue", 
     main = "Truck Network (Directed)")

# Plot 2: Pedestrian Network (Green without arrows)
plot(ped_graph, 
     layout = city_layout, 
     vertex.size = 0.5, 
     vertex.label = NA, 
     edge.color = "forestgreen", 
     main = "Pedestrian Network (Undirected)")

# Reset plotting window back to normal
par(mfrow = c(1, 1))

# ==============================================================================
# STEP 4: SUB-SQUARE GRID CREATION & ROUTING TO CLOSEST SERVICE POINT
# ==============================================================================

# Load the data
sp <- read_excel("data_Maastricht_2025.xlsx", 
                 sheet = "Service Point Locations")

cbs <- read_excel("data_Maastricht_2025.xlsx", 
                  sheet = "CBS Squares")


# Clean CBS data: Remove missing coordinates
cbs <- cbs[!is.na(cbs$X) & !is.na(cbs$Y), ]

# Defining Save nodes
valid_ped_edges <- edges[edges$ped_allowed == TRUE, ]
valid_ped_nodes <- unique(c(valid_ped_edges$V1, valid_ped_edges$V2))
safe_nodes      <- nodes[nodes[["NODE ID"]] %in% valid_ped_nodes, ]

get_nearest_node <- function(target_x, target_y, nodes_df) {
  sq_dist <- (nodes_df$X - target_x)^2 + (nodes_df$Y - target_y)^2
  return(nodes_df[["NODE ID"]][which.min(sq_dist)])
}

sp$Node_ID <- sapply(1:nrow(sp), function(i) get_nearest_node(sp$X[i], sp$Y[i], safe_nodes))

# --- 4.1  Define exact sub-square tile dimensions ----------------------------
main_w <- 7194.917          # Full CBS square width  (coordinate units)
main_h <- 4619.165          # Full CBS square height (coordinate units)

sub_w  <- main_w / 10       # 719.4917 units — width  of one sub-square tile
sub_h  <- main_h / 10       # 461.9165 units — height of one sub-square tile

cat(sprintf(
  "Sub-square dimensions: %.4f (W) x %.4f (H) units\n", sub_w, sub_h))

# --- Define once, before the loop ---

# Count columns: split evenly across 100 sub-squares
count_cols <- c(
  "Population", "Male", "Female",
  "Age0-14", "Age15-24", "Age25-44", "Age45-64", "Age65+",
  "Households", "Single-person households",
  "Multi-person households w/o kids", "Single parent households",
  "Two-parent households", "Houses", "Vacant houses")

# Rate/average columns: same value for every sub-square (don't divide)
rate_cols <- c(
  "Home ownership %", "Rental %", "Social housing %",
  "Avg. home value k€", "OAD", "Median household income",
  "Percentage low income households", "Percentage high income households",
  "Distance nearest supermarket in km")

# --- 4.2  Generate 100 sub-square tiles for every CBS square -----------------
# We pre-allocate a list (one slot per tile) and bind once at the end;
# this is far faster than growing a data frame row-by-row inside a loop.

n_cbs      <- nrow(cbs)
sub_list   <- vector("list", n_cbs * 100L)
slot       <- 1L

for (s in seq_len(n_cbs)) {
  
  # South-west corner of the parent CBS square
  x0 <- cbs$X[s] - main_w / 2
  y0 <- cbs$Y[s] - main_h / 2
  
  for (row_i in 0:9) {        # 10 rows  (south → north)
    for (col_j in 0:9) {      # 10 cols  (west  → east )
      
      xmin <- x0 + col_j * sub_w
      xmax <- xmin + sub_w
      ymin <- y0 + row_i * sub_h
      ymax <- ymin + sub_h
      
      sub_list[[slot]] <- data.frame(
        Parent_Square = cbs$Square[s],
        Sub_ID        = paste0(cbs$Square[s], "_r", row_i, "_c", col_j),
        xmin          = xmin,
        xmax          = xmax,
        ymin          = ymin,
        ymax          = ymax,
        Center_X      = xmin + sub_w / 2,
        Center_Y      = ymin + sub_h / 2,
        
        # Count data: divide by 100 so totals are preserved when summing sub-squares
        setNames(
          as.data.frame(t(
            sapply(count_cols, function(col) {
              val <- cbs[[col]][s]
              ifelse(is.na(val), 0, val) / 100
            })
          )),
          paste0("Sub_", gsub(" ", "_", count_cols))
        ),
        
        # Rate/average data: inherited directly from the parent square
        setNames(
          as.data.frame(t(
            sapply(rate_cols, function(col) cbs[[col]][s])
          )),
          gsub(" ", "_", rate_cols)),
        
        stringsAsFactors = FALSE
      )
      slot <- slot + 1L
    }
  }
}

sub_df <- do.call(rbind, sub_list)
cat(sprintf("Total sub-squares created: %d\n", nrow(sub_df)))


# Snap each sub-square centre to the nearest pedestrian-safe node 
sub_df$Node_ID <- sapply(seq_len(nrow(sub_df)), function(i) {
  get_nearest_node(sub_df$Center_X[i], sub_df$Center_Y[i], safe_nodes)})

# Shortest-path routing: sub-square nodes → service point
# De-duplicate node IDs first — many tiles share the same nearest intersection,
# so we only compute each path once and look the result up afterwards.

unique_nodes <- unique(sub_df$Node_ID)
cat(sprintf(
  "Running shortest-path calculations for %d unique nodes...\n",
  length(unique_nodes)))

walk_dist_matrix <- distances(
  graph   = ped_graph,
  v       = as.character(unique_nodes),
  to      = as.character(sp$Node_ID),
  mode    = "all",
  weights = E(ped_graph)$weight)

# Summarise the matrix into one row per unique node
node_lookup <- data.frame(
  Node_ID          = unique_nodes,
  Walking_Distance = apply(walk_dist_matrix, 1, min),
  SP_Index         = apply(walk_dist_matrix, 1, which.min),
  stringsAsFactors = FALSE)

# Join results back onto every sub-square via its snapped node
sub_df <- merge(sub_df, node_lookup, by = "Node_ID", all.x = TRUE)

# Translate the SP column index into the human-readable Location ID
sub_df$Closest_SP_ID <- sp[["Location ID"]][sub_df$SP_Index]
sub_df$Closest_SP_ID[is.infinite(sub_df$Walking_Distance)] <- NA

# Diagnostic summary
cat("\nWalking-distance summary across all sub-squares (metres):\n")
print(summary(sub_df$Walking_Distance[!is.infinite(sub_df$Walking_Distance)]))
cat(sprintf("Unreachable sub-squares (Inf distance): %d\n",
            sum(is.infinite(sub_df$Walking_Distance))))

# Force all columns to numeric except the identifier and text columns
skip_cols <- c("Node_ID", "Parent_Square", "Sub_ID", 
               "Median_household_income", "Closest_SP_ID")

sub_df <- sub_df %>%
  mutate(across(
    !all_of(skip_cols),
    ~ suppressWarnings(replace_na(as.numeric(as.character(.x)), 0))))

# ==============================================================================
# STEP 5: HIGH-RESOLUTION CATCHMENT MAP  (geom_rect)
# ==============================================================================

# Convert IDs to factors for consistent, reproducible colour mapping
sub_df$Closest_SP_ID <- as.factor(sub_df$Closest_SP_ID)
sp$`Location ID`     <- as.factor(sp$`Location ID`)

# Only plot tiles that were successfully routed
sub_routed <- sub_df[!is.na(sub_df$Closest_SP_ID), ]

# A map of all CBS squares we have
ggplot() +
  geom_rect(
    data = cbs,
    aes(
      xmin = X - main_w / 2,
      xmax = X + main_w / 2,
      ymin = Y - main_h / 2,
      ymax = Y + main_h / 2,
      fill = ifelse(is.na(Population) | Population == 0, "Unpopulated", "Populated")),
    colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Populated" = "steelblue", "Unpopulated" = "grey80")) +
  labs(
    title    = "Maastricht: All CBS Squares (Original Size)",
    subtitle = paste0(nrow(cbs), " squares — ", 
                      sum(!is.na(cbs$Population) & cbs$Population > 0), " populated, ",
                      sum(is.na(cbs$Population)  | cbs$Population == 0), " unpopulated"),
    x = "X Coordinate", y = "Y Coordinate", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14))


catchment_map <- ggplot() +
  
  # Layer A: Filled rectangles — one tile per sub-square
  # xmin/xmax/ymin/ymax are the exact tile boundaries computed in Step 4.2
  geom_rect(
    data = sub_routed,
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = Closest_SP_ID),
    colour = NA,      # No tile borders — creates a seamless grid
    alpha  = 0.80) +
  
  # Layer B: Service Point markers drawn on top
  geom_point(
    data  = sp,
    aes(x = X, y = Y, fill = `Location ID`),
    size  = 5, shape = 21, colour = "black", stroke = 1.5) +
  
  # Layer C: Labels & theme
  labs(
    title    = "Maastricht: High-Resolution Service Point Catchment Areas",
    subtitle = paste0(
      "Each CBS square subdivided into a 10\u00d710 grid of ",
      formatC(nrow(sub_routed), format = "d", big.mark = ","),
      " tiles — coloured by nearest walking hub"),
    x    = "X Coordinate",
    y    = "Y Coordinate",
    fill = "Service Point") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 9, colour = "grey45"),
    axis.text        = element_text(size = 7),
    panel.grid.major = element_line(colour = "grey92"),
    panel.grid.minor = element_blank())

print(catchment_map)

# ==============================================================================
# STEP 6: FIT NON-LINEAR PICKUP RATE CURVE
# ==============================================================================

activity <- read_excel("data_Maastricht_2025.xlsx", 
                       sheet = "Daily Activity")

# Compute observed pickup rate per service point 
pickup_rate_obs <- activity %>%
  group_by(`Location ID`) %>%
  summarise(
    Total_Pickups    = sum(Pickups,    na.rm = TRUE),
    Total_Deliveries = sum(Deliveries, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    `Location ID`        = as.character(`Location ID`),
    Observed_Pickup_Rate = Total_Pickups / (Total_Pickups + Total_Deliveries))

print(pickup_rate_obs)

# Prepare sub-square distances and populations per SP 

sub_for_fit <- sub_df %>%
  filter(
    !is.na(Closest_SP_ID),
    Closest_SP_ID != 0,
    Sub_Population  > 0,
    Walking_Distance > 0) %>%
  mutate(Closest_SP_ID = as.character(Closest_SP_ID)) %>%
  dplyr::select(Closest_SP_ID, Sub_Population, Walking_Distance)

# Define the curve and objective function
# Decay formula:  f(d) = a * exp(-b * d)
# Boundary constraints used for initialisation:
#   f(200)  ≈ 0.80  →  a * exp(-200b)  = 0.80
#   f(2000) ≈ 0.01  →  a * exp(-2000b) = 0.01
pickup_curve <- function(d, a, b) {
  a * exp(-b * d)}

# Objective: sum of squared errors between observed and
# population-weighted predicted pickup rate across all SPs
objective <- function(params) {
  a <- params[1]
  b <- params[2]
  
  sp_ids <- unique(sub_for_fit$Closest_SP_ID)
  sse    <- 0
  
  for (sp_id in sp_ids) {
    
    obs_row <- pickup_rate_obs[pickup_rate_obs$`Location ID` == sp_id, ]
    if (nrow(obs_row) == 0) next
    
    sp_subs <- sub_for_fit[sub_for_fit$Closest_SP_ID == sp_id, ]
    
    # Apply f(d) per sub-square FIRST, then weight by population
    predicted_rate <- sum(sp_subs$Sub_Population *
                            pickup_curve(sp_subs$Walking_Distance, a, b)) /
      sum(sp_subs$Sub_Population)
    
    sse <- sse + (predicted_rate - obs_row$Observed_Pickup_Rate)^2
  }
  
  return(sse)
}

# Fit the curve 
fit <- optim(
  par    = c(a = 1.30, b = 0.00243),
  fn     = objective,
  method = "L-BFGS-B",
  lower  = c(0.01, 0.00001),
  upper  = c(2.00, 0.01000))

a_fit <- fit$par["a"]
b_fit <- fit$par["b"]

cat(sprintf("Fitted parameters : a = %.4f,  b = %.6f\n",  a_fit, b_fit))
cat(sprintf("Boundary check    : f(0) = %.1f%%,  f(200) = %.1f%%,  f(2000) = %.1f%%\n",
            100 * pickup_curve(0,    a_fit, b_fit),
            100 * pickup_curve(200,  a_fit, b_fit),
            100 * pickup_curve(2000, a_fit, b_fit)))

# Store fitted pickup rate on every sub-square 
sub_df$Pickup_Rate <- pmin(pmax(
  pickup_curve(sub_df$Walking_Distance, a_fit, b_fit), 0), 1)

# Compute predicted vs observed per SP for the plot
sp_plot_df <- sub_for_fit %>%
  group_by(Closest_SP_ID) %>%
  summarise(
    # Population-weighted mean distance: approximate x-position only for plotting
    Plot_Distance  = sum(Sub_Population * Walking_Distance) / sum(Sub_Population),
    Predicted_Rate = sum(Sub_Population * pickup_curve(Walking_Distance, a_fit, b_fit)) /
      sum(Sub_Population),
    .groups = "drop") %>%
  rename(`Location ID` = Closest_SP_ID) %>%
  inner_join(
    pickup_rate_obs[, c("Location ID", "Observed_Pickup_Rate")],
    by = "Location ID")

# Bar plot: actual pickups vs deliveries + predicted pickup rate

# Reshape to long format for stacked bars
bar_df <- pickup_rate_obs %>%
  dplyr::inner_join(
    sp_plot_df[, c("Location ID", "Predicted_Rate")], by = "Location ID") %>%
  tidyr::pivot_longer(
    cols      = c(Total_Pickups, Total_Deliveries),
    names_to  = "Type",
    values_to = "Count") %>%
  mutate(
    Type        = ifelse(Type == "Total_Pickups", "Pickups", "Deliveries"),
    `Location ID` = factor(
      `Location ID`,
      levels = pickup_rate_obs$`Location ID`[
        order(pickup_rate_obs$Total_Pickups + 
                pickup_rate_obs$Total_Deliveries, decreasing = TRUE)]))

# Compute total demand per SP for scaling the predicted rate line
totals <- pickup_rate_obs %>%
  mutate(
    Total_Demand  = Total_Pickups + Total_Deliveries,
    `Location ID` = as.character(`Location ID`)) %>%
  dplyr::inner_join(
    sp_plot_df[, c("Location ID", "Predicted_Rate")], by = "Location ID") %>%
  mutate(
    Predicted_Count = Predicted_Rate * Total_Demand,
    `Location ID`   = factor(`Location ID`, levels = levels(bar_df$`Location ID`)))

bar_plot <- ggplot() +
  
  # Stacked bars: Deliveries (bottom) + Pickups (top)
  geom_col(
    data     = bar_df,
    aes(x = `Location ID`, y = Count, fill = Type),
    position = "stack",
    width    = 0.7) +
  
  # Predicted pickup count as a horizontal line per SP
  geom_segment(
    data = totals,
    aes(
      x    = as.numeric(`Location ID`) - 0.4,
      xend = as.numeric(`Location ID`) + 0.4,
      y    = Predicted_Count,
      yend = Predicted_Count),
    colour    = "black",
    linewidth = 1.0,
    linetype  = "dashed") +
  
  # Observed pickup rate label above each bar
  geom_text(
    data  = totals,
    aes(
      x     = `Location ID`,
      y     = Total_Demand + max(totals$Total_Demand) * 0.02,
      label = paste0(round(Predicted_Rate * 100, 1), "%")),
    size    = 3,
    colour  = "grey30") +
  
  scale_fill_manual(values = c("Pickups" = "steelblue", "Deliveries" = "darkorange")) +
  
  labs(
    title    = "Actual Pickups vs Deliveries per Service Point",
    subtitle = "Dashed line = predicted pickups from fitted decay curve  |  % label = predicted pickup rate",
    x        = "Service Point ID",
    y        = "Total Volume (all days)",
    fill     = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 9, colour = "grey40"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top")

print(bar_plot)

# Print the formula aswell

ggplot() +
  stat_function(
    fun      = function(d) pickup_curve(d, a_fit, b_fit) * 100,
    xlim     = c(0, 3000),
    colour   = "steelblue",
    linewidth = 1.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = sprintf("f(d) = %.4f × e^(−%.5f × d)", a_fit, b_fit),
    x     = "Walking Distance (m)",
    y     = "Pickup Rate") +
  theme_minimal()

# ==============================================================================
# STEP 7: RESIDUAL CORRELATION ANALYSIS
# ==============================================================================

# Base model residuals per SP 
residuals_df <- sp_plot_df %>%
  mutate(Residual = Observed_Pickup_Rate - Predicted_Rate)

# Aggregate sub-squares to SP level as percentages 
# Sum raw counts across all sub-squares per SP, then convert to % of total.

sp_pct <- sub_df %>%
  filter(!is.na(Closest_SP_ID), Sub_Population > 0) %>%
  mutate(Closest_SP_ID = as.character(Closest_SP_ID)) %>%
  group_by(Closest_SP_ID) %>%
  summarise(
    # Raw count totals
    Total_Pop          = sum(Sub_Population,                      na.rm = TRUE),
    Total_Age0_14      = sum(Sub_Age0.14,                         na.rm = TRUE),
    Total_Age15_24     = sum(Sub_Age15.24,                        na.rm = TRUE),
    Total_Age25_44     = sum(Sub_Age25.44,                        na.rm = TRUE),
    Total_Age45_64     = sum(Sub_Age45.64,                        na.rm = TRUE),
    Total_Age65plus    = sum(Sub_Age65.,                          na.rm = TRUE),
    Total_HH           = sum(Sub_Households,                      na.rm = TRUE),
    Total_Single_HH    = sum(Sub_Single.person_households,        na.rm = TRUE),
    Total_Multi_HH     = sum(Sub_Multi.person_households_w.o_kids, na.rm = TRUE),
    Total_SingleParent = sum(Sub_Single_parent_households,        na.rm = TRUE),
    Total_TwoParent_HH = sum(Sub_Two.parent_households,           na.rm = TRUE),
    
    # Rate variables: population-weighted mean (already a %, don't sum)
    Avg_HomeOwnership  = weighted.mean(Home_ownership_.,                   Sub_Population, na.rm = TRUE),
    Avg_Rental         = weighted.mean(Rental_.,                           Sub_Population, na.rm = TRUE),
    Avg_SocialHousing  = weighted.mean(Social_housing_.,                   Sub_Population, na.rm = TRUE),
    Avg_HomeValue      = weighted.mean(Avg._home_value_k.,                 Sub_Population, na.rm = TRUE),
    Avg_LowIncome      = weighted.mean(Percentage_low_income_households,   Sub_Population, na.rm = TRUE),
    Avg_HighIncome     = weighted.mean(Percentage_high_income_households,  Sub_Population, na.rm = TRUE),
    Avg_DistSuper      = weighted.mean(Distance_nearest_supermarket_in_km, Sub_Population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Convert raw counts → share of total catchment population / households
  mutate(
    Pct_Age0_14       = Total_Age0_14      / Total_Pop,
    Pct_Age15_24      = Total_Age15_24     / Total_Pop,
    Pct_Age25_44      = Total_Age25_44     / Total_Pop,
    Pct_Age45_64      = Total_Age45_64     / Total_Pop,
    Pct_Age65plus     = Total_Age65plus    / Total_Pop,
    Pct_Single_HH     = Total_Single_HH   / Total_HH,
    Pct_Multi_HH      = Total_Multi_HH    / Total_HH,
    Pct_SingleParent  = Total_SingleParent / Total_HH,
    Pct_TwoParent_HH  = Total_TwoParent_HH / Total_HH
  ) %>%
  # Drop intermediate count columns — keep only final % and rate columns
  dplyr::select(Closest_SP_ID, starts_with("Pct_"), starts_with("Avg_")) %>%
  rename(`Location ID` = Closest_SP_ID)

cat(sprintf("SP composition table: %d SPs × %d variables\n",
            nrow(sp_pct), ncol(sp_pct) - 1))

# Join residuals with SP composition table
residuals_cbs <- residuals_df %>%
  inner_join(sp_pct, by = "Location ID")

# Correlate each composition variable with the base model residual
# Positive correlation → SPs with more of this group are under-predicted
# Negative correlation → SPs with more of this group are over-predicted

predictor_vars <- names(sp_pct)[names(sp_pct) != "Location ID"]

cor_plot_df <- data.frame(
  Variable    = predictor_vars,
  Correlation = sapply(predictor_vars, function(v) {
    cor(residuals_cbs$Residual, residuals_cbs[[v]], use = "complete.obs")})) %>%
  arrange(Correlation) %>%
  mutate(
    Variable  = factor(Variable, levels = Variable),
    Direction = ifelse(Correlation > 0, "Positive", "Negative"))

# plot
ggplot(cor_plot_df, aes(x = Variable, y = Correlation, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = c(-0.15, 0.15),
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey20") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "firebrick")) +
  labs(
    title    = "What explains the base model residuals?",
    subtitle = "Catchment composition (% of population) vs prediction error\nDashed lines = ±0.15 inclusion threshold for Stage 2",
    x = NULL, y = "Correlation with Residual", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "top")


# ==============================================================================
# STEP 8: TWO-STAGE MODEL  (base decay  +  demographic correction)
# ==============================================================================
# Stage 1 (Step 6): f(d) = a * exp(-b*d)  →  base predicted pickup rate
# Stage 2 (here):   residual ~ catchment demographics
# Final:            Stage1 prediction + Stage2 correction, clamped to [0,1]

# Select variables with |correlation| > 0.15
sig_vars <- cor_plot_df %>%
  filter(abs(Correlation) > 0.15) %>%
  pull(Variable) %>%
  as.character()

cat("Variables entering Stage 2 regression:\n"); print(sig_vars)

# Guard: if no variables pass the threshold, warn and stop
if (length(sig_vars) == 0) {
  stop("No variables exceeded the 0.15 correlation threshold. 
        Try lowering the filter or check your CBS data.")
}

# Fit Stage 2 OLS regression
# Predictors are standardised so coefficients are on a comparable scale.
# Response is the base model residual — not the raw pickup rate.

stage2_df <- residuals_cbs %>%
  dplyr::select(Residual, all_of(sig_vars)) %>%
  mutate(across(all_of(sig_vars), scale))

lm_fit <- lm(
  as.formula(paste("Residual ~", paste(sig_vars, collapse = " + "))),
  data = stage2_df)

cat("\nStage 2 regression summary:\n")
print(summary(lm_fit))

# Add Stage 2 correction to base predictions 
residuals_cbs$Stage2_Correction <- predict(lm_fit, newdata = stage2_df)
residuals_cbs$Pred_TwoStage     <- pmin(pmax(
  residuals_cbs$Predicted_Rate + residuals_cbs$Stage2_Correction, 0), 1)

# SSE comparison: base vs two-stage 
sse_base     <- sum((residuals_cbs$Observed_Pickup_Rate - residuals_cbs$Predicted_Rate)^2)
sse_twostage <- sum((residuals_cbs$Observed_Pickup_Rate - residuals_cbs$Pred_TwoStage)^2)

cat(sprintf("\nSSE base model       : %.6f", sse_base))
cat(sprintf("\nSSE two-stage model  : %.6f", sse_twostage))
cat(sprintf("\nImprovement          : %.1f%%\n",
            100 * (sse_base - sse_twostage) / sse_base))

#Plot: Stage 2 coefficients 
# Which demographics push predicted pickup rate up or down?

broom::tidy(lm_fit, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(Direction   = ifelse(estimate > 0, "Positive", "Negative"),
         Significant = p.value < 0.05) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, colour = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), linewidth = 0.8) +
  geom_point(aes(shape = Significant), size = 3) +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05")) +
  scale_colour_manual(values = c("Positive" = "steelblue", "Negative" = "firebrick")) +
  coord_flip() +
  labs(title    = "Stage 2: Demographic Coefficients on Residuals",
       subtitle = "Standardised predictors — filled dot = significant at p < 0.05",
       x = NULL, y = "Effect on Pickup Rate", colour = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 14),
        plot.subtitle   = element_text(size = 9, colour = "grey40"),
        legend.position = "top")

#Plot: absolute residuals side-by-side, ordered small → large
avg_base     <- mean(abs(residuals_cbs$Observed_Pickup_Rate - residuals_cbs$Predicted_Rate))
avg_twostage <- mean(abs(residuals_cbs$Observed_Pickup_Rate - residuals_cbs$Pred_TwoStage))

residuals_cbs %>%
  mutate(
    Residual_Base     = abs(Observed_Pickup_Rate - Predicted_Rate),
    Residual_TwoStage = abs(Observed_Pickup_Rate - Pred_TwoStage),
    `Location ID`     = factor(`Location ID`,
                               levels = `Location ID`[order(Residual_Base)])) %>%
  dplyr::select(`Location ID`, Residual_Base, Residual_TwoStage) %>%
  tidyr::pivot_longer(-`Location ID`, names_to = "Model", values_to = "Residual") %>%
  mutate(Model = recode(Model,
                        "Residual_Base"     = "Base model",
                        "Residual_TwoStage" = "Two-stage model")) %>%
  ggplot(aes(x = `Location ID`, y = Residual, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = avg_base,
             colour = "firebrick", linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = avg_twostage,
             colour = "steelblue",  linewidth = 0.8, linetype = "dashed") +
  scale_fill_manual(values = c("Base model"      = "firebrick",
                               "Two-stage model" = "steelblue")) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     breaks = seq(0, 1, 0.05)) +
  coord_cartesian(clip = "off") +
  labs(title    = "Absolute Residuals: Base vs Two-Stage Model",
       subtitle = "Ordered smallest → largest  |  Dashed = average absolute residual",
       x = "Service Point", y = "Absolute Residual", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title         = element_text(face = "bold", size = 14),
        plot.subtitle      = element_text(size = 9, colour = "grey40"),
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
        legend.position    = "top",
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(5, 60, 5, 5))

# ==============================================================================
# STEP 9: SIMPLIFIED THEORETICAL NODE COSTS
# ==============================================================================

# Create a dataset containing ONLY CBS squares that have valid housing data
cbs_valid <- cbs[!is.na(cbs[["Avg. home value k€"]]), ]

# Helper function to find the index of the closest VALID CBS square
get_closest_valid_idx <- function(target_x, target_y, cbs_valid_df) {
  sq_dist <- (cbs_valid_df$X - target_x)^2 + (cbs_valid_df$Y - target_y)^2
  return(which.min(sq_dist))
}

# CALCULATE THE WEIGHT FROM EXISTING SERVICE POINTS 

# Map existing Service Points to the closest CBS square that HAS data
sp$Valid_CBS_Idx <- sapply(1:nrow(sp), function(i) {
  get_closest_valid_idx(sp$X[i], sp$Y[i], cbs_valid)})

# Extract the raw housing prices for the existing layout
sp$Raw_House_Price_k <- cbs_valid[["Avg. home value k€"]][sp$Valid_CBS_Idx]

# Find the scaling multiplier so the existing service points average exactly 50,000
target_avg     <- 50000
current_sp_avg <- mean(sp$Raw_House_Price_k)
scaling_weight <- target_avg / current_sp_avg


# COMPUTE POSSIBLE COSTS FOR EVERY SINGLE NODE 

# Map ALL candidate nodes to the closest CBS square that HAS data
safe_nodes$Valid_CBS_Idx <- sapply(1:nrow(safe_nodes), function(i) {
  get_closest_valid_idx(safe_nodes$X[i], safe_nodes$Y[i], cbs_valid)})

# Extract the raw house price for each node
safe_nodes$Raw_House_Price_k <- cbs_valid[["Avg. home value k€"]][safe_nodes$Valid_CBS_Idx]

# Apply the weight to calculate the final possible setup cost for every node
safe_nodes$Node_Cost <- safe_nodes$Raw_House_Price_k * scaling_weight


cat("\n--- STEP 9: RESULTS SUMMARY ---\n")
cat(sprintf("Existing Service Points Raw Avg: %.2f k€\n", current_sp_avg))
cat(sprintf("Calculated Scaling Weight      : %.4f\n", scaling_weight))
cat(sprintf("Existing Layout Cost (Checked) : €%.2f\n", mean(sp$Raw_House_Price_k * scaling_weight)))
cat(sprintf("All Nodes Projected Avg Cost   : €%.2f\n", mean(safe_nodes$Node_Cost)))

ggplot() +
  # 1. Plot all ~8,000 candidate nodes with high transparency to handle overlapping
  geom_point(data = safe_nodes, aes(x = X, y = Y, color = Node_Cost), 
             alpha = 0.4, size = 1) +
  
  # 2. Simple, clear color gradient: Blue (Affordable) to Red (Expensive)
  scale_color_gradient(
    low = "dodgerblue3", 
    high = "firebrick", 
    labels = dollar_format(prefix = "€")) +
  
  # 3. Overlay existing service points as bold black dots with a bright green center
  geom_point(data = sp, aes(x = X, y = Y), 
             color = "black", fill = "chartreuse2", shape = 21, size = 3.5, stroke = 1.2) +
  
  # 4. Map documentation labels
  labs(
    title = "Theoretical Service Point Setup Costs in Maastricht",
    subtitle = "Nodes colored by scaled CBS housing prices (Existing SPs highlighted in green)",
    x = "X Coordinate",
    y = "Y Coordinate",
    color = "Setup Cost") +
  
  # 5. Crisp white background theme with no distracting borders or grids
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(), # Removes the background grid lines
    plot.title       = element_text(face = "bold", size = 14, margin = margin(b=4)),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b=15)),
    axis.text        = element_blank(), # Hides raw coordinate values for a clean map look
    axis.ticks       = element_blank(),
    legend.position  = "right")


# ============================================================================== 
# STEP 10: SUB-SQUARE DEMAND EXTRACTOR & PREPARATION FOR ILP
# ============================================================================== 

# Get day by day total demand and split for each sub-square

daily_demand = activity %>% 
  filter(!(`Day Index` %in% c(359,360)), `Day of Week` != "Sun") %>% 
  mutate(Total_Demand = Pickups + Deliveries, `Location ID` = as.factor(`Location ID`)) %>% 
  dplyr::select(-c(`Day of Week`, Date)) 

# Build SP-level correction table 
sp_correction = residuals_cbs %>% 
  dplyr::select(`Location ID`, Stage2_Correction) 

# Apply two-stage model to every sub-square
sub_df = sub_df %>% 
  left_join( sp_correction %>% 
               rename(Closest_SP_ID = `Location ID`), by = "Closest_SP_ID") %>% 
  mutate(Pickup_Rate_TwoStage = pmin(pmax(Pickup_Rate + coalesce(Stage2_Correction, 0), 0), 1), 
         Demand_Weight = coalesce(Sub_Population, 0) * coalesce(Pickup_Rate_TwoStage, 0)) %>% 
  group_by(Closest_SP_ID) %>% 
  mutate(
    Group_Weight = sum(Demand_Weight, na.rm = TRUE), 
    Demand_Share = case_when(is.na(Closest_SP_ID) ~ 0, Group_Weight <= 0 ~ 0, TRUE ~ Demand_Weight / Group_Weight)) %>% 
  ungroup() %>% 
  dplyr::select(-Group_Weight) 

# Now we have estimated demands for each sub square for each day 
daily_sub_demand = daily_demand %>% 
  left_join(sub_df %>% 
              dplyr::select(Sub_ID, Closest_SP_ID, Demand_Share), by = c("Location ID" = "Closest_SP_ID"),
            relationship = "many-to-many") %>%
  mutate(Demand_Share = coalesce(Demand_Share, 0),
         Sub_Pickups = Pickups * Demand_Share,
         Sub_Deliveries = Deliveries * Demand_Share, 
         Demand = Total_Demand * Demand_Share)

## If we want to compress peak days and add them up, we need to do an analysis of how peak days overlap / do not

# Identify the maximum peak day for every single individual sub-square
sub_square_peaks <- daily_sub_demand %>%
  group_by(Sub_ID) %>%
  slice_max(order_by = Demand, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Sub_ID, `Day Index`, Peak_Demand = Demand) %>%
  # CRITICAL FILTER: Drop the 8,029 sub-squares that have 0 demand all year round
  filter(Peak_Demand > 0)

# Capture the total count of active, package-generating neighborhoods
total_active_subs <- nrow(sub_square_peaks)

# Calculate how system peaks concentrate across specific dates
busy_days_concentration <- sub_square_peaks %>%
  group_by(`Day Index`) %>%
  summarise(
    Sub_Squares_Peaking_Here = n(),
    Total_Peak_Volume_Contribution = sum(Peak_Demand, na.rm = TRUE),
    # Calculate the exact proportion of the active network peaking on this date
    Percentage_Of_Active_Network = (n() / total_active_subs) * 100) %>%
  arrange(desc(Sub_Squares_Peaking_Here)) %>%
  # Create a running tally to see how few days govern the network
  mutate(Cumulative_Percentage = cumsum(Percentage_Of_Active_Network))

# This plot shows how one day dominates across almost all cells, meaning we can assume peak days overlap within our system
plot_concentration <- busy_days_concentration %>%
  slice_max(order_by = Sub_Squares_Peaking_Here, n = 10, with_ties = FALSE) %>%
  mutate(
    Is_Dominant = ifelse(`Day Index` == 336, "Extreme Dominance (Day 336)", "Secondary Overlap Days"),
    `Day Label` = reorder(factor(paste0("Day ", `Day Index`)), Sub_Squares_Peaking_Here)) %>%
  ggplot(aes(x = `Day Label`, y = Percentage_Of_Active_Network, fill = Is_Dominant)) +
  geom_col(width = 0.7, color = "black", size = 0.2) +
  # Label each bar with its exact percentage share of the network
  geom_text(aes(label = paste0(round(Percentage_Of_Active_Network, 1), "%")), 
            hjust = -0.15, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Extreme Dominance (Day 336)" = "firebrick3", "Secondary Overlap Days" = "dodgerblue4")) +
  labs(
    title = "Maastricht Network Coincidence: Proof of Peak Day Overlap",
    subtitle = "Percentage of active city sub-squares hitting their annual maximum on the exact same date",
    x = "Critical Dates (Day Index)",
    y = "% of Total Active Network Cells",
    fill = "Context Classification") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom")

# Output the plot to your workspace
print(plot_concentration)

# ILP aggregation of model parameters
ilp_demand_parameters <- daily_sub_demand %>%
  group_by(Sub_ID) %>%
  summarise(
    # Option A: Peak Demand (Best if you want to ensure lockers NEVER overflow)
    Peak_Daily_Demand = max(Demand, na.rm = TRUE),
    # Yearly Volume
    Total_Yearly_Demand = sum(Demand, na.rm = TRUE)) %>%
  ungroup()


# ==============================================================================
# STEP 11: NODE CLUSTERING — CANDIDATE LOCATION SELECTION FOR ILP
# ==============================================================================

# Helper functions

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]}

safe_weighted_mean <- function(x, w) {
  x <- to_numeric_safe(x)
  w <- to_numeric_safe(w)
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok], na.rm = TRUE)}

# Build node-level demand
node_demand <- sub_df %>%
  filter(!is.na(Node_ID)) %>%
  mutate(Node_ID = as.character(Node_ID)) %>%
  group_by(Node_ID) %>%
  summarise(
    Node_Total_Demand   = sum(Demand_Weight,  na.rm = TRUE),
    Node_Total_Pop      = sum(Sub_Population, na.rm = TRUE),
    .groups = "drop")

# Attach demand to safe_nodes

nodes_enriched <- safe_nodes %>%
  mutate(Node_ID = as.character(`NODE ID`)) %>%
  left_join(node_demand, by = "Node_ID") %>%
  mutate(
    Node_Total_Demand = coalesce(Node_Total_Demand, 0),
    Node_Total_Pop    = coalesce(Node_Total_Pop,    0))

# Build and scale feature matrix

coord_weight_nodes  <- 3.0
demand_weight_nodes <- 1

n_candidates <- 100L

feat_raw <- nodes_enriched %>%
  dplyr::select(X, Y, Node_Total_Demand) %>%
  as.data.frame() %>%
  lapply(to_numeric_safe) %>%     # reuse helper from Step 9A
  as.data.frame()

# Impute any NA/Inf 
for (nm in names(feat_raw)) {
  v <- feat_raw[[nm]]
  if (all(is.na(v))) {
    v <- rep(0, length(v))
  } else {
    v[is.na(v)]       <- median(v, na.rm = TRUE)
    v[!is.finite(v)]  <- 0
  }
  feat_raw[[nm]] <- v
}

feat_scaled <- scale(feat_raw)
feat_scaled[is.na(feat_scaled)] <- 0

# Apply weights
feat_scaled[, c("X", "Y")]          <- feat_scaled[, c("X", "Y")]          * coord_weight_nodes
feat_scaled[, "Node_Total_Demand"]  <- feat_scaled[, "Node_Total_Demand"]  * demand_weight_nodes

# K-means clustering on nodes
k_nodes <- as.integer(n_candidates)
k_nodes <- max(1L, min(k_nodes, nrow(unique(as.data.frame(feat_scaled)))))

set.seed(123)
km_nodes <- stats::kmeans(feat_scaled, centers = k_nodes, nstart = 25, iter.max = 200)
nodes_enriched$Candidate_Cluster <- km_nodes$cluster

cat(sprintf("Clusters requested: %d  |  Clusters created: %d\n",
            n_candidates, n_distinct(nodes_enriched$Candidate_Cluster)))

# Pick one candidate node per cluster
# From each cluster we take the node with the HIGHEST demand — not the geometric centroid.

candidate_nodes_clustered <- nodes_enriched %>%
  group_by(Candidate_Cluster) %>%
  arrange(desc(Node_Total_Demand), Node_Cost) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Open_Location_Index = row_number()) %>%
  dplyr::select(
    Open_Location_Index,
    `NODE ID`,
    Node_ID,
    X, Y,
    Node_Total_Demand,
    Node_Total_Pop,
    Node_Cost,
    Candidate_Cluster)

# Force-include existing SP locations as additional candidates
# The ILP will decide whether to open them — this just ensures they are evaluated

sp_candidates <- sp %>%
  mutate(Node_ID = as.character(Node_ID)) %>%
  left_join(
    nodes_enriched %>%
      mutate(Node_ID = as.character(Node_ID)) %>%
      dplyr::select(Node_ID, Node_Total_Demand, Node_Total_Pop, Node_Cost),
    by = "Node_ID") %>%
  mutate(
    Open_Location_Index = nrow(candidate_nodes_clustered) + row_number(),
    `NODE ID`           = Node_ID,
    Node_Total_Demand   = coalesce(Node_Total_Demand, 0),
    Node_Total_Pop      = coalesce(Node_Total_Pop,    0),
    Node_Cost           = coalesce(Node_Cost, mean(safe_nodes$Node_Cost, na.rm = TRUE)),
    Candidate_Cluster   = NA_integer_) %>%
  dplyr::select(Open_Location_Index, `NODE ID`, Node_ID, X, Y,
                Node_Total_Demand, Node_Total_Pop, Node_Cost, Candidate_Cluster)

# Remove any existing SPs already captured by k-means (avoid duplicates)
sp_candidates <- sp_candidates %>%
  filter(!Node_ID %in% candidate_nodes_clustered$Node_ID)

# Combine clustered candidates with forced SP locations
candidate_nodes_clustered <- candidate_nodes_clustered %>%
  mutate(`NODE ID` = as.character(`NODE ID`))

candidate_nodes_clustered <- bind_rows(candidate_nodes_clustered, sp_candidates) %>%
  mutate(Open_Location_Index = row_number())

cat(sprintf("Total candidates after adding existing SPs : %d\n",
            nrow(candidate_nodes_clustered)))

# Plotting

ggplot() +
  geom_point(
    data  = nodes_enriched,
    aes(x = X, y = Y, colour = log1p(Node_Total_Demand)),
    size  = 0.4, alpha = 0.4) +
  scale_colour_gradient(low = "grey85", high = "firebrick",
                        name = "log(demand)") +
  geom_point(
    data   = candidate_nodes_clustered,
    aes(x = X, y = Y),
    shape  = 21, size = 2.5, stroke = 0.8,
    fill   = "black", colour = "white") +
  coord_equal() +
  labs(
    title    = "Node Clustering: Candidate Locations for ILP",
    subtitle = sprintf(
      "%d candidates selected from %d safe nodes  |  demand-weighted k-means (k = %d)",
      nrow(candidate_nodes_clustered), nrow(nodes_enriched), k_nodes),
    x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.position  = "right")

ggplot() +
  geom_point(
    data  = nodes_enriched,
    aes(x = X, y = Y, colour = factor(Candidate_Cluster)),
    size  = 0.8, alpha = 0.5) +
  geom_point(
    data   = candidate_nodes_clustered,
    aes(x = X, y = Y),
    shape  = 21, size = 3, stroke = 1,
    fill   = "black", colour = "white") +
  coord_equal() +
  guides(colour = "none") +   # too many clusters for a readable legend
  labs(
    title    = "Node Clusters with Selected Candidates",
    subtitle = sprintf(
      "%d clusters — background nodes coloured by cluster, white-ringed dots = selected candidate",
      n_distinct(nodes_enriched$Candidate_Cluster)),
    x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, colour = "grey40"))

# ==============================================================================
# STEP 12: ILP PARAMETER PREPARATION & NETWORK SPARSIFICATION
# ==============================================================================

# Static cost parameters
CapacityCost <- 36.50   # €/slot/year  (€0.10/day × 365)
DeliveryCost <- 1.50    # €/km  (truck driving distance node → assigned SP)

# Fixed costs per candidate
candidate_params <- candidate_nodes_clustered %>%
  mutate(
    Candidate_ID = as.character(Node_ID),
    SP_Cost      = Node_Cost,
    APL_Cost     = Node_Cost / 10) %>%
  dplyr::select(Candidate_ID, X, Y, SP_Cost, APL_Cost)

# Aggregate demand to node level 
node_population <- sub_df %>%
  filter(!is.na(Node_ID)) %>%
  group_by(Node_ID = as.character(Node_ID)) %>%
  summarise(Total_Population = sum(Sub_Population, na.rm = TRUE), .groups = "drop")

node_daily_demand <- daily_sub_demand %>%
  left_join(sub_df %>% dplyr::select(Sub_ID, Node_ID), by = "Sub_ID") %>%
  filter(!is.na(Node_ID)) %>%
  group_by(Node_ID = as.character(Node_ID), `Day Index`) %>%
  summarise(Daily_Demand = sum(Demand, na.rm = TRUE), .groups = "drop")

node_demand_params <- node_daily_demand %>%
  group_by(Node_ID) %>%
  summarise(
    Total_Demand = sum(Daily_Demand, na.rm = TRUE),
    Peak_Demand  = max(Daily_Demand, na.rm = TRUE),
    .groups = "drop") %>%
  left_join(node_population, by = "Node_ID") %>%
  filter(Total_Demand > 0)

demand_nodes    <- as.character(node_demand_params$Node_ID)
candidate_nodes <- as.character(candidate_params$Candidate_ID)

cat(sprintf("Aggregated to %d active demand nodes and %d candidates.\n",
            length(demand_nodes), length(candidate_nodes)))

# Walking distances (demand nodes → candidates) — used for pickup pairs
walk_dist_matrix_m <- distances(
  graph   = ped_graph,
  v       = demand_nodes,
  to      = candidate_nodes,
  mode    = "all",
  weights = E(ped_graph)$weight)

# Snap demand nodes and candidates to truck-accessible nodes before

valid_truck_edges <- edges[edges$truck_allowed == TRUE, ]
valid_truck_nodes <- unique(c(valid_truck_edges$V1, valid_truck_edges$V2))
truck_nodes_df    <- nodes[nodes[["NODE ID"]] %in% valid_truck_nodes, ]

# Snap each demand node to its nearest truck-accessible node
demand_truck_snap <- sapply(demand_nodes, function(nid) {
  coords <- safe_nodes[safe_nodes[["NODE ID"]] == nid, c("X", "Y")]
  as.character(get_nearest_node(coords$X, coords$Y, truck_nodes_df))})

# Snap each candidate to its nearest truck-accessible node
candidate_truck_snap <- sapply(candidate_nodes, function(nid) {
  coords <- candidate_params[candidate_params$Candidate_ID == nid, c("X", "Y")]
  as.character(get_nearest_node(coords$X, coords$Y, truck_nodes_df))})

cat(sprintf("Truck-snapped demand nodes    : %d unique truck nodes\n", n_distinct(demand_truck_snap)))
cat(sprintf("Truck-snapped candidate nodes : %d unique truck nodes\n", n_distinct(candidate_truck_snap)))

# Truck driving distances using snapped nodes
truck_dist_matrix_raw <- distances(
  graph   = truck_graph,
  v       = unique(demand_truck_snap),
  to      = unique(candidate_truck_snap),
  mode    = "out",          # directed — trucks follow one-way rules
  weights = E(truck_graph)$weight)

# Map results back to original demand_nodes and candidate_nodes ordering
truck_dist_matrix_m <- truck_dist_matrix_raw[
  demand_truck_snap,
  candidate_truck_snap,
  drop = FALSE]
rownames(truck_dist_matrix_m) <- demand_nodes
colnames(truck_dist_matrix_m) <- candidate_nodes

# Node-level Stage 2 demographic correction
node_corrections <- sub_df %>%
  filter(!is.na(Node_ID)) %>%
  group_by(Node_ID = as.character(Node_ID)) %>%
  summarise(Stage2_Correction = mean(Stage2_Correction, na.rm = TRUE), .groups = "drop") %>%
  mutate(Stage2_Correction = coalesce(Stage2_Correction, 0))

# PICKUP PAIRS — sparse, walk ≤ 2 km

pickup_pairs <- as.data.frame(as.table(walk_dist_matrix_m)) %>%
  rename(Node_ID = Var1, Candidate_ID = Var2, Walk_Dist_m = Freq) %>%
  mutate(
    Node_ID      = as.character(Node_ID),
    Candidate_ID = as.character(Candidate_ID)) %>%
  filter(Walk_Dist_m <= 2000) %>%
  left_join(node_corrections, by = "Node_ID") %>%
  mutate(
    Stage2_Correction = coalesce(Stage2_Correction, 0),
    Base_Rate         = a_fit * exp(-b_fit * Walk_Dist_m),
    Pickup_Rate       = pmin(pmax(Base_Rate + Stage2_Correction, 0), 1)) %>%
  # Distance rank within each demand node (1 = closest candidate)
  group_by(Node_ID) %>%
  arrange(Walk_Dist_m) %>%
  mutate(rank_pickup = row_number()) %>%
  ungroup() %>%
  dplyr::select(Node_ID, Candidate_ID, Walk_Dist_m, Pickup_Rate, rank_pickup)

cat(sprintf("Pickup pairs (walk ≤ 2 km): %d\n", nrow(pickup_pairs)))


# DELIVERY PAIRS — based on truck distances, top 50 per demand node
delivery_pairs_raw <- as.data.frame(as.table(truck_dist_matrix_m)) %>%
  rename(Node_ID = Var1, Candidate_ID = Var2, Truck_Dist_m = Freq) %>%
  mutate(
    Node_ID      = as.character(Node_ID),
    Candidate_ID = as.character(Candidate_ID)) %>%
  filter(is.finite(Truck_Dist_m)) %>%
  group_by(Node_ID) %>%
  arrange(Truck_Dist_m) %>%
  mutate(rank_delivery = row_number()) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  dplyr::select(Node_ID, Candidate_ID, Truck_Dist_m, rank_delivery)

# Build a quick XY lookup for every candidate
candidate_xy <- candidate_params %>%
  dplyr::select(Candidate_ID, X, Y)

# Which demand nodes are currently uncovered?
covered_nodes   <- unique(delivery_pairs_raw$Node_ID)
uncovered_nodes <- setdiff(demand_nodes, covered_nodes)

if (length(uncovered_nodes) > 0) {
  cat(sprintf(
    "Snapping %d uncovered demand node(s) to nearest candidate (Euclidean fallback)...\n",
    length(uncovered_nodes)
  ))
  
  # XY of uncovered demand nodes (from node_demand_params + safe_nodes)
  uncovered_xy <- safe_nodes %>%
    mutate(Node_ID = as.character(`NODE ID`)) %>%
    filter(Node_ID %in% uncovered_nodes) %>%
    dplyr::select(Node_ID, X, Y)
  
  # For each uncovered node, find the single nearest candidate by Euclidean distance
  snap_rows <- lapply(seq_len(nrow(uncovered_xy)), function(i) {
    nx <- uncovered_xy$X[i]
    ny <- uncovered_xy$Y[i]
    sq_dist <- (candidate_xy$X - nx)^2 + (candidate_xy$Y - ny)^2
    best    <- which.min(sq_dist)
    data.frame(
      Node_ID      = uncovered_xy$Node_ID[i],
      Candidate_ID = candidate_xy$Candidate_ID[best],
      Truck_Dist_m = sqrt(sq_dist[best]),   # Euclidean proxy — finite by construction
      rank_delivery = 1L,
      stringsAsFactors = FALSE)
  })
  
  snap_df <- do.call(rbind, snap_rows)
  
  # Append snapped rows to the raw delivery pairs
  delivery_pairs_raw <- bind_rows(delivery_pairs_raw, snap_df)
  
  cat(sprintf("After snap: %d demand nodes now have at least one delivery candidate.\n",
              n_distinct(delivery_pairs_raw$Node_ID)))
}

# Guard: every demand node must have at least one delivery candidate
uncovered <- setdiff(demand_nodes, unique(delivery_pairs_raw$Node_ID))
if (length(uncovered) > 0) {
  stop(sprintf(
    "%d demand node(s) have no reachable candidate in the truck graph. Check graph connectivity.\nFirst 5: %s",
    length(uncovered), paste(head(uncovered, 5), collapse = ", ")))
}

delivery_pairs <- delivery_pairs_raw

cat(sprintf("Delivery pairs (top 20 per node): %d\n", nrow(delivery_pairs)))

# --- Sparsity summary --------------------------------------------------------
total_possible <- length(demand_nodes) * length(candidate_nodes)
cat(sprintf("\n--- ILP PARAMETER DIMENSIONS ---\n"))
cat(sprintf("Candidates (J)            : %d\n",   nrow(candidate_params)))
cat(sprintf("Demand nodes (N)          : %d\n",   nrow(node_demand_params)))
cat(sprintf("Pickup pairs  (p_nj)      : %d  (%.1f%% of N×J)\n",
            nrow(pickup_pairs),   100 * nrow(pickup_pairs)   / total_possible))
cat(sprintf("Delivery pairs (d_nj)     : %d  (%.1f%% of N×J)\n",
            nrow(delivery_pairs), 100 * nrow(delivery_pairs) / total_possible))


# ==============================================================================
# STEP 13: ILP — NEAREST-FACILITY SINGLE-ASSIGNMENT FORMULATION
# ==============================================================================
#
# DECISION VARIABLES
# ------------------
#   x_j   ∈ {0,1}   open a Service Point at candidate j          [idx_x]
#   y_j   ∈ {0,1}   open an APL at candidate j                   [idx_y]
#   p_nj  ∈ {0,1}   node n assigned to j as pickup facility      [idx_p]
#   d_nj  ∈ {0,1}   node n assigned to j as delivery SP          [idx_d]
#   Cap_j ≥ 0       package slot capacity at facility j           [idx_cap]
#
# CONSTRAINTS
# -----------
#   C1  Mutual exclusivity          x_j + y_j ≤ 1
#   C2  Pickup assignment coverage  Σ_j p_nj ≤ 1          ∀ n
#   C3  Pickup only if open         p_nj ≤ x_j + y_j      ∀ (n,j) in pickup_pairs
#   C5  Delivery assignment cover   Σ_j d_nj = 1          ∀ n
#   C6  Delivery only from SP       d_nj ≤ x_j            ∀ (n,j) in delivery_pairs
#   C8  Per-location capacity       Cap_j ≥ Σ_n p_nj*Peak_n  (pickup storage)
#                                   + Σ_n d_nj*(1−PR_nj)*Peak_n  (delivery holding, SP only)
#                                   allowing 2 % per-location overflow →
#                                   Cap_j ≥ 1/1.02 * [assigned peak]
#   C9  Capacity slack              Max demand - capacity
#   C10  City-wide overflow         Network overflow  < 1 - 1/1.01
#
# OBJECTIVE  (minimise)
# ---------------------
#   Setup costs   : Σ_j SP_Cost_j*x_j + APL_Cost_j*y_j
#   Capacity costs: Σ_j Cap_j * CapacityCost
#   Delivery costs: Σ_{n,j} d_nj * (1−PR_nj) * Total_n * DeliveryCost * Truck_km_nj
#                   where Truck_km_nj is the truck driving distance node n → candidate j
# ==============================================================================

J <- nrow(candidate_params)
N <- nrow(node_demand_params)

cand_idx <- setNames(seq_len(J), candidate_params$Candidate_ID)
node_idx <- setNames(seq_len(N), node_demand_params$Node_ID)

# Attach integer indices to pair tables
pickup_pairs <- pickup_pairs %>%
  mutate(n = node_idx[Node_ID], j = cand_idx[Candidate_ID]) %>%
  filter(!is.na(n), !is.na(j)) %>%
  arrange(n, j)

delivery_pairs <- delivery_pairs %>%
  mutate(n = node_idx[Node_ID], j = cand_idx[Candidate_ID]) %>%
  filter(!is.na(n), !is.na(j)) %>%
  arrange(n, j)

NP <- nrow(pickup_pairs)
ND <- nrow(delivery_pairs)

# Variable blocks

idx_x   <- seq_len(J)
idx_y   <- J + seq_len(J)
idx_p   <- 2*J + seq_len(NP)
idx_d   <- 2*J + NP + seq_len(ND)
idx_cap <- 2*J + NP + ND + seq_len(J)
idx_s  <- 2*J + NP + ND + J + seq_len(J)  

n_vars <- 2*J + NP + ND + J + J

peak_n  <- node_demand_params$Peak_Demand
total_n <- node_demand_params$Total_Demand

# Delivery cost coefficient for d_nj
# Pickup rate for delivery pairs comes from pickup_pairs where available
pickup_lookup <- pickup_pairs %>%
  dplyr::select(Node_ID, Candidate_ID, Pickup_Rate) %>%
  rename(PR_from_pickup = Pickup_Rate)

delivery_pairs <- delivery_pairs %>%
  left_join(pickup_lookup, by = c("Node_ID", "Candidate_ID")) %>%
  mutate(
    Pickup_Rate_eff = coalesce(PR_from_pickup, 0),
    Truck_km_eff    = Truck_Dist_m / 1000,
    DeliveryCost_nj = total_n[n] * DeliveryCost * Truck_km_eff
  )

pickup_pairs <- pickup_pairs %>%
  left_join(
    delivery_pairs %>% dplyr::select(Node_ID, Candidate_ID, Truck_Dist_m),
    by = c("Node_ID", "Candidate_ID")
  ) %>%
  mutate(
    # If no truck route exists for this pair, saving is 0 (truck wasn't going there anyway)
    Truck_km_eff = coalesce(Truck_Dist_m / 1000, 0)
  )

# After building delivery_pairs, find nodes with no pickup pair
uncovered_pickup_nodes <- setdiff(node_demand_params$Node_ID, unique(pickup_pairs$Node_ID))

# For those nodes, ensure Pickup_Rate_eff = 0 in delivery_pairs (100% delivery)
delivery_pairs <- delivery_pairs %>%
  mutate(Pickup_Rate_eff = if_else(Node_ID %in% uncovered_pickup_nodes, 0, Pickup_Rate_eff))

# Objective vector
obj <- numeric(n_vars)
obj[idx_x]   <- candidate_params$SP_Cost
obj[idx_y]   <- candidate_params$APL_Cost
obj[idx_cap] <- CapacityCost
obj[idx_d]   <- delivery_pairs$DeliveryCost_nj


# Negative cost on p_nj — when pickup is assigned it saves delivery cost
pickup_delivery_km <- delivery_pairs %>%
  group_by(Node_ID) %>%
  slice_min(Truck_Dist_m, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Node_ID, Best_Truck_km = Truck_km_eff)

pickup_pairs <- pickup_pairs %>%
  left_join(pickup_delivery_km, by = "Node_ID") %>%
  mutate(Best_Truck_km = coalesce(Best_Truck_km, Truck_km_eff))

obj[idx_p] <- -pickup_pairs$Pickup_Rate *
  total_n[pickup_pairs$n] *
  DeliveryCost *
  pickup_pairs$Best_Truck_km

# Variable bounds & types 
var_lb   <- numeric(n_vars)                                        
var_ub   <- c(rep(1, 2*J + NP + ND), rep(Inf, J), rep(Inf, J))   
var_type <- c(rep("B", 2*J + NP + ND), rep("C", J), rep("C", J))

# Pre-index pair row numbers
pickup_pairs$p_idx   <- seq_len(NP)
delivery_pairs$d_idx <- seq_len(ND)

# Accumulators — we build all triplets first, then bind ONCE at the end
all_rows <- vector("list", 11)
all_cols <- vector("list", 11)
all_vals <- vector("list", 11)
all_rhs  <- vector("list", 11)
all_dirs <- vector("list", 11)

current_row <- 0L

# ------------------------------------------------------------------------------
# C1: Mutual exclusivity — x_j + y_j ≤ 1   (J constraints)
# ------------------------------------------------------------------------------
c1_n   <- J
c1_r   <- rep(seq_len(J), 2)                    # row for x_j  then row for y_j
c1_c   <- c(idx_x, idx_y)
c1_v   <- rep(1, 2 * J)

all_rows[[1]] <- c1_r
all_cols[[1]] <- c1_c
all_vals[[1]] <- c1_v
all_rhs [[1]] <- rep(1,  J)
all_dirs[[1]] <- rep("<=", J)

current_row <- current_row + c1_n
cat(sprintf("After C1: %d constraints\n", current_row))

# ------------------------------------------------------------------------------
# C2: Pickup coverage — Σ_j p_nj ≤ 1   ∀ n with at least one pickup pair
# ------------------------------------------------------------------------------
covered_pickup_nodes <- sort(unique(pickup_pairs$n))
node_to_c2_row       <- setNames(
  current_row + seq_along(covered_pickup_nodes),
  covered_pickup_nodes
)
c2_n <- length(covered_pickup_nodes)

all_rows[[2]] <- node_to_c2_row[as.character(pickup_pairs$n)]
all_cols[[2]] <- idx_p[pickup_pairs$p_idx]
all_vals[[2]] <- rep(1, NP)
all_rhs [[2]] <- rep(1,  c2_n)
all_dirs[[2]] <- rep("<=", c2_n)

current_row <- current_row + c2_n
cat(sprintf("After C2: %d constraints\n", current_row))

# ------------------------------------------------------------------------------
# C3: Pickup only if open — p_nj ≤ x_j + y_j   (NP constraints, 3 terms each)
# ------------------------------------------------------------------------------
c3_n   <- NP
c3_base <- current_row + seq_len(NP)            # one constraint per pickup pair

all_rows[[3]] <- c(c3_base, c3_base, c3_base)
all_cols[[3]] <- c(idx_p[pickup_pairs$p_idx],
                   idx_x[pickup_pairs$j],
                   idx_y[pickup_pairs$j])
all_vals[[3]] <- c(rep( 1, NP), rep(-1, NP), rep(-1, NP))
all_rhs [[3]] <- rep(0, NP)
all_dirs[[3]] <- rep("<=", NP)

current_row <- current_row + c3_n
cat(sprintf("After C3: %d constraints\n", current_row))

# C4 removed — ordered exclusion for pickup was over-constraining the LP
# relaxation. The objective naturally prefers nearest-facility assignment.
all_rows[[4]] <- integer(0)
all_cols[[4]] <- integer(0)
all_vals[[4]] <- numeric(0)
all_rhs [[4]] <- numeric(0)
all_dirs[[4]] <- character(0)

# ------------------------------------------------------------------------------
# C5: Delivery coverage — Σ_j d_nj = 1   ∀ n
# ------------------------------------------------------------------------------
covered_delivery_nodes <- sort(unique(delivery_pairs$n))
node_to_c5_row         <- setNames(
  current_row + seq_along(covered_delivery_nodes),
  covered_delivery_nodes
)
c5_n <- length(covered_delivery_nodes)

all_rows[[5]] <- node_to_c5_row[as.character(delivery_pairs$n)]
all_cols[[5]] <- idx_d[delivery_pairs$d_idx]
all_vals[[5]] <- rep(1, ND)
all_rhs [[5]] <- rep(1,  c5_n)
all_dirs[[5]] <- rep("==", c5_n)

current_row <- current_row + c5_n
cat(sprintf("After C5: %d constraints\n", current_row))


# ------------------------------------------------------------------------------
# C6: Delivery only from SP — d_nj ≤ x_j   (ND constraints, 2 terms each)
# ------------------------------------------------------------------------------
c6_n    <- ND
c6_base <- current_row + seq_len(ND)

all_rows[[6]] <- c(c6_base, c6_base)
all_cols[[6]] <- c(idx_d[delivery_pairs$d_idx],
                   idx_x[delivery_pairs$j])
all_vals[[6]] <- c(rep( 1, ND), rep(-1, ND))
all_rhs [[6]] <- rep(0, ND)
all_dirs[[6]] <- rep("<=", ND)

current_row <- current_row + c6_n
cat(sprintf("After C6: %d constraints\n", current_row))

# C7 removed — ordered exclusion for delivery was over-constraining the LP
# relaxation. C6 already ensures delivery only goes to open SPs, and the
# objective naturally prefers nearest-SP assignment via the truck cost term.
all_rows[[7]] <- integer(0)
all_cols[[7]] <- integer(0)
all_vals[[7]] <- numeric(0)
all_rhs [[7]] <- numeric(0)
all_dirs[[7]] <- character(0)

# ------------------------------------------------------------------------------
# C8: Per-location capacity — Cap_j ≥ 1/1.02 × [pickup peak + delivery holding]
#     Only for candidates that have at least one assigned pair.
# ------------------------------------------------------------------------------
overflow_factor <- 1/1.02

active_j <- sort(union(unique(pickup_pairs$j), unique(delivery_pairs$j)))
c8_n     <- length(active_j)

# Map candidate j → its C8 constraint row
j_to_c8_row <- setNames(current_row + seq_along(active_j), active_j)

# Capacity variable term:  +1 * Cap_j
c8_cap_rows <- j_to_c8_row[as.character(active_j)]
c8_cap_cols <- idx_cap[active_j]
c8_cap_vals <- rep(1, c8_n)

# Pickup term:  - 1/1.02 * Peak_n * p_nj
c8_p_rows <- j_to_c8_row[as.character(pickup_pairs$j)]
c8_p_cols <- idx_p[pickup_pairs$p_idx]
c8_p_vals <- -overflow_factor * peak_n[pickup_pairs$n]

# Delivery holding term:  - 1/1.02 * (1 - PR_nj) * Peak_n * d_nj
c8_d_rows <- j_to_c8_row[as.character(delivery_pairs$j)]
c8_d_cols <- idx_d[delivery_pairs$d_idx]
c8_d_vals <- -overflow_factor * (1 - delivery_pairs$Pickup_Rate_eff) * peak_n[delivery_pairs$n]

all_rows[[8]] <- c(c8_cap_rows, c8_p_rows, c8_d_rows)
all_cols[[8]] <- c(c8_cap_cols, c8_p_cols, c8_d_cols)
all_vals[[8]] <- c(c8_cap_vals, c8_p_vals, c8_d_vals)
all_rhs [[8]] <- rep(0, c8_n)
all_dirs[[8]] <- rep(">=", c8_n)

current_row <- current_row + c8_n
cat(sprintf("After C8: %d constraints\n", current_row))

# ------------------------------------------------------------------------------
# C9: Per-facility slack — s_j ≥ AssignedPeak_j - Cap_j   ∀ j
#
# s_j captures the actual overflow at facility j (if any).
# Rearranged for ILP: s_j + Cap_j - AssignedPeak_j ≥ 0
# i.e.: s_j + Cap_j - Σ_n p_nj*Peak_n - Σ_n d_nj*(1-PR_nj)*Peak_n ≥ 0
# ------------------------------------------------------------------------------
c9_n        <- length(active_j)
c9_base     <- current_row + seq_along(active_j)
j_to_c9_row <- setNames(c9_base, active_j)
current_row <- current_row + c9_n

all_rows[[9]] <- c(
  j_to_c9_row[as.character(active_j)],
  j_to_c9_row[as.character(active_j)],
  j_to_c9_row[as.character(pickup_pairs$j)],
  j_to_c9_row[as.character(delivery_pairs$j)]
)
all_cols[[9]] <- c(
  idx_s[active_j],
  idx_cap[active_j],
  idx_p[pickup_pairs$p_idx],
  idx_d[delivery_pairs$d_idx]
)
all_vals[[9]] <- c(
  rep(1, c9_n),
  rep(1, c9_n),
  -peak_n[pickup_pairs$n],
  -(1 - delivery_pairs$Pickup_Rate_eff) * peak_n[delivery_pairs$n]
)
all_rhs [[9]] <- rep(0, c9_n)
all_dirs[[9]] <- rep(">=", c9_n)

# ------------------------------------------------------------------------------
# C10: Average facility overflow ≤ 1%
#      (1/J) * Σ_j s_j ≤ (1-1/1.01) * (1/J) * Σ_j AssignedPeak_j
#      1/J cancels → Σ_j s_j ≤ (1-1/1.01) * Σ_j AssignedPeak_j
#      Rearranged: Σ_j s_j - (1-1/1.01) * Σ_j AssignedPeak_j ≤ 0
#      i.e.: Σ_j s_j - (1-1/1.01)*Σ_{n,j} p_nj*Peak_n
#                    - (1-1/1.01)*Σ_{n,j} d_nj*(1-PR_nj)*Peak_n ≤ 0
# ------------------------------------------------------------------------------
current_row <- current_row + 1L

all_rows[[10]] <- c(
  rep(current_row, length(active_j)),                          # s_j terms
  rep(current_row, nrow(pickup_pairs)),         # p_nj terms
  rep(current_row, nrow(delivery_pairs))        # d_nj terms
)
all_cols[[10]] <- c(
  idx_s[active_j],                                        # +s_j
  idx_p[pickup_pairs$p_idx],                    # p_nj columns
  idx_d[delivery_pairs$d_idx]                   # d_nj columns
)
all_vals[[10]] <- c(
  rep(1, length(active_j)),                                                                    # +s_j
  -(1-1/1.01) * peak_n[pickup_pairs$n],                                                # -(1-1/1.01)*Peak_n*p_nj
  -(1-1/1.01) * (1 - delivery_pairs$Pickup_Rate_eff) * peak_n[delivery_pairs$n]        # -(1-1/1.01)*(1-PR)*Peak_n*d_nj
)
all_rhs [[10]] <- 0
all_dirs[[10]] <- "<="

cat(sprintf("After C10: %d constraints\n", current_row))

# ------------------------------------------------------------------------------
# C11: placeholder (empty — reserved for future constraints)
# ------------------------------------------------------------------------------
all_rows[[11]] <- integer(0)
all_cols[[11]] <- integer(0)
all_vals[[11]] <- numeric(0)
all_rhs [[11]] <- numeric(0)
all_dirs[[11]] <- character(0)

# ------------------------------------------------------------------------------
# Assemble — single bind of all triplets
# ------------------------------------------------------------------------------
n_constraints <- current_row
rows <- unlist(all_rows)
cols <- unlist(all_cols)
vals <- unlist(all_vals)
rhs  <- unlist(all_rhs)
dirs <- unlist(all_dirs)

cat(sprintf("Total constraints: %d\n", n_constraints))
cat(sprintf("Total variables  : %d\n", n_vars))

A <- sparseMatrix(i = rows, j = cols, x = vals,
                  dims = c(n_constraints, n_vars))

# Solve
cat("Solving ILP...\n")
result <- Rglpk_solve_LP(
  obj     = obj,
  mat     = A,
  dir     = dirs,
  rhs     = rhs,
  types   = var_type,
  bounds  = list(lower = list(ind = seq_len(n_vars), val = var_lb),
                 upper = list(ind = seq_len(n_vars), val = var_ub)),
  control = list(
    verbose    = TRUE,
    presolve   = TRUE))

cat(sprintf("Solver status: %d  (%s)\n", result$status,
            if (result$status == 0) "OPTIMAL" else "NON-OPTIMAL / INFEASIBLE"))

# --- Extract solution --------------------------------------------------------
x_sol   <- round(result$solution[idx_x])
y_sol   <- round(result$solution[idx_y])
p_sol   <- round(result$solution[idx_p])
d_sol   <- round(result$solution[idx_d])
cap_sol <- ceiling(result$solution[idx_cap])
optimum_sol <- round(result$optimum, 2)

# --- Facility summary --------------------------------------------------------
open_facilities <- candidate_params %>%
  mutate(
    j             = seq_len(J),
    Open_SP       = x_sol,
    Open_APL      = y_sol,
    Capacity      = cap_sol,
    Type          = case_when(
      Open_SP  == 1 ~ "Service Point",
      Open_APL == 1 ~ "APL",
      TRUE          ~ "Closed"),
    Facility_Cost   = case_when(Open_SP == 1 ~ SP_Cost, Open_APL == 1 ~ APL_Cost, TRUE ~ 0),
    Annual_Cap_Cost = Capacity * CapacityCost) %>%
  filter(Type != "Closed")

# --- Demand routing summary --------------------------------------------------
pickup_assignment <- pickup_pairs %>%
  mutate(Assigned = p_sol) %>%
  filter(Assigned == 1) %>%
  left_join(node_demand_params %>% dplyr::select(Node_ID, Total_Demand, Peak_Demand),
            by = "Node_ID") %>%
  mutate(
    Pickup_Volume   = Pickup_Rate * Total_Demand,
    Delivery_Volume = (1 - Pickup_Rate) * Total_Demand)

delivery_assignment <- delivery_pairs %>%
  mutate(Assigned = d_sol) %>%
  filter(Assigned == 1) %>%
  left_join(node_demand_params %>% dplyr::select(Node_ID, Total_Demand, Peak_Demand),
            by = "Node_ID") %>%
  mutate(
    Delivery_Volume = (1 - Pickup_Rate_eff) * Total_Demand,
    Delivery_Cost   = Delivery_Volume * DeliveryCost * (Truck_Dist_m / 1000))

# --- Solution report ---------------------------------------------------------
true_fixed    <- sum(obj[idx_x]   * x_sol)
true_apl      <- sum(obj[idx_y]   * y_sol)
true_cap      <- sum(obj[idx_cap] * cap_sol)
true_delivery <- sum(obj[idx_d]   * d_sol)
true_savings  <- sum(obj[idx_p]   * p_sol)  # negative
net_delivery  <- true_delivery + true_savings

node_volumes <- node_demand_params %>%
  mutate(Node_ID = as.character(Node_ID)) %>%
  left_join(
    pickup_assignment %>%
      dplyr::select(Node_ID, Pickup_Volume),
    by = "Node_ID") %>%
  mutate(
    Pickup_Volume   = coalesce(Pickup_Volume, 0),
    Delivery_Volume = Total_Demand - Pickup_Volume)

total_pickup_volume   <- sum(node_volumes$Pickup_Volume)
total_delivery_volume <- sum(node_volumes$Delivery_Volume)

cat(sprintf("\n=== SOLUTION SUMMARY ===\n"))
cat(sprintf("Service Points opened : %d\n", sum(open_facilities$Type == "Service Point")))
cat(sprintf("APLs opened           : %d\n", sum(open_facilities$Type == "APL")))
cat(sprintf("\nFixed SP costs        : €%s\n", format(round(true_fixed),       big.mark = ",")))
cat(sprintf("Fixed APL costs       : €%s\n", format(round(true_apl),         big.mark = ",")))
cat(sprintf("Fixed setup costs     : €%s\n", format(round(true_fixed + true_apl), big.mark = ",")))
cat(sprintf("Annual capacity costs : €%s\n", format(round(true_cap),         big.mark = ",")))
cat(sprintf("Net delivery costs    : €%s\n", format(round(net_delivery),     big.mark = ",")))

cat(sprintf("Total objective value : €%s\n", format(round(result$optimum),   big.mark = ",")))
cat(sprintf("\nTotal annual demand   : %.0f packages\n", sum(total_n)))
cat(sprintf("Picked up (annual)    : %.0f packages (%.1f%%)\n",
            total_pickup_volume,  100 * total_pickup_volume  / sum(total_n)))
cat(sprintf("Delivered (annual)    : %.0f packages (%.1f%%)\n",
            total_delivery_volume, 100 * total_delivery_volume / sum(total_n)))

# --- Per facility detail -----------------------------------------------------
cat(sprintf("\n--- Open Facilities ---\n"))
print(open_facilities %>%
        dplyr::select(Candidate_ID, Type, Facility_Cost, Capacity, Annual_Cap_Cost) %>%
        arrange(Type, desc(Capacity)))

# --- Capacity overflow check -------------------------------------------------
cap_check <- open_facilities %>%
  mutate(j = match(Candidate_ID, candidate_params$Candidate_ID)) %>%
  rowwise() %>%
  mutate(
    Assigned_Peak_Pickup   = sum(p_sol[pickup_pairs$j == j] *
                                   peak_n[pickup_pairs$n[pickup_pairs$j == j]], na.rm = TRUE),
    Assigned_Peak_Delivery = sum(d_sol[delivery_pairs$j == j] *
                                   (1 - delivery_pairs$Pickup_Rate_eff[delivery_pairs$j == j]) *
                                   peak_n[delivery_pairs$n[delivery_pairs$j == j]], na.rm = TRUE),
    Assigned_Peak_Total    = Assigned_Peak_Pickup + if_else(Type == "Service Point", Assigned_Peak_Delivery, 0),
    Overflow_Pct           = 100 * pmax(Assigned_Peak_Total / Capacity - 1, 0)) %>%
  ungroup()


global_assigned_peak = sum(cap_check$Assigned_Peak_Total)
global_capacity = sum(cap_check$Capacity)
global_overflow_pct = 100 * pmax(global_assigned_peak - global_capacity, 0) / global_capacity

cat(sprintf("\n--- Capacity Overflow Check ---\n"))
cat(sprintf("Max global overflow : %.2f%%\n", global_overflow_pct))
cat(sprintf("City-wide capacity            : %.0f slots\n",   sum(cap_check$Capacity)))
cat(sprintf("Avg facility overflow         : %.2f%%\n",
            mean(cap_check$Overflow_Pct, na.rm = TRUE)))

if (sum(d_sol) == N) {
  cat("All Demand Served - ILP solution correct")}

cat(sprintf("\n--- Per Facility Overflow Detail ---\n"))
print(cap_check %>%
        dplyr::select(Candidate_ID, Type, Capacity,
                      Assigned_Peak_Pickup, Assigned_Peak_Delivery,
                      Assigned_Peak_Total, Overflow_Pct) %>%
        arrange(desc(Overflow_Pct)))

# --- Demand node routing detail ----------------------------------------------
cat(sprintf("\n--- Pickup Assignment Detail ---\n"))
cat(sprintf("Nodes with pickup assignment  : %d / %d\n",
            n_distinct(pickup_assignment$Node_ID), N))
print(pickup_assignment %>%
        dplyr::select(Node_ID, Candidate_ID, Walk_Dist_m, Pickup_Rate,
                      Pickup_Volume, Delivery_Volume) %>%
        arrange(desc(Pickup_Volume)))

cat(sprintf("\n--- Delivery Assignment Detail ---\n"))
cat(sprintf("Nodes with delivery assignment: %d / %d\n",
            n_distinct(delivery_assignment$Node_ID), N))
print(delivery_assignment %>%
        dplyr::select(Node_ID, Candidate_ID, Truck_Dist_m, Pickup_Rate_eff,
                      Delivery_Volume, Delivery_Cost) %>%
        arrange(desc(Delivery_Cost)))


# ==============================================================================
# SECTION 14: ILP SOLUTION vs CURRENT SERVICE POINT NETWORK
# ==============================================================================

pal_current <- "#E07B39"   # warm orange  = current network
pal_ilp     <- "#2C7BB6"   # steel blue   = ILP solution
pal_kept    <- "#4DAC26"   # green        = locations kept in ILP
pal_closed  <- "#D7191C"   # red          = current locations NOT in ILP
pal_new     <- "#2C7BB6"   # blue         = new ILP-only locations

theme_section14 <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, colour = "grey10"),
    plot.subtitle    = element_text(size = 9,  colour = "grey45", margin = margin(b = 8)),
    plot.caption     = element_text(size = 8,  colour = "grey55", hjust = 0),
    axis.text        = element_text(size = 9,  colour = "grey30"),
    axis.title       = element_text(size = 10, colour = "grey20"),
    panel.grid.major = element_line(colour = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position  = "top",
    legend.text      = element_text(size = 9),
    strip.text       = element_text(face = "bold", size = 10))


# ------------------------------------------------------------------------------
# 14.1  Classify every candidate location: Kept | Closed | New  ---------------
# ------------------------------------------------------------------------------

# Current SP Node IDs (from Step 4)
current_sp_nodes <- as.character(sp$Node_ID)

# ILP-selected SP nodes
ilp_sp_nodes <- as.character(
  open_facilities %>%
    filter(Type == "Service Point") %>%
    pull(Candidate_ID))

# ILP-selected APL nodes
ilp_apl_nodes <- as.character(
  open_facilities %>%
    filter(Type == "APL") %>%
    pull(Candidate_ID))

ilp_any_open <- union(ilp_sp_nodes, ilp_apl_nodes)

# Classification
kept_nodes   <- intersect(current_sp_nodes, ilp_any_open)   # current AND in ILP
closed_nodes <- setdiff(current_sp_nodes,   ilp_any_open)   # current but NOT in ILP
new_nodes    <- setdiff(ilp_any_open,       current_sp_nodes) # ILP-only

cat(sprintf(
  "\n=== 14.1 Location classification ===\n  Kept   : %d\n  Closed : %d\n  New    : %d\n",
  length(kept_nodes), length(closed_nodes), length(new_nodes)))

# Build a unified map data frame
map_current <- sp %>%
  mutate(
    Node_ID = as.character(Node_ID),
    Status  = case_when(
      Node_ID %in% kept_nodes   ~ "Kept",
      Node_ID %in% closed_nodes ~ "Closed",
      TRUE                      ~ "Current Only"),
    Network = "Current")

map_ilp <- open_facilities %>%
  mutate(
    Node_ID = as.character(Candidate_ID),
    Status  = case_when(
      Node_ID %in% kept_nodes ~ "Kept",
      TRUE                    ~ "New"), 
    Network = "ILP")

# ------ PLOT 14.1A : Side-by-side map of current vs ILP locations -------------

p_map_current <- ggplot() +
  # CBS square outlines for context
  geom_rect(
    data = cbs %>% filter(!is.na(Population) & Population > 0),
    aes(
      xmin = X - main_w / 2, xmax = X + main_w / 2,
      ymin = Y - main_h / 2, ymax = Y + main_h / 2),
    fill = "grey96", colour = "grey85", linewidth = 0.2) +
  geom_point(
    data  = map_current,
    aes(x = X, y = Y, colour = Status),
    size  = 3.5, shape = 16) +
  scale_colour_manual(
    values = c("Kept" = pal_kept, "Closed" = pal_closed, "Current Only" = pal_current),
    name   = NULL) +
  coord_equal() +
  labs(
    title    = "Current Service Points",
    subtitle = sprintf(
      "%d locations  |  %d kept in ILP  |  %d closed by ILP",
      nrow(map_current), length(kept_nodes), length(closed_nodes))) +
  theme_section14 +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

p_map_ilp <- ggplot() +
  geom_rect(
    data = cbs %>% filter(!is.na(Population) & Population > 0),
    aes(
      xmin = X - main_w / 2, xmax = X + main_w / 2,
      ymin = Y - main_h / 2, ymax = Y + main_h / 2),
    fill = "grey96", colour = "grey85", linewidth = 0.2) +
  geom_point(
    data  = map_ilp,
    aes(x = X, y = Y, colour = Status, shape = Type),
    size  = 3.5) +
  scale_colour_manual(
    values = c("Kept" = pal_kept, "New" = pal_new),
    name   = NULL) +
  scale_shape_manual(
    values = c("Service Point" = 16, "APL" = 17),
    name   = NULL) +
  coord_equal() +
  labs(
    title    = "ILP-Optimised Network",
    subtitle = sprintf(
      "%d open locations  |  %d Service Points  |  %d APLs  |  %d new sites",
      nrow(map_ilp),
      sum(open_facilities$Type == "Service Point"),
      sum(open_facilities$Type == "APL"),
      length(new_nodes))) +
  theme_section14 +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

# Combined map panel
grid.arrange(p_map_current, p_map_ilp, ncol = 2,
             top = grid::textGrob(
               "Section 14.1 — Service Point Locations: Current vs ILP",
               gp = grid::gpar(fontface = "bold", fontsize = 15)))


# ------------------------------------------------------------------------------
#  14.2 Capacity Analysis & Visualization  -------------------------------------
# ------------------------------------------------------------------------------

# Make baseline data table for capacities and demands
baseline_daily <- activity %>%
  filter(!(`Day Index` %in% c(359,360)), `Day of Week` != "Sun") %>%
  mutate(
    Daily_Demand = Pickups + Deliveries,
    Location_ID = as.character(`Location ID`)) %>%
  group_by(Location_ID) %>%
  summarise(
    Mean_Daily_Demand  = mean(Daily_Demand),
    Peak_Daily_Demand  = max(Daily_Demand),
    Total_Annual_Demand = sum(Daily_Demand),
    Capacity           = Peak_Daily_Demand / 1.02,
    Mean_Utilisation = 100 * Mean_Daily_Demand / Capacity,
    Peak_Utilisation = 100 * Peak_Daily_Demand / Capacity)

# Make same table for ilp
ilp_daily_facility <- node_daily_demand %>%
  left_join(
    delivery_assignment %>%
      select(Node_ID, Candidate_ID) %>%
      distinct(),
    by = "Node_ID") %>%
  group_by(Candidate_ID, `Day Index`) %>%
  summarise(Daily_Demand = sum(Daily_Demand), .groups = "drop") %>%
  group_by(Candidate_ID) %>%
  summarise(
    Mean_Daily_Demand = mean(Daily_Demand),
    Peak_Daily_Demand = max(Daily_Demand),
    .groups = "drop") %>%
  right_join(
    cap_check %>%
      select(
        Candidate_ID,
        Type,
        Capacity,
        Assigned_Peak_Total),
    by = "Candidate_ID") %>%
  mutate(
    Mean_Utilisation   = 100 * Mean_Daily_Demand / Capacity,
    Peak_Utilisation   = 100 * Peak_Daily_Demand / Capacity,
    Design_Utilisation = 100 * Assigned_Peak_Total / Capacity)

# Capacities
cap_plot1 <- cap_check %>%
  arrange(desc(Capacity)) %>%
  mutate(Candidate_ID = factor(Candidate_ID, levels = Candidate_ID)) %>%
  ggplot(
    aes(Candidate_ID, Capacity, fill = Type)) +
  geom_col(alpha = 0.85) +
  geom_hline(
    yintercept = mean(cap_check$Capacity),
    linetype = "dashed",
    colour = "firebrick") +
  scale_fill_manual(
    values = c("Service Point" = "#2C7BB6", "APL" = "#E69F00")) +
  labs(
    title = "Facility Capacities",
    subtitle = "Dashed line = Average",
    x = "Facility",
    y = "Capacity") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
print(cap_plot1)

# Utilization spread
cap_plot2 <- ilp_daily_facility %>%
  ggplot(aes(Design_Utilisation, fill = Type)) +
  geom_histogram(
    bins = 15,
    alpha = 0.75,
    colour = "white") +
  scale_fill_manual(
    values = c(
      "Service Point" = "#2C7BB6",
      "APL" = "#E69F00")) +
  labs(
    title = "Facility Peak Utilisation",
    x = "Utilisation (%)",
    y = "Number of Facilities") +
  theme_minimal()
print(cap_plot2)

# Capacity distribution
capacity_compare <- bind_rows(
  baseline_daily %>%
    transmute(
      Network = "Baseline",
      Capacity),
  ilp_daily_facility %>%
    transmute(
      Network = "ILP",
      Capacity))

cap_plot3 = ggplot(capacity_compare, aes(x = Capacity, fill = Network)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c(
    "Baseline" = "#E69F00",
    "ILP" = "#0072B2")) +
  labs(
    title = "Capacity Distribution Comparison",
    x = "Capacity",
    y = "Density") +
  theme_minimal()
print(cap_plot3)

# Equality Graph
baseline_lorenz <- baseline_daily %>%
  arrange(Capacity) %>%
  mutate(
    CumCap = cumsum(Capacity) / sum(Capacity),
    CumFac = row_number() / n())

ilp_lorenz <- ilp_daily_facility %>%
  arrange(Capacity) %>%
  mutate(
    CumCap = cumsum(Capacity) / sum(Capacity),
    CumFac = row_number() / n())

cap_plot4 = ggplot() +
  geom_line(
    data = baseline_lorenz,
    aes(CumFac, CumCap, colour = "Baseline"),
    linewidth = 1.2) +
  geom_line(
    data = ilp_lorenz,
    aes(CumFac, CumCap, colour = "ILP"),
    linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Model", values = c(
    "Baseline" = "#E69F00",
    "ILP" = "#0072B2")) +
  labs(
    title = "Capacity Equality Across Facilities",
    x = "Cumulative Share of Facilities",
    y = "Cumulative Share of Capacity") +
  theme_minimal()
print(cap_plot4)

# Map of capacity locations
cap_plot5 <- ggplot() +
  geom_rect(
    data = cbs %>%
      filter(
        !is.na(Population),
        Population > 0),
    aes(
      xmin = X - main_w/2,
      xmax = X + main_w/2,
      ymin = Y - main_h/2,
      ymax = Y + main_h/2),
    fill = "grey96",
    colour = "grey85",
    linewidth = 0.2) +
  geom_point(
    data = map_ilp,
    aes(
      x = X,
      y = Y,
      size = Capacity,
      colour = Type,
      shape = Type),
    alpha = 0.85) +
  scale_colour_manual(
    values = c(
      "Service Point" = "#2C7BB6",
      "APL" = "#E69F00")) +
  scale_shape_manual(
    values = c(
      "Service Point" = 16,
      "APL" = 17)) +
  scale_size_continuous(range = c(4, 18)) +
  coord_equal() +
  labs(
    title = "ILP Network Capacity",
    colour = NULL,
    shape = NULL,
    size = "Capacity") +
  theme_section14 +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom")
print(cap_plot5)

# --- Build a summary table for all existing SPs ---
sp_cost_lookup <- sp %>%
  transmute(
    Location_ID = as.character(`Location ID`),
    SP_Cost     = candidate_params$SP_Cost[
      match(as.character(Node_ID), candidate_params$Candidate_ID)])

closed_analysis <- baseline_daily %>%
  mutate(Location_ID = as.character(Location_ID)) %>%
  left_join(
    pickup_rate_obs %>%
      transmute(
        Location_ID          = as.character(`Location ID`),
        Observed_Pickup_Rate),
    by = "Location_ID") %>%
  left_join(sp_cost_lookup, by = "Location_ID") %>%
  mutate(
    Status = case_when(
      Location_ID %in% as.character(sp$`Location ID`[sp$Node_ID %in% kept_nodes])   ~ "Kept",
      Location_ID %in% as.character(sp$`Location ID`[sp$Node_ID %in% closed_nodes]) ~ "Closed",
      TRUE ~ "Other"),
    Cost_Per_Package = SP_Cost / Total_Annual_Demand) %>%
  filter(Status %in% c("Kept", "Closed"))

# Helper to get a sensible upper trim per variable
trim_upper <- function(x, p = 0.95) quantile(x, p, na.rm = TRUE)

p14_3a <- ggplot(closed_analysis, aes(x = Status, y = Mean_Daily_Demand, fill = Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("Kept" = pal_kept, "Closed" = pal_closed)) +
  coord_cartesian(ylim = c(0, trim_upper(closed_analysis$Mean_Daily_Demand))) +
  labs(
    title    = "Daily Demand: Closed vs Kept Service Points",
    subtitle = "Closed locations consistently handle fewer parcels per day",
    x        = NULL,
    y        = "Mean Daily Demand (packages)") +
  theme_section14 + theme(legend.position = "none")

p14_3b <- ggplot(closed_analysis, aes(x = Status, y = Observed_Pickup_Rate, fill = Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("Kept" = pal_kept, "Closed" = pal_closed)) +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, trim_upper(closed_analysis$Observed_Pickup_Rate))) +
  labs(
    title    = "Pickup Rate: Closed vs Kept Service Points",
    subtitle = "Closed locations show lower customer pickup activity",
    x        = NULL,
    y        = "Observed Pickup Rate") +
  theme_section14 + theme(legend.position = "none")

p14_3c <- ggplot(closed_analysis, aes(x = Status, y = Cost_Per_Package, fill = Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("Kept" = pal_kept, "Closed" = pal_closed)) +
  scale_y_continuous(labels = scales::dollar_format(prefix = "€")) +
  coord_cartesian(ylim = c(0, trim_upper(closed_analysis$Cost_Per_Package))) +
  labs(
    title    = "Fixed Cost per Package: Closed vs Kept Service Points",
    subtitle = "Closed locations carry a disproportionately high cost burden",
    x        = NULL,
    y        = "Annual Fixed Cost per Package") +
  theme_section14 + theme(legend.position = "none")

print(p14_3a)
print(p14_3b)
print(p14_3c)

# --- Print summary stats ---
cat("\n=== 14.3 Closed vs Kept Summary ===\n")
closed_analysis %>%
  group_by(Status) %>%
  summarise(
    N                    = n(),
    Avg_Daily_Demand     = round(mean(Mean_Daily_Demand), 1),
    Avg_Pickup_Rate      = scales::percent(mean(Observed_Pickup_Rate, na.rm = TRUE), 0.1),
    Avg_Cost_Per_Package = round(mean(Cost_Per_Package, na.rm = TRUE), 2),
    .groups = "drop") %>%
  print()

# ==============================================================================
# SECTION 15: COST DISTRIBUTION ANALYSIS 
# ==============================================================================

# --- Setup cost: match each existing SP node to its candidate_params entry ---
current_setup_costs <- sp %>%
  mutate(Node_ID = as.character(Node_ID)) %>%
  left_join(
    candidate_params %>%
      transmute(Node_ID = as.character(Candidate_ID), SP_Cost),
    by = "Node_ID") %>%
  mutate(SP_Cost = coalesce(SP_Cost, mean(candidate_params$SP_Cost, na.rm = TRUE)))

# --- Capacity cost: based on peak demand / 1.02 rule, same as baseline_daily ---
# baseline_daily already has Capacity = Peak_Daily_Demand / 1.02 per SP
current_cap_cost <- baseline_daily %>%
  mutate(
    Location_ID    = as.character(Location_ID),
    Annual_Cap_Cost = Capacity * CapacityCost)

# --- Delivery cost: reconstruct for the CURRENT network
# Each demand node is currently served by its closest existing SP (Closest_SP_ID from sub_df).
# Delivery cost = (1 - pickup_rate) × total_demand × DeliveryCost × truck_km to that SP

current_delivery_cost <- node_demand_params %>%
  mutate(Node_ID = as.character(Node_ID)) %>%
  left_join(
    sub_df %>%
      filter(!is.na(Node_ID), !is.na(Closest_SP_ID)) %>%
      mutate(Node_ID = as.character(Node_ID)) %>%
      group_by(Node_ID) %>%
      slice(1) %>%
      ungroup() %>%
      select(Node_ID, Closest_SP_ID),
    by = "Node_ID") %>%
  left_join(
    node_corrections %>% mutate(Node_ID = as.character(Node_ID)),
    by = "Node_ID") %>%
  mutate(
    Stage2_Correction = coalesce(Stage2_Correction, 0),
    # Walking distance to current SP (already in sub_df via Node_ID → Closest_SP_ID)
    Base_Pickup_Rate  = a_fit * exp(-b_fit * coalesce(
      sub_df$Walking_Distance[match(Node_ID, as.character(sub_df$Node_ID))], 0)),
    Effective_PR      = pmin(pmax(Base_Pickup_Rate + Stage2_Correction, 0), 1),
    # Truck distance: use the delivery_pairs table to find truck km to current SP
    Truck_km = coalesce(
      delivery_pairs$Truck_Dist_m[
        match(paste(Node_ID, Closest_SP_ID),
              paste(as.character(delivery_pairs$Node_ID),
                    as.character(delivery_pairs$Candidate_ID)))] / 1000,
      0),
    Current_Delivery_Cost = (1 - Effective_PR) * Total_Demand * DeliveryCost * Truck_km)

# --- Combine into one per-SP cost table for the current network ---
current_costs_per_sp <- current_setup_costs %>%
  select(Node_ID, SP_Cost) %>%
  left_join(
    current_cap_cost %>% select(Location_ID, Annual_Cap_Cost),
    by = c("Node_ID" = "Location_ID")) %>%
  left_join(
    current_delivery_cost %>%
      group_by(Closest_SP_ID) %>%
      summarise(Annual_Delivery_Cost = sum(Current_Delivery_Cost, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(Closest_SP_ID = as.character(Closest_SP_ID)),
    by = c("Node_ID" = "Closest_SP_ID")) %>%
  mutate(
    Annual_Cap_Cost      = coalesce(Annual_Cap_Cost, 0),
    Annual_Delivery_Cost = coalesce(Annual_Delivery_Cost, 0),
    Total_Annual_Cost    = SP_Cost + Annual_Cap_Cost + Annual_Delivery_Cost,
    Network              = "Current",
    Type                 = "Service Point")

ilp_delivery_per_facility <- delivery_assignment %>%
  group_by(Candidate_ID) %>%
  summarise(Annual_Delivery_Cost = sum(Delivery_Cost, na.rm = TRUE), .groups = "drop") %>%
  mutate(Candidate_ID = as.character(Candidate_ID))

# Pickup cost savings: negative delivery cost credited back per facility
ilp_pickup_savings_per_facility <- pickup_assignment %>%
  mutate(
    Pickup_Saving = Pickup_Rate * Total_Demand * DeliveryCost *
      coalesce(Truck_km_eff, 0)) %>%
  group_by(Candidate_ID) %>%
  summarise(Annual_Pickup_Saving = sum(Pickup_Saving, na.rm = TRUE), .groups = "drop") %>%
  mutate(Candidate_ID = as.character(Candidate_ID))

ilp_costs_per_facility <- open_facilities %>%
  mutate(Candidate_ID = as.character(Candidate_ID)) %>%
  left_join(ilp_delivery_per_facility,    by = "Candidate_ID") %>%
  left_join(ilp_pickup_savings_per_facility, by = "Candidate_ID") %>%
  mutate(
    SP_Cost              = Facility_Cost,
    Annual_Delivery_Cost = coalesce(Annual_Delivery_Cost, 0) -
      coalesce(Annual_Pickup_Saving, 0),
    Total_Annual_Cost    = SP_Cost + Annual_Cap_Cost + Annual_Delivery_Cost,
    Network              = "ILP",
    Node_ID              = Candidate_ID)

cost_compare <- bind_rows(
  current_costs_per_sp %>%
    select(Node_ID, Type, Network, SP_Cost,
           Annual_Cap_Cost, Annual_Delivery_Cost, Total_Annual_Cost),
  ilp_costs_per_facility %>%
    select(Node_ID, Type, Network, SP_Cost,
           Annual_Cap_Cost, Annual_Delivery_Cost, Total_Annual_Cost))

# Summary table
cost_summary <- cost_compare %>%
  group_by(Network, Type) %>%
  summarise(
    N_Facilities         = n(),
    Total_Setup_Cost     = sum(SP_Cost,              na.rm = TRUE),
    Total_Cap_Cost       = sum(Annual_Cap_Cost,      na.rm = TRUE),
    Total_Delivery_Cost  = sum(Annual_Delivery_Cost, na.rm = TRUE),
    Total_Annual_Cost    = sum(Total_Annual_Cost,    na.rm = TRUE),
    Avg_Cost_Per_Facility = mean(Total_Annual_Cost,  na.rm = TRUE),
    .groups = "drop")

# Extra graphs for the repport + ------------------------

cat("\n=== SECTION 15: COST DISTRIBUTION SUMMARY ===\n")
print(cost_summary)

theme_15 <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey45"),
    legend.position  = "top",
    panel.grid.minor = element_blank())

pal_network <- c("Current" = "#E07B39", "ILP" = "#2C7BB6")

# --- Plot 15.1: Setup / fixed cost distribution ---
p15_setup <- ggplot(cost_compare,
                    aes(x = SP_Cost, fill = Network)) +
  geom_density(alpha = 0.55, colour = NA) +
  geom_rug(aes(colour = Network), alpha = 0.5, linewidth = 0.4) +
  scale_fill_manual(values  = pal_network, name = NULL) +
  scale_colour_manual(values = pal_network, name = NULL) +
  scale_x_continuous(labels = scales::dollar_format(prefix = "€")) +
  labs(
    title    = "Distribution of Annual Setup Costs per Facility",
    subtitle = "ILP shifts weight toward lower-cost locations; APLs pull the tail inward",
    x        = "Annual Setup Cost (€)",
    y        = "Density") +
  theme_15

print(p15_setup)

cost_stacked <- cost_compare %>%
  group_by(Network) %>%
  summarise(
    `Setup Cost`    = sum(SP_Cost,              na.rm = TRUE),
    `Capacity Cost` = sum(Annual_Cap_Cost,      na.rm = TRUE),
    `Delivery Cost` = sum(Annual_Delivery_Cost, na.rm = TRUE),
    .groups = "drop") %>%
  tidyr::pivot_longer(
    cols      = c(`Setup Cost`, `Capacity Cost`, `Delivery Cost`),
    names_to  = "Cost_Component",
    values_to = "Amount") %>%
  mutate(Cost_Component = factor(Cost_Component,
                                 levels = c("Setup Cost", "Capacity Cost", "Delivery Cost")))

# NOTE: These ILP cost components are hardcoded from the solver output of Step 13
# (setup = fixed SP + APL costs, capacity, net delivery after pickup savings).
# They are not recomputed here. Re-running Step 13 with different parameters
# requires updating these three values, or deriving them from true_fixed + true_apl,
# true_cap and net_delivery.

ilp_cost_components <- tibble(
  Network        = "ILP",
  Cost_Component = c("Setup Cost", "Capacity Cost", "Delivery Cost"),
  Amount         = c(1015383, 319740, 817814)
)

# Current network totals aggregated from per-SP reconstruction
current_cost_components <- cost_compare %>%
  filter(Network == "Current") %>%
  summarise(
    `Setup Cost`    = sum(SP_Cost,              na.rm = TRUE),
    `Capacity Cost` = sum(Annual_Cap_Cost,      na.rm = TRUE),
    `Delivery Cost` = sum(Annual_Delivery_Cost, na.rm = TRUE)) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = "Cost_Component",
                      values_to = "Amount") %>%
  mutate(Network = "Current")

cost_stacked <- bind_rows(current_cost_components, ilp_cost_components) %>%
  mutate(
    Network        = factor(Network, levels = c("Current", "ILP")),
    Cost_Component = factor(Cost_Component,
                            levels = c("Setup Cost", "Capacity Cost", "Delivery Cost")))

# Network totals for annotation
network_totals <- cost_stacked %>%
  group_by(Network) %>%
  summarise(Total = sum(Amount), .groups = "drop")

p15_stacked <- ggplot(cost_stacked,
                      aes(x = Network, y = Amount, fill = Cost_Component)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  # Per-component labels inside bars
  geom_text(
    aes(label = paste0("€", scales::comma(round(Amount / 1e3)), "k")),
    position  = position_stack(vjust = 0.5),
    size      = 3.2, colour = "white", fontface = "bold") +
  # Total cost label above each bar
  geom_text(
    data     = network_totals,
    aes(x = Network, y = Total,
        label = paste0("Total: €", scales::comma(round(Total / 1e3)), "k")),
    inherit.aes = FALSE,
    vjust    = -0.5, size = 3.5, fontface = "bold", colour = "grey20") +
  scale_fill_manual(
    values = c(
      "Setup Cost"    = "#1B4F8A",
      "Capacity Cost" = "#4A90D9",
      "Delivery Cost" = "#A8C8F0"),
    name = NULL) +
  scale_y_continuous(
    labels = scales::dollar_format(prefix = "€"),
    expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Total Network Cost Breakdown: Current vs ILP",
    subtitle = "Stacked by cost component — ILP values match solver output exactly",
    x        = NULL,
    y        = "Total Annual Cost (€)") +
  theme_15

print(p15_stacked)

p15_individual <- ggplot(cost_stacked %>% filter(Network == "Current"),
                         aes(x = Cost_Component, y = Amount, fill = Cost_Component)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = paste0("€", scales::comma(round(Amount / 1e3)), "k")),
    vjust  = -0.5,
    size   = 3.2, colour = "grey20", fontface = "bold") +
  scale_fill_manual(
    values = c(
      "Setup Cost"    = "#1B4F8A",
      "Capacity Cost" = "#4A90D9",
      "Delivery Cost" = "#A8C8F0"),
    name = NULL) +
  scale_y_continuous(
    labels = scales::dollar_format(prefix = "€"),
    expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Current Network Cost Breakdown",
    subtitle = "Individual cost components for the existing network",
    x        = NULL,
    y        = "Total Annual Cost (€)") +
  theme_15 +
  theme(legend.position = "none")

print(p15_individual)


p_utilisation <- baseline_daily %>%
  arrange(Mean_Utilisation) %>%
  mutate(Location_ID = factor(Location_ID, levels = Location_ID)) %>%
  ggplot(aes(x = Location_ID)) +
  geom_segment(
    aes(xend = Location_ID, y = Mean_Utilisation, yend = Peak_Utilisation),
    colour = "grey70", linewidth = 0.6
  ) +
  geom_point(
    aes(y = Peak_Utilisation),
    colour = "#1B4F8A", size = 2.8, shape = 16
  ) +
  geom_point(
    aes(y = Mean_Utilisation),
    colour = "#A8C8F0", size = 2.8, shape = 16
  ) +
  geom_hline(
    yintercept = mean(baseline_daily$Mean_Utilisation),
    linetype = "dashed", colour = "#A8C8F0", linewidth = 0.6
  ) +
  geom_hline(
    yintercept = mean(baseline_daily$Peak_Utilisation),
    linetype = "dashed", colour = "#1B4F8A", linewidth = 0.6
  ) +
  annotate("text",
           x = 35.5, y = mean(baseline_daily$Mean_Utilisation),
           label = paste0("avg. mean ", round(mean(baseline_daily$Mean_Utilisation)), "%"),
           hjust = 0, size = 2.8, colour = "#A8C8F0"
  ) +
  annotate("text",
           x = 35.5, y = mean(baseline_daily$Peak_Utilisation),
           label = paste0("avg. peak ", round(mean(baseline_daily$Peak_Utilisation)), "%"),
           hjust = 0, size = 2.8, colour = "#1B4F8A"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title    = "Service Point Capacity Utilisation",
    subtitle = "Light blue = mean daily utilisation  ·  Dark blue = peak daily utilisation  ·  Ordered by mean",
    x        = "Service Point",
    y        = "Capacity Utilisation (%)"
  ) +
  theme_15 +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank()
  )

print(p_utilisation)

p_demand_dist <- baseline_daily %>%
  arrange(desc(Mean_Daily_Demand)) %>%
  mutate(Location_ID = factor(Location_ID, levels = Location_ID)) %>%
  ggplot(aes(x = Location_ID, y = Mean_Daily_Demand)) +
  geom_col(fill = "#1B4F8A", width = 0.7, alpha = 0.85) +
  geom_hline(
    yintercept = mean(baseline_daily$Mean_Daily_Demand),
    linetype = "dashed", colour = "#A8C8F0", linewidth = 0.7
  ) +
  annotate("text",
           x = 35.5, y = mean(baseline_daily$Mean_Daily_Demand),
           label = paste0("avg. ", round(mean(baseline_daily$Mean_Daily_Demand)), " pkg/day"),
           hjust = 0, size = 2.8, colour = "#A8C8F0"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = scales::comma
  ) +
  labs(
    title    = "Daily Demand per Service Point — Current Network",
    subtitle = "Average daily packages handled per location  ·  Ordered highest to lowest  ·  Dashed = network average",
    x        = "Service Point",
    y        = "Mean Daily Demand (packages)"
  ) +
  theme_15 +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank()
  )

print(p_demand_dist)

# The APL sensitivity analysis in the report (APL setup cost, adoption rate,
# home delivery cost, demand intensity) is not included in this script.
# It was run by re-solving the Step 13 ILP with modified parameters.
# See the report for the method and results.