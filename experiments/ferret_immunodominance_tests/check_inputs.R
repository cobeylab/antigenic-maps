## Visualize 2D inputs

even = read_rds('inputs/2D_even_immunodominance_inputs.rds')
skewed = read_rds('inputs/2D_skewed_immunodominance_inputs.rds')
E1 = read_rds('inputs/2D_E1_immunodominance_inputs.rds')
E2 = read_rds('inputs/2D_E2_immunodominance_inputs.rds')
E3 = read_rds('inputs/2D_E3_immunodominance_inputs.rds')

cowplot::plot_grid(
even$ag_ab_coords %>%
  ggplot() +
  geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .2) +
  geom_point(aes(x = c1_Ag, y = c2_Ag)) +
  ggtitle('even'),
skewed$ag_ab_coords %>%
  ggplot() +
  geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .2) +
  geom_point(aes(x = c1_Ag, y = c2_Ag)) +
  ggtitle('skewed')
)

cowplot::plot_grid(
E1$ag_ab_coords %>%
  ggplot() +
  geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .5) +
  geom_point(aes(x = c1_Ag, y = c2_Ag)) +
  ggtitle('E1'),
E2$ag_ab_coords %>%
  ggplot() +
  geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .5) +
  geom_point(aes(x = c1_Ag, y = c2_Ag)) +
  ggtitle('E2'),
E3$ag_ab_coords %>%
  ggplot() +
  geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .5) +
  geom_point(aes(x = c1_Ag, y = c2_Ag)) +
  ggtitle('E3')
)



## Visualize 3D inputs
even = read_rds('inputs/3D_even_immunodominance_inputs.rds')
skewed = read_rds('inputs/3D_skewed_immunodominance_inputs.rds')
E1 = read_rds('inputs/3D_E1_immunodominance_inputs.rds')
E2 = read_rds('inputs/3D_E2_immunodominance_inputs.rds')
E3 = read_rds('inputs/3D_E3_immunodominance_inputs.rds')

cowplot::plot_grid(
  even$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .2) +
    geom_point(aes(x = c1_Ag, y = c2_Ag)) +
    ggtitle('even'),
  skewed$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .2) +
    geom_point(aes(x = c1_Ag, y = c2_Ag)) +
    ggtitle('skewed'),
  even$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c2_Ab, y = c3_Ab, color = epitope), alpha = .2) +
    geom_point(aes(x = c2_Ag, y = c3_Ag)) +
    ggtitle('even'),
  skewed$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c2_Ab, y = c3_Ab, color = epitope), alpha = .2) +
    geom_point(aes(x = c2_Ag, y = c3_Ag)) +
    ggtitle('skewed')
)


cowplot::plot_grid(
  E1$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .5) +
    geom_point(aes(x = c1_Ag, y = c2_Ag)) +
    ggtitle('E1'),
  E2$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .5) +
    geom_point(aes(x = c1_Ag, y = c2_Ag)) +
    ggtitle('E2'),
  E3$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), alpha = .5) +
    geom_point(aes(x = c1_Ag, y = c2_Ag)) +
    ggtitle('E3'),
  
  E1$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c3_Ab, y = c2_Ab, color = epitope), alpha = .5) +
    geom_point(aes(x = c3_Ag, y = c2_Ag)) +
    ggtitle('E1'),
  E2$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c3_Ab, y = c2_Ab, color = epitope), alpha = .5) +
    geom_point(aes(x = c3_Ag, y = c2_Ag)) +
    ggtitle('E2'),
  E3$ag_ab_coords %>%
    ggplot() +
    geom_point(aes(x = c3_Ab, y = c2_Ab, color = epitope), alpha = .5) +
    geom_point(aes(x = c3_Ag, y = c2_Ag)) +
    ggtitle('E3')
  
)
