# User-defined functions for HI_Coastal project

mapplot <- function(x) {
  ggplot(mapdat) +
    geom_sf(aes(geometry = geometry)) + 
    aes(fill = as.factor(.data[[x]]))+ 
    theme_void() + 
    theme(title = element_text(size = 14),
          plot.subtitle = element_text(4),
          plot.caption = element_text(size = 12, hjust = 1),
          legend.text = element_text(size = 10),
          legend.position = c(0.9, 0.8),
          text = element_text(size = 10)) +
    scale_fill_manual(breaks = c(0,1), values = c("white", "#2AA8AA"), 
                      labels = c("Not selected", "Selected"), na.value = "light grey")+
    labs(fill = "LSOA group",
         title = "Location of 'coastal' LSOAs",
         subtitle = paste0("Definition: ",
                           definitions$definition[definitions$var == {{x}}] )) + 
    annotate("text", x=-5, y=53, label = paste0("Total selected LSOAs:  ",
                                                scales::comma(sum(mapdat[[x]][mapdat[[x]] == 1])),
                                                "\nSelected population:  ",
                                                scales::comma(sum(mapdat$total[mapdat[[x]] == 1])))) +
    coord_sf()
  ggsave(paste0(Sys.getenv("out"), "Coastal_map_", x, ".jpg"), dpi = 600, width = 210, height = 297, units = "mm")
}

plot_fun <- function(mod, by_out) {
  rs <- by_out |>
    select(paste0("resids_",mod)) |>
    unnest_legacy()
  
  pred <- by_out |>
    select(paste0("pred_", mod)) |>
    unnest_legacy()
  
  f <- inner_join(rs, pred, by=c("outcome","LSOA11", "RGN"))
  f |>
    ggplot(aes(x = pred, y = resid)) +
    geom_point() + 
    facet_wrap(~ outcome) + 
    geom_hline(yintercept = 0, colour = "blue")
  ggsave(paste0(Sys.getenv("interim"),"/residsplot_",mod,".png"))
}

glance_fun <- function(gmod, by_out) {
  by_out |>
    mutate(glance = map(.data[[gmod]], broom::glance)) |>
    unnest(glance) |>
    select(outcome, p.value, BIC, adj.r.squared) |>
    mutate(model = gmod) |>
    rename(model_p_value = p.value)
}



addmodels <- function(col, by_out){
  by_out <- by_out |>
    mutate({{col}} := map(data,  models[modelsets== col]))
}

addmodels2 <- function(col, by_out){
  by_out <- by_out |>
    mutate({{col}} := map(data,  models[modelsets== col]))
}