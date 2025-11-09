#!/usr/bin/env Rscript
library(scales)
library(tidyverse)
library(vroom)
library(ggforce)


# 1 ── Read & basic preprocessing -------------------------------------------------

df0 <- read_tsv("6-size_speed-HQ661k-long_fmt.tsv")

df1 <- df0 %>%
    rename(cpu  = `User time (s) - XZ-decompression and search`, size = `Batch size (Gb)`, index = `Candidate index`) %>%
    mutate(cpu = cpu / 60)

# 2 ── One point per index ----------------------

single_idx <- df1 %>%
    group_by(index) %>%
    summarise(
        cpu  = sum(cpu),
        # total minutes
        size = sum(size),
        # total GB
        .groups = "drop"
    ) %>%
    ungroup()


# 3 ── Non-dominated rows for each batch variant  ---------------

df2 <- df1 %>%
    group_by(experiment) %>%
    arrange(size, cpu, .by_group = TRUE) %>%
    mutate(dominated = cummin(cpu) < cpu) %>%
    filter(!dominated) %>%             # keep only non-dominated
    select(experiment, cpu, size) %>%
    group_split()


# 4 ── Global Pareto frontier by combining indices ------------------------------

frontier <- tibble(cpu = 0, size = 0)     # start at origin

for (bv in df2) {
    # --------------------- patched block ------------------- #
    frontier <- frontier            %>%
        rename(
            cpu_f  = cpu,
            # <- temporary tag
            size_f = size
        )         %>%
        crossing(rename(
            bv,
            cpu_b  = cpu,
            # <- tag on batch variant
            size_b = size
        )) %>%
        transmute(cpu  = cpu_f  + cpu_b, size = size_f + size_b)     # aggregate
    # ------------------------------------------------------ #

    frontier <- frontier            %>%
        distinct()                     %>%
        arrange(size, cpu)             %>%
        mutate(dominated = cummin(cpu) < cpu) %>%
        filter(!dominated)
}

# 5 ── Plot ----------------------------------------------------------------------

#shape and color for all indexes
shape_vals <- c(0, 1, 2, 3, 4, 5, 6)
col_vals <- hue_pal()(6)[c(4,1,3,2,5,6,4)]

p <- ggplot(frontier |> arrange(size), aes(size, cpu)) +
    geom_line() +
    geom_point(
        # cumulative single-index points
        data  = single_idx,
        aes(size, cpu, color = index, shape = index),
        size   = 3
    ) +
    scale_shape_manual(values = shape_vals) +
    scale_color_manual(values = col_vals) +
    labs(x = "Disk [GB]", y = "CPU time [min]")  +
    ylim(0, NA) +
    xlim(0, NA) +
    theme_bw()

p +
    facet_zoom(
        x = size <= 150,
        y = cpu  <= 3110,
        zoom.size = 1,     # inset same height as main panel
        show.area = FALSE   # grey rectangle marks zoomed region
    )

h <- 7
w <- 17
u <- "cm"

ggsave(
  "9-pareto-optimal-decompression_included.pdf",
  width = w,
  height = h,
  unit = u
)

