library(dplyr)
library(stringr)

tissue_list <- c(
  'brain', 'thyroid', 'thymus', 'heart', 'lung', 'liver', 'spleen', 'pancreas', 'stomach', 'small', 'large', 'kidney', 'adrenal', 'muscle', 'adipose', 
  'bone', 'testis', 'sperm'
)

complete_tissue_list <- c(
  'Brain', 'Thyroid', 'Thymus', 'Heart', 'Lung', 'Liver', 'Spleen', 'Pancreas', 'Stomach', 'Small Intestine', 'Large Intestine', 'Kidney', 'Adrenal Gland', 'Muscle', 'Adipose Tissue', 
  'Bone Marrow', 'Testis', 'Sperm'
)

generate_colorbar <- function(vmax, cmap, revert = FALSE, height = 400, width = 100) {
  library(ggplot2)
  library(RColorBrewer)
  
  colors <- RColorBrewer::brewer.pal(9, cmap)
  if(revert) colors <- rev(colors)
  col_gradient <- colorRampPalette(colors)(100)
  
  df <- data.frame(
    ratio = seq(0, vmax, length.out = 100),
    x = factor(1)
  )
  
  safe_vmax <- gsub("\\.", "_", format(round(vmax, 2), nsmall = 2))
  filename <- sprintf("colorbar_%s_%s.png", cmap, safe_vmax)
  if(revert) filename <- paste0("reverted_", filename)
  output_path <- file.path("user", filename)
  
  dir.create("user", showWarnings = FALSE, recursive = TRUE)
  
  p <- ggplot(df, aes(x = x, y = ratio, fill = ratio)) +
    geom_tile(width = 0.2, height = 1) + 
    scale_fill_gradientn(colors = col_gradient, limits = c(0, vmax)) +
    scale_y_reverse(
      breaks = seq(0, vmax, length.out = 5),
      labels = function(x) sprintf("%.2f", x),
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 0, color = "black"),
      axis.ticks.y = element_line(color = "#00000000"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 15, 2, 5)
    )
  
  ggsave(
    output_path,
    plot = p,
    width = width/100,
    height = height/100,
    dpi = 100,
    units = "in"
  )
  
  return(output_path)
}

plot_mRNA_tissues <- function(mRNA, data_info, data_array, colname='gene_info', bin_length=10, vmax=0.12, cmap='Reds', title='default', min_length=10, max_length=300, y_tick=TRUE,color_revert=FALSE) {
  # # print(mRNA, colname)
  # print(data_info)
  corr_df <- data_info %>% filter(!!as.name(colname) == mRNA)

  # corr_df <- data_info %>% filter(!!as.name(colname) == mRNA)
  corr_idx <- corr_df$raw_idx 
  
  corr_tissue <- corr_df$sample_info
  
  result_df <- list()
  
  for (tissue in tissue_list) {
    if (tissue %in% corr_tissue) {
      pos_idx <- which(corr_tissue == tissue)
      polya_list <- c()
      
      for (pi in pos_idx) {
        pos <- corr_idx[pi] 
        tmp_list <- data_array[[pos+1]]
        polya_list <- c(polya_list, tmp_list[tmp_list >= min_length & tmp_list < max_length])
      }
      
      ratio <- get_one_sample_dfs(polya_list, bin_length)$count
      result_df[[tissue]] <- ratio
    } else {
      ratio <- get_one_sample_dfs(c(), bin_length)$count
      result_df[[tissue]] <- ratio
    }
  }
  
  result_df <- as.data.frame(do.call(rbind, result_df))
  rownames(result_df) <- tissue_list
  # print(result_df)

  # Generate color palette and reverse it if needed
  color_palette <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = cmap))(100)
  if (color_revert) {
    color_palette <- rev(color_palette)  # Reverse the colormap if the checkbox is checked
  }

  col_labels <- rep("", ncol(result_df))
  valid_breaks <- unique(c(min_length, seq(50, max_length, by = 50)))
  valid_breaks <- valid_breaks[valid_breaks >= min_length & valid_breaks <= max_length]
  col_indices <- (valid_breaks - min_length) / bin_length + 1 
  col_labels[col_indices] <- paste0("\n", as.character(valid_breaks))

  
  ph <- pheatmap(
    result_df,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    breaks = seq(0, vmax, length.out = 101),
    fontsize_row = 15,
    display_numbers = FALSE,
    border_color = NA,
    legend = FALSE,
    show_colnames = TRUE,
    labels_col = col_labels,
    labels_row = complete_tissue_list,
    angle_col = 0,
    fontsize_col = 10,
    show_rownames = y_tick,
    main = mRNA,
    silent = TRUE,
    margins = c(12, 12)
  )

  # axis
  offset         <- 0.8 * 5 / bin_length
  tick_positions <- offset + (valid_breaks - min_length) / bin_length * 0.655
  # print(tick_positions)
  
  ph_grob <- ph[[4]]

  upper_boundary <- -11.9 
  
  axis_grob <- grid::grobTree(
    grid::segmentsGrob(
      x0 = unit(tick_positions / ncol(result_df), "npc"),
      x1 = unit(tick_positions / ncol(result_df), "npc"),
      y0 = unit(upper_boundary, "npc"),
      y1 = unit(-12.1, "npc"),
      gp = grid::gpar(lwd = 1, col = "black")
    ),
    grid::segmentsGrob(
      x0 = unit(tick_positions[[1]] / ncol(result_df), "npc"),
      x1 = unit(tick_positions[[length(tick_positions)]] / ncol(result_df), "npc"),
      y0 = unit(upper_boundary, "npc"),
      y1 = unit(upper_boundary, "npc"),
      gp = grid::gpar(lwd = 1, col = "black")
    ),
    # grid::textGrob(
    #   label = valid_breaks,
    #   x = unit(tick_positions / ncol(result_df), "npc"),
    #   y = unit(-11.6, "cm"),
    #   just = c("center", "top"),
    #   gp = grid::gpar(fontsize = 8)
    # ),
    vp = grid::viewport(y = 0.93, height = 0.07) 
  )
  
  combined_grob <- grid::grobTree(
    ph_grob,
    axis_grob
  )
  
  return(combined_grob)
}

get_one_sample_dfs <- function(polya_list, bin_length=10, max_length=300) {
  d1 <- data.frame(length=seq(from=ceiling(10/bin_length), to=ceiling(max_length/bin_length) - 1, by=1))

  polya_list <- na.omit(polya_list)
  
  if (length(polya_list) == 0) {
    d <- d1
    d$count <- 0
  } else {
    x <- table(cut(polya_list, breaks=seq(1, 500, by=bin_length), include.lowest=TRUE))
    d <- data.frame(length=seq_along(x), count=as.numeric(x))
    d <- merge(d1, d, by='length', all.x=TRUE)
    d$count[is.na(d$count)] <- 0
    d$count <- d$count / sum(d$count)
  }
  return(d)
}

save_plot_as_png <- function(plot_obj, filename, width = 700, height = 800) {
  png(filename, width = width, height = height, res = 150)
  grid::grid.newpage()
  grid::grid.draw(plot_obj)
  dev.off()
}
