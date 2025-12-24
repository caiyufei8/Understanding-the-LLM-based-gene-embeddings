library(this.path, lib = '~/myRlibs')
script_dir <- this.dir()
setwd(script_dir)
library(purrr, lib = '~/myRlibs')
# install.packages("BiocManager",lib="~/myRlibs",repos='https://cran.r-project.org')
# library(BiocManager,lib='~/myRlibs')
# BiocManager::install("msigdb",lib="~/myRlibs")
# BiocManager::install("GSEABase",lib="~/myRlibs")
# BiocManager::install("ExperimentHub",lib="~/myRlibs")
# BiocManager::install("qusage",lib="~/myRlibs")
library(qusage, lib = '~/myRlibs')
library(msigdb, lib = '~/myRlibs')
library(GSEABase, lib = '~/myRlibs')
library(ExperimentHub, lib = '~/myRlibs')
# install.packages("farver", lib = '~/myRlibs')
# install.packages("RColorBrewer", lib = '~/myRlibs')
# install.packages("ggplot2", lib = '~/myRlibs')
# install.packages("labeling", lib = '~/myRlibs')
library(farver, lib = '~/myRlibs')
library(RColorBrewer, lib = '~/myRlibs')
library(ggplot2, lib = '~/myRlibs')
library(labeling, lib = '~/myRlibs')
# install.packages("ggrepel", lib = '~/myRlibs')
library(ggrepel, lib = '~/myRlibs')
# install.packages("cowplot", lib = '~/myRlibs')
library(cowplot, lib = '~/myRlibs')
# install.packages("RColorBrewer", lib = '~/myRlibs')
library(RColorBrewer, lib = '~/myRlibs')
# install.packages("plotrix", lib = '~/myRlibs')
library(plotrix, lib = '~/myRlibs')
library(dplyr, lib = '~/myRlibs')
library(tidyr, lib = '~/myRlibs')
# install.packages("tibble", lib = '~/myRlibs')
library(tibble, lib = '~/myRlibs')
# install.packages("ggtext", lib = '~/myRlibs')
library(ggtext, lib = '~/myRlibs')
# install.packages("ggnewscale", lib = '~/myRlibs')
library(ggnewscale, lib = '~/myRlibs')
plot_h_barplot<-function(bar_names,bar_values){
  df <- data.frame(
    name = factor(bar_names, levels = rev(bar_names)), # Reverse levels for correct order in plot
    value = bar_values
  )
  
  df <- df %>%
    arrange(desc(value))
  df$name <- factor(df$name, levels = rev(df$name))
  
  # Create the bar plot
  p<-ggplot(df, aes(y = name, x = value)) +
    geom_col(fill = "skyblue") + # Bars for bar_values
    geom_text(aes(label = value), hjust = -0.5, color = "darkred", size=5) + # Value on the bar
     labs(
      # title = "Bar Plot with Additional Values",
      x = "Number of times detected as significant",
      y = "Hallmark pathways"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(face = "plain", size = 18), # X 軸和 Y 軸標籤加粗
      axis.text.y = element_text(face = "plain", size = 12), # Bold names on y-axis
      panel.grid.major.y = element_blank(), # Remove horizontal grid lines
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )+
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  return(p)
}
plot_h_barplot <- function(bar_names, bar_values) {
  df <- data.frame(
    name = factor(bar_names, levels = rev(bar_names)), # Reverse levels for correct order in plot
    value = bar_values
  )
  
  df <- df %>%
    arrange(desc(value))
  df$name <- factor(df$name, levels = rev(df$name))
  
  # Create the bar plot
  p <- ggplot(df, aes(y = name, x = value)) +
    geom_col(fill = "skyblue", width = 0.7) + # Slightly narrower bars for more label space
    geom_text(aes(label = value), hjust = -0.5, color = "darkred", size = 30 / .pt) + # Value labels, size in points
    theme_classic() +
    theme(
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_blank(), # Remove y-axis title
      axis.text.y = element_text(face = "plain", size = 30, margin = margin(r = 10), lineheight = 0.8), # Add margin and adjust lineheight for y-axis text
      axis.line.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.major.y = element_blank(), # Remove horizontal grid lines
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
    coord_cartesian(xlim = c(0, max(bar_values) * 1.2)) # Double the x-axis range
  return(p)
}
plot_h_barplot<-function(bar_names,bar_values){
  df <- tibble::tibble(
    pathway = bar_names,
    score = bar_values
  )
  
  ## 讓 pathway 按照數值從大到小排列（保持原圖的排序）
  df <- df %>%
    arrange(desc(score)) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  ## 2. 設定 x 軸刻度 (用來畫「在橫條後面的格線」) ----------------------------
  x_breaks <- pretty(df$score)  # 例如 0, 50, 100, 150, ...
  
  ## 3. 繪圖 ---------------------------------------------------------------
  p <- ggplot(df, aes(x = pathway, y = score)) +
    ## 先畫在「橫條後面」的格線（用 yintercept；等會再 coord_flip）
    geom_hline(yintercept = x_breaks,
               linetype = "dashed",
               linewidth = 0.3,
               colour = "#DDDDDD") +
    ## 再畫橫條（實際上先是直條，最後 flip）
    geom_col(width = 0.7,
             fill = "#77AADD") +  # 柔和一點的藍色，比較精緻
    ## 在橫條右側標出數值
    geom_text(aes(label = score),
              hjust = -0.15,            # 往右挪一點，跑到 bar 外面
              size  = 10,              # 字稍微大一點
              family = "Helvetica",
              color = "#A63A2A") +   # 字體可改成你系統有的，如 "Arial"
    ## 轉成橫條圖
    coord_flip() +
    ## x 軸（原本的 score）刻度與留白
    scale_y_continuous(
      breaks = x_breaks,
      expand = expansion(mult = c(0, 0.12))  # 右側留一點空間給數字
    ) +
    ## 主體風格：簡潔期刊風
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    theme(
      axis.title.x  = element_blank(),
      axis.title.y  = element_blank(),
      axis.text.y   = element_text(size = 24, hjust = 1, colour = "black"),
      axis.text.x   = element_text(size = 24, hjust = 1, colour = "black"),
      panel.grid.major = element_blank(),    # 關掉預設格線，避免跟自訂的重疊
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 24, 5.5, 12) # 右邊多一點空間，數字不會被切掉
    )
  return(p)
}
nature_barplot <- function(bin_names, bin_values) {
  df <- data.frame(
    names = factor(bin_names, levels = bin_names),
    values = bin_values
  )
  ggplot(df, aes(x = names, y = values)) +
    geom_bar(stat = "identity", fill = "lightblue", width = 0.7) +
    geom_text(aes(label = values), vjust = -0.5, color = "#C0392B", size = 40) +
    theme_classic() +  # 基礎字體大小14pt，無襯線
    theme(
      axis.title.y = element_blank(),  # 隱藏y軸標題
      axis.text.y = element_blank(),   # 隱藏y軸刻度標籤
      axis.ticks.y = element_blank(),  # 隱藏y軸刻度線
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),  # 隱藏x軸標題
      axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),    # 移除網格線
      plot.title = element_blank(),     # 移除圖表標題
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}
nature_barplot <- function(bin_names, bin_values) {
  df <- data.frame(
    bin   = bin_names,
    idx   = 1:9,
    count = bin_values
  )
  
  # 用 bar 的位置 (idx) 和高度 (count) 擬合 kernel density
  x_rep   <- rep(df$idx, df$count)        # 每個 idx 重複 count 次
  dens    <- density(x_rep, from = 1, to = 9, adjust = 2.5)
  max_cnt <- max(df$count)
  
  # 把 density 線縮放到「柱子上方、紅字下方」的高度區間
  line_base   <- max_cnt * 1.05          # 線的基準高度（略高於最高的 bar）
  line_height <- max_cnt * 0.20          # 線的振幅
  dens_df <- data.frame(
    x = dens$x,
    y = dens$y * sum(df$count)   # ← 關鍵只有這一行
  )
  
  # 繪圖
  p <- ggplot(df, aes(x = idx, y = count)) +
    # 淡藍色柱子
    geom_col(fill = "#77AADD", width = 0.6) +
    # 上方 kernel density 線（在 bar 上方、紅字下方）
    geom_line(
      data = dens_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = 3,
      color = "#3E6F8E"
    ) +
    # 紅色數字標註
    geom_text(
      aes(label = count, y = count + 40),
      color = "#A63A2A",
      size = 35
    ) +
    # X 軸使用類別標籤
    scale_x_continuous(
      breaks = df$idx,
      labels = df$bin
    ) +
    # labs(x = NULL, y = "Count") +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x  = element_text(angle = 60, hjust = 1, color = "black", size =120),
      axis.text.y  = element_text(color = "black", size = 120),
      # axis.title.y = element_text(color = "black", size = 60),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      # 垂直方向的淡灰網格線
      panel.grid.major.y = element_line(color = "#DDDDDD"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank()
    )
  return(p)
}
draw_combined_pie_chart<-function(labels,A_values,B_values,A_subtitle,B_subtitle){
  n_colors <- length(labels)
  colors <- brewer.pal(n = min(n_colors, 8), name = "Set1")  # Set1 提供高辨識度顏色
  if (n_colors > 8) {
    colors <- colorRampPalette(colors)(n_colors)  # 若類別數超過 8，自動插值生成更多顏色
  }
  data_A <- data.frame(
    category = labels,
    value = A_values
  ) %>%
    mutate(percentage = value / sum(value) * 100,
           label = paste0(value, "\n", round(percentage, 1), "%"),
           pie = A_subtitle)
  
  data_B <- data.frame(
    category = labels,
    value = B_values
  ) %>%
    mutate(percentage = value / sum(value) * 100,
           label = paste0(value, "\n", round(percentage, 1), "%"),
           pie = B_subtitle)
  
  # 合併資料
  data_combined <- rbind(data_A, data_B) %>%
    mutate(pie = factor(pie, levels = c(A_subtitle, B_subtitle)))
  
  # 繪製並列餅狀圖
  ggplot(data_combined, aes(x = 1, y = value, fill = category)) +
    geom_bar(stat = "identity", position = "fill") +  # 使用 position = "fill" 確保填滿圓
    coord_polar("y") +
    geom_text(aes(label = label, x=1.34), 
              position = position_fill(vjust = 0.5),
              size = 7) +
    facet_wrap(~ pie, ncol = 2) +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 32, face = "bold", margin = margin(b = 15)),
      strip.background = element_rect(fill = "transparent", color = NA),  # 移除背景框
      legend.text = element_text(size = 32),
      legend.title = element_blank(),
      legend.key.size = unit(3, "lines"),  # 放大圖例圖標
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
}
draw_pie_chart <- function(subset_names, integer_values) {
  df <- data.frame(
    group = factor(subset_names, levels = subset_names, ordered = T),
    value = integer_values
  )
  df <- df %>%
    dplyr::mutate(
      pos = sum(value) - cumsum(value) + 0.5 * value,
      label_text = paste0(group, ": ", value)
    )
  p <- ggplot(df, aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    # geom_text_repel(aes(y = pos, label = group),#label_text),
                    ## angle = 45,       # 將文字旋轉 45 度
                    ## vjust = -0.5,
                    ## hjust = -0.1,
                    # size = 15,                 # 標籤字體大小
                    # nudge_x = 0.9,              # 將標籤向外推移
                    # segment.size = 0.5,       # 連接線段的粗細
                    # segment.color = NA,
                    # show.legend = F) +
    geom_text_repel(aes(y = pos, label = value),#label_text),
                    size = 40,                 # 標籤字體大小
                    nudge_x = 0.3,              # 將標籤向外推移
                    segment.size = 0.5,       # 連接線段的粗細
                    segment.color = NA,
                    show.legend = F) +
    scale_fill_brewer(palette = "Set2") +
    theme_void() +
    theme(# legend.position = "none",
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 100),
          legend.key.size = unit(5, "lines"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  return(p)
}
draw_pie_chart <- function(subset_names, integer_values) {
  df <- data.frame(
    group = factor(subset_names, levels = subset_names, ordered = TRUE),
    value = integer_values
  )
  df <- df %>%
    dplyr::mutate(
      pos = sum(value) - cumsum(value) + 0.5 * value,
      label_text = paste0(group, ": ", value)
    )
  
  ## 一組和前兩張圖統一的冷灰藍色調
  pie_cols <- c(
    "#99DDFF",  # 最淺藍灰，和橫條圖的 bar 顏色接近
    "#77AADD",
    "#44BB99",
    "#EE8866",
    "#FFAABB",
    "#DDDDDD"   # 最深一點，作為對比
  )
  
  p <- ggplot(df, aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text_repel(
      aes(y = pos, label = value),
      size = 40,
      nudge_x = 0.3,
      segment.size = 0.5,
      segment.color = NA,
      show.legend = FALSE
    ) +
    ## 這裡用手動配色替代原來的 Set2
    scale_fill_manual(values = pie_cols[seq_along(subset_names)]) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.title = element_blank(),
      legend.text  = element_text(size = 100),
      legend.key.size = unit(5, "lines"),
      plot.margin = margin(t = 0, r = 0, b = 10, l = 0)
    )
  
  return(p)
}
draw_percentage_bar_chart <- function(subset_names, percent, percent_lables) {
  df <- data.frame(
    names = factor(subset_names, levels = subset_names),
    values = percent,
    labels = percent_lables
  )
  ggplot(df, aes(x = names, y = values)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = labels), vjust = -0.5, color = "#A63A2A", size = 40) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +  # 基礎字體大小14pt，無襯線
    theme(
      axis.title.y = element_blank(),  # 隱藏y軸標題
      axis.text.y = element_blank(),   # 隱藏y軸刻度標籤
      axis.ticks.y = element_blank(),  # 隱藏y軸刻度線
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),  # 隱藏x軸標題
      axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),    # 移除網格線
      plot.title = element_blank(),     # 移除圖表標題
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}
draw_percentage_bar_chart <- function(subset_names, percent, percent_lables) {
  df <- data.frame(
    names  = factor(subset_names, levels = subset_names),
    values = percent,
    labels = percent_lables
  )
  
  ggplot(df, aes(x = names, y = values, fill = names)) +   # ★ 把 fill 映射到 names
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = labels), vjust = -0.5,
              color = "#A63A2A", size = 40) +
    
    ## ★ 用和餅圖相同的一組色號（按你的實際色號改這行）
    scale_fill_manual(values = c(
      BioCarta       = "#99DDFF",
      KEGG_MEDICUS   = "#77AADD",
      PID            = "#44BB99",
      Reactome       = "#EE8866",
      WikiPathways   = "#FFAABB",
      KEGG_LEGACY    = "#DDDDDD"
    )) +
    
    theme_classic() +
    # ylab("Percentage") +                                   # ★ 加 y 軸標題
    theme(
      #axis.title.y = element_text(size = 80, color = "black"),  # ★ 顯示 y 軸標題
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 120, color = "black"),  # ★ 顯示 y 軸刻度
      axis.ticks.y = element_line(),                            # ★ 顯示 y 軸刻度線
      axis.line.y  = element_line(),
      
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 120, color = "black",
                                  angle = 45, hjust = 1, vjust = 1),
      panel.grid   = element_blank(),
      plot.title   = element_blank(),
      plot.margin  = margin(t = 0, r = 0, b = 0, l = 0),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "#DDDDDD"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x * 100, "%"),
      expand = expansion(mult = c(0, 0.1))
    )
}
yes_or_no<-function(bool, a, b)
{
  if (bool) {
    return(a)
  } else {
    return(b)
  }
}
draw_three_barplots <- function(labels, vector_a, vector_b, vector_c, a_lab, b_lab, c_lab, need_legend, label_size){
  data <- data.frame(
    Label = rep(labels, 3),
    Value = c(vector_a, vector_b, vector_c),
    Group = rep(c(a_lab, b_lab, c_lab), each = length(labels))
  )
  ggplot(data, aes(x = factor(Label, levels = unique(labels)), y = Value,
                   fill = factor(Group, levels = c(a_lab, b_lab, c_lab)),
                   group = factor(Group, levels = c(a_lab, b_lab, c_lab)))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    geom_text(aes(label = as.integer(Value)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = label_size, color='red') +
    # scale_fill_manual(values = setNames(c("blue", "green", "red"), c(a_lab, b_lab, c_lab))) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    # labs(x = x_lab, y = y_lab, fill = "Embedding: ") +
    theme(legend.position = yes_or_no(need_legend,"top","none"),
          # axis.title = element_text(face = "plain", size = 100), # X 軸和 Y 軸標籤加粗
          axis.text.x = element_text(face='plain',angle = 45, hjust = 1, vjust = 1, size=100),
          legend.text = yes_or_no(need_legend,element_text(size = 100),element_blank()),
          legend.title = element_blank(),
          legend.key.size = yes_or_no(need_legend,unit(5, "lines"),element_blank()),  # 放大圖例圖標
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          axis.title.y = element_blank(),  # 隱藏y軸標題
          axis.text.y = element_blank(),   # 隱藏y軸刻度標籤
          axis.ticks.y = element_blank(),  # 隱藏y軸刻度線
          axis.line.y = element_blank(),
          axis.title.x = element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}
draw_three_barplots <- function(labels, vector_a, vector_b, vector_c,
                                a_lab, b_lab, c_lab,
                                need_legend, label_size){
  ## ---- 位置：加大每組之間的間距 ----
  x_pos <- (seq_along(labels) - 1) * 1.5   # 原本大約是 1，改成 1.5 拉開組距
  
  data <- data.frame(
    Label = rep(labels, 3),
    x     = rep(x_pos, 3),
    Value = c(vector_a, vector_b, vector_c),
    Group = rep(c(a_lab, b_lab, c_lab), each = length(labels))
  )
  
  ## ---- 顏色：與前面圖形和諧的三個色號 ----
  fill_cols  <- c(
    a_lab = "#A2C8C2",  # full-description：柔和青綠
    b_lab = "#F0C9A9",  # half-description：柔和杏色
    c_lab = "#B3C2E0"   # gene-symbol：柔和藍紫
  )
  line_cols <- c(
    a_lab = "#4C9085",  # density 線條用同色系稍深
    b_lab = "#C4763E",
    c_lab = "#4F6AA3"
  )
  names(fill_cols) <- c(a_lab, b_lab, c_lab)
  names(line_cols) <- c(a_lab, b_lab, c_lab)
  
  ## ---- 為每個 group 做 kernel density，並把 y 放大到和柱高同一尺度 ----
  max_y <- max(data$Value)
  
  make_density_df <- function(vals, x_vals, group_name){
    # 為了簡單且可控，複製 x，次數=該 bin 的值
    x_rep <- rep(x_vals, times = vals)
    # 如果某組全是 0，直接回空 data.frame
    if (length(x_rep) == 0) {
      return(data.frame())
    }
    dens <- density(x_rep,
                    bw  = 0.8,                    # 稍大一點的帶寬，讓曲線不要太彎
                    from = min(x_pos) - 0.5,
                    to   = max(x_pos) + 0.5,
                    adjust = 2.5)
    data.frame(
      x     = dens$x,
      y     = dens$y / max(dens$y) * max_y * 0.9,  # 縮放到接近最高柱
      Group = group_name
    )
  }
  
  dens_a <- make_density_df(vector_a, x_pos, a_lab)
  dens_b <- make_density_df(vector_b, x_pos, b_lab)
  dens_c <- make_density_df(vector_c, x_pos, c_lab)
  dens_df <- bind_rows(dens_a, dens_b, dens_c)
  
  ## ---- 繪圖 ----
  ggplot(data, aes(x = x, y = Value,
                   fill  = factor(Group, levels = c(a_lab, b_lab, c_lab)),
                   group = factor(Group, levels = c(a_lab, b_lab, c_lab)))) +
    # 柱子：加大 dodge 寬度，配合前面設定的 x_pos 拉開組距
    geom_col(position = position_dodge(width = 0.9), width = 0.9) +
    # 每個 group 的 density 曲線
    geom_line(data = dens_df,
              aes(x = x, y = y,
                  colour = factor(Group, levels = c(a_lab, b_lab, c_lab))),
              inherit.aes = FALSE,
              linewidth = 3) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = line_cols, guide = "none") +
    # 數字標籤
    geom_text(aes(label = as.integer(Value)),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = label_size, color = "#A63A2A") +
    # X 軸轉回原來的文字標籤
    scale_x_continuous(breaks = x_pos, labels = labels) +
    # Y 軸顯示 & 上方留一點空間
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_classic() +
    theme(
      legend.position = yes_or_no(need_legend, "top", "none"),
      # axis.title.y    = element_text(size = 80),
      axis.title.y    = element_blank(),
      axis.text.y     = element_text(size = 60, color = "black"),
      axis.ticks.y    = element_line(),
      axis.line.y     = element_line(),
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(face = "plain",
                                     angle = 45, hjust = 1, vjust = 1,
                                     size  = 120, color = "black"),
      legend.text     = yes_or_no(need_legend, element_text(size = 120),
                                  element_blank()),
      legend.title    = element_blank(),
      legend.key.size = yes_or_no(need_legend, unit(4, "lines"),
                                  unit(0, "lines")),
      plot.margin     = margin(t = 0, r = 0, b = 0, l = 0),
      # 加上柔和的水平格線
      panel.grid.major.y = element_line(color = "#DDDDDD", linewidth = 0.4),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank()
    )# +
    # labs(y = "Count", fill = "Embedding")
}
draw_three_barplots <- function(labels, vector_a, vector_b, vector_c,
                                a_lab, b_lab, c_lab,
                                need_legend, label_size){
  ## ---- 位置：加大每組之間的間距 ----
  x_pos <- (seq_along(labels) - 1) * 1.5   # 原本大約是 1，改成 1.5 拉開組距
  
  data <- data.frame(
    Label = rep(labels, 3),
    x     = rep(x_pos, 3),
    Value = c(vector_a, vector_b, vector_c),
    Group = rep(c(a_lab, b_lab, c_lab), each = length(labels))
  )
  
  ## ---- 顏色：與前面圖形和諧的三個色號 ----
  fill_cols  <- c(
    a_lab = "#EE8866",  # full-description：柔和青綠
    b_lab = "#44BB99",  # half-description：柔和杏色
    c_lab = "#AAAA00"   # gene-symbol：柔和藍紫
  )
  line_cols <- c(
    a_lab = "#E96235",  # density 線條用同色系稍深
    b_lab = "#36967B",
    c_lab = "#7A7A00"
  )
  dens_fill_cols <- c(
    a_lab = "#F5B8A3",  # full-description 的填充（例）
    b_lab = "#78CEB6",  # half-description 的填充（例）
    c_lab = "#F5F500"   # gene-symbol 的填充（例）
  )
  names(dens_fill_cols) <- c(a_lab, b_lab, c_lab)
  names(fill_cols) <- c(a_lab, b_lab, c_lab)
  names(line_cols) <- c(a_lab, b_lab, c_lab)
  fill_cols <- setNames(c("#EE8866", "#44BB99", "#AAAA00"), c(a_lab, b_lab, c_lab))
  dens_fill_cols <- setNames(c("#F5B8A3", "#36967B", "#7A7A00"), c(a_lab, b_lab, c_lab))
  line_cols <- setNames(c("#F5B8A3", "#78CEB6", "#F5F500"), c(a_lab, b_lab, c_lab))
  
  ## ---- 為每個 group 做 kernel density，並把 y 放大到和「自己那組柱高」同一尺度 ----
  make_density_df <- function(vals, x_vals, group_name){
    # 為了簡單且可控，複製 x，次數=該 bin 的值
    x_rep <- rep(x_vals, times = vals)
    # 如果某組全是 0，直接回空 data.frame
    if (length(x_rep) == 0 || all(vals == 0)) {
      return(data.frame())
    }
    dens <- density(x_rep,
                    bw  = 0.8,                    # 稍大一點的帶寬，讓曲線不要太彎
                    from = min(x_pos) - 0.5,
                    to   = max(x_pos) + 0.5,
                    adjust = 1)
    # group_max <- max(vals)  # 該 group 自己的最高柱子
    n_total <- sum(vals)                       # 该组总计数
    bin_w   <- median(diff(sort(unique(x_vals))))  # bin 宽度（你这里约等于 1.5）
    scale_k <- n_total * bin_w
    data.frame(
      x     = dens$x,
      # y     = dens$y / max(dens$y) * group_max * 1.05,  # 線高 ≈ 自己最高柱子的 1.05 倍
      y     = dens$y * scale_k,
      Group = group_name
    )
  }
  
  dens_a <- make_density_df(vector_a, x_pos, a_lab)
  dens_b <- make_density_df(vector_b, x_pos, b_lab)
  dens_c <- make_density_df(vector_c, x_pos, c_lab)
  dens_df <- bind_rows(dens_a, dens_b, dens_c)
  data$Group   <- factor(data$Group,   levels = c(a_lab, b_lab, c_lab))
  dens_df$Group <- factor(dens_df$Group, levels = c(a_lab, b_lab, c_lab))
  
  ## ---- 繪圖 ----
  ggplot(data, aes(x = x, y = Value,
                   # fill  = factor(Group, levels = c(a_lab, b_lab, c_lab)),
                   # fill = Group,
                   # group = factor(Group, levels = c(a_lab, b_lab, c_lab))
                   group = Group)) +
    # 柱子：加大 dodge 寬度，配合前面設定的 x_pos 拉開組距
    geom_area(
      data = dens_df,
      aes(x = x, y = y,
          # fill = factor(Group, levels = c(a_lab, b_lab, c_lab))),
          fill = Group),
      inherit.aes = FALSE,
      alpha = 0.25,
      position = "identity"
    ) +
    scale_fill_manual(values = dens_fill_cols) +
    ggnewscale::new_scale_fill() +
    geom_col(aes(fill = Group),
             position = position_dodge(width = 0.9), width = 0.9,
             color = "black", linewidth = 0.6) +
    scale_fill_manual(values = fill_cols) +
    # 每個 group 的 density 曲線
    geom_line(data = dens_df,
              aes(x = x, y = y,
                  colour = factor(Group, levels = c(a_lab, b_lab, c_lab))),
              inherit.aes = FALSE,
              linewidth = 4) +
    scale_colour_manual(values = line_cols, guide = "none") +
    # 數字標籤
    geom_text(aes(label = as.integer(Value)),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = label_size, color = "#A63A2A") +
    # X 軸轉回原來的文字標籤
    scale_x_continuous(breaks = x_pos, labels = labels) +
    # Y 軸顯示 & 上方留一點空間
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_classic() +
    theme(
      legend.position = yes_or_no(need_legend, "top", "none"),
      # axis.title.y    = element_text(size = 80),
      axis.title.y    = element_blank(),
      axis.text.y     = element_text(size = 120, color = "black"),
      axis.ticks.y    = element_line(),
      axis.line.y     = element_line(),
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(face = "plain",
                                     angle = 45, hjust = 1, vjust = 1,
                                     size  = 120, color = "black"),
      legend.text     = yes_or_no(need_legend, element_text(size = 120),
                                  element_blank()),
      legend.title    = element_blank(),
      legend.key.size = yes_or_no(need_legend, unit(4, "lines"),
                                  unit(0, "lines")),
      plot.margin     = margin(t = 0, r = 0, b = 0, l = 0),
      # 加上柔和的水平格線
      panel.grid.major.y = element_line(color = "#DDDDDD", linewidth = 0.4),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank()
    )# +
  # labs(y = "Count", fill = "Embedding")
}
p<-0.05
human_protein_genes<-read.csv("../embeddings/original_text-embedding-3-small_filtered.csv")$X
h_pathways<-readRDS(paste0("./data/pathways/h_test.rds"))
c2_pathways<-readRDS(paste0("./data/pathways/c2_test.rds"))

# Original Hallmark
original_h<-read.csv('original_h.csv')
original_h<-original_h[which(original_h$padj<p/1536),]
# fdr_pvals <- p.adjust(original_h$pval, method = "BH")
# original_h<-original_h[which(fdr_pvals<p),]
pathway_counts<-table(original_h$pathway)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(h_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bar_names <- sub("HALLMARK_", "", bar_names)
ploted<-plot_h_barplot(bar_names, bar_values)
ggsave("./final_plots/h_pathways_barplot.eps",plot=ploted,device='eps', width = 15, height = 30, units = "in")
dimension_counts<-c()
hp_counts<-c()
for (key in names(h_pathways)) {
  genes <- h_pathways[[key]]
  print(paste("Key:", key, "| human protein counts:", sum(genes %in% human_protein_genes)))
  dimension_counts<-c(dimension_counts,sum(original_h$pathway == key))
  hp_counts<-c(hp_counts,sum(genes %in% human_protein_genes))
}
df<-data.frame(
  pathway = names(h_pathways),
  dimension_counts = dimension_counts,
  hp_counts = hp_counts
)
write.csv(df, file = 'h_pathways_counts.csv', row.names = FALSE)

# Original C2
original_c2<-read.csv('original_c2.csv')
original_c2<-original_c2[which(original_c2$padj<p/1536),]
# fdr_pvals <- p.adjust(original_c2$pval, method = "BH")
# original_c2<-original_c2[which(fdr_pvals<p),]
pathway_counts<-table(original_c2$pathway)
sorted_indices <- order(pathway_counts)
sorted_pathway_counts <- pathway_counts[sorted_indices]
only_one_number<-sorted_pathway_counts[which(sorted_pathway_counts==1)]
length(only_one_number)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(c2_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bin_names<-c('0',
             '1',
             '2~5',
             '6~10',
             '11~20',
             '21~50',
             '51~100',
             '101~200',
             '>=201')
original_c2_bin_values<-c(sum(bar_values==0),
              sum(bar_values==1),
              sum(bar_values>=2 & bar_values<=5),
              sum(bar_values>=6 & bar_values<=10),
              sum(bar_values>=11 & bar_values<=20),
              sum(bar_values>=21 & bar_values<=50),
              sum(bar_values>=51 & bar_values<=100),
              sum(bar_values>=101 & bar_values<=200),
              sum(bar_values>=201))
ploted<-nature_barplot(bin_names, original_c2_bin_values)
ggsave("./final_plots/c2_pathways_barplot.eps",plot=ploted,device='eps', width = 35, height = 30, units = "in")
dimension_counts<-c()
hp_counts<-c()
for (key in names(c2_pathways)) {
  genes <- c2_pathways[[key]]
  print(paste("Key:", key, "| human protein counts:", sum(genes %in% human_protein_genes)))
  dimension_counts<-c(dimension_counts,sum(original_c2$pathway == key))
  hp_counts<-c(hp_counts,sum(genes %in% human_protein_genes))
}
useful_indices<- which(dimension_counts > 0)
df<-data.frame(
  pathway = names(c2_pathways)[useful_indices],
  dimension_counts = dimension_counts[useful_indices],
  hp_counts = hp_counts[useful_indices]
)
write.csv(df, file = 'c2_pathways_counts.csv', row.names = FALSE)

# C2 subsets
unique_detected_c2_pathways<-unique(original_c2$pathway)
subset_names <- c('BioCarta',
                  'KEGG_MEDICUS',
                  'PID',
                  'Reactome',
                  'WikiPathways',
                  'KEGG_LEGACY')
file_names <- c(
  'c2.cp.biocarta.v2024.1.Hs.symbols.gmt',
  'c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt',
  'c2.cp.pid.v2024.1.Hs.symbols.gmt',
  'c2.cp.reactome.v2024.1.Hs.symbols.gmt',
  'c2.cp.wikipathways.v2024.1.Hs.symbols.gmt',
  'c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt'
)
existing_values<-c()
missing_values<-c()
for (i in 1:length(subset_names))
{
  gene_set <- read.gmt(paste('./c2_subsets/', file_names[i], sep = ''))
  for(gene_set_name in names(gene_set)){
    if(sum(gene_set[[gene_set_name]] %in% human_protein_genes) < 10){
      gene_set<-gene_set[names(gene_set) != gene_set_name]
    }
  }
  existing_values<-c(existing_values, sum(names(gene_set) %in% unique_detected_c2_pathways))
  missing_values<-c(missing_values, sum(!(names(gene_set) %in% unique_detected_c2_pathways)))
}
# ploted<-draw_combined_pie_chart(subset_names,existing_values,missing_values,
#                         'Significant: 2744, 88.69%','Insignificant: 350, 11.31%')
# ggsave("./final_plots/c2_pathways_pie_chart.eps",plot=ploted,device='eps', width = 20, height = 10, units = "in")
ploted<-draw_pie_chart(subset_names,existing_values+missing_values)
ggsave("./final_plots/c2_pathways_pie_chart.eps",plot=ploted,device='eps', width = 35, height = 30, units = "in")
missing_percents<-missing_values/(existing_values+missing_values)
ploted<-draw_percentage_bar_chart(subset_names,missing_percents,
                                  paste0(round(missing_percents*100,2),"%"))
ggsave("./final_plots/c2_pathways_percentage_bar_chart.eps",plot=ploted,device='eps', width = 40, height = 30, units = "in")

# name-half C2
gene_names_c2<-read.csv('gene_names_c2.csv')
gene_names_c2<-gene_names_c2[which(gene_names_c2$padj<p/1536),]
# fdr_pvals <- p.adjust(gene_names_c2$pval, method = "BH")
# gene_names_c2<-gene_names_c2[which(fdr_pvals<p),]
pathway_counts<-table(gene_names_c2$pathway)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(c2_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bin_names<-c('0',
             '1',
             '2~5',
             '6~10',
             '11~20',
             '21~50',
             '51~100',
             '101~200',
             '>=201')
gene_names_c2_bin_values<-c(sum(bar_values==0),
              sum(bar_values==1),
              sum(bar_values>=2 & bar_values<=5),
              sum(bar_values>=6 & bar_values<=10),
              sum(bar_values>=11 & bar_values<=20),
              sum(bar_values>=21 & bar_values<=50),
              sum(bar_values>=51 & bar_values<=100),
              sum(bar_values>=101 & bar_values<=200),
              sum(bar_values>=201))

half_length_c2<-read.csv('half_length_c2.csv')
half_length_c2<-half_length_c2[which(half_length_c2$padj<p/1536),]
# fdr_pvals <- p.adjust(half_length_c2$pval, method = "BH")
# half_length_c2<-half_length_c2[which(fdr_pvals<p),]
pathway_counts<-table(half_length_c2$pathway)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(c2_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bin_names<-c('0',
             '1',
             '2~5',
             '6~10',
             '11~20',
             '21~50',
             '51~100',
             '101~200',
             '>=201')
half_length_c2_bin_values<-c(sum(bar_values==0),
                            sum(bar_values==1),
                            sum(bar_values>=2 & bar_values<=5),
                            sum(bar_values>=6 & bar_values<=10),
                            sum(bar_values>=11 & bar_values<=20),
                            sum(bar_values>=21 & bar_values<=50),
                            sum(bar_values>=51 & bar_values<=100),
                            sum(bar_values>=101 & bar_values<=200),
                            sum(bar_values>=201))

plotted<-draw_three_barplots(bin_names,
                             original_c2_bin_values,
                             half_length_c2_bin_values,
                             gene_names_c2_bin_values,
                             'full-description',
                             'half-description',
                             'gene-symbol',
                             F, 30)
ggsave("./final_plots/c2_pathways_three_barplots.pdf",plot=plotted,device=cairo_pdf, width = 45, height = 30, units = "in")

# name-half H
original_h<-read.csv('original_h.csv')
original_h<-original_h[which(original_h$padj<p/1536),]
# fdr_pvals <- p.adjust(original_h$pval, method = "BH")
# original_h<-original_h[which(fdr_pvals<p),]
pathway_counts<-table(original_h$pathway)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(h_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bin_names<-c('0',
             '1',
             '2~5',
             '6~10',
             '11~20',
             '21~50',
             '51~100',
             '101~200',
             '>=201')
original_h_bin_values<-c(sum(bar_values==0),
                           sum(bar_values==1),
                           sum(bar_values>=2 & bar_values<=5),
                           sum(bar_values>=6 & bar_values<=10),
                           sum(bar_values>=11 & bar_values<=20),
                           sum(bar_values>=21 & bar_values<=50),
                           sum(bar_values>=51 & bar_values<=100),
                           sum(bar_values>=101 & bar_values<=200),
                           sum(bar_values>=201))

gene_names_h<-read.csv('gene_names_h.csv')
gene_names_h<-gene_names_h[which(gene_names_h$padj<p/1536),]
# fdr_pvals <- p.adjust(gene_names_h$pval, method = "BH")
# gene_names_h<-gene_names_h[which(fdr_pvals<p),]
pathway_counts<-table(gene_names_h$pathway)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(h_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bin_names<-c('0',
             '1',
             '2~5',
             '6~10',
             '11~20',
             '21~50',
             '51~100',
             '101~200',
             '>=201')
gene_names_h_bin_values<-c(sum(bar_values==0),
                            sum(bar_values==1),
                            sum(bar_values>=2 & bar_values<=5),
                            sum(bar_values>=6 & bar_values<=10),
                            sum(bar_values>=11 & bar_values<=20),
                            sum(bar_values>=21 & bar_values<=50),
                            sum(bar_values>=51 & bar_values<=100),
                            sum(bar_values>=101 & bar_values<=200),
                            sum(bar_values>=201))

half_length_h<-read.csv('half_length_h.csv')
half_length_h<-half_length_h[which(half_length_h$padj<p/1536),]
# fdr_pvals <- p.adjust(half_length_h$pval, method = "BH")
# half_length_h<-half_length_h[which(fdr_pvals<p),]
pathway_counts<-table(half_length_h$pathway)
bar_names<- names(pathway_counts)
bar_values<-as.numeric(pathway_counts)
missing_names<-setdiff(names(h_pathways),bar_names)
bar_names<-c(bar_names, missing_names)
bar_values<-c(bar_values, rep(0, length(missing_names)))
bin_names<-c('0',
             '1',
             '2~5',
             '6~10',
             '11~20',
             '21~50',
             '51~100',
             '101~200',
             '>=201')
half_length_h_bin_values<-c(sum(bar_values==0),
                             sum(bar_values==1),
                             sum(bar_values>=2 & bar_values<=5),
                             sum(bar_values>=6 & bar_values<=10),
                             sum(bar_values>=11 & bar_values<=20),
                             sum(bar_values>=21 & bar_values<=50),
                             sum(bar_values>=51 & bar_values<=100),
                             sum(bar_values>=101 & bar_values<=200),
                             sum(bar_values>=201))

plotted<-draw_three_barplots(bin_names,
                             original_h_bin_values,
                             half_length_h_bin_values,
                             gene_names_h_bin_values,
                             'full-description',
                             'half-description',
                             'gene-symbol',
                             T,30)
ggsave("./final_plots/h_pathways_three_barplots.pdf",plot=plotted,device=cairo_pdf, width = 45, height = 30, units = "in")

draw_heatmap <- function(dat) {
  long_dat <- dat %>%
    pivot_longer(-model, names_to = "var", values_to = "score") %>%
    mutate(
      # 路徑集（Hallmark / C2）
      pathway = if_else(grepl("^H_", var),
                        "Hallmark pathways", "C2 pathways"),
      # Full / Half / Symbol
      type = case_when(
        grepl("Full",   var) ~ "Full",
        grepl("Half",   var) ~ "Half",
        grepl("Symbol", var) ~ "Symbol"
      ),
      # x 軸因子：把 pathway 與 type 組在一起方便排序
      x = factor(
        paste(pathway, type, sep = "_"),
        levels = c("Hallmark pathways_Full",
                   "Hallmark pathways_Half",
                   "Hallmark pathways_Symbol",
                   "C2 pathways_Full",
                   "C2 pathways_Half",
                   "C2 pathways_Symbol")
      ),
      # y 軸模型順序與圖中一致（第一列在最上面）
      model = factor(model, levels = rev(dat$model)),
      # 轉百分比
      pct = score * 100
    )
  
  ## 顏色設定 --------------------------------------------------------------
  # 上面小條：Full / Half / Symbol
  type_cols <- c(
    "Full"   = "#2FB47C",  # 深藍
    "Half"   = "#61CBF4",  # 淺藍
    "Symbol" = "#ffbf00"   # 黃
  )
  
  # 底下大條：Hallmark / C2 pathways
  pathway_cols <- c(
    "Hallmark pathways" = "#9b59b6",  # 紫
    "C2 pathways"       = "#1f78b4"   # 藍
  )
  
  # 熱圖主色階：藍(低) – 白(中) – 紅(高)
  heat_cols <- c("#313695", "#ffffff", "#a50026")
  
  ## 1. 上方 Full / Half / Symbol 色條 -----------------------------------
  type_bar_dat <- long_dat %>%
    distinct(x, type)
  
  top_bar <- ggplot(type_bar_dat, aes(x = x, y = 1, fill = type)) +
    geom_tile() +
    scale_x_discrete(
      labels = c("Full", "Half", "Symbol",
                 "Full", "Half", "Symbol")
    ) +
    scale_fill_manual(values = type_cols) +
    guides(fill = "none") +
    theme_void() +
    theme(
      plot.margin = margin(5.5, 5.5, 0, 5.5),
      axis.text.x = element_blank()
    )
  
  ## 2. 中間主熱圖 --------------------------------------------------------
  heatmap <- ggplot(long_dat, aes(x = x, y = model, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", pct)), size = 3) +
    scale_x_discrete(
      labels = c("Full", "Half", "Symbol",
                 "Full", "Half", "Symbol")
    ) +
    scale_fill_gradientn(
      colours = heat_cols,
      limits  = c(40, 100),
      name    = "Percentage"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 10),
      plot.margin  = margin(0, 5.5, 0, 5.5)
    )
  
  ## 3. 下方 Hallmark / C2 pathways 色條 -----------------------------------
  pathway_bar_dat <- long_dat %>%
    distinct(x, pathway)
  
  bottom_bar <- ggplot(pathway_bar_dat, aes(x = x, y = 1, fill = pathway)) +
    geom_tile() +
    scale_fill_manual(values = pathway_cols) +
    scale_x_discrete(labels = rep("", 6)) +
    guides(fill = "none") +
    theme_void() +
    theme(
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    ) +
    # 在底色條上加上文字標籤
    annotate("text", x = 2, y = 1, label = "Hallmark pathways",
             vjust = 1.8, size = 4, colour = "#4a1486") +
    annotate("text", x = 5, y = 1, label = "C2 pathways",
             vjust = 1.8, size = 4, colour = "#08306b")
  
  ## 4. 使用 patchwork 合併三張圖 ----------------------------------------
  final_plot <- top_bar / heatmap / bottom_bar +
    plot_layout(heights = c(0.15, 1, 0.22))
  
  final_plot
}
draw_heatmap <- function(dat) {
  # 保留原始順序，之後再反轉
  dat <- dat %>% mutate(model = factor(model, levels = model))
  
  long_dat <- dat %>%
    pivot_longer(-model, names_to = "var", values_to = "score") %>%
    mutate(
      pathway = if_else(grepl("^H_", var),
                        "Hallmark pathways", "C2 pathways"),
      type = case_when(
        grepl("Full",   var) ~ "Full",
        grepl("Half",   var) ~ "Half",
        grepl("Symbol", var) ~ "Symbol",
        TRUE ~ NA_character_
      ),
      col_id = as.numeric(factor(
        paste(pathway, type, sep = "_"),
        levels = c("Hallmark pathways_Full",
                   "Hallmark pathways_Half",
                   "Hallmark pathways_Symbol",
                   "C2 pathways_Full",
                   "C2 pathways_Half",
                   "C2 pathways_Symbol")
      )),
      model = factor(model, levels = rev(levels(model))),
      pct   = score * 100
    )
  
  type_cols <- c(
    "Full"   = "#EE8866",
    "Half"   = "#44BB99",
    "Symbol" = "#BBCC33"
  )
  pathway_cols <- c(
    "Hallmark pathways" = "#99DDFF",
    "C2 pathways"       = "#AAAA00"
  )
  heat_cols <- c("#77AADD", "#ffffff", "#EE8866")
  
  x_breaks  <- 1:6
  x_labels  <- c("Full", "Half", "Symbol", "Full", "Half", "Symbol")
  # ✅ 關鍵 1：整體 x 範圍從 0 開始，不再開到負數
  xlim_all  <- c(0, 6.5)
  
  ## 1. x 軸文字 --------------------------------------------------------
  type_bar_dat <- long_dat %>%
    distinct(col_id) %>%
    arrange(col_id)
  type_bar_dat$type_lab <- x_labels
  
  xlab_plot <- ggplot(type_bar_dat, aes(x = col_id, y = 0)) +
    geom_blank() +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels
    ) +
    scale_y_continuous(limits = c(-0.5, 0.5),
                       expand = expansion(mult = 0)) +
    coord_cartesian(xlim = xlim_all, clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      panel.background = element_blank(),
      axis.title   = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_text(size = 18, margin = margin(b = 3), colour = "black"),
      plot.margin  = margin(t = 1, r = 4, b = 0, l = 4)
    )
  
  ## 2. 上方色帶 --------------------------------------------------------
  type_bar2 <- long_dat %>%
    distinct(col_id, type)
  
  type_bar <- ggplot(type_bar2, aes(x = col_id, y = 0, fill = type)) +
    geom_tile(width = 1, height = 0.35) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = rep("", 6)
    ) +
    scale_y_continuous(limits = c(-0.5, 0.5),
                       expand = expansion(mult = 0)) +
    scale_fill_manual(values = type_cols) +
    coord_cartesian(xlim = xlim_all, clip = "off") +
    guides(fill = "none") +
    theme_void() +
    theme(
      plot.margin = margin(t = 0, r = 4, b = 2, l = 4)
    )
  
  ## 3. 主熱圖 + 行標籤 -------------------------------------------------
  label_dat <- tibble(
    model = levels(long_dat$model),
    # ✅ 關鍵 2：標籤 anchor 放在 0.3，介於 panel 左邊 (0) 和第一欄中心 (1) 之間
    x     = 0.4,
    col   = case_when(
      model == "biobert-base-cased-v1.1" ~ "#e41a1c",
      model == "MedEmbed-small-v0.1"     ~ "#1a9850",
      TRUE                               ~ "grey20"
    ),
    face  = ifelse(model %in% c("biobert-base-cased-v1.1",
                                "MedEmbed-small-v0.1"),
                   "bold", "plain")
  )
  
  heatmap <- ggplot(long_dat, aes(x = col_id, y = model, fill = pct)) +
    geom_tile(width = 1, height = 1, color = "white") +
    geom_text(aes(label = sprintf("%.2f", pct)), size = 5) +
    geom_text(
      data = label_dat,
      aes(x = x, y = model, label = model, colour = col, fontface = face),
      inherit.aes = FALSE,
      hjust = 1,
      size  = 5
    ) +
    scale_colour_identity(guide = "none") +
    scale_y_discrete(labels = NULL) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = rep("", 6)
    ) +
    scale_fill_gradientn(
      colours = heat_cols,
      limits  = c(25, 100),
      name    = "Percentage"
    ) +
    coord_cartesian(xlim = xlim_all, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 13),
      plot.margin  = margin(t = 0, r = 4, b = 0, l = 4)
    )
  
  ## 4. 下方 Hallmark / C2 色帶 ----------------------------------------
  pathway_bar_dat <- long_dat %>%
    distinct(col_id, pathway)
  
  bottom_bar <- ggplot(pathway_bar_dat, aes(x = col_id, y = 0, fill = pathway)) +
    geom_tile(width = 1, height = 0.25) +
    scale_fill_manual(values = pathway_cols) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = rep("", 6)
    ) +
    scale_y_continuous(limits = c(-0.9, 0.5),
                       expand = expansion(mult = 0)) +
    coord_cartesian(xlim = xlim_all, clip = "off") +
    guides(fill = "none") +
    theme_void() +
    theme(
      plot.margin = margin(t = 0, r = 4, b = 30, l = 150),
    ) +
    annotate("text", x = 2, y = -0.9, label = "Hallmark pathways",
             size = 6, colour = "#4a1486") +
    annotate("text", x = 5,   y = -0.9, label = "C2 pathways",
             size = 6, colour = "#08306b")
  
  ## 5. 垂直堆疊 ---------------------------------------------------------
  final_plot <- xlab_plot / type_bar / heatmap / bottom_bar +
    plot_layout(heights = c(0.06, 0.04, 1, 0.08))
  
  final_plot
}
summary <- read.csv('summary.csv', row.names = 1)
dataset_names <- c('biobert-base-cased-v1.1',
                   'original',
                   'stella-base-en-v2',
                   'gte-tiny',
                   'NoInstruct-small-Embedding-v0',
                   'e5-small',
                   'GIST-all-MiniLM-L6-v2',
                   'GIST-small-Embedding-v0',
                   'e5-small-v2',
                   'MedEmbed-small-v0.1',
                   'bge-small-en-v1.5',
                   'gte-small')
display_names <- c('biobert-base-cased-v1.1',
                   'OpenAI',
                   'stella-base-en-v2',
                   'gte-tiny',
                   'NoInstruct-small-Embedding-v0',
                   'e5-small',
                   'GIST-all-MiniLM-L6-v2',
                   'GIST-small-Embedding-v0',
                   'e5-small-v2',
                   'MedEmbed-small-v0.1',
                   'bge-small-en-v1.5',
                   'gte-small')
column_names <- c("H_Full", "H_Half", "H_Symbol", "C2_Full", "C2_Half", "C2_Symbol")
df <- data.frame(matrix(0, nrow = length(display_names), ncol = length(column_names),
                        dimnames = list(display_names, column_names)))
for (i in 1:length(dataset_names)) {
  df[display_names[i], 1] <- summary[dataset_names[i], 4]
  if ('original' == dataset_names[i])
  {
    df[display_names[i], 2] <- summary['half_length', 4]
    df[display_names[i], 3] <- summary['gene_names', 4]
  }
  else
  {
    df[display_names[i], 2] <- summary[paste0('half_length_', dataset_names[i]), 4]
    df[display_names[i], 3] <- summary[paste0('gene_names_', dataset_names[i]), 4]
  }
  df[display_names[i], 4] <- summary[dataset_names[i], 8]
  if ('original' == dataset_names[i])
  {
    df[display_names[i], 5] <- summary['half_length', 8]
    df[display_names[i], 6] <- summary['gene_names', 8]
  }
  else
  {
    df[display_names[i], 5] <- summary[paste0('half_length_', dataset_names[i]), 8]
    df[display_names[i], 6] <- summary[paste0('gene_names_', dataset_names[i]), 8]
  }
}
write.csv(df, file = 'final_heatmap_data.csv', row.names = TRUE)
# install.packages("markdown", lib = '~/myRlibs')
# install.packages("litedown", lib = '~/myRlibs')
library(markdown, lib = '~/myRlibs')
library(litedown, lib = '~/myRlibs')
# install.packages("readr", lib = '~/myRlibs')
library(readr, lib = '~/myRlibs')
# install.packages("forcats", lib = '~/myRlibs')
library(forcats, lib = '~/myRlibs')
# install.packages("patchwork", lib = '~/myRlibs')
library(patchwork, lib = '~/myRlibs')
library(grid, lib = '~/myRlibs')
library(cowplot, lib = '~/myRlibs')
dat <- read_csv("final_heatmap_data.csv", show_col_types = FALSE)
colnames(dat)[1] <- "model"
plotted <- draw_heatmap(dat)
ggsave("./final_plots/slm.eps",plot=plotted,device=cairo_ps, width = 10,
       height = 8, units = "in", dpi = 300)
