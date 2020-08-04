
# 第4章 法則と語句の重みおよび特徴語抽出 ----------------------------------------------------

# テキストの前処理 -------------------------------------------------------------------

# 利用パッケージ
library(RMeCab)
library(dplyr)
library(tidyr)
library(ggplot2)


# 抽出する品詞を指定
PoS <- c("名詞")


# ファイルパスを指定
dir_path <- "text_data_cp932/kobushi"

# 形態素解析
res_mecab <- docDF(dir_path, type = 1, pos = PoS)

# 単語文書行列を作成
term_doc_mat <- res_mecab[, -(1:3)] %>% 
  as.matrix()

# 文書単語行列を作成
doc_term_mat <- res_mecab[, -(1:3)] %>% 
  t()


# ファイルパスを指定
dir_path_x <- "text_data_cp932/kobushi"
dir_path_y <- "text_data_cp932/tsubaki"

# 形態素解析
res_mecab_x <- docDF(dir_path_x, type = 1, pos = PoS)
res_mecab_y <- docDF(dir_path_y, type = 1, pos = PoS)

# 頻度を集計
freq_df_x <- tidyr::tibble(
  term = res_mecab_x[["TERM"]], 
  freq_x = rowSums(res_mecab_x[, -c(1:3)])
)
freq_df_y <- tidyr::tibble(
  term = res_mecab_y[["TERM"]], 
  freq_y = rowSums(res_mecab_y[, -c(1:3)])
)

# 単語文書行列を作成
freq_df_xy <- dplyr::full_join(
  freq_df_x, freq_df_y, by = "term"
) %>% 
  dplyr::mutate(
    freq_x = replace_na(freq_x, 0), 
    freq_y = replace_na(freq_y, 0)
  ) %>% 
  dplyr::mutate(freq_xy = freq_x + freq_y)

# 頻度による順位付け
freq_rank_df_x <- freq_df_x %>% 
  dplyr::arrange(dplyr::desc(freq_x)) %>% 
  dplyr::mutate(rank = dplyr::min_rank(-freq_x))

# 総単語数
N_x <- sum(freq_df_x[["freq_x"]])


# 4.1 ジップの法則 --------------------------------------------------------------

# 相対頻度
f_j <- freq_rank_df_x[["freq_x"]] / N_x

# ランク
r_j <- freq_rank_df_x[["rank"]]

# ジップの法則
c_j <- f_j * r_j
summary(c_j)

# 作図
Zipf_df <- tibble(
  rank = r_j, 
  relative_freq = f_j
)
ggplot(Zipf_df, aes(x = rank, y = relative_freq)) + 
  geom_point() + 
  labs(title = "Zipf's law", x = "rank", y = "relative freq")
ggplot(freq_rank_df_x, aes(x = rank, y = freq_x / sum(freq_x))) + 
  geom_point() + 
  labs(title = "Zipf's law", x = "rank", y = "relative freq")


# 拡張版ジップの法則？


# 相対頻度の対数
log_f_j <- log(f_j)

# ランクの対数
log_r_j <- log(r_j)

# Zipf-Mandelbrot法則:線形回帰
res_lm <- lm(log_f_j ~ log_r_j)
summary(res_lm)

# 推定パラメータ
log_c <- res_lm$coefficients[["(Intercept)"]]
a <- res_lm$coefficients[["log_r_j"]]

# 作図
ZM_df <- tibble(
  log_f = log_f_j, 
  log_r = log_r_j, 
  hat_log_f = log_c + a * log_r
)
ggplot(ZM_df) + 
  geom_point(mapping = aes(x = log_r, y = log_f), position = "jitter") + 
  geom_line(mapping = aes(x = log_r, y = hat_log_f)) + 
  labs(title = "Zipf-Mandelbrot law", x = "log rank", y = "log relative freq")


# 4.3.3 TF-IDF重み付け --------------------------------------------------------

# テキストにおける語句tの頻度
tf_ij <- doc_term_mat

# テキストの総数
N <- nrow(doc_term_mat)

# 語句tを含むテキストの数
df_j <- colSums(doc_term_mat > 0)

# IDF
IDF_j <- log(N / df_j)

# TF-IDF
TF_IDF <- t(t(tf_ij) * IDF_j)

# tf_ij-IDF重み付け
TFIDF_w_ij <- t(log(t(tf_ij) + 1) * IDF_j)

# TCF重み
TCF_w_ij <- t((t(tf_ij) * IDF_j) / sqrt(colSums(t(tf_ij) * IDF_j)))

# ITC重み
ITC_w_ij <- t((t(log(tf_ij + 1)) * IDF_j) / sqrt(colSums((t(log(tf_ij + 1) * IDF_j))^2)))


# 4.3.4 エントロピー重み付け --------------------------------------------------------

# エントロピー重み付け
Entropy_w_ij <- log(tf_ij + 1) * (1 + colSums(t(tf_ij) / df_j * log(t(tf_ij) / df_j + 1e-7)) / log(N))


# 4.3.5 相互情報量による共起頻度の重み付け -------------------------------------------------

# 語句iの頻度:N_x_i > 0
N_x_i <- freq_df_x[["freq_x"]]
N_y_k <- freq_df_y[["freq_y"]]

# 語句iの出現確率:最尤法
p_x_i <- N_x_i / sum(N_x_i)
p_y_k <- N_y_k / sum(N_y_k)

# シャノンのエントロピー
H_x <- - sum(p_x_i * log(p_x_i), na.rm = TRUE)
H_y <- - sum(p_y_k * log(p_y_k), na.rm = TRUE)

# 条件チェック
0 <= H_x | H_x <= log(sum(N_x_i))
0 <= H_y | H_y <= log(sum(N_y_k))

# 頻度表:
freq_df_xy2 <- tidyr::tibble(
  term_x = rep(freq_df_xy[["term"]], times = nrow(freq_df_xy)), 
  freq_x = rep(freq_df_xy[["freq_x"]], times = nrow(freq_df_xy)), 
  term_y = rep(freq_df_xy[["term"]], each = nrow(freq_df_xy)), 
  freq_y = rep(freq_df_xy[["freq_y"]], each = nrow(freq_df_xy))
) %>% 
  dplyr::mutate(freq_xy = freq_x + freq_y)

# 語句jの頻度:N_x_j >= 0
N_x_j <- freq_df_xy2[["freq_x"]]
N_y_j <- freq_df_xy2[["freq_y"]]

# 語句jの出現確率
p_x_j <- N_x_j / sum(N_x_j)
p_y_j <- N_y_j / sum(N_y_j)

# 結合エントロピー？
H_xy = - sum(p_x_j * p_y_j * log(p_x_j * p_y_j), na.rm = TRUE)

# 条件チェック
H_xy <= H_x + H_y


# 条件付きエントロピー


# 語句jの頻度:N_x_j >= 0
N_x_j <- freq_df_xy[["freq_x"]]
N_y_j <- freq_df_xy[["freq_y"]]

# 語句jの出現確率
p_x_j <- N_x_j / sum(N_x_j)
p_y_j <- N_y_j / sum(N_y_j)

# カルバックライブラーダイバージェンス
KLD_x.y <- sum(p_x_j * log(p_x_j / p_y_j) * (p_y_j / p_y_j), na.rm = TRUE)
KLD_y.x <- sum(p_y_j * log(p_y_j / p_x_j) * (p_x_j / p_x_j), na.rm = TRUE)

# 条件チェック
KLD_x.y >= 0
KLD_y.x >= 0
KLD_x.y == KLD_y.x # 通常はFALSE

# ユークリッド距離
ED_xy <- sqrt(sum((p_x_j - p_y_j)^2))
ED_yx <- sqrt(sum((p_y_j - p_x_j)^2))

# 条件チェック
ED_xy == ED_yx


## 相互情報量



