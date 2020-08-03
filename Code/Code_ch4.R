
# 第4章 法則と語句の重みおよび特徴語抽出 ----------------------------------------------------

# テキストの前処理 -------------------------------------------------------------------

# 利用パッケージ
library(RMeCab)
library(dplyr)
library(tidyr)

# ファイルパスを指定
dir_path <- "text_data_cp932/kobushi"

# 抽出する品詞を指定
PoS <- c("名詞")

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
  TERM = res_mecab_x[["TERM"]], 
  x = rowSums(res_mecab_x[, -c(1:3)])
)
freq_df_y <- tidyr::tibble(
  TERM = res_mecab_y[["TERM"]], 
  y = rowSums(res_mecab_y[, -c(1:3)])
)

# 単語文書行列を作成
freq_df_xy <- dplyr::full_join(
  freq_df_x, freq_df_y, by = "TERM"
) %>% 
  dplyr::mutate(
    x = replace_na(x, 0), 
    y = replace_na(y, 0)
  ) %>% 
  dplyr::mutate(xy = x + y)


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
N_x_i <- freq_df_x[["x"]]
N_y_k <- freq_df_y[["y"]]

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
  TERM_x = rep(freq_df_xy[["TERM"]], times = nrow(freq_df_xy)), 
  x = rep(freq_df_xy[["x"]], times = nrow(freq_df_xy)), 
  TERM_y = rep(freq_df_xy[["TERM"]], each = nrow(freq_df_xy)), 
  y = rep(freq_df_xy[["y"]], each = nrow(freq_df_xy))
) %>% 
  dplyr::mutate(xy = x + y)

# 語句jの頻度:N_x_j >= 0
N_x_j <- freq_df_xy2[["x"]]
N_y_j <- freq_df_xy2[["y"]]

# 語句jの出現確率
p_x_j <- N_x_j / sum(N_x_j)
p_y_j <- N_y_j / sum(N_y_j)

# 結合エントロピー？
H_xy = - sum(p_x_j * p_y_j * log(p_x_j * p_y_j), na.rm = TRUE)

# 条件チェック
H_xy <= H_x + H_y


# 条件付きエントロピー


# 語句jの頻度:N_x_j >= 0
N_x_j <- freq_df_xy[["x"]]
N_y_j <- freq_df_xy[["y"]]

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



