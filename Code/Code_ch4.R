
# 第4章 法則と語句の重みおよび特徴語抽出 ----------------------------------------------------

# テキストの前処理 -------------------------------------------------------------------

# 利用パッケージ
library(RMeCab)
library(dplyr)
library(tidyr)
library(ggplot2)


# 抽出する品詞を指定
PoS <- c("名詞")

# 削除する語を指定
stop_words <- "[\\(\\)（）!?！？％,\\.…']"

# 1グループの各テキストを比較 -------------------------------------------------------------

# ファイルパスを指定
dir_path <- "text_data_cp932/kobushi"

# 形態素解析
res_mecab <- RMeCab::docDF(dir_path, type = 1, pos = PoS)

# 単語文書行列を作成
term_doc_mat <- res_mecab %>% 
  dplyr::filter(!grepl(stop_words, TERM)) %>% 
  .[, -(1:3)] %>% 
  as.matrix()
colnames(term_doc_mat) <- NULL

# 文書単語行列を作成
doc_term_mat <- res_mecab %>% 
  dplyr::filter(!grepl(stop_words, TERM)) %>% 
  .[, -(1:3)] %>% 
  t()
rownames(doc_term_mat) <- NULL


# 2グループを比較 ----------------------------------------------------------------

# ファイルパスを指定
dir_path_x <- "text_data_cp932/kobushi"
dir_path_y <- "text_data_cp932/tsubaki"

# 形態素解析
res_mecab_x <- RMeCab::docDF(dir_path_x, type = 1, pos = PoS) %>% 
  dplyr::filter(!grepl(stop_words, TERM))
res_mecab_y <- RMeCab::docDF(dir_path_y, type = 1, pos = PoS) %>% 
  dplyr::filter(!grepl(stop_words, TERM))

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


# 延べ語数
N_x <- sum(freq_df_x[["freq_x"]])
N_y <- sum(freq_df_y[["freq_y"]])

# 異なり語数
V_x <- nrow(freq_df_x)
V_y <- nrow(freq_df_y)


# 4.1 ジップの法則 --------------------------------------------------------------

# 頻度による順位付け
freq_rank_df_x <- freq_df_x %>% 
  dplyr::arrange(dplyr::desc(freq_x)) %>% 
  dplyr::mutate(rank = dplyr::min_rank(-freq_x))
freq_rank_df_y <- freq_df_y %>% 
  dplyr::arrange(dplyr::desc(freq_y)) %>% 
  dplyr::mutate(rank = dplyr::min_rank(-freq_y))


# ランク
r <- freq_rank_df_x[["rank"]]

# 相対頻度
f_r <- freq_rank_df_x[["freq_x"]] / N_x

# ジップの法則
c_r <- f_r * r
summary(c_r)

# 作図
Zipf_df <- tibble(
  rank = r, 
  relative_freq = f_r
)
ggplot(Zipf_df, aes(x = rank, y = relative_freq)) + 
  geom_point() + 
  labs(title = "Zipf's law", x = "rank", y = "relative freq")


# 作図
ggplot(freq_rank_df_y, aes(x = rank, y = freq_y / N_y)) + 
  geom_point() + 
  labs(title = "Zipf's law", x = "rank", y = "relative freq")


# ランクの対数
log_r <- log(r)

# 相対頻度の対数
log_f_r <- log(f_r)

# Zipf-Mandelbrot法則:線形回帰
res_lm <- lm(log_f_r ~ log_r)
summary(res_lm)

# 推定パラメータ
a <- res_lm$coefficients[["log_r"]]
log_c <- res_lm$coefficients[["(Intercept)"]]

# 作図
Zipf_df <- tibble(
  log_f_r = log_f_r, 
  log_r = log_r, 
  hat_log_f_r = log_c + a * log_r
)
ggplot(Zipf_df) + 
  geom_point(mapping = aes(x = log_r, y = log_f_r), position = "jitter") + # 散布図
  geom_line(mapping = aes(x = log_r, y = hat_log_f_r)) + # 回帰直線
  labs(title = "Zipf law", x = "log rank", y = "log relative freq")


# 4.2.1 延べ語数と異なり語数を用いた指標 --------------------------------------------------------------

# トークン比
TTR_x <- V_x / N_x
TTR_y <- V_y / N_y
TTR_x; TTR_y

# GuiraudのR
R_x <- V_x / sqrt(N_x)
R_y <- V_y / sqrt(N_y)
R_x; R_y

# HerdanのC
C_x <- log(V_x) / log(N_x)
C_y <- log(V_y) / log(N_y)
C_x; C_y

# Somersのs
s_x <- log(log(V_x)) / log(log(N_x))
s_y <- log(log(V_y)) / log(log(N_y))
s_x; s_y

# Maasのa^2
a2_x <- (log(N_x) - log(V_x)) / log(log(N_x))
a2_y <- (log(N_y) - log(V_y)) / log(log(N_y))
a2_x; a2_y

# TudavaのLN
LN_x <- (1 - V_x^2) / (V_x^2 * log(N_x))
LN_y <- (1 - V_y^2) / (V_y^2 * log(N_y))
LN_x; LN_y

# Duastのk
k_x <- log(V_x) / log(log(N_x))
k_y <- log(V_y) / log(log(N_y))
k_x; k_y

# DugastのU
U_x <- log(log(N_x)) / (log(N_x) - log(V_x))
U_y <- log(log(N_y)) / (log(N_y) - log(V_y))
U_x; U_y


# 4.2.2 頻度スペクトルを用いた指標 -----------------------------------------------------

# 出現頻度ごとの単語数(異なり語数)表
V_mN_df_x <- freq_df_x %>% 
  dplyr::count(freq_x)
V_mN_df_y <- freq_df_y %>% 
  dplyr::count(freq_y)

# 出現頻度ベクトル
m_x <- V_mN_df_x[["freq_x"]]
m_y <- V_mN_df_y[["freq_y"]]

# m回出現した単語数(異なり語数)ベクトル
V_mN_x <- V_mN_df_x[["n"]]
V_mN_y <- V_mN_df_y[["n"]]

# 処理の検証
sum(m_x * V_mN_x) == N_x; sum(m_y * V_mN_y) == N_y # 延べ語数
sum(V_mN_x) == V_x; sum(V_mN_y) == V_y # 異なり語数


# YuleのK
K_x <- 10^4 * (sum(m_x^2 * V_mN_x) - N_x) / N_x^2
K_y <- 10^4 * (sum(m_y^2 * V_mN_y) - N_y) / N_y^2
K_x; K_y

# SimpsonのD
D_x <- sum(V_mN_x * m_x / N_x * (m_x - 1) / (N_x - 1))
D_y <- sum(V_mN_y * m_y / N_y * (m_y - 1) / (N_y - 1))
K_x; K_y

# 出現頻度が2の単語数(異なり語数)
V_2N_x <- V_mN_df_x %>% 
  dplyr::filter(freq_x == 2) %>% 
  .[["n"]]
V_2N_y <- V_mN_df_y %>% 
  dplyr::filter(freq_y == 2) %>% 
  .[["n"]]

# SichelのS
S_x <- V_2N_x / V_x
S_y <- V_2N_y / V_y
S_x; S_y

# 出現頻度が1の単語数(異なり語数)
V_1N_x <- V_mN_df_x %>% 
  dplyr::filter(freq_x == 1) %>% 
  .[["n"]]
V_1N_y <- V_mN_df_y %>% 
  dplyr::filter(freq_y == 1) %>% 
  .[["n"]]

# HonoreのH
H_x <- 100 * log(N_x) / (1 - V_1N_x / V_x)
H_y <- 100 * log(N_y) / (1 - V_1N_y / V_y)
H_x; H_y


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
numer_ij <- t((t(tf_ij) * IDF_j))
TCF_w_ij <- numer_ij / sqrt(rowSums(numer_ij))

# ITC重み
numer_ij <- t(t(log(tf_ij + 1)) * IDF_j)
ITC_w_ij <- numer_ij / sqrt(rowSums(numer_ij)^2)


# 4.3.4 エントロピー重み付け --------------------------------------------------------

# エントロピー重み付け
tmp_ij <- t(t(tf_ij) / df_j)
Entropy_w_ij <- t(t(log(tf_ij + 1)) * (1 + colSums(tmp_ij * log(tmp_ij + 1e-7)) / log(N)))


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



# 4.4.1 カイ二乗統計量 -----------------------------------------------------------

# 単語文書行列
n_ij <- term_doc_mat

# 語句ごとの出現頻度
n_i <- rowSums(n_ij)

# 文書ごとの単語数
n_j <- colSums(n_ij)

# 総単語数
n <- sum(n_ij)

# 期待度数
e_ij <- rep(n_i, times = length(n_j)) * rep(n_j, each = length(n_i)) / n %>% 
  matrix(nrow = length(n_i), ncol = length(n_j))


# カイ二乗統計量
chi2 <- sum((n_ij - e_ij)^2 / e_ij)


# 尤度比統計量
G2 <- 2 * sum(n_ij * log(n_ij / e_ij), na.rm = TRUE)


# 語句を指定
i <- which.max(freq_df_xy[["freq_xy"]])

# 各文書の語句iの出現頻度
n_1j <- freq_df_xy[i, c("freq_x", "freq_y")] %>% 
  as.vector()

# 各文書の語句i以外の出現頻度
n_2j <- colSums(freq_df_xy[-i, c("freq_x", "freq_y")]) %>% 
  as.vector()

# 語句iの総出現頻度
n_1 <- sum(n_1j)

# 語句i以外の総出現頻度
n_2 <- sum(n_2j)

# 文書番号を2つ指定
d1 <- 1
d2 <- 2

# Fisherの正確確率
log_p_numer <- lgamma(n_1 + 1) + lgamma(n_2 + 1) + lgamma(n_j[d1] + 1) + lgamma(n_j[d2] + 1)
log_p_denom <- lgamma(n_1j[d1] + 1) + lgamma(n_1j[d2] + 1) + lgamma(n_2j[d1] + 1) + lgamma(n_2j[d2] + 1) + lgamma(n_1 + n_2 + 1)
p_n <- exp(log_p_numer - log_p_denom)
p_n


