
# 第5章 テキストの特徴分析 -----------------------------------------------------------


# 5.1 特徴分析のデータの形式 ---------------------------------------------------------

# 利用パッケージ
library(RMeCab)
library(dplyr)


# 抽出する品詞を指定
PoS <- c("名詞")

# ファイルパスを指定
dir_path <- "text_data_cp932/kobushi"

# 形態素解析
res_mecab <- docDF(dir_path, type = 1, pos = PoS)


# 文書単語行列を作成
doc_term_mat <- res_mecab %>% 
  dplyr::mutate(n = rowSums(.[, -(1:3)])) %>% # 全文書での頻度を集計
  dplyr::mutate(
    TERM = dplyr::case_when(
      n <= 1 ~ "(OTHERS)", 
      TRUE ~ TERM
    )
  ) %>% # 指定した頻度以下の単語を"(OTHERS)"に統一
  dplyr::select(-c("POS1", "POS2")) %>% # 不要な列を除く
  dplyr::group_by(TERM) %>% # グループ化
  dplyr::summarise_all(sum) %>% # グループごとに頻度を合計
  #dplyr::filter(TERM != "(OTHERS)") %>% # "(OTHERS)"を除く
  dplyr::select(-c("TERM", n)) %>% # 不要な列を除く
  t() # 転置


# 5.2 特異値分解 ---------------------------------------------------------------

# 特異値分解
res_svd <- svd(doc_term_mat)

# 特異値を取得
d <- res_svd[["d"]]

# Uを取得
U <- res_svd[["u"]]

# Vを取得
V <- res_svd[["v"]]

# 検証
dplyr::near(doc_term_mat, U %*% diag(d) %*% t(V))


# 行列を作成
mat <- matrix(1:15, nrow = 3, ncol = 5)
mat

# 特異値分解
res <- svd(mat)
d <- res[["d"]]
U <- res[["u"]]
V <- res[["v"]]
U; diag(d); V

# 検証
U %*% diag(d) %*% t(V)


# 5.3 主成分分析 ---------------------------------------------------------------



# 5.4 対応分析 ----------------------------------------------------------------



