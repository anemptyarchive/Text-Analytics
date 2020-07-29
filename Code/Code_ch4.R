
# 第4章 法則と語句の重みおよび特徴語抽出 ----------------------------------------------------

# 形態素解析 -------------------------------------------------------------------

# 利用パッケージ
library(RMeCab)
library(dplyr)

# ファイルパスを指定
dir_path <- "text_data_cp932/kobushi"

# 抽出する品詞を指定
PoS <- c("名詞")

# 形態素解析
res_mecab <- docDF(dir_path, type = 1, pos = PoS)

# 単語文書行列を作成
term_doc_df <- res_mecab[, -(1:3)]
term_doc_mat <- res_mecab[, -(1:3)] %>% 
  as.matrix()

# 文書単語行列を作成
doc_term_mat <- res_mecab[, -(1:3)] %>% 
  t()


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






