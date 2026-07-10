library(splines)
library(parallel)
library(pracma)
library(testthat)

src_dir <- "D:/文档/Factor model/多研究函数型因子模型实践/Rcode/multiFSYNC/R"
for (f in list.files(src_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

test_dir("D:/文档/Factor model/多研究函数型因子模型实践/Rcode/multiFSYNC/tests/testthat")
