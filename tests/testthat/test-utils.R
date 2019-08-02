context("Utils - bpgmm toolbox")

test_that("Vector is correctly transformed to adjacency matrix", {
  ZOneDim1 <- c(1,2,3,1,2,2)
  ZOneDim2 <- c(1,1,3,1,2,2)

  Zmat1 <- matrix(c(1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0),3,6)
  Zmat2 <- matrix(c(1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0),3,6)

  expect_equal(getZmat(ZOneDim1,3,6),Zmat1)
  expect_equal(getZmat(ZOneDim2,3,6),Zmat2)

  expect_equal(get_Z_mat(ZOneDim1,3,6),Zmat1)
  expect_equal(get_Z_mat(ZOneDim2,3,6),Zmat2)
})


test_that("Log scale ratio calculation is done correctly", {
  expect_equal(calculateRatio(log(2),log(3)),2/3)
  expect_equal(calculateRatio(log(3),log(2)),3/2)
  expect_equal(calculateRatio(log(14),log(13)),14/13)
  expect_equal(calculateRatio(log(32),log(2)),32/2)
  expect_equal(calculateRatio(log(23),log(2)),23/2)
})


test_that("Convert list of string to vector of string correctly", {
  stringList = list()
  stringList[[1]] = c("a", "bb")
  stringList[[2]] = c("aww", "bqqb","sdf")
  expect_equal(listToStrVec(stringList) ,c("(abb)", "(awwbqqbsdf)"))
})
