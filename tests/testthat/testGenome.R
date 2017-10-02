library(testthat)
library(AnaCoDa)

context("Genome")

test_that("general genome functions", {
  expect_equal(testGenome("UnitTestingData"), 0)
})

g <- new(Genome)

test_that("check Index", {
  expect_equal(g$checkIndex(2, 1, 10), TRUE)
  expect_equal(g$checkIndex(5, 1, 10), TRUE)
  expect_equal(g$checkIndex(20, 1, 30), TRUE)
  expect_equal(g$checkIndex(5, 4, 6), TRUE)
  
  #Checking invalid cases
  expect_equal(g$checkIndex(11, 3, 10), FALSE)
  expect_equal(g$checkIndex(5, 6, 11), FALSE)
})

g1 <- new(Gene, "ATGGCCACTATTGGGTCTTAG", "TEST001", "TEST001 Test Gene")
g2 <- new(Gene, "TGGGATTACCAA", "TEST002", "TEST002 Test Gene")
g3 <- new(Gene, "TTGGAAACCACA", "TEST003", "TEST003 Test Gene");
g4 <- new(Gene, "TGGGATTACCCC", "TEST004", "TEST004 Test Gene");
s1 <- new(Gene, "TGGGATTACCAA", "TEST011", "TEST011 Test Gene");
#simulated gene

g$addGene(g1, FALSE)
g$addGene(g2, FALSE)
g$addGene(g3, FALSE)
g$addGene(g4, FALSE)
g$addGene(s1, TRUE)

test_that("get Gene By Index", {
  expect_equal(g$getGeneByIndex(1, FALSE), g1)
  expect_equal(g$getGeneByIndex(2, FALSE), g2)
  expect_equal(g$getGeneByIndex(3, FALSE), g3)
  expect_equal(g$getGeneByIndex(4, FALSE), g4)
  expect_equal(g$getGeneByIndex(1, TRUE), s1)
  
  #Checking invalid cases: invalid indices, therefore returning first gene
  #First check is on simulated genes
  expect_equal(g$getGeneByIndex(2, TRUE), s1)
  expect_equal(g$getGeneByIndex(6, FALSE), g1)
})

test_that("get Gene By Id", {
  expect_equal(g$getGeneById("TEST001", FALSE), g1)
  expect_equal(g$getGeneById("TEST002", FALSE), g2)
  expect_equal(g$getGeneById("TEST003", FALSE), g3)
  expect_equal(g$getGeneById("TEST004", FALSE), g4)
  expect_equal(g$getGeneById("TEST011", TRUE), s1)
})
  
# Create two slices of genomes:
# t1, containing indices from 0 to 3
# t2, containing indices from 1 to 2 with one simulated gene
t1 <- new(Genome)
t2 <- new(Genome)

t1$addGene(g1, FALSE)
t1$addGene(g2, FALSE)
t1$addGene(g3, FALSE)

t2$addGene(g2, FALSE)
t2$addGene(g3, FALSE)
t2$addGene(s1, TRUE)

test_that("get Genome For Gene Indices R", {
  expect_equal(g$getGenomeForGeneIndices(c(0,1,2), FALSE), t1)
  expect_equal(g$getGenomeForGeneIndices(c(1,2), FALSE), t2)
  expect_equal(g$getGenomeForGeneIndices(c(0), TRUE), t2)
})
