context("convertToNCBIGeneID")

#
# Case I : Not unique rowID
#
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "A", "B", "C")
LefttoRight <- rbind(
  c("A", "1"),
  c("B", "2"),
  c("C", "3")
)
input <- convertToNCBIGeneID(input, rowID, LefttoRight)
expect_true(identical(dim(input), c(3L, 5L)))

#
# Case II : Not unique LefttoRight
#
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "B", "C", "D")
LefttoRight <- rbind(
  c("A", "1"),
  c("B", "2"),
  c("B", "4"),
  c("C", "3"),
  c("D", "5")
)
input <- convertToNCBIGeneID(input, rowID, LefttoRight)
expect_true(identical(dim(input), c(4L, 5L)))

#
# Case III : Not unique rowID and not unique LefttoRight
#
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "B", "B", "D")
LefttoRight <- rbind(
  c("A", "1"),
  c("B", "2"),
  c("B", "4"),
  c("D", "5")
)
input <- convertToNCBIGeneID(input, rowID, LefttoRight)
expect_true(identical(dim(input), c(3L, 5L)))


#
# Case IV : Missing/not unique rowID and missing/not unique LefttoRight
#
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", NA, "B", "B")
LefttoRight <- rbind(
  c("A", "1"),
  c("B", "2"),
  c("B", "4"),
  c("D", NA)
)
input <- convertToNCBIGeneID(input, rowID, LefttoRight)
expect_true(identical(dim(input), c(2L, 5L)))
