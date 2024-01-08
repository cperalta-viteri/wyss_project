







V <- torch_rand(5,6,7)


O <- torch_zeros(2L,6L,7L)
O$index_reduce(1L,c(1L,1L,1L,2L,2L),V,reduce = "mean")


O$scatter_reduce(1L,array(1L,c(5,6,7)),V,reduce="sum")







