## compare partition probabilities

cp <- readRDS("CPP_partitionProbs.rds")
icrp <- t(readRDS("iCRP_partitionProbs.rds"))
lsp <- readRDS("LSP_partitionProbs.rds")

cp["11222", ]
plot(lsp["11222", ], icrp["11222", ])
plot(cp["11222", 1:11 ], icrp["11222", ])

cp<- cp[!(rownames(cp) == "11222"), 1:10]
lsp<- lsp[!(rownames(lsp) == "11222"), 1:10]
icrp<- icrp[!(rownames(icrp) == "11222"), 1:10]

values <- matrix(NA, nrow = 9, ncol = 10)

values[1, ] <- apply(icrp, 2, function(x) names(sort(x, TRUE)[1]))
values[2, ] <- apply(cp, 2, function(x) names(sort(x, TRUE)[1]))
values[3, ] <- apply(lsp, 2, function(x) names(sort(x, TRUE)[1]))

values[4, ] <- apply(icrp, 2, function(x) names(sort(x, TRUE)[2]))
values[5, ] <- apply(cp, 2, function(x) names(sort(x, TRUE)[2]))
values[6, ] <- apply(lsp, 2, function(x) names(sort(x, TRUE)[2]))

values[7, ] <- apply(icrp, 2, function(x) names(sort(x, TRUE)[3]))
values[8, ] <- apply(cp, 2, function(x) names(sort(x, TRUE)[3]))
values[9, ] <- apply(lsp, 2, function(x) names(sort(x, TRUE)[3]))

rownames(values) <- rep(c("iCRP", "CP", "LSP"), 3)
xtable::xtable(values)


# 
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrcccccccccc}
# \hline
# iCRP & $\alpha$ & 0 & 0.1 & 0.2 & 0.3 & 0.4 & 0.5 & 0.6 & 0.7 & 0.8 & 0.8 \\ 
# \hline
# &1st & 11111 & 11111 & 11111 & 11111 & 11111 & 11111 & 12222 & 12222 & 12222 & 12222 \\ 
# &2nd & 11112 & 12222 & 12222 & 12222 & 12222 & 12222 & 12111 & 12111 & 12111 & 12111 \\ 
# &3 & 11121 & 12111 & 12111 & 12111 & 12111 & 12111 & 11111 & 12333 & 12333 & 12333 \\ 
# \hline
# CP & $\psi$ & 0 &  1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9  \\ 
# \hline
# &1st & 11111 & 11111 & 11111 & 11111 & 11111 & 12333 & 12333 & 12333 & 12333 & 12333 \\ 
# &2nd & 11112 & 12222 & 12333 & 12333 & 12333 & 11111 & 11111 & 11111 & 11232 & 11232 \\ 
# &3rd & 11121 & 12111 & 12222 & 12222 & 11232 & 11232 & 11232 & 11232 & 11233 & 11233 \\ 
# \hline
# LSP & $\tau$ & 3 & 2.7 & 2.4 & 2.1 & 1.8 & 1.5 & 1.2 & 0.9 & 0.6 & 0.3  \\ 
# \hline
# &1st  & 11111 & 11111 & 11111 & 11111 & 11111 & 11111 & 11111 & 11111 & 11223 & 11223 \\ 
# &2nd & 11112 & 11112 & 11112 & 11112 & 11223 & 11223 & 11223 & 11223 & 11111 & 11111 \\ 
# &3rd & 11122 & 11223 & 11223 & 11223 & 11112 & 11221 & 11221 & 11212 & 11221 & 11232 \\ 
# \hline
# \end{tabular}
# \end{table}

