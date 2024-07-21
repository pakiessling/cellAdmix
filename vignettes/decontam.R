knn_adjacency_matrix = function(df, k, weighted=FALSE) {
    nns = df[,c("x", "y", "z")] %>% FNN::get.knn(k=k)
    xn = nns$nn.index

    i = rep(1:nrow(xn), each=ncol(xn))
    j = t(xn) %>% as.vector()
    x = if(weighted) {t(nns$nn.dist) %>% as.vector()} else rep(1, length(i))
    adj_mat = sparseMatrix(i=i, j=j, x=x)
    return(adj_mat)
}

#' This is a much faster implementation than `knn_count_mat`
knn_count_matrix = function(df, k, all.genes=NULL) {
    adj_mat = knn_adjacency_matrix(df, k)
    j = if (is.null(all.genes)) factor(df$gene) else factor(df$gene, levels=all.genes)
    ncm = adj_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
    colnames(ncm) = levels(j)
    rownames(ncm) = rownames(df)

    return(ncm)
}

#' Get KNN counts for each molecule
#' @param df Data frame with molecule coordinates
#' @param k Neighborhood size
knn_count_mat = function(df_cell, k, include_i = FALSE) {

    k = as.integer(k)

    xn = df_cell %>% select(x, y, z) %>% FNN::get.knn(k = k) %>% .$nn.index

    gene_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$gene, .$id)}
    id_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$mol_id, .$id)}

    nmat = apply(xn, 2, function(i){gene_dict[i]})

    if (include_i) {
        nmat = cbind(df_cell$gene, nmat)
    }

    counts = nmat %>%
        reshape2::melt() %>%
        select(id = Var1, gene_j = value) %>%
        mutate(
          gene_i = gene_dict[id],
          id = id_dict[id]
        )

    counts = counts %>%
        {data.table::as.data.table(.)} %>%
        count(id, gene_i, gene_j, .drop = T) %>%
        {data.table::dcast(., id + gene_i ~ gene_j, value.var = 'n', fill = 0)}

    return(counts)

}

#' Get KNN counts for each molecule
#' @param df Data frame with molecule coordinates
#' @param h Neighborhood size
#' @param ncores Number of cores to use
#' @return Data frame with KNN counts for each molecule
get_knn_counts = function(df, h, include_i = FALSE, ncores = 1) {

    h = as.integer(h)

    cells = df %>% count(cell) %>% filter(n > h) %>% pull(cell)

    N = df %>%
        filter(cell %in% cells) %>%
        split(.$cell) %>%
        mclapply(
            mc.cores = ncores,
            function(df_cell) {

                cell = unique(df_cell$cell)

                counts = knn_count_mat(df_cell, include_i = include_i, k = h) %>%
                    mutate(cell = cell, .before = 1)

                return(counts)

        }) %>%
        bind_rows()

    N[is.na(N)] = 0

    N = as.data.frame(N)

    return(N)
}

##### Factorization #####

run_vrnmf = function(X, n.comp, ncores=1) {
    cvr = suppressMessages({
        vrnmf::vol_preprocess(X) %>% vrnmf::volnmf_main(
            n.comp=k, wvol = 2e-2, n.iter = 3e+3, vol.iter = 1e+1, c.iter = 1e+1,
            extrapolate=TRUE, accelerate=TRUE, verbose=FALSE
        )
    })

    cvr$D = vrnmf::infer_intensities(cvr$C, t(X), n.cores=ncores) %>%
        t() %>% set_rownames(rownames(X))

    return(cvr)
}

run_nmf_raw = function(X, k, method = 'weighted', nrun = 1, ncores = 1, ct_dict = NULL) {

    if (method == 'weighted') {
        if (!is.null(ct_dict)) {
            # perform per-celltype weighting
            ct_dict = ct_dict[names(ct_dict) %in% rownames(X)]
            ct_dict = sort(ct_dict)

            Z = lapply(
                unique(ct_dict),
                function(ct) {
                    X_ct = X[names(ct_dict[ct_dict == ct]),]
                    Z_ct = matrix(rep(1/colSums(X_ct), nrow(X_ct)), nrow = nrow(X_ct), byrow = TRUE)
                    Z_ct[is.infinite(Z_ct)] = 1
                    rownames(Z_ct) = rownames(X_ct)
                    colnames(Z_ct) = colnames(X_ct)
                    return(Z_ct)
                }
            ) %>%
            Reduce('rbind', .)

            Z = Z[names(ct_dict),]
            X = X[names(ct_dict),]
        } else {
            Z = matrix(rep(1/colSums(X), nrow(X)), nrow = nrow(X), byrow = TRUE)
        }

        res = nmf(X, rank = k, method = 'ls-nmf', weight = Z, nrun = nrun, .opt=paste0('vp', ncores), seed=0)
    } else {
        res = nmf(X, rank = k, nrun = nrun, .opt=paste0('vp', ncores), seed=0)
    }

    return(res)
}

loadings_to_long_df = function(loadings) {
    loadings = t(loadings) %>%
        as.data.frame() %>%
        setNames(1:ncol(.)) %>%
        tibble::rownames_to_column('gene') %>%
        reshape2::melt(id.var = 'gene', variable.name = 'factor', value.name = 'loading') %>%
        group_by(factor) %>%
        mutate(frac = loading/sum(loading)) %>%
        ungroup()

    return(loadings)
}

#' @param ct_dict Cell type dictionary for each molecule id
run_nmf = function(N, k, method = 'weighted', ct_dict = NULL, nrun = 30, ncores = NULL) {

    if (is.null(ncores)) {
        ncores = nrun
    }

    if (!is.null(ct_dict)) {
        if (!all(N$id %in% names(ct_dict))) {
            stop('ct_dict does not contain all ids in N!')
        }
        ct_dict = ct_dict[N$id]
    }

    X = N %>%
        select(-any_of(c('gene_i', 'cell'))) %>%
        tibble::column_to_rownames('id') %>%
        as.matrix

    X = X[rowSums(X) > 0,colSums(X) > 0]

    if (method == "vrnmf") {
        res = run_vrnmf(X, k)
        loadings = res$C %>% t() %>% set_colnames(colnames(X))
        scores = res$D
    } else {
        res = run_nmf_raw(X, k, method=method, ct_dict = ct_dict, nrun=nrun, ncores=ncores)
        loadings=res@fit@H
        scores = res@fit@W
    }

    H = loadings_to_long_df(loadings) %>% mutate(method = method)

    res = list(
        loadings = H,
        scores = scores
    )

    return(res)
}


scores_to_admixture_labels = function(scores, mol.df, cell.type) {

    main_comp = apply(scores, 1, which.max)
    is_admix_pred = main_comp %>% {. != which.max(table(.))}
    is_admix_real = mol.df %>% {.[match(names(is_admix_pred), .$mol_id),]} %>%
        {.$compartment != cell.type}

    return(data.frame(main_component=main_comp, is_admix_pred=is_admix_pred, is_admix_real=is_admix_real))
}

#### MRF segmentation #####

#' Get graph from cell
#' @param df_cell Data frame with molecule coordinates
#' @param h Neighborhood size
#' @param distinct_edges Whether to keep only distinct edges
#' @return Graph object
get_graph = function(df_cell, h = 10, distinct_edges = FALSE) {

    h = as.integer(h)

    if (h >= nrow(df_cell)) {
        h = min(h, nrow(df_cell) - 1)
        message('h is too large, setting to ', h, '...')
    }

    if (nrow(df_cell) <= h) {
        return(NULL)
    }

    xn = df_cell %>% select(x, y, z) %>% FNN::get.knn(k = h) %>% .$nn.index

    gene_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$gene, .$mol_id)}
    id_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$mol_id, .$id)}

    xn = apply(xn, 2, function(i){id_dict[i]})
    rownames(xn) = df_cell$mol_id

    edges = xn %>% as.data.frame() %>%
        tibble::rownames_to_column('x') %>%
        reshape2::melt(id.var = 'x', value.name = 'y') %>%
        select(x, y)

    G = graph_from_data_frame(edges, directed = FALSE)

    # only keep distinct edges
    if (distinct_edges) {
        G = simplify(G)
    }

    V(G)$gene = gene_dict[names(V(G))]

    V(G)$id = names(V(G))

    # V(G)$centrality = estimate_betweenness(G, directed = F, cutoff = -1)

    return(G)
}

#' Run MRF
#' @param df_cell Data frame with molecule coordinates
#' @param expr Data frame with gene expression
#' @param k Number of factors
#' @param h Neighborhood size
#' @param t Transition probability
#' @return Data frame with MRF segmentation
run_crf = function(df_cell, expr, k, h = 10, t = 0.05) {

    G = df_cell %>% get_graph(h = h)

    V = igraph::as_data_frame(G, 'vertices') %>%
        left_join(
            expr %>% reshape2::dcast(gene ~ factor, value.var = 'frac', fill = 0),
            by = join_by(gene)
        )

    adj = as_adjacency_matrix(G)

    crf = make.crf(adj, k)

    crf$node.pot = V %>%
        select(as.character(1:k)) %>%
        as.matrix

    A = matrix(rep(t, k^2), nrow = k)
    diag(A) = 1-(k-1)*t

    for (i in 1:length(crf$edge.pot)) {
        crf$edge.pot[[i]] = A
    }

    out <- decode.lbp(crf, verbose = T)

    df_out = data.frame(id = V$name, factor = out) %>%
        mutate(cell = unique(df_cell$cell))

    return(df_out)
}

