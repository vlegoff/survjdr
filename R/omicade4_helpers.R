predict_mcia = function(x_old, x_new, mcia_obj) {

    recalculer = function(tab, scorcol) {
        for (k in 1:nblock) {
            soustabk = tab[, indicablo==veclev[k]]
            uk = scorcol[indicablo==veclev[k]]
            #uk = uk/norm(uk, "2")
            soustabk.hat = t(apply(soustabk, 1, function(x) sum(x*uk)*uk))
            soustabk = soustabk - soustabk.hat
            tab[, indicablo==veclev[k]] = soustabk
        }
        return(tab)
    }
    normaliserparbloc = function(scorcol) {
        for (k in 1:nblock) {
            w1 = scorcol[indicablo==veclev[k]]
            w2 = sqrt(sum(w1*w1))
            w1 = w1/w2
            scorcol[indicablo==veclev[k]] = w1
        }
        return(scorcol)
    }

    # x_old should be ktab object
    indicablo_col = mcia_obj$TC[,1]
    indicablo_row = mcia_obj$TL[,1]
    indicablo = indicablo_col
    veclev = unique(indicablo)
    nblock = length(veclev)
    cw = x_old$cw
    lw = rep(x_old$lw[1], nrow(x_new))
    lw_old = x_old$lw
    
    #normalizing per block
    option = mcia_obj$call$option
    if (is.null(option)) option = "inertia"
    if (option == "internal") {
        if (is.null(x_old$tabw)) {
            warning("internal weights not found (prediction); uniform weights are used")
            option = "uniform"
        }
    }
    Xsepan = ade4::sepan(x_old, nf=4)
    rank.fac = factor(rep(1:nblock, Xsepan$rank))
    tabw = vector("numeric", nblock)
    func_map = list(
        "lambda1" = \(i) return(1/Xsepan$Eig[rank.fac == i][1]),
        "inertia" = \(i) return(1/sum(Xsepan$Eig[rank.fac==i])),
        "uniform" = \(i) return(1),
        "internal" = \(i) return(X$tabw[i])
    )
    for (i in 1:nblock) {
        tabw[i] = func_map[[option]](i)
        x_old[[i]] = sqrt(tabw[i]) * x_old[[i]]
        x_new[, indicablo_col==veclev[i]] = x_new[, indicablo_col==veclev[i]] *
            sqrt(tabw[i])
    }
    
    mat_old = as.matrix(Reduce(cbind, lapply(1:nblock, \(i) x_old[[i]])))
    mat_old = mat_old * sqrt(lw_old)
    mat_old = t(t(mat_old) * sqrt(cw))
    x_new = x_new * sqrt(lw)
    x_new = t(t(x_new) * sqrt(cw))
    y_sinvar = matrix(0, nrow(x_new), mcia_obj$nf)
    
    for (i in 1:mcia_obj$nf) {
        
        a_synvar = t(t(mcia_obj$SynVar[,i]) %*% mat_old)
        y_sinvar[,i] = x_new %*% a_synvar
        y_sinvar[,i] = y_sinvar[,i]/norm(sqrt(lw)*y_sinvar[,i], "2")
        
        a_synvar = normaliserparbloc(a_synvar)
        x_new = recalculer(x_new, a_synvar)
        mat_old = recalculer(mat_old, a_synvar)
    }
    return(list(SynVar = y_sinvar))
}

nsc_transform = function(df) {
    df = as.matrix(df)
    col = ncol(df)
    N = sum(df)
    row.w = apply(df, 1, sum)/N
    col.w = apply(df, 2, sum)/N
    df = t(apply(df, 1, function(x) if (sum(x)==0) col.w else x/sum(x)))
    df = sweep(df, 2, col.w)
    df = ncol(df)*df
    return(list(tab=data.frame(t(df)), row.w=rep(1, col)/col, col.w=row.w))
}

predict_omicade4 = function(x_new, mcia_omicade4) {
    x_old = ade4::ktab.list.dudi(lapply(mcia_omicade4$coa, ade4:::t.dudi))
    x_new = lapply(x_new, omicade4:::array2ade4, pos=TRUE)
    x_new = lapply(x_new, \(x) nsc_transform(x)$tab)
    x_new = as.matrix(Reduce(cbind, x_new))
    return(predict_mcia(x_old, x_new, mcia_omicade4$mcoa))
}
