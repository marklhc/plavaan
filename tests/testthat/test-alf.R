test_that("alf is computed correctly", {
    x <- matrix(c(-2, -1, 0.5, 1, 2), ncol = 1)
    result <- alf(x, eps = 1e-6)
    expected <- sqrt(abs(x))
    expect_equal(result, expected, tolerance = 1e-6)
})

test_that("composite_pair_loss computes correct sum", {
    x1 <- rbind(
        c(1, 1.2, 1.2),
        c(1, 0.6, 0.6),
        c(1, 1.2, 0.9)
    )
    result <- composite_pair_loss(x1, fun = alf, eps = 1e-16)
    expected <- sum(
        sqrt(abs(x1[1, ] - x1[2, ])),
        sqrt(abs(x1[1, ] - x1[3, ])),
        sqrt(abs(x1[2, ] - x1[3, ]))
    )
    expect_equal(result, expected * 2 / 3, tolerance = 1e-3)
    res1 <- composite_pair_loss(x1[, 1], fun = alf)
    res2 <- composite_pair_loss(x1[, 2], fun = alf)
    res3 <- composite_pair_loss(x1[, 3], fun = alf)
    expect_true(res3 > res2 & res2 > res1)
})

test_that("composite_pair_loss handles missing data properly", {
    x2 <- rbind(
        c(1, NA, 1.2),
        c(1, 0.6, NA),
        c(1, 1.2, 0.9)
    )
    result <- composite_pair_loss(x2, fun = alf, eps = 1e-16)
    expected <- sum(
        sqrt(abs(x2[1, 1] - x2[2, 1])),
        sqrt(abs(x2[1, 3] - x2[3, 3])),
        sqrt(abs(x2[2, 2] - x2[3, 2]))
    )
    expect_equal(result, expected * 2 / 3, tolerance = 1e-3)
    g1 <- gr_cpl(x2, gr_l0a)
    g2 <- numDeriv::grad(
        function(x) {
            mm <- matrix(NA, nrow = nrow(x2), ncol = ncol(x2))
            mm[!is.na(x2)] <- x
            composite_pair_loss(mm, fun = l0a)
        },
        na.omit(as.vector(x2))
    )
    expect_equal(na.omit(g1), g2, ignore_attr = TRUE)
})

test_that("gr_cpl() computes the right gradient", {
    ld_mat <- structure(
        c(
            2.13324692100505,
            2.98948361240368,
            2.26448876139142,
            2.92469178270292,
            1.95718901005903,
            2.6885288670966,
            2.69142987862754,
            2.81722345109902
        ),
        dim = c(4L, 2L)
    )
    g1 <- gr_cpl(t(ld_mat), gr_alf)
    g2 <- numDeriv::grad(
        function(x) composite_pair_loss(matrix(x, nrow = 2), fun = alf),
        as.vector(t(ld_mat))
    )
    expect_equal(g1, g2)
    g3 <- gr_cpl(t(ld_mat), gr_l0a)
    g4 <- numDeriv::grad(
        function(x) composite_pair_loss(matrix(x, nrow = 2), fun = l0a),
        as.vector(t(ld_mat))
    )
    expect_equal(g3, g4)
    ld_mat2 <- structure(
        c(
            1.71976901143716,
            2.37695170791305,
            2.17073464443798,
            2.49999081186403,
            1.71881091139299,
            2.37464212814812,
            2.1917965406107,
            2.48314297929977
        ),
        dim = c(4L, 2L)
    )
    g5 <- gr_cpl(t(ld_mat2), gr_l0a, trans = log, gr_trans = function(x) 1 / x)
    g6 <- numDeriv::grad(
        function(x) composite_pair_loss(log(x), fun = l0a),
        t(ld_mat2)
    )
})

test_that("gr_cpl() computes the right gradient with multiple groups", {
    ld_mat <- structure(
        c(
            0.949629979618972,
            1.03082275557671,
            0.994691979624166,
            1.09084353768524,
            1.03546694498492,
            1.07190071490556,
            1.01064097248655,
            1.14899911389366,
            1.044765508461,
            1.12567318028583,
            1.01011416906684,
            1.10052636664253,
            1.02154340252047,
            1.17087864385868
        ),
        dim = c(2L, 7L)
    )
    g1 <- gr_cpl(t(ld_mat), gr_alf)
    g2 <- numDeriv::grad(
        function(x) composite_pair_loss(matrix(x, nrow = 7), fun = alf),
        as.vector(t(ld_mat))
    )
    expect_equal(g1, g2)
    g3 <- gr_cpl(t(ld_mat), gr_l0a)
    g4 <- numDeriv::grad(
        function(x) composite_pair_loss(matrix(x, nrow = 7), fun = l0a),
        as.vector(t(ld_mat))
    )
    expect_equal(g3, g4)
})
