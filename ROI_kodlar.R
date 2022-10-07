#ROI ile ilgili tüm paketler 

library(ROI)
library(ROI.models.globalOptTests)
library(ROI.models.miplib)
library(ROI.models.netlib)
library(ROI.plugin.alabama)
library(ROI.plugin.clp)
library(ROI.plugin.deoptim)
library(ROI.plugin.ecos)
library(ROI.plugin.glpk)
library(ROI.plugin.ipop)
library(ROI.plugin.lpsolve)
library(ROI.plugin.mosek)
library(ROI.plugin.msbinlp)
library(ROI.plugin.neos)
library(ROI.plugin.nloptr)
library(ROI.plugin.optimx)
library(ROI.plugin.osqp)
library(ROI.plugin.qpoases)
library(ROI.plugin.quadprog)
library(ROI.plugin.scs)
library(ROI.plugin.symphony)
library(ompr.roi)


con<-rbind(c(5,7,2),
           c(3,2,-9),
           c(1,3,1))

con_dir<-c("<=","<=","<=")

rhs<-c(61,35,31)

lp<-OP(objective = L_objective(c(3,7,-12)),
       constraints = L_constraint(con,
                                  dir = con_dir,
                                  rhs = rhs),
       bounds = V_bound(li = 3,
                        ui = 3,
                        lb = -10,
                        ub = 10,
                        nobj = 3),
       maximum = TRUE)
lp

#Alternatif kodlama

lp<-OP()

objective(lp)<-L_objective(c(3,7,-12))

constraints(lp)<-L_constraint(con,
                              dir = con_dir,
                              rhs = rhs)

bounds(lp)<-V_bound(li = 3,
                    ui = 3,
                    lb = -10,
                    ub = 10,
                    nobj = 3)

maximum(lp)<-TRUE

lp

param<-rep.int(1, length(objective(lp)))

objective(lp)(param)

terms(objective(lp))

(ROI_solve(lp))

# ROI Çözücüleri

ROI, birden çok optimizasyon çözücüsü için tek bir arayüz oluþturan bir R paketidir. Optimizasyon problemleri için farklý çözücüleri kolayca deðiþtirmek isteyen R kullanýcýlarýna büyük kolaylýk saðlar. 

ROI_installed_solvers()

ROI_applicable_solvers(lp)

(lp_sol<-ROI_solve(lp, solver = "alabama", start=double(3)))

(lp_sol<-ROI_solve(lp, solver = "glpk"))

solution(lp_sol, type = "aux")

solution(lp_sol, type = "msg")

# Ýkili programlama ile ilgili örnek
  
blp<-OP(objective = L_objective(c(-1,-1,-1,-1,-99)),
        constraints = L_constraint(rbind(c(1,1,0,0,0), c(0,0,1,1,0), c(0,0,0,1,1)),
                                   dir = c("<=","<=","<="),
                                   rhs = rep.int(1,3)),
        types = rep("B",5L))

blp_sol<-ROI_solve(blp,
                   solver = "glpk",
                   method = "msbinlp",
                   nsol_max = 32)

blp_sol

solution(blp_sol)

lp_inaccurate_sol<-ROI_solve(lp,
                             solver = "glpk",
                             tol = 1e-32)

solution(lp_inaccurate_sol, force = TRUE)

solution(lp_sol, type = "dual")

solution(lp_sol, type = "aux")

solution(lp_sol, type = "msg")

solution(lp_sol, type = "objval")

solution(lp_sol, type = "status")

solution(lp_sol, type = "status_code")

# Optimizasyon Problemlerini Okuma ve Yazma

lp_file<-tempfile()

ROI_write(lp, lp_file, "lp_lpsolve")

writeLines(readLines(lp_file))

ROI_read(lp_file, "lp_lpsolve")


tmpfile <- tempfile()
con <- gzcon(url("http://www.zib.de/koch/perplex/data/netlib/mps/capri.mps.gz"))
writeLines(readLines(con), tmpfile)
close(con)
(capri <- ROI_read(tmpfile, type="mps_fixed", "lpsolve"))

(sol <- ROI_solve(capri))

# Test Koleksiyonlarý
netlib()

netlib("boeing1")

problem_name<-"boeing1"

model<-netlib(problem_name)

ROI_solve(model)

#Ýkinci dereceden programlama ile ilgili örnek
qp<-OP(Q_objective(Q = diag(2), L = c(-1,0)),
       L_constraint(c(4,6), ">=", 10))

qp_sol<-ROI_solve(qp, "qpoases")

qp_sol

solution(qp_sol)

#Ek kýsýt eklememiz gerekirse

#7x_1 + 11x_2 - x_1^2 - x_2^2 ??? 40
  
qcqp<-qp

constraints(qcqp)<-rbind(constraints(qp),
                         Q_constraint(-diag(c(2,2)),
                                      L = c(7,11),
                                      ">=", 40))

ROI_applicable_solvers(qcqp)

qcqp_sol<-ROI_solve(qcqp, "alabama", start = c(5,5))

solution(qcqp_sol, type = "msg")

#Doðrusal olmayan programlama ile ilgili örnek

nlp1<-OP(maximum = TRUE, 
         bounds = V_bound(ud = 42, nobj = 3L))

objective(nlp1)<-F_objective(F = function(x) prod(x),
                             n =3,
                             G = function(x) c(prod(x[-1]),
                                               prod(x[-2]),
                                               prod(x[-3])))

constraints(nlp1)<-F_constraint(F = function(x) x[1]+2*x[2]+2*x[3], 
                                dir = "<=", 
                                rhs = 72,
                                J= function(x)c(1,2,2))

nlp1

nlp1_sol1<-ROI_solve(nlp1, "alabama", 
                     start = c(10,10,10))

nlp1_sol1

solution(nlp1_sol1)

nlp1_sol2<-ROI_solve(nlp1, "alabama",
                     start = double(3))
nlp1_sol2

solution(nlp1_sol2, type = "msg")

#Karýþýk tamsayýlý programlama ile ilgili örnek

A<-rbind(c(5,3),c(7,1))

milp<-OP(c(5,7),
         constraints = L_constraint(A, c(">=",">="),
                                    c(7,9)),
         types = c("B", "I"))

milp_sol<-ROI_solve(milp, solver = "glpk")

milp_sol

solution(milp_sol)

# Doðrusal Regresyon Problemi
create_L1_problem <- function(x, y) {
  p <- ncol(x) + 1L
  m <- 2 * nrow(x)
  L <- cbind(1, x, diag(nrow(x)), -diag(nrow(x)))
  bnds <- V_bound(li = seq_len(p), lb = rep(-Inf, p), nobj = p+m)
  OP(objective = L_objective(c(rep(0, p), rep(1, m))),
     constraints = L_constraint(L, dir = eq(nrow(x)), rhs = y),
     bounds = bnds)
}

data("stackloss", package = "datasets")

l1p<-create_L1_problem(x = as.matrix(stackloss)[,-4],
                       y = stackloss[,4])


L1_res<-ROI_solve(l1p, solver = "glpk")

roi_sol<-solution(L1_res)[1:ncol(stackloss)]

# lmtest paketiyle regresyonu çözüp model katsayýlarýný karþýlaþtýralým

attach(stackloss)
lmodel<-lm(stack.loss~Air.Flow+Water.Temp+Acid.Conc., data = stackloss)
lm_sol<-lmodel$coefficients

difference<-(roi_sol-lm_sol)

cbind(roi_sol,lm_sol, difference)

# Gezgin Satýcý Problemi

library(knitr)
library(dplyr)
library(ggplot2)
n<-10

max_x<-500
max_y<-500

set.seed(123456)
cities <- data.frame(id = 1:n, x = runif(n, max = max_x), y = runif(n, max = max_y))

ggplot(cities, aes(x,y)) + geom_point()

distance <- as.matrix(stats::dist(select(cities, x, y), diag = TRUE, upper = TRUE))
dist_fun <- function(i, j) {
  vapply(seq_along(i), function(k) distance[i[k], j[k]], numeric(1L))
}

library(ompr)

model <- MIPModel() %>%
  # we create a variable that is 1 iff we travel from city i to j
  add_variable(x[i, j], i = 1:n, j = 1:n, 
               type = "integer", lb = 0, ub = 1) %>%
  
  # a helper variable for the MTZ formulation of the tsp
  add_variable(u[i], i = 1:n, lb = 1, ub = n) %>% 
  
  # minimize travel distance
  set_objective(sum_expr(dist_fun(i, j) * x[i, j], i = 1:n, j = 1:n), "min") %>%
  
  # you cannot go to the same city
  set_bounds(x[i, i], ub = 0, i = 1:n) %>%
  
  # leave each city
  add_constraint(sum_expr(x[i, j], j = 1:n) == 1, i = 1:n) %>%
  #
  # visit each city
  add_constraint(sum_expr(x[i, j], i = 1:n) == 1, j = 1:n) %>%
  
  # ensure no subtours (arc constraints)
  add_constraint(u[i] >= 2, i = 2:n) %>% 
  add_constraint(u[i] - u[j] + 1 <= (n - 1) * (1 - x[i, j]), i = 2:n, j = 2:n)

model

library(ompr.roi)
library(ROI.plugin.glpk)

result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))

result
head(result$solution)

solution <- get_solution(result, x[i, j]) %>% 
  filter(value > 0) 
kable(head(solution, 3))

paths <- select(solution, i, j) %>% 
  rename(from = i, to = j) %>% 
  mutate(trip_id = row_number()) %>% 
  tidyr::gather(property, idx_val, from:to) %>% 
  mutate(idx_val = as.integer(idx_val)) %>% 
  inner_join(cities, by = c("idx_val" = "id"))

kable(head(arrange(paths, trip_id), 4))

ggplot(cities, aes(x, y)) + 
  geom_point() + 
  geom_line(data = paths, aes(group = trip_id)) + 
  ggtitle(paste0("Optimal route with cost: ", round(objective_value(result), 2)))

# Sudoku yaratma ve çözme


library(sudokuAlt)

sudoku_puzzle <- makeGame()
plot(sudoku_puzzle)

library(slam)
library(ROI)
library(ROI.plugin.msbinlp)
library(ROI.plugin.glpk)

as.matrix.sudoku <- function(x) matrix(as.numeric(x), 9, 9)

to_col_index <- function(i, j, v) {
    (v - 1) * 81 + (i - 1) * 9 + j
}

index_to_triplet <- function(idx) {
    .index_to_triplet <- function(idx) {
        v <- (idx - 1) %/% 81 + 1
        idx <- idx - (v - 1) * 81
        i <- (idx - 1)%/% 9 + 1
        idx <- idx - (i - 1) * 9
        c(i = i, j = idx, v = v)
    }
    t(sapply(idx, .index_to_triplet))
}

solve_sudoku <- function(M, solver, solve = TRUE, ...) {
    stm <- simple_triplet_matrix
    seq9_9 <- rep(1:9, 9)
    seq9_9e <- rep(1:9, each = 9)
    seq81_9e <- rep(seq_len(81), each = 9)
    ones729 <- rep.int(1, 9^3)
    M <- as.matrix(M)
    M[is.na(M)] <- 0
    M <- as.simple_triplet_matrix(M)
    ## setup OP
    op <- OP(double(9^3))
    ## basic setting (fixed coefficients)
    j <- mapply(to_col_index, M$i, M$j, M$v)
    nfv <- length(M$i) ## number of fixed variables
    A0 <- stm(i = seq_len(nfv), j = j, v = rep.int(1, nfv), ncol = 9^3)
    LC0 <- L_constraint(A0, eq(nfv), rep.int(1, nfv))
    ## sum_{i=1:n} x_ijk = 1,  j = 1:n, k = 1:n
    only_one_k_in_each_column <- function(j, k) {
        sapply(1:9, function(i) to_col_index(i, j, k))
    }
    j <- unlist(mapply(only_one_k_in_each_column, seq9_9e, seq9_9, SIMPLIFY = FALSE))
    A1 <- stm(i = seq81_9e, j = j, v = ones729, nrow = 81, ncol = 9^3)
    ## sum_{j=1:n} x_ijk = 1
    only_one_k_in_each_row <- function(i, k) {
        sapply(1:9, function(j) to_col_index(i, j, k))
    }
    j <- unlist(mapply(only_one_k_in_each_row, seq9_9e, seq9_9, SIMPLIFY = FALSE))
    A2 <- stm(i = seq81_9e, j = j, v = ones729, nrow = 81, ncol = 9^3)
    only_one_k_in_each_submatrix <- function(blocki, blockj, k) { 
        i <- (blocki - 1) * 3 + 1:3
        j <- (blockj - 1) * 3 + 1:3
        coo <- expand.grid(i = i, j = j, v = k)
        mapply(to_col_index, i = coo$i, j = coo$j, v = coo$v)
    }
    coo <- expand.grid(i = 1:3, j = 1:3, k = 1:9)
    j <- unlist(mapply(only_one_k_in_each_submatrix, 
                       blocki = coo$i, blockj = coo$j, k  = coo$k, SIMPLIFY = FALSE))
    A3 <- stm(i = seq81_9e, j = j, v = ones729, ncol = 9^3)
    ## at every position in the matrix must be one value
    fill_matrix <- function(i, j) {
        sapply(1:9, function(k) to_col_index(i, j, k))
    }
    j <- unlist(mapply(fill_matrix, i = seq9_9e, j = seq9_9, SIMPLIFY = FALSE))
    A4 <- stm(i = seq81_9e, j = j, v = ones729, ncol = 9^3)
    A <- rbind(A1, A2, A3, A4)
    LC1 <- L_constraint(A, eq(nrow(A)), rep.int(1, nrow(A)))
    constraints(op) <- rbind(LC0, LC1)
    types(op) <- rep.int("B", 9^3)
    if (!solve) return(op)
    s <- ROI_solve(op, solver = solver, ...)
    sol <- solution(s)
    to_sudoku_solution <- function(sol) {
        coo <- index_to_triplet(which(as.logical(sol)))
        sudoku_solution <- as.matrix(stm(coo[,1], coo[,2], coo[,3]), nrow = 9, ncol = 9)
        structure(sudoku_solution, class = c("sudoku", "matrix"))
    }
    
    if ( any(lengths(sol) > 1L) & length(sol) > 1L ) {
        lapply(solution(s), to_sudoku_solution)
    } else {
        if ( length(sol) == 1L ) {
            to_sudoku_solution(sol[[1L]])
        } else {
            to_sudoku_solution(sol)
        }
    }
}

sudoku_is_valid_solution <- function(x) {
    .sudoku_is_valid_solution <- function(x) {
        stopifnot(inherits(x, "sudoku"))
        seq19 <- seq_len(9)
        for (i in seq19) {
            if ( any(sort(as.vector(x[i, ])) != seq19) ) return(FALSE)
            if ( any(sort(as.vector(x[, i])) != seq19) ) return(FALSE)
        }
        for (i in 1:3) {
            for (j in 1:3) {
                block <- x[(i-1) * 3 + 1:3, (j-1) * 3 + 1:3]
                if ( any(sort(as.vector(block)) != seq19) ) return(FALSE)
            }
        }
        return(TRUE)
    }
    if ( is.list(x) ) {
        sapply(x, .sudoku_is_valid_solution)
    } else {
        .sudoku_is_valid_solution(x)
    }
}

sudoku_solution <- solve_sudoku(sudoku_puzzle, solver = "glpk")
sudoku_is_valid_solution(sudoku_solution)

sudoku_solution

has_unique_solution <- inherits(sudoku_solution, "sudoku")

sudoku_solution <- solve_sudoku(sudoku_puzzle, solver = "msbinlp", nsol_max = 10L)






































