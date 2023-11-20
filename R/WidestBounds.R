#' Widest identification bounds
#'
#' This function computes the widest identification bounds of the conditional expectation \eqn{E(Y\mid X)} when some values of \eqn{Y} are not observed.
#' We use the partial identifiability approach introduced by Manski (1989). Let \eqn{Z} be a binary random variable such that \eqn{Z=1} if \eqn{Y} is observed, and 0 if not.
#'
#' By the law of total probability, we know that
#' \deqn{E(Y\mid X)=E(Y\mid X,Z=0)P(Z=0\mid X)+E(Y\mid X,Z=1)P(Z=1\mid X)}
#' However, \eqn{E(Y\mid X,Z=0)} is not observed because it depends on the non observed values of \eqn{Y}.
#'
#' If \eqn{Y} is bounded by \eqn{Y0} and \eqn{Y1} (i.e., \eqn{Y0\leq Y\leq Y1}), then
#'
#' \deqn{Y0 P(Z=0\mid X)+E(Y\mid X,Z=1)P(Z=1\mid X)\leq E(Y\mid X)\leq Y1 P(Z=0\mid X)+E(Y\mid X,Z=1)P(Z=1\mid X)}
#'
#' In this function, the identification bounds for specific values for Y0 and Y1 are computed.
#' @param database data frame.
#' @param X vector with the number of the columns in database that are used as covariates X. Must be of length 1 or 2.
#' @param Y number of the column in database indicating the response variable.
#' @param Z number of the column in database that indicates the binary variable Z in data.
#' @param Y0 The minimum possible value in the range of Y.
#' @param Y1 The maximum possible value in the range of Y.
#' @param v dataframe indicating the values to make the prediction of Y. The number of columns must be the same as the length of X.
#' @param meanmodel model used to estimate the regression model \eqn{E(Y\mid X,Z=1)}. By default is a linear model (lm).
#' @param linkprob model used to estimate the probability \eqn{P(Z=1\mid X)}.
#'
#' @return
#' \code{data.bounds} A data frame that contains the observed values of X and Y. Additionally it contains predicted values of Y and the lower and upper bounds.
#'
#' \code{plot.bounds} A plot with the widest identification bounds.
#'
#' \code{widthBound} If v is provided, the width of the bound is shown for each row of v.
#' @export
#'
#'
#' @references Manski, C. (1989). Anatomy of the selection problem. The Journal of Human Resources, 24(3), 343–360.
#'
#' @examples
#' data("DataSelection")
#'
#' WidestBounds(DataSelection,X=c(6,5),Y=2,Z=7,Y0=1,Y1=7,meanmodel = "lm",linkprob = "probit", v=data.frame(c(650,700),c(650,600)))
#'
#' WidestBounds(DataSelection,X=6,Y=2,Z=7,Y0=1,Y1=7,meanmodel = "lm",linkprob = "probit", v=data.frame(c(650,700)))


WidestBounds <- function(database,X,Y,Z,Y0,Y1,
                        v = NULL,meanmodel = "lm",linkprob = "logit")
{
  if (is.character(Y0) | is.character(Y1))
    stop("Both Y0 and Y1 must be numeric")
  if (is.na(Y0) | is.na(Y1))
    stop("Both Y0 and Y1 must be numeric")
  if(Y0>=Y1)
    stop("Y0 must be lower than Y1")

  if(length(Y)!= 1 )
    stop("The length of Y must be 1")
  if(any(is.character(Y)))
    stop("Y must be numeric")
  if(!(Y %in% 1:dim(database)[2]))
    stop("Y must be a column number of database")

  nX = length(X)
  if(!(nX %in% 1:2) )
    stop("The length of X must be 1 or 2")
  if(any(is.character(X)))
    stop("X must be numeric")
  if(!(sum(X %in% 1:dim(database)[2])==nX))
    stop("X must be a column number of database")

  if(length(Z)!= 1 )
    stop("The length of Z must be 1")
  if(any(is.character(Z)))
    stop("Z must be numeric")
  if(!(Z %in% 1:dim(database)[2]))
    stop("Z must be a column number of database")

  if(any(as.numeric(table(c(X,Y,Z)))!=1))
    stop("X, Y and Z must be different numbers")

  database_names <- names(database)
  names(database)[Y] <- "Y"
  names(database)[Z] <- "Z"
  names(database)[X] <- paste0("X",1:nX)

  if(!is.numeric(database$Z))
    stop("Z must be numeric column")
  if( !any(names(table(database$Z)) %in% as.character(c(0,1))))
    stop("Z must be dichotomic (0-1) column of database")

  if(!is.null(v)) if(!is.data.frame(v)) stop("v must be a data frame")

  if(nX==1){
    if(!is.null(v)){
      if(dim(v)[2] !=nX)
        stop("The number of columns of v must be 1")
      names(v) <- "X1"
      if (any(is.character(v$X1)))
        stop("The vector v must be a numeric vector")
      if (any(is.na(v$X1)))
        stop("Some values of v are NA")
    }

    database = database %>%
      dplyr::arrange(dplyr::desc(X1))

    data.obs = database %>%
      dplyr::filter(Z == 1)

    data.XZ = database %>%
      dplyr::select(X1, Z)


    #predicted Y using observed data
    if(meanmodel == "lm") mean.Y = stats::lm(Y ~ X1, data = data.obs)
    predY = stats::predict(mean.Y)

    #model for p(z \mid x)
    if(linkprob == "logit") LogitZ = stats::glm(Z ~ X1,family=binomial(link='logit') ,data = data.XZ)
    if(linkprob == "probit") LogitZ = stats::glm(Z ~ X1,family=binomial(link='probit') ,data = data.XZ)
    fittedZ = LogitZ$fitted.values

    prob.Z1 = data.XZ %>%
      dplyr::mutate(fittedZ) %>%
      dplyr::filter(Z == 1) %>%
      dplyr::rename("prob.Z1"="fittedZ") %>%
      dplyr::select(prob.Z1)

    data.bounds = data.obs %>%
      dplyr::select(X1,Y) %>%
      dplyr::mutate(predY,prob.Z1 = prob.Z1$prob.Z1,
                    LB = round(predY * prob.Z1 + Y0 * (1 - prob.Z1),3),
                    UB = round(predY * prob.Z1 + Y1 * (1 - prob.Z1),3)) %>%
      dplyr::select(-prob.Z1) %>%
      dplyr::arrange(X1)


    plot.bounds = ggplot2::ggplot() +
      ggplot2::geom_point(data = data.bounds,ggplot2::aes(x = X1, y = Y, shape = "Observed Y")) +
      ggplot2::ylim(Y0,Y1) +
      ggplot2::ylab(database_names[Y]) +
      ggplot2::xlab(database_names[X]) +
      ggplot2::geom_line(data = data.bounds,ggplot2::aes(x = X1, y = predY, color = "Predicted Y in observed data", linetype = "Predicted Y in observed data"),linewidth = 1.1) +
      ggplot2::geom_line(data = data.bounds, ggplot2::aes(x = X1, y = LB, color = "Lower bound for predicted Y",linetype = "Lower bound for predicted Y"), linewidth = 1.1) +
      ggplot2::geom_line(data = data.bounds, ggplot2::aes(x = X1, y = UB, color = "Upper bound for predicted Y", linetype = "Upper bound for predicted Y"), linewidth = 1.1) +
      ggplot2::scale_shape(name = "", breaks = "Observed Y") +
      ggplot2::scale_color_manual(name = "", breaks = c("Upper bound for predicted Y",
                                                        "Predicted Y in observed data", "Lower bound for predicted Y"),
                                  values = c(`Upper bound for predicted Y` = "darkgreen",
                                             `Predicted Y in observed data` = "blue",
                                             `Lower bound for predicted Y` = "darkgreen")) +
      ggplot2::scale_linetype_manual(name = "", breaks = c("Upper bound for predicted Y",
                                                           "Predicted Y in observed data", "Lower bound for predicted Y"),
                                     values = c(`Upper bound for predicted Y` = "dotted",
                                                `Predicted Y in observed data` = "solid",
                                                `Lower bound for predicted Y` = "dotdash")) +
      ggplot2::theme_bw() + ggplot2::theme(legend.direction = "vertical",
                                           legend.position = "right", legend.box = "vertical", legend.key = ggplot2::element_blank(),
                                           legend.background = ggplot2::element_rect(fill = "white",
                                                                                     colour = "white"))


    if(!is.null(v)){
      if (any(v$X1 >= max(data.obs$X1) | v$X1 <= min(data.obs$X1)))
        stop(paste("Some values of v are out of the range of covariable X"))

      predictV = stats::predict(LogitZ, list(X1 = v$X1), type = "response")
      predictM = stats::predict(mean.Y, newdata = list(X1 = v$X1))
      UpperBound = predictM * predictV + Y1 * (1 - predictV)
      LowerBound = predictM * predictV + Y0 * (1 - predictV)
      width = UpperBound - LowerBound
      df_unc = data.frame(v, LowerBound, UpperBound, width) %>%
        dplyr::arrange(v)


      names(df_unc)[1] <- database_names[X]
    }

    names(data.bounds)[1:2] <- database_names[c(X,Y)]

    list.return = list(data.bounds %>% `rownames<-`(NULL), plot.bounds)
    names(list.return) = c("Data bounds", "Plot bounds")

    if(!is.null(v)){
      list.return = list(data.bounds %>% `rownames<-`(NULL), plot.bounds,df_unc)
      names(list.return) = c("Data bounds", "Plot bounds","widthBound")
    }
  }

  if(nX==2){
    if(!is.null(v)){
      if(dim(v)[2] != nX)
        stop("The number of columns of v must be 2")
      names(v) <- c("X1","X2")
      if (any(is.character(v$X1)))
        stop("The data frame v must contain numeric elements")
      if (any(is.character(v$X2)))
        stop("The data frame v must contain numeric elements")
      if (any(is.na(v)))
        stop("Some values of v are NA")
    }

    database = database %>%
      dplyr::arrange(dplyr::desc(X1),dplyr::desc(X2))

    data.obs = database %>%
      dplyr::filter(Z == 1)

    data.XZ = database %>%
      dplyr::select(X1,X2, Z)


    #predicted Y using observed data
    if(meanmodel == "lm") mean.Y = stats::lm(Y ~ X1+X2, data = data.obs)
    predY = stats::predict(mean.Y)

    #model for p(z \mid x)
    if(linkprob == "logit") LogitZ = stats::glm(Z ~ X1+X2,family=binomial(link='logit') ,data = data.XZ)
    if(linkprob == "probit") LogitZ = stats::glm(Z ~ X1+X2,family=binomial(link='probit') ,data = data.XZ)
    fittedZ = LogitZ$fitted.values

    #grid values plot
    grid.lines = 20
    x1seq <- seq(min(data.obs$X1), max(data.obs$X1), length.out = grid.lines)
    x2seq <- seq(min(data.obs$X2), max(data.obs$X2), length.out = grid.lines)
    x1x2 <- expand.grid( X1 = x1seq, X2 = x2seq)


    predY.x1x2 <- matrix(stats::predict(mean.Y, newdata = x1x2),
                         nrow = grid.lines, ncol = grid.lines)

    prob.Z1.x1x2 = matrix(stats::predict(LogitZ, newdata = list( X1= x1x2$X1,X2= x1x2$X2), type="response"),
                          nrow = grid.lines, ncol = grid.lines)

    LB=predY.x1x2*prob.Z1.x1x2+Y0*(1-prob.Z1.x1x2)
    UB=predY.x1x2*prob.Z1.x1x2+Y1*(1-prob.Z1.x1x2)


    prob.Z1 = data.XZ %>%
      dplyr::mutate(fittedZ) %>%
      dplyr::filter(Z == 1) %>%
      dplyr::rename("prob.Z1"="fittedZ") %>%
      dplyr::select(prob.Z1)

    data.bounds = data.obs %>%
      dplyr::select(X1,X2,Y) %>%
      dplyr::mutate(predY,prob.Z1 = prob.Z1$prob.Z1,
                    LB = round(predY * prob.Z1 + Y0 * (1 - prob.Z1),3),
                    UB = round(predY * prob.Z1 + Y1 * (1 - prob.Z1),3)) %>%
      dplyr::select(-prob.Z1) %>%
      dplyr::arrange(X1,X2)


    #3D plot
    plot_and_list <- function() {
      plotlist <- list()

      par(mar = c(5, 4, 4, 2) - 1)
      plot3D::scatter3D(data.obs$X1, data.obs$X2, data.obs$Y,
                pch = 19, cex = 0.5, colvar = NULL, col = "black",
                theta = 20, phi = 20, bty = "b", zlim = c(1, 7),
                xlab = database_names[X[1]], ylab = database_names[X[2]], zlab = database_names[Y],
                surf = list(x = x1seq, y = x2seq, z = predY.x1x2,
                            facets = TRUE,
                            col = plot3D::ramp.col(col = c("blue", "blue"), alpha = 0.3)),
                main = paste("Predicted", database_names[Y])
      )
      par(mar = c(5, 4, 4, 2) + 0.1)
      plot3D::persp3D(x = x1seq, y = x2seq, z = LB,
              colvar = NULL, add = TRUE, col = "darkgreen", alpha = 0.5)
      plot3D::persp3D(x = x1seq, y = x2seq, z = UB,
              add = TRUE, colvar = NULL, col = "darkgreen", alpha = 0.5)
      legend("bottomleft",
             legend=c("Identification bounds under",
                      "weakly-informative assumption",
                      paste("Predicted", database_names[Y], "in selected")),
             col=c("darkgreen","white","blue","white"), pch=15,
             cex=0.65,bg = "transparent")

      return(plotlist)
    }

    names(data.bounds)[1:3] <- database_names[c(X,Y)]
    list.return <- list(data.bounds %>% `rownames<-`(NULL),
                        PlotBound = plot_and_list())

    names(list.return) = c("Data bounds", "Plot bounds")

    # width of bounds for specific X1 and X2 values
    if(!is.null(v)){

      if (any(v$X1 >= max(data.obs$X1) | v$X1 <= min(data.obs$X1)))
        stop(paste("Some values of v are out of the range of covariable X"))

      if (any(v$X2 >= max(data.obs$X2) | v$X2 <= min(data.obs$X2)))
        stop(paste("Some values of v are out of the range of covariable X"))

      predictV = stats::predict(LogitZ, newdata = list(X1=v$X1, X2=v$X2), type="response")
      predictM = stats::predict(mean.Y,newdata = list(X1=v$X1, X2=v$X2))


      UpperBound = predictM * predictV + Y1 * (1 - predictV)
      LowerBound = predictM * predictV + Y0 * (1 - predictV)
      width = UpperBound - LowerBound

      df_unc = data.frame(v, width) %>%
        dplyr::arrange(v)

      names(df_unc)[1:2] <- database_names[X]

      list.return <- list(data.bounds,PlotBound = plot_and_list(),`Width prediction` = df_unc)
      names(list.return) = c("Data bounds", "Plot bounds","widthBound")
    }
  }


  return(list.return)
}




