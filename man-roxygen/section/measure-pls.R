#' @section measure-pls: 
#' Two measures of accuracy are available: Correlation (\code{cor}), as well as
#' the Residual Sum of Squares (\code{RSS}). For \code{cor}, the parameters
#' which would maximise the correlation between the predicted and the actual
#' components are chosen. The \code{RSS} measure tries to predict the held-out
#' data by matrix reconstruction and seeks to minimise the error between actual
#' and predicted values. For \code{mode='canonical'}, The X matrix is used to
#' calculate the \code{RSS}, while for others modes the \code{Y} matrix is used.
#' This measure gives more weight to any large errors and is thus sensitive to
#' outliers. It also intrinsically selects less number of features on the
#' \code{Y} block compared to \code{measure='cor'}.
##Four measures of accuracy are available: Mean Absolute Error (\code{MAE}),
##Mean Square Error(\code{MSE}), \code{Bias} and \code{R2}. Both MAE and MSE
##average the model prediction error. MAE measures the average magnitude of
##the errors without considering their direction. It is the average over the
##fold test samples of the absolute differences between the Y predictions and
##the actual Y observations. The MSE also measures the average magnitude of
##the error. Since the errors are squared before they are averaged, the MSE
##tends to give a relatively high weight to large errors. The Bias is the
##average of the differences between the Y predictions and the actual Y
##observations and the R2 is the correlation between the predictions and the
##observations. All those measures are averaged across all Y variables in the
##PLS2 case.