#' @title Fuzzy Cognitive Maps (FCMs) Inference
#' @description Estimates the Inference of a Fuzzy Cognitive Map.  FCM is proven to be capable of causal inference and is applicable to complex decision problems where numerous interlinked dependent variables influence one another.
#'
#' @param activation_vec 1 x m dataframe which contains the initial concept values. A concept is turned on or activated by making its vector element 1 or 0 or in [0, 1].
#' @param weight_mat m x m dataframe which stores the weights assigned to the pairs of concepts. The weights are usually normalized to the interval [0, 1] or [−1, +1].
#' @param iter The required number of iterations in order to reach the FCM convergence. Defaults to 20.
#' @param infer Select an Inference Rule ('k' Kosko, 'mk' modified Kosko, 'r' Rescale,'kc' Kosko-clamped, 'mkc' modified Kosko-clamped or 'rc' Rescale-clamped). Default value is set to 'k'
#' @param transform Contains the Transformation functions ('b' Bivalent,  'tr' Trivalent,  's' Sigmoid or 't' Hyperbolic tangent). The transformation function is used to reduce unbounded weighted sum to a certain range, which hinders quantitative analysis, but allows for qualitative comparisons between concepts. Default value is set equal to 's'.
#' @param lambda A parameter that determines the steepness of the sigmoid and hyperbolic tangent function at values around 0. Different lambda value may perform more appropriate for different problems.
#' @param e Epsilon (e) is a residual, describing the minimum error difference among the subsequent concepts. Its value depends on the application type. Defaults to to 0.001.
#'
#' @return Returns iter x m dataframe which contains the concepts' values of each iteration after the the transformation function.
#' @export
#' @author Zoumpolia Dikopoulou <dikopoulia@gmail.com>, <zoumpolia.dikopoulou@uhasselt.be>
#' @author Elpiniki Papageorgiou <epapageorgiou@teiste.gr>, <e.i.papageorgiou75@gmail.com>
#'
#' @references B. Kosko, "Fuzzy cognitive maps", International Journal of Man-Machine Studies 24, p.p. 65-75, 1986.
#' @references Groumpos, P.P, Stylios, C.D.; "Modelling supervisory control systems using fuzzy cognitive maps", Chaos, Solitons & Fractals, Volume 11, Issues 1–3, p.p. 329–336, 2000.
#' @references Papageorgiou E.I., "Fuzzy Cognitive Maps for Applied Sciences and Engineering From Fundamentals to Extensions and Learning Algorithms", Intelligent Systems Reference Library, Volume 54, 2014.
#' @references Papageorgiou E.I., Stylios C.D., GroumposP.P. , "Unsupervised learning techniques for finetuning fuzzy cognitive map causal links.", Int. J. Human Comput. Stud. Vol. 64, pp. 727–743, 2006.
#' @examples \dontrun{
#' # Example for the FCM inference
#'
#'
#' ### Input data
#'
#' act.vec <- data.frame(1, 1, 1, 0, 0, 0)    # Create the activation vector
#' colnames(act.vec) <- c("C1", "C2", "C3", "C4", "C5", "C6")  # Change the column names
#'
#' C1 = c(0.0, 0.0, 0.6, 0.9, 0.0, 0.0)
#' C2 = c(0.1, 0.0, 0.0, 0.0, 0.0, 0.0)
#' C3 = c(0.0, 0.7, 0.0, 0.0, 0.9, 0.0)
#' C4 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.9)
#' C5 = c(0.0, 0.0, 0.0, 0.0, 0.0, -0.9)
#' C6 = c(-0.3, 0.0, 0.0, 0.0, 0.0, 0.0)
#' C7 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.8)
#'
#' # Create the weight matrix
#' w.mat <- matrix(c(C1, C2, C3, C4, C5, C6, C7), nrow = 6, ncol = 6, byrow = TRUE)
#' colnames(w.mat) <- c("C1", "C2", "C3", "C4", "C5", "C6")
#' w.mat <- as.data.frame(w.mat)
#'
#'
#'
#' ### Select the arguments for the fcm.infer function
#'
#' output <- fcm.infer(act.vec, w.mat, 25, "r", "s")
#' View(output$values)          # View the concept values for each iteration
#'
#'
#' ### Visualize the concepts' values for each state
#'
#' library (reshape2)
#' library (ggplot2)
#' iterations <- as.numeric(rownames(output$values))  # create a numeric vector named "iterations"
#' df <- data.frame(iterations, output$values)   # add "iterations" in the "output$values" dataframe
#' df2 <- melt(df, id="iterations")   #transform df into long formats
#' ggplot(data=df2,
#'        aes(x=iterations, y=value, group=variable, colour=variable)) +
#'        theme_bw() + geom_line(size=0.7) + geom_point(size = 3)
#'
#' }



  fcm.infer <- function (activation_vec, weight_mat, iter = 20, infer = 'k', transform = 's', lambda = 1, e = 0.001) {



    # ------------------------------------------ checks on function input ------------------------------------------------------------------------------------ #

    # Check the values of the activation vector
    if (length(which(activation_vec > 1)) & length(which(activation_vec > -1))) {
      stop ("Please check the concepts' values of the activation vector. They must be in the range -1 and 1.")
    }


    # Check the weights of the matrix
    if (length(which(weight_mat > 1)) & length(which(weight_mat > -1)) ) {
      stop ("Please check the weights of the matrix. They must be in the range -1 and 1.")
    }


    # Check for missing values
    if (sum(is.na(activation_vec)) > 0) {
      stop ("Please check the activation vector for missing values.")
    }


    if (sum(is.na(weight_mat)) > 0) {
      stop ("Please check the weight matrix for missing values.")
    }


    # Check the variable of the transformation function
    if(iter <= 0 ) stop ("The iterations must be higher than zero.")


    # Check the variable of the Inference Rule
    if(sum(!(infer %in% c('k', 'mk', 'r', 'kc', 'mkc', 'rc'))) > 0) stop ("For the Inference Rule only Kosko 'k', modified Kosko 'mk',  Rescale 'r', Kosko-clamped 'kc', modified Kosko-clamped 'mkc' or Rescale-clamped 'rc' variables are allowed.")


    # Check the variable of the transformation function
    if(sum(!(transform %in% c('b', 'tr', 's', 't'))) > 0)
      stop ("For the transformation functions only Bivalent 'b', Trivalent 'tr', Sigmoid 's' or
            Hyperbolic tangent 't' variables are allowed.")


    # Check the variable of the lambda value
    if((lambda <= 0) || (lambda >= 10)) stop ("Lambda value should be in the range 1 to 10.")


    # Check the variable of e parameter
    if(sum(!(e %in% c(0.01, 0.001, 0.0001, 0.00001, 0.000001))) > 0)
      stop ("Select one of the possible e values: 0.01, 0.001, 0.0001, 0.00001 or 0.000001.")


    # ------------------------------------------ Input values ------------------------------------------------------------------------------------ #


    m <- ncol(weight_mat)


    # ------------------------------------------ Inference Rules  ------------------------------------------------------------------------------------ #


    mylist <- list()
    for(i in 1:(iter-1)) {

      if(i == 1) {
        if (infer == "k" || infer == "kc"){
          initial_vec <- colSums(t(activation_vec) * weight_mat)
        } else if  (infer == "mk" || infer == "mkc"){
          initial_vec <- activation_vec + colSums(t(activation_vec) * weight_mat)
        } else if (infer == "r" || infer == "rc"){
          initial_vec <- (2 * activation_vec - 1) + colSums(t((2 * activation_vec) - 1) * weight_mat)
        }

        if (transform == "s") {
          initial_vec <- 1/(1+exp(- lambda * initial_vec)) }
        if (transform == "t") {
          initial_vec <- tanh(lambda * initial_vec)
        }

      } else {
        # calculates the new vector (for the second until the last iteration or time step)
        if (infer == "k" || infer == "kc"){
          initial_vec <- colSums(t(initial_vec) * weight_mat)
        } else if  (infer == "mk" || infer == "mkc"){
          initial_vec <- initial_vec + colSums(t(initial_vec) * weight_mat)
        } else if (infer == "r" || infer == "rc"){
          initial_vec <- (2 * initial_vec - 1) + colSums(t((2 * initial_vec) - 1) * weight_mat)
        }

        if (transform == "s") {
          initial_vec <- 1/(1+exp(- lambda * initial_vec)) }
        if (transform == "t") {
          initial_vec <- tanh(lambda * initial_vec)
        }
      }

      if (transform == "b") {
        for(j in 1:m) {
          if (initial_vec[j] > 0){
            initial_vec[j] <- 1
          } else if (initial_vec[j] <= 0){
            initial_vec[j] <- 0
          }
        }
      }

      if (transform == "tr") {
        for(j in 1:m) {
          if (initial_vec[j] > 0){
            initial_vec[j] <- 1
          } else if (initial_vec[j] < 0){
            initial_vec[j] <- - 1
          } else initial_vec[j] <- 0
        }
      }

      if (infer == "kc" || infer == "mkc" || infer == "rc"){
        for(k in 1:m) {
          if(activation_vec[k] == 1) {
            initial_vec[k] <- (initial_vec[k] = 1)
          }
        }
      }
      mylist[[i]] <- initial_vec     # insert each produced stabilized vector in the list

   }



    steps_t <- (as.data.frame (do.call("rbind",mylist)))   # transform the produced stabilized vectors into a data frame
    step_1 <- as.numeric(activation_vec)

    # Insert the activation vector in the first row of the dataframe that contains the stabilized vectors of all time steps
    A <- (rbind(step_1, steps_t))
    last_conv <- as.double(A[iter,] - A[(iter-1),])   # check if the steady state has been reached of the last two iterations
    Res_e <- (length(last_conv[last_conv <= e]))    # Set the residual value (epsillon "e") equal to 0.001


    if ( Res_e < m)  {
      cat("\n WARNING: More iterations are required to reach the convergence.\n \n")
    } else {

      mylist1 <- list()
      for(i in 2:(iter)){
        subst <- abs(apply(A, 2, function(x) x[i] - x[i-1]))   # subtraction of "ith" - "(i-1)th" state
        mylist1[[i]] <- subst     # Save all subtraction vectors in a list
      }
      subst.mat <- do.call("rbind",mylist1)


      w <- as.data.frame(matrix(e, (iter - 1), m))    # Create a dataframe [(iterations - 1), m)] of values = epsillon


      mylist3 <- list()
      for(i in 1:(iter-1)){
        if(all(subst.mat[i,] < w[i,]))    # Check for the converged state
        {
          cv <- 1      # The concepts' value (cv) is converged
        }
        else {
          cv <- 2      # The concepts' value is NOT converged
        }
        mylist3[[i]] <- cv
      }
      cv.mat<-do.call("rbind",mylist3)


      conv_state <- min(which(cv.mat == 1))
      cat(sprintf("\n The concepts' values are converged in %ith state (e <= %f) \n", conv_state + 1, e))
      cat("\n")
      print(A[(conv_state + 1),], row.names = FALSE)
      cat("\n")
    }


     outlist <- list('values'= A)     # the concepts values in each state
    return (outlist)

 }

