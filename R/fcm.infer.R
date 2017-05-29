# Fuzzy Cognitive Maps (FCMs) Inference
# Provides a selection of 6 different inference rules and 4 threshold functions in order to obtain the inference of the FCM 
# (Fuzzy Cognitive Map). Moreover, the 'fcm' package returns a data frame of the concepts' values of each state after the inference 
# procedure. Fuzzy cognitive maps were introduced by Kosko (1986) providing ideal causal cognition tools for modeling and simulating 
# dynamic systems. 

# activation_vec: 1 x m data frame which contains the initial concept values. A concept is turned on or activated by making its vector element 1 or 0 or in [0, 1].
# weight_mat: m x m data frame which stores the weights assigned to the pairs of concepts. The weights are usually normalized to the interval [0, 1] or [−1, +1].
# iter: The required number of iterations in order to reach the FCM convergence. Defaults to 20.
# infer: Select an Inference Rule ('k' Kosko, 'mk' modified Kosko, 'r' Rescale,'kc' Kosko-clamped, 'mkc' modified Kosko-clamped or 'rc' Rescale-clamped). Default value is set to 'k'
# transform: Contains the Transformation functions ('b' Bivalent,  'tr' Trivalent,  's' Sigmoid or 't' Hyperbolic tangent). The transformation function is used to reduce unbounded weighted sum to a certain range, which hinders quantitative analysis, but allows for qualitative comparisons between concepts. Default value is set equal to 's'.
# lambda: A parameter that determines the steepness of the sigmoid and hyperbolic tangent function at values around 0. Different lambda value may perform more appropriate for different problems.
# e: Epsilon (e) is a residual, describing the minimum error difference among the subsequent concepts. Its value depends on the application type. Defaults to to 0.001.

## Returns a iter x m data frame which contains the concepts' values of each iteration after the the transformation function.

#####################################################################################################################################################################################################################
# author Zoumpolia Dikopoulou <dikopoulia@gmail.com>, <zoumpolia.dikopoulou@uhasselt.be>
# author Elpiniki Papageorgiou <epapageorgiou@teiste.gr>, <e.i.papageorgiou75@gmail.com>

# References 
# B. Kosko, "Fuzzy cognitive maps", International Journal of Man-Machine Studies 24, p.p. 65-75, 1986.
# Groumpos, P.P, Stylios, C.D.; "Modelling supervisory control systems using fuzzy cognitive maps", Chaos, Solitons & Fractals, Volume 11, Issues 1–3, p.p. 329–336, 2000.
# Papageorgiou E.I., "Fuzzy Cognitive Maps for Applied Sciences and Engineering From Fundamentals to Extensions and Learning Algorithms", Intelligent Systems Reference Library, Volume 54, 2014.
# Papageorgiou E.I., Stylios C.D., GroumposP.P. , "Unsupervised learning techniques for finetuning fuzzy cognitive map causal links.", Int. J. Human Comput. Stud. Vol. 64, pp. 727–743, 2006.
######################################################################################################################################################################################################################

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

