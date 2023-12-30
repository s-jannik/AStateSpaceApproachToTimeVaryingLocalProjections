# Author: Jannik Steenbergen

############## Various functions for data cleaning and manipulation ############

# fct. for data cleaning
preprocess_data <- function(data, name){
  # Arguments
    # name: name of dataset. Either "fred_qd" or "eurozone"
  
  if (name == "eurozone"){
    
    # Remove a set of variables from the dataset 
    data <- data[,colnames(data) %in% c("X", "HICP", "UR", "YC3M")]   # where "X" is dates
    
    # create quarterly sequence
    q_seq <- seq(as.Date("1990-01-01"), by="quarter", length.out = 132)
    
    # convert quarterly sequence to preferred formatting
    q_seq <- paste(as.POSIXlt(q_seq)$year+1900, quarters(q_seq))
    
    # Remove column with dates and add rownames
    data <- data[,-1]
    row.names(data) <- q_seq
  }
  
  if (name == "fred_qd"){
    
    # Remove header
    data <- data[-c(1,2), ]
    
    # Remove last obs (NA)
    data <- data[-nrow(data), ]
    
    # Extract data and rename columns 
    data <- data[,c("sasdate", "CPIAUCSL", "FEDFUNDS", "UNRATE")]
    colnames(data) <- c("sasdate", "CPI", "IR", "UR")    # IR: interest rate
    
    # create quarterly sequence
    q_seq <- seq(as.Date("1959-03-01"), by="quarter", length.out = 259)
    
    # convert quarterly sequence to preferred formatting
    q_seq <- paste(as.POSIXlt(q_seq)$year+1900, quarters(q_seq))
    
    # Remove column with dates and add rownames
    data <- data[,-1]
    row.names(data) <- q_seq
    
  }
  
  return(data)
}

# A data cleaning function:
get_dataset <- function(data, maxlags){
  # Description: 
  # takes the 3-dimensional dataset and a vector maxlags = c(maxlags_cpi, maxlags_ir, maxlags_ur)
  
  T <- nrow(data)
  maxlag <- max(maxlags)
  
  new_data <- matrix(NA, nrow=(T-maxlag), ncol=(3+sum(maxlags)))
  
  new_data <- data[(maxlag+1):T, c("CPI", "IR")]
  new_data <- cbind(new_data, rep(1, nrow(new_data)))
  colnames(new_data) = c("CPI", "IR", "ones")
  
  print(new_data)
  
  for (var in 1:length(maxlags)){
    
    if (maxlags[var] != 0){
      for (i in 1:maxlags[var]){
        if (i != maxlag) lag <- data[(maxlag+1-i):(T-i),var]
        if (i == maxlag) lag <- data[1:(T-i),var]
        
        colnames <- c(colnames(new_data), paste(colnames(data)[var], "l", i, sep=""))
        new_data <- cbind(new_data, lag)
        colnames(new_data) <- colnames
        
      }
    }
  }
  return(new_data) 
}


lag_func <- function(y, lag_order){
  
  y <- c(rep(NA, lag_order), y[1:(length(y)-lag_order)])
  
  return(y)
}

get_prediction_data <- function(data, max_lag, h){
  # description:
  # returns prediction data
  # inputs:
  # data: data as dataframe
  # max_lag: maximum lag order
  # h: forecasting horizon  
  
  pred_data <- data.frame(row.names=rownames(data))
  
  # define dependent variable (first variable in data per default)
  pred_data[,1] <-data[,1]
  
  for (lag in 1:max_lag){
    
    for (v in colnames(data)){
      
      v_name <- paste(v, "_l", lag, sep="")
      
      pred_data[,v_name] <- lag_func(y=data[,v], lag_order=(lag+(h-1)))
      
    }
  }
  
  return(pred_data)
}


################### Function for stationarity transformations ##################

# fct performing stationarity transformations 
transform_data_tcode <- function(data, tcode, use_which_data="all_available", as_ts_object){
  # Arguments:
  # data: dataframe with qseq as rownames 
  # tcode: a named vector of tcodes for each variable
  # use_which_data: either "all_available" or "largest_balanced"
  # "all_available" uses the largest dataset available for each variable. The size differs across variables
  # "largest_balanced" uses the length of the shortest TRANSFORMED variable for all TRANSFORMED variables.
  # this means that inflation uses more observations than other series since it is year-on-year difference, 
  # but transformed inflation is just as long as the shortest transformed variable. 
  # Output
  # transformed data as ts object 
  
  
  # transform data to ts object
  start <- as.numeric(strsplit(gsub("\\D", " ", rownames(data)[1]), "  ")[[1]])
  start_year <- start[1]
  start_quarter <- start[2]
  
  ts_data <- ts(data, frequency=4, start=c(start_year,start_quarter))
  
  # initialize empty ts object for transformed data
  ts_t_data <- ts(frequency=4, start=start(ts_data))
  
  for (var in colnames(ts_data)){
    
    if (tcode[var] == 1){
      
      t_var <- ts_data[,var]                # transformed variable
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 2){
      
      t_var <- diff(ts_data[,var])
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 3){
      
      t_var <- diff(ts_data[,var], differences=2)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 4){
      
      t_var <- log(ts_data[,var])
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 5){
      
      t_var <- diff(log(ts_data[,var]))
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 6){
      
      t_var <- diff(log(ts_data[,var]), differences=2)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 7){
      
      t_var <- diff(ts_data[,var]/lag(ts_data[,var]) - 1)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 8){
      
      t_var <- diff(log(ts_data[,var]))*400
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    
    if (tcode[var] == 9){
      
      t_var <- diff(log(ts_data[,var]), lag=4)*100
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
    if (tcode[var] == 10){
      
      t_var <- diff(diff(log(ts_data[,var]), lag=4)*100)
      ts_t_data <- cbind(ts_t_data, t_var)
      
    }
  }
  
  # remove support column
  ts_t_data <- ts_t_data[,-1]
  
  colnames(ts_t_data) <- colnames(ts_data)
  
  # return
  if(as_ts_object == TRUE) {
    return(ts_t_data)
  }
  
  if(as_ts_object == FALSE){
    
    t_data <- as.data.frame(ts_t_data)
    rownames(t_data) <- rownames(data)
    return(t_data)
  }
  
  
}
