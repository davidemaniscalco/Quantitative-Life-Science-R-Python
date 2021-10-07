source('./Config.R')

# *******************************************************************************************************
# Builds ONE item of NN elements with probability ff that each element is 1
# *******************************************************************************************************
Build_Item <- function(NN,ff){
    return(sample(0:1,NN,replace=TRUE,prob=c(1-ff,ff)))
}

# *******************************************************************************************************
# Builds a (LL x NN) matrix containing LL items of length NN. Still ff is the prob
# *******************************************************************************************************
Build_Item_Matrix <- function(NN,LL,ff){
    return(matrix(sample(0:1,NN*LL,replace=TRUE,prob=c(1-ff,ff)),nrow=LL,ncol=NN,byrow=TRUE))
}

# *****************************************************************************************************************
# Calculates the Overlap between two items. It is a scalar product
# ****************************************************************************************************************
Overlap <- function(x,y){
    return(x %*% y)
}

# ******************************************************************************************************************
# Returns a (LL x LL) symmetric matrix, whose elements represents overlaps.
#     first it generates LL items of length NN with prob ff. Then it calculates all the overlaps 
#     between its elements and put them in the matrix. If diagonal=FALSE the self-overlaps are
#     put to 0.
# ******************************************************************************************************************
Build_Symmetric_Similarity_Matrix <- function(NN,LL,ff,diagonal){
    x <- Build_Item_Matrix(NN,LL,ff)
    SM <- x %*% t(x)
    if(!diagonal){
        diag(SM) <- 0
    }
    return(SM)
}

# *******************************************************************************************************************
# Returns a (LL x LL) symmetric matrix, whose elements represents overlaps.
#     first it generates TWO*LL items of length NN with prob ff, to be put in TWO different item matrices.
#     Then it calculates all the overlaps between the elements of these 2 matrices and put them in the matrix. 
#     If diagonal=FALSE the self-overlaps are put to 0.
# *******************************************************************************************************************
Build_Asymmetric_Similarity_Matrix <- function(NN,LL,ff,diagonal){
    x <- Build_Item_Matrix(NN,LL,ff)
    y <- Build_Item_Matrix(NN,LL,ff)
    ASM <- x %*% t(y)
    if(!diagonal){
        diag(ASM) <- 0
    }
    return(ASM)
}

# *********************************************************************************************************************
# Returns a (LL x LL) symmetric matrix with 0 elements on the diagonal. The other elements are integer 
# sampled from start to end
# *********************************************************************************************************************
old_SAMPLE_Simmetric_Matrix <- function(start,end,LL){
    if(start==0){
        warning('WARNING: You have chosen start=0')
    }
    vec = sample(start:end,LL*(LL-1)/2,replace=TRUE)
    mat <- matrix(0,nrow=LL,ncol=LL)
    mat[lower.tri(mat)] <- vec
    mat <- t(mat)
    mat[lower.tri(mat)] <- vec
    return(mat)
}

# **********************************************************************************************************************
SAMPLE_Simmetric_Matrix <- function(LL){
    vec = runif(LL*(LL-1)/2)
    mat <- matrix(0,nrow=LL,ncol=LL)
    mat[lower.tri(mat)] <- vec
    mat <- t(mat)
    mat[lower.tri(mat)] <- vec
    return(mat)
}

# **********************************************************************************************************************
SAMPLE_Asymmetric_Matrix <- function(LL){
    mat <- matrix(runif(LL**2),nrow=LL,ncol=LL)
    diag(mat) <- 0
    return(mat)
}

## **********************************************************************************************************************
# Inputs: a similarity matrix, either symmetric or asymmetric, and the starting item.
# Starts with starts_item. The next recalled item will be the one with the maximum overlap, i.e. the item
# corrsponding to the maximum value of the start_item row of that matrix. All the reaclled items are
# saved in the Recalled_Items vector. The process is repeated until an element is recalled twice. Then a vector
# containing all the recalled items, each one once, is returned
## *********************************************************************************************************************
Simple_Recall_Process <- function(Similarity_Matrix,start_item){
    diag(Similarity_Matrix) <- 0                                       #fundamental to not have self-overlaps
    max_vector = apply(Similarity_Matrix,FUN='which.max',MARGIN=1)
    
    Recalled_Items <- start_item
    next_item <- max_vector[tail(Recalled_Items,1)]
    Recalled_Items <- c(Recalled_Items,next_item)

    while(!any((tail(Recalled_Items,1)==Recalled_Items[1:(length(Recalled_Items)-1)]))){
        next_item <- max_vector[tail(Recalled_Items,1)]
        Recalled_Items <- c(Recalled_Items,next_item)
    }

    Recalled_Items <- Recalled_Items[-length(Recalled_Items)]
    
    return(Recalled_Items)
}

## ***********************************************************************************************************************
# Starts as the Simple_Recall_Process. Then are calculated the max_vector and the second_max_vector, containing
# respectively the element with the maximum overlap and the element with the second max overlap for each row.
# Then it is build the (2 x LL) happened_transitions matrix, full of zeros. When a transition of the
# 'first type' happens TO the element ii, then the element ii of the FIRST row of this matrix is set to 1.
# Exactly the same with the second row, with the transitions of the SECOND type.
# Loop. Next item is calculated using max_vector, but if it is exactly the previous item of the recalled list,
# then I calculated next item using max_second_vector. For each happened transition, both of the 1st or of the
# 2nd type, I put +1 in the happened_transition matrix. As soon as the happened_transition matrix contains a 2,
# then this means that I have entered in a loop, and the process stops.
## ***********************************************************************************************************************
Double_Recall_Process <- function(Similarity_Matrix,start_item){
    
    diag(Similarity_Matrix) <- 0                                       #fundamental to not have self-overlaps
    
    # build max_vector and also max_second_vector
    max_vector = apply(Similarity_Matrix,FUN='which.max',MARGIN=1)
    copy = t(Similarity_Matrix)
    copy[max_vector+seq(0,nrow(copy)**2-nrow(copy),nrow(copy))] = 0    #I leave all the maximum
    copy = t(copy)
    max_second_vector = apply(copy,FUN='which.max',MARGIN=1)
    
    # keep tracked, by putting zeros, happened transitions
    happened_transitions = matrix(0,nrow=2,ncol=ncol(Similarity_Matrix))
    
    # first iteration
    Recalled_Items <- start_item
    next_item <- max_vector[tail(Recalled_Items,1)]
    happened_transitions[1,tail(Recalled_Items,1)] <- happened_transitions[1,tail(Recalled_Items,1)] + 1
    Recalled_Items <- c(Recalled_Items,next_item)
   
    # A transition has happened twice <-> one element of this matrix is =2
    while(all(happened_transitions!=2)){
        next_item <- max_vector[tail(Recalled_Items,1)]
        if(next_item==Recalled_Items[length(Recalled_Items)-1]){    #if I had to come back, I go to 2nd max
            next_item <- max_second_vector[tail(Recalled_Items,1)]
            happened_transitions[2,tail(Recalled_Items,1)] <- happened_transitions[2,tail(Recalled_Items,1)] + 1
            #ii = ii + 1
        }else{
            happened_transitions[1,tail(Recalled_Items,1)] <- happened_transitions[1,tail(Recalled_Items,1)] + 1
        }
        Recalled_Items <- c(Recalled_Items,next_item)
    }
    
    return(unique(Recalled_Items)) #unique is needed, the extra transitions are contained in the vec andmustbe left
}


###########################################################################################################
# COMPUTE RECALLED MATRIX
###########################################################################################################
Compute_Recalled_Matrix <- function(Items_List, Repetitions, NN, ff, symmetric, single, verbose, each,with_ff){

    Recalled_Matrix <- matrix(0, nrow=Repetitions, ncol=length(Items_List))
    start_time <- Sys.time()
    if(with_ff){
        if(!symmetric){
            if(!single){
                cat('ASYMMETRIC \n DOUBLE \n ff = ',ff,'\n')
                flush.console()
                for(jj in seq(Repetitions)){
                    ii=1
                    if(verbose){
                        if(jj %% each == 0){
                            cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                            flush.console()
                        }
                    }
                    for(LL in Items_List){
                        ASM <- Build_Asymmetric_Similarity_Matrix(NN,LL,ff,FALSE)
                        Recalled_Matrix[jj,ii] <- length(unique(Double_Recall_Process(ASM,starting_item)))
                        ii <- ii + 1
                    }
                }

            }else{
                cat('ASYMMETRIC \n SIMPLE \n ff = ',ff,'\n')
                flush.console()
                for(jj in seq(Repetitions)){
                    ii=1
                    if(verbose){
                        if(jj %% each == 0){
                            cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                            flush.console()
                        }
                    }
                    for(LL in Items_List){
                        ASM <- Build_Asymmetric_Similarity_Matrix(NN,LL,ff,FALSE) 
                        Recalled_Matrix[jj,ii] <- length(unique(Simple_Recall_Process(ASM,starting_item)))
                        ii <- ii + 1
                    }
                }
            }
        }else{
            if(!single){
                cat('SYMMETRIC \n DOUBLE \n ff = ',ff,'\n')
                flush.console()
                for(jj in seq(Repetitions)){
                    ii=1
                    if(verbose){
                        if(jj %% each == 0){
                            cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                            flush.console()
                        }
                    }
                    for(LL in Items_List){
                        ASM <- Build_Symmetric_Similarity_Matrix(NN,LL,ff,FALSE)
                        Recalled_Matrix[jj,ii] <- length(unique(Double_Recall_Process(ASM,starting_item)))
                        ii <- ii + 1
                    }
                }

            }else{
                cat('SYMMETRIC \n SIMPLE \n ff = ',ff,'\n')
                flush.console()
                for(jj in seq(Repetitions)){
                    ii=1
                    if(verbose){
                        if(jj %% each == 0){
                            cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                            flush.console()
                        }
                    }
                    for(LL in Items_List){
                        ASM <- Build_Symmetric_Similarity_Matrix(NN,LL,ff,FALSE) 
                        Recalled_Matrix[jj,ii] <- length(unique(Simple_Recall_Process(ASM,starting_item)))
                        ii <- ii + 1
                    }
                }
            }
        }
    }else{
        if(symmetric){
            cat('RANDOM SAMPLED SYMMETRIC MATRIX \n')
            flush.console()
            for(jj in seq(Repetitions)){
                ii=1
                if(verbose){
                    if(jj %% each == 0){
                        cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                        flush.console()
                    }
                }
                for(LL in Items_List){
                    mat <- SAMPLE_Simmetric_Matrix(LL)
                    Recalled_Matrix[jj,ii] <- length(unique(Double_Recall_Process(mat,starting_item)))
                    ii <- ii + 1
                }
            }
        }else{
            if(!single){
                cat('RANDOM SAMPLED ASYMMETRIC MATRIX, DOUBLE \n')
                flush.console()
                for(jj in seq(Repetitions)){
                    ii=1
                    if(verbose){
                        if(jj %% each == 0){
                            cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                            flush.console()
                        }
                    }
                    for(LL in Items_List){
                        mat <- SAMPLE_Asymmetric_Matrix(LL)
                        Recalled_Matrix[jj,ii] <- length(unique(Double_Recall_Process(mat,starting_item)))
                        ii <- ii + 1
                    }
                }
            }else{
                cat('RANDOM SAMPLED ASYMMETRIC MATRIX, SIMPLE \n')
                flush.console()
                for(jj in seq(Repetitions)){
                    ii=1
                    if(verbose){
                        if(jj %% each == 0){
                            cat('Repetition:',jj,'  Time:',Sys.time()-start_time,'\n')
                            flush.console()
                        }
                    }
                    for(LL in Items_List){
                        mat <- SAMPLE_Asymmetric_Matrix(LL)
                        Recalled_Matrix[jj,ii] <- length(unique(Simple_Recall_Process(mat,starting_item)))
                        ii <- ii + 1
                    }
                }
            }
            
        }
        
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    
    return(Recalled_Matrix)
}

#############################################################################################################
#############################################################################################################
#                                                 MAIN
#############################################################################################################

Recalled_Matrix <- Compute_Recalled_Matrix(Items_List, Repetitions, NN, ff, symmetric, single, verbose, each, with_ff)
write.table(Recalled_Matrix,filename,row.names=FALSE,col.names=FALSE)
warnings()
