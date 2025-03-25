################
# Handling zero
################

handling_zero_forecast <- function(DJI_return_density_train_trans)
{
    DJI_return_density_train_transformation = DJI_return_density_train_trans * (10^6)
    n_1 = ncol(DJI_return_density_train_transformation)
    epsilon = sapply(1:n_1, function(X) max(DJI_return_density_train_transformation[,X] - round(DJI_return_density_train_transformation[,X], 2)))
    
    DJI_CoDa_mat = matrix(NA, m, n_1)
    for(ik in 1:n_1)
    { 
        index = which(round(DJI_return_density_train_transformation[,ik], 2) == 0)
        if(length(index) == 0)
        {
            DJI_CoDa_mat[,ik] =  DJI_return_density_train_transformation[,ik]/(10^6)
        }
        else
        {
            DJI_CoDa_mat[,ik] = replace(DJI_return_density_train_transformation[,ik], index, epsilon[ik])
            DJI_CoDa_mat[-index,ik] = (DJI_return_density_train_transformation[-index,ik] * (1 - (length(index) * epsilon[ik])/(10^6)))/(10^6)
        }
    }    
    return(DJI_CoDa_mat)
}
