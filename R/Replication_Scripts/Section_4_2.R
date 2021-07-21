# SECTION 4.2: Medical Data Classification with Structural Information

# ======================================== #
#        Set the working directory         #
# ======================================== #

if (file.exists("Reproduction_data.rda")){
    load('Reproduction_data.rda')
} else {
    print("Set the working directory to the folder that contains the reproduction data")
        print("setwd(PATH-TO-CURR-DIR)")
}    
# ======================================== #
#       Load the necessary packages        #
# ======================================== #

library(SQUIC)

if (!require("Matrix", quietly = TRUE)) {
    install.packages("Matrix")
}
library(Matrix)

if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
library(devtools)

if (!require("datamicroarray", quietly = TRUE)) {
    devtools::install_github('ramhiser/datamicroarray')
}
library(datamicroarray)

if (!require("caret", quietly = TRUE)) {
    install.packages("caret")
}
library(caret)

if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
library(ggplot2)
# ======================================== #

# ======================================== #
#             Linear DA Function           #
# ======================================== #
Linear_DA <- function(Data_Train, labels_train, 
                            Data_Test, labels_test, Theta_est)
{
        
        # Initialize  
        N_train      = length(labels_train);
        N_test       = length(labels_test);
        N_classes    = max(labels_train);
        p            = dim(Data_Train)[1];
        rho_k        = matrix(0, ncol = N_classes, nrow = N_test)
        
        for(k in 1:N_classes){
                # prior probability for class k
                prior      = length(which(labels_train == k))/N_train
                
                # sample mean for class k
                mu         = rowMeans(Data_Train[,which(labels_train == k)])
                mu_mat     = kronecker(matrix(1,1,ncol(Data_Test)), mu);                
                
                # x - \mu
                x_minus_mu = Data_Test - 0.5 * mu_mat
                
                # Linear Discriminant function
                rho        = diag(t(x_minus_mu) %*% Theta_est %*% mu_mat + log(prior))
                rho_k[,k]  = rho
        }
        
        # find the predicted classes
        predicted = c();
        for(i in 1:N_test){
                predicted[i]=which.max(rho_k[i,])
        }

        # Evaluate accuracy of labelling assignment
        Results = confusionMatrix(factor(predicted),factor(labels_test))
        
        return(list(Results = Results, predicted = predicted))
}

# datasets under consideration
datasets  = c('burczynski','yeoh','shipp','alon','sorlie')
# optimal scalar tuning parameters
lambdas   = c(0.7, 0.8, 0.8, 0.8, 0.5)
# parameters for the estimation of the matrix tuning parameter
eta       = 0.1; lambda_bias    = 0.95;
ACC_Bias  = vector(); ACC_No_Bias = vector(); time_Bias = vector(); time_No_Bias = vector();

# Load Data (n x p)
N_data = length(datasets);
data(list=datasets, package = 'datamicroarray')

for (i in 1:N_data){
set.seed(1991)
curr_dataset  = datasets[i];
Data          = eval(as.symbol(curr_dataset))$x;
Data          = data.matrix(Data, rownames.force = NA);
labels        = as.numeric(eval(as.symbol(curr_dataset))$y);

# Split data in test & train
trainIndex    = createDataPartition(labels, p = .7, list = FALSE)
trainIndex    = as.numeric(trainIndex)
Data_Train    = t(Data[ trainIndex,])
Data_Test     = t(Data[-trainIndex,])
labels_train  = labels[trainIndex]
labels_test   = labels[-trainIndex]

# Normalizing data
N_train = dim(Data_Train)[2];
a       = (apply(Data_Train,1,var))*(N_train-1)/N_train;
a       = 1/sqrt(a);
D       = Matrix::Diagonal(x=a);
iD      = Matrix::Diagonal(x=1/a);
Data_Train_scale = D%*%Data_Train;

# No Graphical Bias 
lambda_no_bias   = lambdas[i];
res_Null         = SQUIC(Y = Data_Train_scale, lambda = lambda_no_bias, M = NULL, verbose = 1)
time_No_Bias[i]  = res_Null$info_time_total;
Theta_Null       = D %*% res_Null$X %*% D;

# Linear Discriminant Analysis
LDA_Null       = Linear_DA(Data_Train, labels_train,
                        Data_Test, labels_test, Theta_Null)
ACC_No_Bias[i] = LDA_Null$Results$overall[1];

# Load the MST graphical bias
A_tree = eval(as.name(paste(as.symbol(curr_dataset), '_tree', sep="")));

# Tree Bias
M_tree        = eta*A_tree;
res_tree      = SQUIC(Y = Data_Train_scale, lambda = lambda_bias, M = M_tree, verbose = 1)
time_Bias[i]  = res_tree$info_time_total;
Theta_tree    = D %*% res_tree$X %*% D;

# Linear Discriminant Analysis
LDA_tree   = Linear_DA(Data_Train, labels_train,
                        Data_Test, labels_test, Theta_tree)
ACC_Bias[i] = LDA_tree$Results$overall[1];
}
# Plot time and ACC for all datasets
times_all = cbind(time_No_Bias,time_Bias);
ACC_all   = cbind(ACC_No_Bias,ACC_Bias);
rownames(times_all) <- datasets; rownames(ACC_all) <- datasets;

# Open the pdf file
pdf("Figure4a.pdf") 
barplot(t(times_all), beside=TRUE, font.axis=2,
        # main = "Figure 4a: Timings for the LDA experiments",
        xlab="Dataset", ylab="Time (sec)",
        # cex.axis=1.5, cex.names = 1.5, cex.main = 1.5, cex.lab=1.5,
        col=c("chartreuse4","brown4"),
        ylim=c(0,max(times_all)+20), legend=c("Scalar Tuning Parameter", "Matrix Tuning Parameter"),
        args.legend = list(x = "top", inset=c(-0.0, -0.00),
                           bty = "n", ncol = 2)
);
# Close the pdf file
dev.off()

# Open the pdf file
pdf("Figure4b.pdf") 
barplot(t(ACC_all),beside=TRUE, col=c("chartreuse4","brown4"),
        # main = "Figure 4b: Accuracy for the LDA experiments",
        font.axis=2, ylab="ACC", xlab="Dataset",
        legend=c("Scalar Sparsity Parameter", "Matrix Sparsity Parameter"),
        ylim=c(0.5,1), xpd = FALSE,
        args.legend = list(x = "top", inset=c(-0.0, -0.00),
                           bty = "n", ncol = 2)
        );
dev.off()


print('===================================')
print(paste('Sec 4.2: LDA with SQUIC'))
print('===================================')
LDA_Table = data.frame(datasets, lambdas, 
                      times_all, ACC_all);
print(LDA_Table, digits = 3)
