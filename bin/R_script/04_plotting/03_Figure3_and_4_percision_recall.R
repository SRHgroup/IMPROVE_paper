
load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")

# --------------------------------------------
#               Figure 3  percicion recall curves 
# --------------------------------
## random forrest 


pred_nnalign <- prediction(Pred_Modelling$Prediction,Pred_Modelling$response)

# get_pr_dat <- function(model = "improve", predictions = Pred_Modelling$prediction_rf, labels  = Pred_Modelling$response) {
# pred <- prediction(predictions, labels)
# perf <- performance(pred, "prec", "rec")
# recall <- as.data.frame(perf@x.values)
# percision <- as.data.frame(perf@y.values)
# pr_dat <- cbind(recall,percision)
# colnames(pr_dat) <- c("recall","percision")
# pr_dat$model <- model
# return(pr_dat)
# }

improve_tme_dat <- get_pr_dat(model = "improve_tme", predictions = Pred_Modelling$prediction_rf_tme, labels  = Pred_Modelling$response)
improve_dat <- get_pr_dat(model = "improve", predictions = Pred_Modelling$prediction_rf, labels  = Pred_Modelling$response)
NNalign_dat <- get_pr_dat(model = "NNalign", predictions = Pred_Modelling$Prediction, labels  = Pred_Modelling$response)
rankel_dat <- get_pr_dat(model = "Rankel", predictions = Pred_Modelling$RankEL_minus, labels  = Pred_Modelling$response)
simple_dat <- get_pr_dat(model = "simple", predictions = pred_df_simple$prediction_rf, labels  = pred_df_simple$response)
auc(improve_dat$recall, improve_dat$precision)

pr_all_dat <- bind_rows(improve_dat,NNalign_dat) %>% bind_rows(.,rankel_dat ) %>% bind_rows(., improve_tme_dat) %>% bind_rows(.,simple_dat)
pr_all_dat <- na.omit(pr_all_dat)

pr_curve_all  <- ggplot(pr_all_dat , aes(x =recall, y = percision )) +
  geom_point(aes(color = model))  +
  geom_line(aes(color = model)) +
  scale_x_log10() + 
  scale_y_continuous(limits = c(0,0.3)) + 
  theme_bw()
ggsave(pr_curve_all, file = "results/PaperPlots/other/pr_curveall.png")

library(MLmetrics)
PRAUC(Pred_Modelling$prediction_rf_tme, Pred_Modelling$response)
PRAUC(Pred_Modelling$prediction_rf, Pred_Modelling$response)
PRAUC(Pred_Modelling$Prediction, Pred_Modelling$response)
PRAUC(Pred_Modelling$RankEL_minus, Pred_Modelling$response)
PRAUC(pred_df_simple$prediction_rf, pred_df_simple$response)





# --------------- old ------------------



perf_nnalign <- performance(pred_nnalign, "prec", "rec")
png(file="results/PaperPlots/other/pr_improve.png")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "... Precision/Recall graphs ...")
dev.off()
png(file="results/PaperPlots/other/pr_nnalign.png")
plot(perf_nnalign,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "... Precision/Recall graphs ...")
dev.off()
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)


library(yardstick)
truth <- rbinom(200, 1, 0.9)
pred_1 <- rbinom(200, 1, 0.8)
pred_2 <- rbinom(200, 1, 0.7)

dat <- data.frame(truth = as.factor(truth), pred_1 = as.numeric(pred_1), pred_2 = as.numeric(pred_2))

pr_curve(dat, truth, pred_1)


pr.curve(scores.class0=fg,scores.class1=bg)


library(plotly)
library(tidymodels)
library(ggplot2)

set.seed(0)
X <- matrix(rnorm(10000),nrow=500)
y <- sample(0:1, 500, replace=TRUE)
data <- data.frame(X,y)
data$y <- as.factor(data$y)
X <- subset(data,select = -c(y))
logistic_glm <-
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") %>%
  fit(y ~ ., data = data)

y_scores <- logistic_glm %>%
  predict(X, type = 'prob')

y_score <- y_scores$.pred_1
db <- data.frame(data$y, y_score)

z <- roc_curve(data = db, 'data.y', 'y_score')
z$specificity <- 1 - z$specificity
colnames(z) <- c('threshold', 'tpr', 'fpr')

fig1 <- plot_ly(x= y_score, color = data$y, colors = c('blue', 'red'), type = 'histogram', alpha = 0.5, nbinsx = 50) %>%
  layout(barmode = "overlay")
fig1


fig2 <- plot_ly(data = z, x = ~threshold) %>%
  add_trace(y = ~fpr, mode = 'lines', name = 'False Positive Rate', type = 'scatter')%>%
  add_trace(y = ~tpr, mode = 'lines', name = 'True Positive Rate', type = 'scatter')%>%
  layout(title = 'TPR and FPR at every threshold')
fig2 <- fig2 %>% layout(legend=list(title=list(text='<b> Rate </b>')))
fig2




library(plotly)
library(tidymodels)
library(fastDummies)



Pred_Modelling$response <- as.factor(Pred_Modelling$response)
pr_pred <- pr_curve(data = Pred_Modelling, response, prediction_rf)
aps_pred <- mean(Pred_Modelling$prediction_rf)
pred <- paste('pred (AP =',toString(round(aps_pred,2)),')',sep = '')

Pred_Modelling$response <- as.factor(Pred_Modelling$response)
pr_pred_tme <- pr_curve(data = Pred_Modelling, response, prediction_rf_tme)
aps_pred_tme <- mean(Pred_Modelling$prediction_rf_tme)
pred_tme <- paste('pred_tme (AP =',toString(round(aps_pred_tme,2)),')',sep = '')

max(Pred_Modelling$prediction_rf)
Pred_Modelling$response <- as.factor(Pred_Modelling$response)
pr_RankEL <- pr_curve(data = Pred_Modelling, response, RankEL_minus)
aps_pred_rankel <- mean(Pred_Modelling$RankEL_minus)
rankEL <- paste('rankel (AP =',toString(round(aps_pred_rankel,2)),')',sep = '')


max(Pred_Modelling$Prediction)
pr_NNalign <- pr_curve(data = Pred_Modelling, response, Prediction)
aps_pred_nnalign <- mean(Pred_Modelling$RankEL_minus)
rankEL <- paste('rankel (AP =',toString(round(aps_pred_nnalign,2)),')',sep = '')


pr_pred$model <- "IMPROVE"
pr_pred_tme$model <- "IMPROVE_TME"
pr_RankEL$model <- "RankEL"
pr_NNalign$model <- "NNalign"

percision_recall_all <- bind_rows(pr_pred,pr_pred_tme) %>% 
  bind_rows(.,pr_RankEL ) %>% 
  bind_rows(.,pr_NNalign )


ggplot(percision_recall_all , aes(x =recall, y = precision )) +
  geom_point(aes(color = model))  +
  geom_line(aes(color = model))  +
  scale_x_log10()



pr_curve(Pred_Modelling, response,prediction_rf)
roc_auc(Pred_Modelling, response,prediction_rf)
roc_auc(Pred_Modelling, response,Prediction)


# Create an empty figure, and add a new line for each class
fig <- plot_ly()%>%
  add_segments(x = 0, xend = 1, y = 1, yend = 0, line = list(dash = "dash", color = 'black'), showlegend = FALSE) %>%
  add_trace(data = pr_pred,x = ~recall, y = ~precision, mode = 'lines', name = pred, type = 'scatter')%>%
  add_trace(data = pr_RankEL,x = ~recall, y = ~precision, mode = 'lines', name = rankEL, type = 'scatter')%>%
  add_trace(data = pr_pred_tme,x = ~recall, y = ~precision, mode = 'lines', name = pred_tme, type = 'scatter')%>%
  layout(xaxis = list(
    title = "Recall"
  ), yaxis = list(
    title = "Precision"
  ),legend = list(x = 100, y = 0.5))
fig

data(two_class_example)

library(ggplot2)
library(dplyr)
pr_curve(Pred_Modelling, truth, Class1) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  coord_equal() +
  theme_bw()

