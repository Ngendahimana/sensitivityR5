
data(toy)
psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family = binomial(), 
    data = toy)
toy$ps <- psmodel$fitted
toy$linps <- psmodel$linear.predictors

X <- toy$linps  ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
Y <- toy$out3.time
match1 <- Match(Y = Y, Tr = Tr, X = X, M = 1, replace = FALSE, ties = FALSE)

match.it <- matchit(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, data = toy, 
    method = "nearest", ratio = 1)



sensSurv(data = toy, exp = "treated", outcome = "out2", failtime = "out3.time", Gamma = 2.4, 
    Gammainterval = 0.01, alpha = 0.05, object_name = match.it, plot_title = "Time To Event Outcome Sensitivity Plot")

# test code
data = toy
exp = "treated"
outcome = "out2"
failtime = "out3.time"
Gamma = 2.4
Gammainterval = 0.01
alpha = 0.05
object_name = match.it
plot_title = "Time To Event Outcome Sensitivity Plot"


sensSurv = function(data = toy, object_name = match1, exp, outcome, failtime, Gamma, alpha, 
    Gammainterval, plot_title = NULL) {
    results = list()
    
    if (class(object_name) == "Match") {
        extractor = (c(object_name$index.treated, object_name$index.control))
        group_id = c(c(1:length(object_name$index.treated)), c(1:length(object_name$index.control)))
        data1 = data[extractor, ]
        data1$match = group_id
        
        
    } else if (class(object_name) == "matchit") {
        
        data$rId = row.names(data)
        data2 = data[, c("rId", exp)]
        # data2$rId = row.names(data2)
        names(data2)[2] = "exp"
        t_id = data2$rId[data2$exp == 1]
        c_id = object_name$match.matrix
        extractor = as.numeric(c(t_id, c_id))
        data1 = data[extractor, ]
        k2 = length(extractor)/2
        match = rep(1:k2, 2)
        data1$match = match
        
        
    } else {
        print("Matching functions not known")
    }
    
    
    data1s = subset(data1, select = c("match", exp, outcome, failtime))
    names(data1s)[c(2:4)] = c("exp", "outcome", "failtime")
    data2s = subset(data1s, exp == 0)
    data3s = subset(data1s, exp == 1)
    data4s = subset(data2s, select = c("match", "exp", "outcome", "failtime"))
    data5s = subset(data3s, select = c("match", "exp", "outcome", "failtime"))
    data6s = full_join(data4s, data5s, by = "match")
    
    names(data6s)[c(4, 7)] = c("failtimeNotexp", "failtimeExp")
    data6s$timediff = data6s$failtimeNotexp - data6s$failtimeExp
    wonpairs = sum(data6s$timediff != 0)
    expoutlive = sum(data6s$timediff < 0)
    
    results$wonpairs = sum(data6s$timediff != 0)
    results$expoutlive = sum(data6s$timediff < 0)
    
    gamVal = seq(1, Gamma, by = Gammainterval)
    pplus = 1/(1 + gamVal)
    pminus = gamVal/(1 + gamVal)
    
    table1 = data.frame(cbind(gamVal, pplus, pminus))
    table1$expTplus = wonpairs * table1$pplus
    table1$expTminus = wonpairs * table1$pminus
    table1$sd_expT = sqrt(wonpairs * table1$pplus * (1 - table1$pplus))
    
    
    
    for (i in 1:length(gamVal)) {
        table1$pupper[i] = round(min(1, 2 * pnorm((expoutlive - table1[i, 5])/table1[i, 6], 
            lower.tail = FALSE)), 4)
    }
    
    for (i in 1:length(gamVal)) {
        table1$plower[i] = round(min(1, 2 * pnorm((expoutlive - table1[i, 4])/table1[i, 6], 
            lower.tail = FALSE)), 4)
    }
    
    
    table1$min = abs(alpha - table1$pupper)
    
    
    
    vrt = table1[table1$min == min(table1$min), ]$gamVal
    hrz = table1[table1$min == min(table1$min), ]$pupper
    vrt1 = round(table1[table1$min == min(table1$min), ]$gamVal, 2)
    hrz1 = round(table1[table1$min == min(table1$min), ]$pupper, 2)
    
    
    plot = ggplot(data = table1, aes(x = gamVal, y = pupper)) + geom_line() + geom_point(aes(x = vrt, 
        y = hrz)) + ylab("p upper bound") + xlab("gamma (Bias)") + theme_bw() + annotate("text", 
        x = vrt + 0.1 * vrt, y = hrz, label = paste0("(", vrt1, ",", hrz1, ")")) + labs(title = plot_title, 
        caption = paste("matching done by", class(object_name), "function"))
    
    results$plot = plot
    results$upperbound_pval = hrz = table1[table1$min == min(table1$min), ]$pupper
    results$Gamma = table1[table1$min == min(table1$min), ]$gamVal
    
    return(results)
    
    
}


