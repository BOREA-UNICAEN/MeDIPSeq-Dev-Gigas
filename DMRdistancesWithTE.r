gen <- readRDS("matrixReloaded.data")

message("Extracting A...")
extractA <- data.frame(Feature = "",
                       DistA28c = 0, stringsAsFactors = FALSE)

# Extract data for A:
for (i in 1:length(gen))
{
 if(!is.null(gen[[names(gen)[i]]]$map$DistA28c[1]))
 {
  if(!is.infinite(gen[[names(gen)[i]]]$map$DistA28c[1]))
  {
   subs <- subset(gen[[names(gen)[i]]]$map, select = c(Feature, DistA28c))
  } else subs <- NULL
 } else subs <- NULL

 if (!is.null(subs)) extractA <- rbind(extractA, subs)
}

message("Extracting B...")
extractB <- data.frame(Feature = "",
                       DistBInt = 0, stringsAsFactors = FALSE)

# Extract data for B:
for (i in 1:length(gen))
{
 if(!is.null(gen[[names(gen)[i]]]$map$DistBInt[1]))
 {
  if(!is.infinite(gen[[names(gen)[i]]]$map$DistBInt[1]))
  {
   subs <- subset(gen[[names(gen)[i]]]$map, select = c(Feature, DistBInt))
  } else subs <- NULL
 } else subs <- NULL

 if (!is.null(subs)) extractB <- rbind(extractB, subs)
}

message("Extracting C...")
extractC <- data.frame(Feature = "",
                       DistCSpa = 0, stringsAsFactors = FALSE)

# Extract data for C:
for (i in 1:length(gen))
{
 if(!is.null(gen[[names(gen)[i]]]$map$DistCSpa[1]))
 {
  if(!is.infinite(gen[[names(gen)[i]]]$map$DistCSpa[1]))
  {
   subs <- subset(gen[[names(gen)[i]]]$map, select = c(Feature, DistCSpa))
  } else subs <- NULL
 } else subs <- NULL

 if (!is.null(subs)) extractC <- rbind(extractC, subs)
}

extractA <- extractA[-1, ]
extractB <- extractB[-1, ]
extractC <- extractC[-1, ]

# Clean up the feature levels:
greppata <- grep("LTR|LINE", extractA$Feature)
extractA$Feature[greppata] <- "TE-I"
greppata <- grep("LTR|LINE", extractB$Feature)
extractB$Feature[greppata] <- "TE-I"
greppata <- grep("LTR|LINE", extractC$Feature)
extractC$Feature[greppata] <- "TE-I"
greppata <- grep("DNA|RC", extractA$Feature)
extractA$Feature[greppata] <- "TE-II"
greppata <- grep("DNA|RC", extractB$Feature)
extractB$Feature[greppata] <- "TE-II"
greppata <- grep("DNA|RC", extractC$Feature)
extractC$Feature[greppata] <- "TE-II"
greppata <- grep("similarity|Repeat", extractA$Feature)
extractA$Feature[greppata] <- "REP"
greppata <- grep("similarity|Repeat", extractB$Feature)
extractB$Feature[greppata] <- "REP"
greppata <- grep("similarity|Repeat", extractC$Feature)
extractC$Feature[greppata] <- "REP"

# Clean intron and promoter names:
greppata <- grep("intrn", extractA$Feature); extractA$Feature[greppata] <- "INT"
greppata <- grep("intrn", extractB$Feature); extractB$Feature[greppata] <- "INT"
greppata <- grep("intrn", extractC$Feature); extractC$Feature[greppata] <- "INT"
greppata <- grep("prom", extractA$Feature); extractA$Feature[greppata] <- "PRO"
greppata <- grep("prom", extractB$Feature); extractB$Feature[greppata] <- "PRO"
greppata <- grep("prom", extractC$Feature); extractC$Feature[greppata] <- "PRO"

# Remove unknowns:
extractA <- extractA[ extractA$Feature != "Unknown" , ]
extractB <- extractB[ extractB$Feature != "Unknown" , ]
extractC <- extractC[ extractC$Feature != "Unknown" , ]

# Remove mRNA because they make chier:
extractA <- extractA[ extractA$Feature != "mRNA" , ]
extractB <- extractB[ extractB$Feature != "mRNA" , ]
extractC <- extractC[ extractC$Feature != "mRNA" , ]

# Prepare a single data frame to improve plotting:
distframe <- data.frame(
  Feature = c(extractA$Feature,extractB$Feature,extractC$Feature),
  Dist = c(extractA$DistA28c,extractB$DistBInt,extractC$DistCSpa),
  Comp = c(rep("A28c",nrow(extractA)),rep("BInt",nrow(extractB)),rep("CSpa",nrow(extractC))))

# Plot:
xp <- c(0, 0.8, 1.6, 2.4, 3.2, 4.0,
        5.6, 6.4, 7.2, 8, 8.8, 9.6,
        11.2, 12, 12.8, 13.6, 14.4, 15.2)
colvec <- rep(c("aquamarine3","cadetblue3","chartreuse3","darkorchid3",
                "deeppink3","deeppink3"), 3)
xnames <- c("CDS","INT","PRO","REP","TE-I","TE-II")

boxw <- boxplot(Dist ~ Feature + Comp, data = distframe, outline = FALSE, at = xp,
                ylab = "Distance to DMR (kbp)", notch = TRUE,
                axes = FALSE, frame.plot = TRUE, main = "",
                col = colvec,
                ylim = c(2.1e+04, 7.8e+05),
                par(mar = par("mar")),
                pars = list(whisklty = 1, whiskcol = "grey40",
                            staplelty = 1, staplecol = "grey40"))
#mtext(c("A","B","C"), side = 1, , line = 2.2, at = c(1.6, 7.2, 12.8))
axis(side = 1, at = xp, labels = rep(xnames, 3), las = 2)
axis(side = 2, las = 1, labels = c(0, 200, 400, 600, 800),
     at = c(0, 2e05, 4e05, 6e05, 8e05))

text(c(1.6, 7.2, 12.8), 770000, labels = c("A","B","C"), cex = 1.3)
#model <- kruskal.test(distances ~ block, data = work)

#pairm <- pairwise.wilcox.test(work$distances, work$block,
#                         p.adjust.method = "bonferroni")
