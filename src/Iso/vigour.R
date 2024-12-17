if (!require("Iso")) {
    install.packages("Iso")
}
library(Iso)

matplot(vigour[, 1], vigour[, 2:6],
        main = "Growth vigour of stands of New Brunswick spruce",
        xlab = "year", ylab = "vigour", type = "b")
