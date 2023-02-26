test.glmBoot <- function() {
    
    reserve <- glmBoot(test.triangle, 1e3)

    checkTrue(!any(is.na(reserve)))
}