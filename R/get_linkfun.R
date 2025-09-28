
get.linkfun <- function () {
  expression({
    if (is.character(link)) {
      linkf <- catch.conditions ({
        get(link, mode = "function", envir = parent.frame())
      })$value

      if (any(class(linkf) %in% c("simpleError", "error",
                                  "condition", "try-error"))) {
        switch (link,
                logit = {
                  link <- logit() # msbreg::logit()
                },
                probit = {
                  link <- probit() # msbreg::probit()
                },
                cloglog = {
                  link <- cloglog() # msbreg::cloglog()
                },
                cauchit = {
                  link <- cauchit() # msbreg::cauchit()
                },
                {
                  print(link)
                  stop("Error: 'link' not recognized", call. = FALSE)
                })
      }
      else {
        link <- linkf
      }
    }

    if (is.function(link)) {
      if (identical(link, boot::logit)) {
        link <- logit() # msbreg::logit()
      }
      else {
        linkf <- catch.conditions ({
          link()
        })$value

        if (any(class(linkf) %in% c("simpleError", "error",
                                    "condition", "try-error"))) {
          print(link)
          stop("Error: 'link' not recognized", call. = FALSE)
        }

        link <- linkf
      }
    }

    if (is.null(link$link)) {
      print(link)
      stop("'link' not recognized")
    }
  })
}
