ctv2htmldebug <- function (x, file = NULL, cran = FALSE, css = NULL, packageURL = NULL, 
          reposname = "CRAN") 
{
  
  if (is.character(x)) 
    x <- read.ctv(x, cran = cran)
  if (is.null(file)) 
    file <- paste0(x$name, ".html")
  if (is.null(css) & cran) 
    css <- "../CRAN_web.css"
  if (is.null(x$url) & cran) 
    x$url <- paste0("https://CRAN.R-project.org/view=", x$name)
  if (is.null(packageURL)) {
    packageURL <- if (cran) 
      "../packages/"
    else "https://CRAN.R-project.org/package=%s"
  }
  ampersSub <- function(x) gsub("&", "&amp;", x)
  obfuscate <- function(x) paste(sprintf("&#x%x;", as.integer(sapply(unlist(strsplit(gsub("@", 
                                                                                          " at ", x), NULL)), charToRaw))), collapse = "")
  for (i in 1:length(x)) if (is.character(x[[i]])) 
    Encoding(x[[i]]) <- "unknown"
  title <- paste0(reposname, " Task View: ", ctv:::htmlify(x$topic))
  htm1 <- c("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">", 
            "<html xmlns=\"http://www.w3.org/1999/xhtml\">", "<head>", 
            paste0("  <title>", title, "</title>"), if (!is.null(css)) paste0("  <link rel=\"stylesheet\" type=\"text/css\" href=\"", 
                                                                              css, "\" />"), "  <meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\" />", 
            sprintf("  <meta name=\"citation_title\" content=\"%s\" />", 
                    title), sprintf("  <meta name=\"citation_author\" content=\"%s\" />", 
                                    ctv:::htmlify(x$maintainer)), sprintf("  <meta name=\"citation_publication_date\" content=\"%s\" />", 
                                                                    x$version), if (!is.null(x$url)) sprintf("  <meta name=\"citation_public_url\" content=\"%s\" />", 
                                                                                                             x$url), sprintf("  <meta name=\"DC.title\" content=\"%s\" />", 
                                                                                                                             title), sprintf("  <meta name=\"DC.creator\" content=\"%s\" />", 
                                                                                                                                             ctv:::htmlify(x$maintainer)), sprintf("  <meta name=\"DC.issued\" content=\"%s\" />", 
                                                                                                                                                                             x$version), if (!is.null(x$url)) sprintf("  <meta name=\"DC.identifier\" content=\"%s\" />", 
                                                                                                                                                                                                                      x$url), "</head>", "", "<body>", paste0("  <h2>", 
                                                                                                                                                                                                                                                              reposname, " Task View: ", ctv:::htmlify(x$topic), "</h2>"), 
            paste0("  <table summary=\"", x$name, " task view information\">"), 
            paste0("    <tr><td valign=\"top\"><b>Maintainer:</b></td><td>", 
                   ctv:::htmlify(x$maintainer), "</td></tr>"), if (!is.null(x$email)) paste0("    <tr><td valign=\"top\"><b>Contact:</b></td><td>", 
                                                                                       obfuscate(x$email), "</td></tr>"), paste0("    <tr><td valign=\"top\"><b>Version:</b></td><td>", 
                                                                                                                                 ctv:::htmlify(x$version), "</td></tr>"), if (!is.null(x$url)) paste0("    <tr><td valign=\"top\"><b>URL:</b></td><td><a href=\"", 
                                                                                                                                                                                                ctv:::htmlify(x$url), "\">", ctv:::htmlify(x$url), "</a></td></tr>"), 
            if (!is.null(x$source)) paste0("    <tr><td valign=\"top\"><b>Source:</b></td><td><a href=\"", 
                                           ctv:::htmlify(x$source), "\">", ctv:::htmlify(x$source), "</a></td></tr>"), 
            "  </table>")
  htm2 <- x$info
  pkg2html <- if (grepl("%s", packageURL, fixed = TRUE)) {
    function(a, b) paste0("    <li><a href=\"", sprintf(packageURL, 
                                                        a), "\">", a, "</a>", if (b) 
                                                          " (core)"
                          else "", "</li>")
  }
  else {
    function(a, b) paste0("    <li><a href=\"", packageURL, 
                          a, "/index.html\">", a, "</a>", if (b) 
                            " (core)"
                          else "", "</li>")
  }
  htm3 <- c(paste0("  <h3>", reposname, " packages:</h3>"), 
            "  <ul>", sapply(1:NROW(x$packagelist), function(i) pkg2html(x$packagelist[i, 
                                                                                       1], x$packagelist[i, 2])), "  </ul>")
  htm4 <- c("  <h3>Related links:</h3>", "  <ul>", sapply(x$links, 
                                                          function(x) paste0("    <li>", x, "</li>")), "  </ul>")
  if (!is.null(x$otherlinks)) {
    htm4 <- c(htm4, "", "  <h3>Other resources:</h3>", "  <ul>", 
              sapply(x$otherlinks, function(x) paste0("    <li>", 
                                                      x, "</li>")), "  </ul>")
  }
  print(class(htm1))
  print(class(htm2))
  print(class(htm3))
  print(class(htm4))
  if(is.list(htm1))
    stop("header section is a list and not a vector")
  if(is.list(htm2))
    stop("body is a list and not a vector")
  if(is.list(htm3))
    stop("package list section is a list and not a vector")
  if(is.list(htm4))
    stop("links section is a list and not a vector")
  
  htm <- c(htm1, "", htm2, "", htm3, "", htm4, "", "</body>", 
           "</html>")
  #print(head(htm))
  htm.len <- sapply(htm, length)
  print(htm[htm.len > 1])
  print(table(htm.class <- sapply(htm, class)))
  stopifnot(all(inherits(htm.class, "character")))
  writeLines(htm, con = file)
  invisible(htm)
}

ctv2htmldebug("Distributions.md")
