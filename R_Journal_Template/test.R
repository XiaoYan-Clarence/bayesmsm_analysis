# Generated by `rjournal_pdf_article()` using `knitr::purl()`: do not edit by hand
# Please edit test.Rmd to modify this file

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(plotly)
library(ggplot2)
library(palmerpenguins)
library(kableExtra)


## ----penguins-alison, out.width = "100%", out.height = "30%", fig.cap = "Artwork by \\@allison\\_horst", fig.alt="A picture of three different penguins with their species: Chinstrap, Gentoo, and Adelie. "----
knitr::include_graphics("figures/penguins.png")


## ----penguins-tab-interactive, eval = knitr::is_html_output(), layout = "l-body-outset"----
#> knitr::kable(head(penguins), format = "html", caption = "A basic table")


## ----penguins-tab-static, eval = knitr::is_latex_output()---------------------
knitr::kable(head(penguins), format = "latex", caption = "A basic table") %>% 
  kableExtra::kable_styling(font_size = 7)


## ----penguins-plotly, echo = TRUE, fig.height = 5, fig.cap="A basic interactive plot made with the plotly package on palmer penguin data. Three species of penguins are plotted with bill depth on the x-axis and bill length on the y-axis. When hovering on a point, a tooltip will show the exact value of the bill depth and length for that point, along with the species name.", include=knitr::is_html_output(), eval=knitr::is_html_output(), fig.alt = "A scatterplot of bill length against bill depth, both measured in millimetre. The three species are shown in different colours and loosely forms three clusters. Adelie has small bill length and large bill depth, Gentoo has small bill depth but large bill length, and Chinstrap has relatively large bill depth and bill length."----
#> p <- penguins %>%
#>   ggplot(aes(x = bill_depth_mm, y = bill_length_mm,
#>              color = species)) +
#>   geom_point()
#> ggplotly(p)


## ----penguins-ggplot, echo = TRUE, fig.height = 5, fig.cap="A basic non-interactive plot made with the ggplot2 package on palmer penguin data. Three species of penguins are plotted with bill depth on the x-axis and bill length on the y-axis. Visit the online article to access the interactive version made with the plotly package.", include=knitr::is_latex_output(), eval=knitr::is_latex_output()----
penguins %>% 
  ggplot(aes(x = bill_depth_mm, y = bill_length_mm, 
             color = species)) + 
  geom_point()

