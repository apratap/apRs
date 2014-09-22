library(ggvis)
install.packages("ggvis")
install.packages("dplyr")

library("ggvis")
library("dplyr")

head(mtcars)

p <- ggvis(mtcars, x=~wt, y=~mpg)
layer_points(p)

mtcars %>% ggvis(x=~wt,y=~mpg) %>% layer_points()

mtcars %>% ggvis(x=~mpg, y=~disp) %>% mutate(disp=disp / 61.0237) %>% layer_points()

mtcars %>% ggvis(~mpg, ~disp, stroke=~vs) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, fill=~vs) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, size=~vs) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, shape=~factor(cyl)) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, fill:="red", stroke:="black") %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, size:=300, opacity:=0.4) %>% layer_points()
mtcars %>% ggvis(~mpg, ~disp, shape:= "cross") %>% layer_points()


mtcars %>% 
  ggvis(~wt, ~mpg, 
        size := input_slider(10, 100),
        opacity := input_slider(0, 1)
  ) %>% 
  layer_points()


mtcars %>% ggvis(~wt, ~mpg) %>% 
  layer_points() %>% 
  add_tooltip(function(df) df$wt)


