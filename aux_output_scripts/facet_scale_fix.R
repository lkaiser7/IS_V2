# override scale for custom ranges or breaks in ggplot2
# modifying "init_scales()" from facet class in ggplot2
# https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/ 

# function to override scale in select facet panel
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

# implement function above to override facet wrap
CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

# contructor function to implement all changes to plot
facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

# # EXAMPLE CODE #
# 
# # test data
# set.seed(123)
# test_large_values_data <- 1:9 %>% set_names() %>%
#   map_dfr(~data.frame(
#     x_variable = runif(1, 100, 20000) + runif(10, -100, 100),
#     y_variable = runif(10, 0, 1)
#   ), .id = "facet_name")
# 
# # store plot
# p_annoying_x_scale <- ggplot(test_large_values_data, aes(x_variable, y_variable)) +
#   geom_point() +
#   facet_wrap(~facet_name, scales = "free", ncol = 4)
# 
# # initial plot and scales
# p_annoying_x_scale
# 
# # adjusted scales
# p_annoying_x_scale +
#   facet_wrap_custom(~facet_name, scales = "free", ncol = 4, scale_overrides = list(
#     scale_override(1, scale_x_continuous(breaks = c(5750, 5900))),
#     scale_override(6, scale_x_continuous(breaks = c(17800, 17900)))
#   ))

### END FUNCTION OVERRIDE ###