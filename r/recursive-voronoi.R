library("deldir")
library("gpclib")
library("surveillance")
library("jsonlite")
library("purrr")
library("erer")
library("rmutil")
library("data.table")

source("voronoi.R")
source("kaleidescope.R")
source("util.R")
source("draw.R")

get_size <- function(element, idx_year) {
  if ('children' %in% names(element)) {
    size <- sum(unlist(lapply(element$children, function(x) get_size(x, idx_year))))
  } else if ('size' %in% names(element)){
    size <- element$size[idx_year]
  } else {
    size <- 0
  }
  if (size == 'NULL') {
    size <- 0
  }
  size
}

recursive_match <- function(categories, treemap, idx_year, level=1) {
  objs <- c()
  for (c_idx in 1:length(categories)) {
    category <- categories[[c_idx]]
    obj <- list()
    category_name <- category$name
    obj$name <- category_name
    obj$size <- get_size(category, idx_year)
    obj$level <- level
    if ('chg' %in% names(category)){
      if (is.null(category$chg[[idx_year]])) {
        obj$monthly_change <- NA
      } else {
        obj$monthly_change <- category$chg[[idx_year]]
      }
    } else {
      obj$monthly_change <- NA
    }
    if ('chg2' %in% names(category)){
      if (is.null(category$chg2[[idx_year]])) {
        obj$yearly_change <- NA
      } else {
        obj$yearly_change <- category$chg2[[idx_year]]
      }
    } else {
      obj$yearly_change <- NA
    }
    idx <- match(category_name, treemap$names)[[1]]

    if (is.na(idx)) {
      obj$target <- NA
      obj$vertices <- NA
      obj$x <- NA
      obj$y <- NA
    } else {
      obj$target <- treemap$t[[idx]]
      obj$x <- treemap$s$x[[idx]]
      obj$y <- treemap$s$y[[idx]]
      if (treemap$a[[idx]] == 0) {
        obj$vertices <- NA
      } else {
        x_coords <- get.pts(treemap$k[[idx]])[[1]]$x
        y_coords <- get.pts(treemap$k[[idx]])[[1]]$y
        obj$vertices <- c()
        for (i in 1:length(x_coords)) {
          aux <- list()
          aux$x <- x_coords[[i]]
          aux$y <- y_coords[[i]]
          obj$vertices <- c(obj$vertices, list(aux))
        }
      }
    }
    obj$leaf <- !('children' %in% names(category))

    objs <- c(objs, list(obj))

    if ('children' %in% names(category)) {
      if(length(category$children) > 0) {
        sub_objs <- recursive_match(category$children, treemap, idx_year, level=level+1)
        objs <- c(objs, sub_objs)
      }
    }
  }
  objs
}

calculate_treemap <- function(categories, center, region, extent, idx_year) {
  labels <- map(categories, function(x) x$name)
  temp <- unlist(map(categories, function(x) get_size(x, idx_year)))

  nCells <- length(labels)
  filename <- Reduce(function(a,b) paste(a,b, sep='-'), labels)
  filepath <- paste("temp/", filename, ".cin")

  treemap <- NULL
  attempt <- 0
  while (is.null(treemap) && attempt <= 200) {
    attempt <- attempt + 1
    print(paste("ATTEMPT ", attempt))
    print(paste(labels))
    if (file.exists(filepath)) {
      data <- fread(filepath)
      x <- data$x
      y <- data$y
      w <- data$w
      if (attempt == 10) {
        unlink(filepath)
      }
    } else {
      w <- runif(nCells, 1, 100)

      # extent <- 10
      dx <- extent * runif(nCells) - extent/2
      dy <- extent * runif(nCells) - extent/2
      x <- center$x + dx
      y <- center$y + dy
    }
    target <- temp/sum(temp)
    tryCatch({
      treemap <- allocate(labels, list(x=x, y=y), w, region, target, debug=TRUE, maxIteration=100)
      for (idx_treemap in 1:length(treemap$t)) {
        # if ((treemap$t[[idx_treemap]] == 0 && treemap$a[[idx_treemap]] != 0) ||
        if (treemap$t[[idx_treemap]] != 0 && treemap$a[[idx_treemap]] == 0) {
              treemap <- NULL
              break
            }
      }
      
    },
    warning = function(warn) {
     print(paste("WARNING: ", warn))
     treemap <- NULL
     unlink(filepath)
    },
    error = function(err) {
     print(paste("ERROR ", filename, ": ", err))
     treemap <- NULL
     unlink(filepath)
    },
    finally = function(f) {
     continue
    })
  }

  if (is.null(treemap)) {
    print(" TREEMAP IS NULL ")
    stop()
  }
  # if (!file.exists(filepath)) {
  df <- data.frame(x=treemap$s$x, y=treemap$s$y, w=treemap$w)
  fwrite(df, file=filepath)
  # }
  treemap
}

merge_treemaps <- function(treemap1, treemap2) {
  namess <- c(treemap1$names, treemap2$names)
  ks <- c(treemap1$k, treemap2$k)
  xs <- c(treemap1$s$x, treemap2$s$x)
  ys <- c(treemap1$s$y, treemap2$s$y)
  ss <- list(x=xs, y=ys)
  ws <- c(treemap1$w, treemap2$w)
  as <- c(treemap1$a, treemap2$a)
  ts <- c(treemap1$t, treemap2$t)

  merged_treemap <- list(names=namess, k=ks, s=ss, w=ws, a=as, t=ts)
  merged_treemap
}

recursive_treemap <- function(categories, center, region, extent, idx_year, level=1) {
  if (length(categories) == 1){
    ss <- list(x=center[[1]], y=center[[2]])
    as <- pi * extent**2
    treemap <- list(names=c(categories[[1]]$name), k=c(region), s=c(ss), w=c(1), a=c(as), t=c(1))
  } else {
    treemap <- calculate_treemap(categories, center, region, extent, idx_year)
  }

  for (idx in 1:length(categories)) {
    category <- categories[[idx]]
    if ('children' %in% names(category)) {
      print(category$name)
      if(length(category$children) > 0) {
        center <- list(x=treemap$s$x[idx], y=treemap$s$y[idx])
        radius <- sqrt(treemap$a[[idx]]/pi)
        sub_treemap <- recursive_treemap(category$children, center, treemap$k[[idx]], radius, idx_year, level=level+1)
        treemap <- merge_treemaps(treemap, sub_treemap)
      }
    }
  }
  treemap
}

result <- fromJSON('input/test3b-age.json', simplifyVector=FALSE)

center <- list(x=500, y=500)
disc <- discpoly(c(500,500), 3/2*500, npoly=64, class = "gpc.poly")
categories <- result[[1]]$children

nMonths <- length(result[[1]]$children[[1]]$children[[1]]$children[[1]]$children[[1]]$children[[1]]$chg)

for (idx_year in 19:21) {
  tryCatch({
    filename <- paste("output/age/treemap", idx_year, ".json", sep="")
    print(filename)
    if (file.exists(filename)) {
      next
    }
    treemap <- recursive_treemap(categories, center, disc, 500, idx_year, level=1)
    categories_match <- recursive_match(categories, treemap, idx_year)
    json_data <- toJSON(categories_match, pretty=TRUE)
    write(json_data, filename)
    folder_name = paste('output/age/temp', idx_year)
    dir.create(folder_name)
    file.copy("temp", folder_name, recursive=TRUE)
  },
  warning = function(warn) {
   print(paste("WARNING: ", warn))
  },
  error = function(err) {
   print(paste("ERROR ", filename, ": ", err))
  },
  finally = function(f) {
   continue
  })

}

# drawRegions(treemap, label=FALSE)
