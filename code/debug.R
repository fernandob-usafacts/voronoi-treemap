
#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

# Code to help debug steps in the Voronoi Treemap algorithm

# Plot a tidied cell
# Use semitrans blue so can see backtracks
# Draw nice big blob at start
# Draw slightly smaller blob at end (so can see if got back to start)

drawTidyCell <- function(...) {
    UseMethod("drawTidyCell")
}

drawTidyCell.default <- function(cell, 
                                 lineCol=rgb(0,0,1,.5), 
                                 lineWidth=1,
                                 startFill=rgb(0,0,1,.5),
                                 endFill=rgb(1,0,0,.5),
                                 newpage=TRUE) {
    if (is.null(cell))
        return()
    if (newpage) {
        grid.newpage()
        pushViewport(viewport(width=.9, height=.9,
                              xscale=c(-2*scale, 2*scale),
                              yscale=c(-2*scale, 2*scale)))
    }
    grid.circle(cell$x[1], cell$y[1], r=unit(2, "mm"),
                default.units="native",
                gp=gpar(col=NA, fill=startFill))
    grid.lines(cell$x, cell$y,  
                default.units="native",
                gp=gpar(col=lineCol, lwd=lineWidth))
    N <- length(cell$x)
    grid.circle(cell$x[N], cell$y[N], r=unit(1, "mm"),
                default.units="native",
                gp=gpar(col=NA, fill=endFill))
}

drawTidyCell.multipleCells <- function(cell, 
                                       lineCol=rgb(0,0,1,.5), 
                                       lineWidth=1,
                                       startFill=rgb(0,0,1,.5),
                                       endFill=rgb(1,0,0,.5),
                                       newpage=TRUE) {
    drawTidyCell(cell[[1]],
                 lineCol, lineWidth, startFill, endFill, newpage)
    for (i in 2:length(cell)) {
        drawTidyCell(cell[[i]],
                     lineCol, lineWidth, startFill, endFill, newpage=FALSE)
    }
}

drawRoughCell <- function(cell,
                         lineCol=rgb(0,0,1,.5), 
                         lineWidth=1,
                         startFill=rgb(0,0,1,.5),
                         endFill=rgb(1,0,0,.5),
                         newpage=TRUE) {
    if (is.null(cell$border))
        return()
    if (newpage) {
        grid.newpage()
        pushViewport(viewport(width=.9, height=.9,
                              xscale=c(-2*scale, 2*scale),
                              yscale=c(-2*scale, 2*scale)))
    }
    grid.segments(cell$border[, 1], cell$border[, 2],
                  cell$border[, 3], cell$border[, 4],
                  default.units="native",
                  gp=gpar(col=lineCol, lwd=lineWidth))
}

drawTidyCells <- function(cells) {
    for (i in 1:length(cells)) {
        grid.newpage()
        pushViewport(viewport(width=.9, height=.9))
        pushViewport(viewport(xscale=c(-2*scale, 2*scale),
                              yscale=c(-2*scale, 2*scale)))
        grid.rect(gp=gpar(col=NA, fill="grey80"))
        for (j in 1:length(cells)) {
            drawTidyCell(cells[[j]], 
                         lineCol="white",
                         lineWidth=3,
                         startFill=NA, endFill=NA, 
                         newpage=FALSE)
        }
        drawTidyCell(cells[[i]], newpage=FALSE)
    }
}
    
drawRoughCells <- function(cells) {
    for (i in 1:length(cells)) {
        grid.newpage()
        pushViewport(viewport(width=.9, height=.9))
        pushViewport(viewport(xscale=c(-2*scale, 2*scale),
                              yscale=c(-2*scale, 2*scale)))
        grid.rect(gp=gpar(col=NA, fill="grey80"))
        for (j in 1:length(cells)) {
            drawRoughCell(cells[[j]], 
                         lineCol="white",
                         lineWidth=3,
                         startFill=NA, endFill=NA, 
                         newpage=FALSE)
        }
        drawRoughCell(cells[[i]], newpage=FALSE)
    }
}

drawTidyCellsOnePage <- function(cells) {
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(3, 4)))
    for (i in 1:length(cells)) {
        pushViewport(viewport(layout.pos.col=(i-1)%%4 + 1,
                              layout.pos.row=(i-1)%/%4 + 1))
        grid.rect()
        drawTidyCell(cells[[i]], newpage=FALSE)
        popViewport()
    }
    popViewport()
}

drawSites <- function(s, w, scol="grey", r=unit(.5, "mm"),
                      weights=TRUE, newpage=TRUE) {
    if (newpage) {
        grid.newpage()
        pushViewport(viewport(xscale=c(-1*scale, 1*scale),
                              yscale=c(-1*scale, 1*scale)))
    }
    grid.circle(s$x, s$y, r=r, 
                gp=gpar(col=NA, fill=scol),
                default.units="native")
    if (weights) {
        col <- ifelse(w < 0, "red", "black")
        grid.circle(s$x, s$y, r=abs(w), 
                    gp=gpar(col=col),
                    default.units="native")
    }
}

drawRegion <- function(region, gp=gpar(), newpage=TRUE) {
    if (newpage) {
        grid.newpage()
        pushViewport(viewport(xscale=c(-2*scale, 2*scale),
                              yscale=c(-2*scale, 2*scale)))
    }
    pts <- getpts(region)
    grid.polygon(pts$x,
                 pts$y,
                 gp=gp,
                 default.units="native")
}
