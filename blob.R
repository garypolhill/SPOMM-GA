#!/usr/bin/env Rscript

# blob.R
#
# This is a script to create a plot of the species distribution generated
# by a SPOMM-GA run.
#
# Alessandro Gimona, 26 November 2010

# Copyright (C) 2010  Macaulay Institute
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the Licence, or
# (at your option) any later version
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

require(lattice)
args <- commandArgs(TRUE)

if(length(args) != 3) {
  print("Usage: blob.R <species distribution file> <levels> <PDF>",
        quote = FALSE)
  q(status = 1)
}

input.file <- args[1]

levels <- as.numeric(args[2])
 
sp <- read.table(input.file, header = TRUE)
 
names(sp)[3] <- "prop"
 
z <- matrix(sp$prop, 50, 50)

x = y = 1:50

pdf(args[3], paper = "a4", width = 7, height = 7)
par(mai = c(1, 1, 2, 0))
 
levelplot(z ~ sp$x * sp$y, contour = FALSE, labels = FALSE,
          at = 0:levels / levels, col.regions = rev(terrain.colors(levels)),
          xlab = "x", ylab = "y", main = args[1])
q(status = 0)
