# define the igraph class as S4
setClass("igraph")

# define the SpatialLinesNetwork class
setClass("SpatialLinesNetwork", representation(sl = "SpatialLines", g = "igraph", 
                                               nb = "list"), validity = function(object) {
                                                 stopifnot(length(object@sl) == length(E(object@g)))
                                                 stopifnot(length(object@nb) == length(V(object@g)))
                                               })

# creates a SpatialLines object, and converts it into a SpatialLinesNetwork
SpatialLinesNetwork = function(sl) {
  stopifnot(is(sl, "SpatialLines"))
  if (!is(sl, "SpatialLinesDataFrame")) 
    sl = new("SpatialLinesDataFrame", sl, data = data.frame(id = 1:length(sl)))
  if (!all(sapply(sl@lines, length) == 1)) 
    stop("SpatialLines is not simple: each Lines element should have only a single Line")
  startEndPoints = function(x) {
    firstLast = function(L) {
      cc = coordinates(L)[[1]]
      rbind(cc[1, ], cc[nrow(cc), ])
    }
    do.call(rbind, lapply(x@lines, firstLast))
  }
  s = startEndPoints(sl)
  zd = zerodist(SpatialPoints(s))
  pts = 1:nrow(s)
  
  # the following can't be done vector-wise, there is a progressive effect:
  if (nrow(zd) > 0) {
    for (i in 1:nrow(zd)) pts[zd[i, 2]] = pts[zd[i, 1]]
  }
  stopifnot(identical(s, s[pts, ]))
  
  # map to 1:length(unique(pts))
  pts0 = match(pts, unique(pts))
  node = rep(1:length(sl), each = 2)
  nb = lapply(1:length(unique(pts)), function(x) node[which(pts0 == x)])
  g = graph(pts0, directed = FALSE)  # edges
  nodes = s[unique(pts), ]
  g$x = nodes[, 1]  # x-coordinate vertex
  g$y = nodes[, 2]  # y-coordinate vertex
  g$n = as.vector(table(pts0))  # nr of edges
  # line lengths:
  sl$length = sapply(sl@lines, function(x) LineLength(x@Lines[[1]]))
  E(g)$weight = sl$length
  # create list with vertices, starting/stopping for each edge?  add for
  # each SpatialLines, the start and stop vertex
  pts2 = matrix(pts0, ncol = 2, byrow = TRUE)
  sl$start = pts2[, 1]
  sl$end = pts2[, 2]
  new("SpatialLinesNetwork", sl = sl, g = g, nb = nb)
}


# The following function converts a set of points into a SpatialPolygons or into a SpatialLines object:
dd <- function(x, ..., to = "polygons") {
  stopifnot(is(x, "Spatial"))
  cc = coordinates(x)
  dd = deldir(list(x = cc[, 1], y = cc[, 2]), ...)
  if (to == "polygons") {
    tm = triMat(dd)
    fn = function(ix) {
      pts = tm[ix, ]
      pts = c(pts, pts[1])
      Polygons(list(Polygon(rbind(cc[pts, ]))), ix)
    }
    SpatialPolygons(lapply(1:nrow(tm), fn), proj4string = x@proj4string)
  } else if (to == "lines") {
    segs = dd$delsgs
    lst = vector("list", nrow(segs))
    for (i in 1:nrow(segs)) lst[i] = Lines(list(Line(cc[c(segs[i, 5], segs[i, 
                                                                           6]), ])), i)
    SpatialLines(lst, proj4string = x@proj4string)
  } else stop("argument to should be polygons or lines")
}