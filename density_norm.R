# compute distance between two densities

density_dist <- function(d1, d2, Grid)
{
    Positive = (d1 > 1e-6) & (d2 > 1e-6)
    Grid = Grid[Positive]
    d1 = d1[Positive]
    d2 = d2[Positive]
    grid.size=length(Grid)
    first_integral=c()
    for(t in 1:grid.size){
      first_integral[t]=trapz(Grid,(log(d1/d1[t])-log(d2/d2[t]))^2)
    }
    return(trapz(Grid,first_integral)/2)
}

# find mode

getmode <- function(v)
{
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
