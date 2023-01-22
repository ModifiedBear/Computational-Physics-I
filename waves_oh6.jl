using LinearAlgebra
using CairoMakie
using BenchmarkTools

function circle_indices(X, Y, center::Tuple, radius)
  return [(x - center[1])^2 + (y - center[2])^2 <= radius^2 for x in X, y in Y]
end

function init(nx, ny)
  # initial values
  L = 5
  #  calculations
  #dx, dy = L/nx, L/ny
  xs = LinRange(-L,L,nx);
  ys = LinRange(-L,L,ny);

  c = 0.3;
  h = 0.5 ;# spatial width
  k = 0.5; # time step width
  
  # remember, higher precision means lower stability
  # choose h and k accordingly

  α = ones(nx, ny) .* ((c*k) / h)^2; # alpha squared ()
  #α[134:138, 1:45] .= 0;
  #α[134:138, 55:95] .= 0;
  #α[134:138, 105:end] .= 0;
  α[circle_indices(xs, ys, (-2,-2), 0.5)] .= 0
  κ = sqrt.(α) .* (k/h) .* 1.; # no tengo idea por qué hay que multiplicar por 8 para que salga bonito

  #α[134:138, 1:45] .= 0;
  #α[130:134,:] .= 0;
  #α[130:134,70:80] .= ((c*k) / h)^2;
  #α[134:138, 1:45] .= 0;
  #α[134:138, 55:95] .= 0;
  #α[134:138, 105:end] .= 0;
  α[100:120, 100:120] .= 0;
  #α[30:50, 20:60] .= 0;
  
  #Zₙ₊₁ = [exp(- 0.5 * (x)^2 - 0.5 * (y)^2) for x in xs, y in ys]
  Zₙ₊₁ = zeros(nx, ny);
  Zₙ   = zeros(nx, ny);
  Zₙ₋₁ = zeros(nx, ny);
  #data = zeros(N, nx+1, ny+1)
  
  return Zₙ, Zₙ₊₁, Zₙ₋₁, α, κ, xs, ys;
end

function oh6(vₙ, vₙ₊₁, α, κ, nx, ny,n)
  # better than the vectorized format
  bndry = 3;
  for ii in 1:n
    # you NEED to use copy(), this isn't python
    ω = 0.1
    #u_new[nx-2:nx-1,:] .= 10 .* sin(ω * ii)
    #M[ii,:,:] = copy(vₙ₊₁)
    #vₙ₊₁[nx-bndry-1:nx-bndry,div(ny,3):ny-div(ny,3)] .= 10 .* sin(ω * ii);
    vₙ₊₁[div(nx,2)-1:div(nx,2)+1,div(ny,2)-1:div(ny,2)+1] .= 10 .* cos(ω * ii);
    vₙ₋₁ = copy(vₙ);
    vₙ = copy(vₙ₊₁);
    for ii in 4:nx-3
      for jj in 4:ny-3
        vₙ₊₁[ii, jj] =α[ii,jj] * (  2  * vₙ[ii,jj-3]
                                  - 27 * vₙ[ii,jj-2]
                                  + 270* vₙ[ii,jj-1]
  
                                  + 2  * vₙ[ii-3,jj]
                                  - 27 * vₙ[ii-2,jj]
                                  + 270* vₙ[ii-1,jj]
                                  - 980* vₙ[ii,  jj]
                                  + 270* vₙ[ii+1,jj]
                                  - 27 * vₙ[ii+2,jj]
                                  + 2  * vₙ[ii+2,jj]
  
                                  + 270* vₙ[ii,jj+1]
                                  - 27 * vₙ[ii,jj+2]
                                  + 2  * vₙ[ii,jj+2]
                                  ) / 180 + 2 * vₙ[ii, jj] - vₙ₋₁[ii, jj];
        for kk in 1:bndry
          vₙ₊₁[kk, jj]     = vₙ[kk+1,   jj]   + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[kk+1,     jj]-vₙ[kk,     jj]); # x = 0 
          vₙ₊₁[nx-3+kk,jj] = vₙ[nx-4+kk , jj] + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[nx-4+kk,  jj]-vₙ[nx-3+kk,jj]); # x = N
          vₙ₊₁[ii, kk]     = vₙ[ii,  kk+1]    + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[ii,     kk+1]-vₙ[ii,     kk]); # y = 0
          vₙ₊₁[ii,ny-3+kk] = vₙ[ii, ny-4+kk]  + (κ[ii,jj]-1)/(κ[ii,jj]+1) * (vₙ₊₁[ii,  ny-4+kk]-vₙ[ii,ny-3+kk]);# y = N
        end
      end
    end
  end
  
  return vₙ₊₁
end  

function get_plot(data,medium, xs, ys)
  set_theme!(theme_black())
  custom_cmap=[RGBAf(0, 0, 0, 0.5) for c in to_colormap(:Greys)]
  fig=Figure()
  ax=Axis(fig[1,1])
  ax.aspect = DataAspect() 
  hm=Makie.heatmap!(ax,xs, ys, data,colormap=:batlow, interpolate=true)
  
  mask = ones(size(medium))
  mask[medium .> minimum(medium)] .= NaN;
  
  bl=Makie.heatmap!(ax, xs, ys, mask; colormap=custom_cmap, interpolate=false)
  Colorbar(fig[1,2],hm)
  rowsize!(fig.layout, 1, ax.scene.px_area[].widths[2]) # set colorbar height
  fig
end

begin
  N = 1000
  nx, ny = 150,150
  U, U_new, U_old, Alpha, Kappa, xs, ys = init(nx,ny);
  data = oh6(U, U_new, Alpha, Kappa,nx, ny, N);
end

get_plot(data, Alpha,xs,ys)